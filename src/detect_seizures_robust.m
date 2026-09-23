function seizure_results = detect_seizures_robust(data, cfg)
% DETECT_SEIZURES_ROBUST  Parallel seizure-detection branch for degraded
% signal (attenuated gain and/or 50Hz line contamination), selected via
% cfg.cases.profiles.<case>.seizure_mode = 'robust'. src/detect_seizures.m
% (the 'legacy' branch) is NEVER called or modified by this file -- see
% README.md "dos ramas" and KNOWN_ISSUES.md.
%
% Motivation (measured on animal 005, see 005-s/ground_truth.csv): the
% legacy detector applies its 15s minimum-duration filter to each
% threshold-crossing segment BEFORE any merging, so a real seizure that
% dips briefly below threshold gets fragmented into sub-15s pieces and
% silently discarded (this is what lost GT1/GT4 in the legacy baseline).
% Also, the amplitude-energy metric alone cannot separate a sustained-
% moderate false positive from a real, irregular, briefly-extreme seizure.
% Hysteresis was tried and explicitly rejected: it chains interictal-spike
% trains into one giant (~1000s) false event with no recall benefit.
%
% STAGE 1 (energy): numerically identical to detect_seizures.m, via the
%   shared seizure_energy_trace.m utility (same normalization, per-block
%   bandpass, per-block Hilbert envelope^power, movmean, global
%   median*factor threshold -- every parameter from cfg.seizure.*).
% STAGE 2 (line length): line_length_trace.m on the same band-passed
%   signal, same trimmed block extents (line-length/line-length-median-
%   global is scale-invariant, see that file).
% STAGE 3 (event construction, NO hysteresis of any kind):
%   mask = energy > threshold; merge crossings separated by a gap
%   <= cfg.seizure_robust.merge_gap_s, bounded so a merge never pushes an
%   event's span past cfg.seizure_robust.max_duration_s (a naturally long
%   single crossing, i.e. not a merge product, is NOT truncated); THEN
%   apply cfg.seizure_robust.min_duration_s. Order matters: merge-then-
%   duration is what recovers a fragmented real seizure that the legacy
%   duration-then-nothing order discards. Events over max_duration_s are
%   KEPT and flagged over_max_duration=true, never silently dropped.
% STAGE 4 (line-length filter, the single most important requirement):
%   ll_ratio = MEDIAN(line length within the event) / ll_median_global.
%   Never max/percentile -- a single spike would dominate either and
%   erase exactly the separation the median gives between real seizures
%   (sustained irregularity) and sustained-moderate false positives.
%   Kept iff ll_ratio >= cfg.seizure_robust.ll_threshold (default 1.75).
%
% Confidence metrics (ll_ratio, peak_energy_ratio, hf_ratio_db,
% envelope_cv) are attached to every surviving candidate for traceability
% and human review. They are INFORMATIVE ONLY -- never used as exclusion
% filters here or anywhere downstream. In particular hf_ratio_db (high-
% frequency power [cfg.seizure_robust.hf_band] relative to
% [cfg.seizure_robust.hf_reference_band]) partly reflects movement
% artifact; gating on it would make the detector blind to non-convulsive
% electrographic seizures, so it is reported, never filtered (see
% KNOWN_ISSUES.md).
%
%   seizure_results = detect_seizures_robust(data, cfg)
%
%   data : struct from load_lfp_txt.m, loaded from a *_clean.txt file
%   cfg  : struct from pipeline_config.m
%
% OUTPUT (struct seizure_results):
%   .seizures  : table of STAGE-4 survivors (id, start_s, end_s,
%                duration_s, start_abs, end_abs, block_id,
%                over_max_duration, ll_ratio, peak_energy_ratio,
%                hf_ratio_db, envelope_cv) -- same core columns as
%                detect_seizures.m's .seizures, plus confidence columns.
%   .candidates: table of EVERY Stage-3 survivor (same columns plus
%                candidate_id, kept), i.e. matched AND line-length-
%                rejected candidates, each tagged with its own ll_ratio.
%                Feed this to evaluate_detections.m's opts.candidates_table
%                to sweep ll_threshold without re-running detection.
%   .metrics   : median_energy, threshold, pct_above (Stage 1, identical
%                definition to detect_seizures.m), plus n_candidates,
%                n_dropped_short, n_dropped_ll, n_kept, n_over_max.
%   .qc        : n_blocks_total, n_blocks_rejected_short, rejected_blocks
%                (identical shape to detect_seizures.m's .qc).
%   .fs, .file

    signal = data.signal(:);
    fs = data.fs;
    t_rel = data.t_rel(:);
    session_start = data.session_start;

    trace = seizure_energy_trace(data, cfg);
    trimmed_blocks = trimmed_blocks_of(trace.blocks);
    [ll_full, ll_median_global] = line_length_trace(trace.bp_full, trimmed_blocks, fs, cfg);

    merge_gap_samples = round(cfg.seizure_robust.merge_gap_s * fs);
    max_duration_samples = round(cfg.seizure_robust.max_duration_s * fs);
    min_duration_s = cfg.seizure_robust.min_duration_s;
    ll_threshold = cfg.seizure_robust.ll_threshold;

    hf_ok = cfg.seizure_robust.hf_band(2) < fs / 2;

    candidate_rows = cell(0, 13);
    n_dropped_short = 0;
    cand_id = 0;

    for b = 1:numel(trace.blocks)
        blk = trace.blocks(b);
        above = blk.energy > trace.threshold;
        raw_segs = mask_to_segments(above);
        merged = merge_bounded(raw_segs, merge_gap_samples, max_duration_samples);

        for k = 1:size(merged, 1)
            ls = merged(k, 1);
            le = merged(k, 2);
            duration_s = (le - ls + 1) / fs;
            if duration_s < min_duration_s
                n_dropped_short = n_dropped_short + 1;
                continue;
            end

            gs = blk.trimmed_start + ls - 1;
            ge = blk.trimmed_start + le - 1;
            over_max = duration_s > cfg.seizure_robust.max_duration_s;

            energy_seg = trace.energy_full(gs:ge);
            ll_seg = ll_full(gs:ge);
            ll_ratio = median(ll_seg, 'omitnan') / ll_median_global;
            peak_energy_ratio = max(energy_seg, [], 'omitnan') / trace.threshold;
            envelope_cv = std(energy_seg, 0, 'omitnan') / mean(energy_seg, 'omitnan');

            if hf_ok
                raw_seg = signal(gs:ge);
                hf_power = bandpower(raw_seg, fs, cfg.seizure_robust.hf_band);
                ref_power = bandpower(raw_seg, fs, cfg.seizure_robust.hf_reference_band);
                hf_ratio_db = 10 * log10(hf_power / ref_power);
            else
                hf_ratio_db = NaN;
            end

            kept = ll_ratio >= ll_threshold;
            cand_id = cand_id + 1;
            start_s = t_rel(gs);
            end_s = t_rel(ge);

            candidate_rows(end+1, :) = {cand_id, start_s, end_s, duration_s, ...
                abs_time_or_nat(start_s, session_start), abs_time_or_nat(end_s, session_start), ...
                blk.block_id, over_max, ll_ratio, peak_energy_ratio, hf_ratio_db, envelope_cv, kept}; %#ok<AGROW>
        end
    end

    if isempty(candidate_rows)
        candidates = empty_candidate_table(session_start);
    else
        candidates = cell2table(candidate_rows, 'VariableNames', ...
            {'candidate_id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', ...
             'block_id', 'over_max_duration', 'll_ratio', 'peak_energy_ratio', 'hf_ratio_db', 'envelope_cv', 'kept'});
    end

    kept_mask = candidates.kept;
    seizures = candidates(kept_mask, {'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', ...
        'block_id', 'over_max_duration', 'll_ratio', 'peak_energy_ratio', 'hf_ratio_db', 'envelope_cv'});
    seizures.id = (1:height(seizures))';
    seizures = seizures(:, [end, 1:end-1]);

    n_seizures = height(seizures);
    fprintf('detect_seizures_robust: %s -> %d seizure(s) (median_energy=%.3e, threshold=%.3e, %d/%d blocks usable, %d candidate(s), %d dropped-short, %d dropped-ll)\n', ...
        data.file, n_seizures, trace.median_energy, trace.threshold, numel(trace.blocks), size(trace.raw_blocks, 1), ...
        height(candidates), n_dropped_short, sum(~kept_mask));

    n_merged_segments = n_dropped_short + height(candidates); % Stage-3 output, before the min_duration_s filter
    seizure_results = struct();
    seizure_results.seizures = seizures;
    seizure_results.candidates = candidates;
    seizure_results.metrics = struct( ...
        'median_energy', trace.median_energy, 'threshold', trace.threshold, 'pct_above', trace.pct_above, ...
        'n_candidates', height(candidates), 'n_dropped_short', n_dropped_short, ...
        'n_dropped_ll', sum(~kept_mask), 'n_kept', n_seizures, 'n_over_max', sum(seizures.over_max_duration), ...
        'n_segments', n_merged_segments, 'n_rejected', n_merged_segments - n_seizures);
        % n_segments/n_rejected: aliases matching detect_seizures.m's field names/semantics
        % (segments found before the duration filter / total not present in the final
        % .seizures table), so build_seizure_summary_row (run_pipeline_edf.m) needs no
        % seizure_mode branch. n_segments here is POST-merge (Stage 3), since "a segment"
        % for this branch is a merged crossing, not a raw pre-merge one.
    seizure_results.qc = struct( ...
        'n_blocks_total', size(trace.raw_blocks, 1), 'n_blocks_rejected_short', height(trace.rejected_blocks), ...
        'rejected_blocks', trace.rejected_blocks);
    seizure_results.fs = fs;
    seizure_results.file = data.file;

    %% ---- output paths --------------------------------------------------
    out_dir = cfg.seizure.output_dir;
    if isempty(out_dir)
        out_dir = pwd;
    end
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    [~, fname, ~] = fileparts(data.file);
    base = strrep(fname, '_clean', '');

    % Remove this channel's own seizure output files from any PREVIOUS run
    % before writing new ones -- otherwise, if an earlier run (legacy or
    % robust, at a different case/threshold) found more seizures than this
    % one, the extra numbered files would be silently orphaned (not
    % overwritten, since fewer files are written this time) and the folder
    % would show stale, higher counts than seizures_events.csv reports.
    % Matches this base name only, so it never touches another channel's
    % files; a wildcard delete() with zero matches is a silent no-op.
    delete(fullfile(out_dir, [base '_seizure*']));

    mat_file = '';
    fig_file = '';
    png_file = '';

    if n_seizures > 0
        mat_file = fullfile(out_dir, [base '_seizures.mat']);
        save(mat_file, 'seizure_results');

        try
            gap_blocks_t = gap_segments_seconds(data.valid_mask(:), t_rel);
            [fig_file, png_file] = save_panorama_figure_robust(out_dir, base, cfg, data.meta, ...
                t_rel, signal, trace.bp_full, trace.energy_full, trace.threshold, ...
                ll_full, ll_median_global, seizures, gap_blocks_t, n_seizures);

            save_zoom_figures_robust(out_dir, base, cfg, t_rel, signal, trace.bp_full, trace.energy_full, ...
                trace.threshold, ll_full, ll_median_global, seizures, gap_blocks_t);
        catch ME
            warning('detect_seizures_robust:FigureSaveFailed', ...
                'Could not generate/save seizure figures for %s: %s. Seizure detection results are unaffected.', ...
                data.file, ME.message);
        end
    end

    seizure_results.mat_file = mat_file;
    seizure_results.fig_file = fig_file;
    seizure_results.png_file = png_file;
end

%% ======================================================================
function trimmed = trimmed_blocks_of(blocks)
    if isempty(blocks)
        trimmed = zeros(0, 2);
    else
        trimmed = [vertcat(blocks.trimmed_start), vertcat(blocks.trimmed_end)];
    end
end

function merged = merge_bounded(raw_segs, merge_gap_samples, max_duration_samples)
% Merge consecutive crossing segments separated by a gap <=
% merge_gap_samples, UNLESS doing so would push the merged span past
% max_duration_samples -- in which case the merge stops there (a new
% event starts) rather than silently exceeding the cap. A single raw
% segment that is already longer than max_duration_samples on its own
% (no merge involved) passes through untouched; it is flagged, never
% truncated, by the caller's over_max_duration check.
    if isempty(raw_segs)
        merged = zeros(0, 2);
        return;
    end
    merged = zeros(0, 2);
    cur_s = raw_segs(1, 1);
    cur_e = raw_segs(1, 2);
    for i = 2:size(raw_segs, 1)
        gap = raw_segs(i, 1) - cur_e - 1;
        candidate_span = raw_segs(i, 2) - cur_s + 1;
        if gap <= merge_gap_samples && candidate_span <= max_duration_samples
            cur_e = raw_segs(i, 2);
        else
            merged(end+1, :) = [cur_s, cur_e]; %#ok<AGROW>
            cur_s = raw_segs(i, 1);
            cur_e = raw_segs(i, 2);
        end
    end
    merged(end+1, :) = [cur_s, cur_e];
end

function abs_t = abs_time_or_nat(t_rel_s, session_start)
    if isnat(session_start)
        abs_t = NaT;
        abs_t.TimeZone = session_start.TimeZone;
    else
        abs_t = rel_to_abs_time(t_rel_s, session_start);
    end
end

function T = empty_candidate_table(session_start)
    T = table('Size', [0 13], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'datetime', 'datetime', ...
                           'double', 'logical', 'double', 'double', 'double', 'double', 'logical'}, ...
        'VariableNames', {'candidate_id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', ...
                           'block_id', 'over_max_duration', 'll_ratio', 'peak_energy_ratio', 'hf_ratio_db', ...
                           'envelope_cv', 'kept'});
    T.start_abs.TimeZone = session_start.TimeZone;
    T.end_abs.TimeZone = session_start.TimeZone;
end

%% ======================================================================
% Figure helpers below are an independent implementation for this branch
% (NOT calls into detect_seizures.m, which is private/local to that file
% and must not be modified or imported from). Visual layout intentionally
% mirrors detect_seizures.m's panorama/zoom figures for reviewer
% continuity, plus a 4th line-length panel, the robust branch's defining
% signal.

function gap_blocks_t = gap_segments_seconds(valid_mask, t_rel)
    gap_idx = mask_to_segments(~valid_mask);
    if isempty(gap_idx)
        gap_blocks_t = zeros(0, 2);
    else
        gap_blocks_t = [t_rel(gap_idx(:, 1)), t_rel(gap_idx(:, 2))];
    end
end

function shade_gaps(gap_blocks_t)
    for i = 1:size(gap_blocks_t, 1)
        xregion(gap_blocks_t(i, 1), gap_blocks_t(i, 2), 'FaceColor', [0.85 0.85 0.85], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    end
end

function shade_seizures_robust(seizures)
    for k = 1:height(seizures)
        if seizures.over_max_duration(k)
            edge_color = [0.6 0 0.6]; % over-max events get a distinct (purple) outline, still shown, never hidden
        else
            edge_color = 'r';
        end
        xregion(seizures.start_s(k), seizures.end_s(k), 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
            'EdgeColor', edge_color, 'LineWidth', 2);
    end
end

function lbl = region_label(meta, fallback)
    if isfield(meta, 'region') && ~isempty(meta.region)
        lbl = meta.region;
    else
        lbl = fallback;
    end
end

function [fig_file, png_file] = save_panorama_figure_robust(out_dir, base, cfg, meta, ...
        t_rel, signal, bp_full, energy_full, threshold, ll_full, ll_median_global, seizures, gap_blocks_t, n_seizures)
    target_points = 20000;

    fig = figure('Position', [50, 50, 1400, 1100], 'Visible', 'off');

    [t1, y1] = decimate_minmax(t_rel, signal, target_points);
    subplot(4, 1, 1);
    shade_gaps(gap_blocks_t); hold on;
    plot(t1, y1, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
    shade_seizures_robust(seizures);
    title(sprintf('%s - Cleaned LFP (%d seizures, robust branch)', region_label(meta, base), n_seizures), 'Interpreter', 'none');
    xlabel('Time (s)'); ylabel('uV'); grid on; hold off;

    [t2, y2] = decimate_minmax(t_rel, bp_full, target_points);
    subplot(4, 1, 2);
    shade_gaps(gap_blocks_t); hold on;
    plot(t2, y2, 'Color', [0.2 0.4 0.7], 'LineWidth', 0.5);
    shade_seizures_robust(seizures);
    title(sprintf('Band-passed [%g-%g Hz]', cfg.seizure.bandpass_band(1), cfg.seizure.bandpass_band(2)));
    xlabel('Time (s)'); ylabel('Normalized'); grid on; hold off;

    [t3, y3] = decimate_minmax(t_rel, energy_full, target_points);
    subplot(4, 1, 3);
    shade_gaps(gap_blocks_t); hold on;
    plot(t3, y3, 'Color', [0.2 0.6 0.3], 'LineWidth', 0.5);
    yline(threshold, 'r--', sprintf('Threshold (median x %.1f)', cfg.seizure.median_factor), 'LineWidth', 2);
    shade_seizures_robust(seizures);
    title('Energy Metric (Stage 1, identical to legacy)'); xlabel('Time (s)'); ylabel('Energy (a.u.)'); grid on; hold off;

    ll_norm = ll_full / ll_median_global;
    [t4, y4] = decimate_minmax(t_rel, ll_norm, target_points);
    subplot(4, 1, 4);
    shade_gaps(gap_blocks_t); hold on;
    plot(t4, y4, 'Color', [0.6 0.3 0.6], 'LineWidth', 0.5);
    yline(cfg.seizure_robust.ll_threshold, 'r--', sprintf('ll\\_threshold (%.2f)', cfg.seizure_robust.ll_threshold), 'LineWidth', 2);
    shade_seizures_robust(seizures);
    title('Line-Length Ratio (Stage 4 filter, median-based)'); xlabel('Time (s)'); ylabel('ll / ll\_median\_global'); grid on; hold off;

    fig_file = fullfile(out_dir, [base '_seizures.fig']);
    png_file = fullfile(out_dir, [base '_seizures.png']);
    savefig(fig, fig_file);
    saveas(fig, png_file, 'png');
    close(fig);
end

function save_zoom_figures_robust(out_dir, base, cfg, t_rel, signal, bp_full, energy_full, ...
        threshold, ll_full, ll_median_global, seizures, gap_blocks_t)

    ll_norm = ll_full / ll_median_global;

    for k = 1:height(seizures)
        zoom_start = max(seizures.start_s(k) - cfg.seizure.zoom_margin_s, t_rel(1));
        zoom_end = min(seizures.end_s(k) + cfg.seizure.zoom_margin_s, t_rel(end));
        idx_zoom = (t_rel >= zoom_start) & (t_rel <= zoom_end);

        fig = figure('Position', [50, 50, 1400, 1100], 'Visible', 'off');

        subplot(4, 1, 1);
        shade_gaps(gap_blocks_t); hold on;
        plot(t_rel(idx_zoom), signal(idx_zoom), 'k', 'LineWidth', 0.8);
        xregion(seizures.start_s(k), seizures.end_s(k), 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 2);
        xline(seizures.start_s(k), 'r--', 'Start', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
        xline(seizures.end_s(k), 'r--', 'End', 'LineWidth', 2, 'LabelHorizontalAlignment', 'right');
        over_max_note = '';
        if seizures.over_max_duration(k)
            over_max_note = ' [OVER MAX DURATION]';
        end
        title(sprintf('Seizure #%d (robust branch)%s', k, over_max_note), 'Interpreter', 'none');
        xlabel('Time (s)'); ylabel('Voltage (uV)'); grid on; xlim([zoom_start, zoom_end]); hold off;

        subplot(4, 1, 2);
        shade_gaps(gap_blocks_t); hold on;
        plot(t_rel(idx_zoom), bp_full(idx_zoom), 'k', 'LineWidth', 0.8);
        xregion(seizures.start_s(k), seizures.end_s(k), 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 2);
        xline(seizures.start_s(k), 'r--', 'LineWidth', 2);
        xline(seizures.end_s(k), 'r--', 'LineWidth', 2);
        title(sprintf('Band-passed [%g-%g Hz]', cfg.seizure.bandpass_band(1), cfg.seizure.bandpass_band(2)), 'Interpreter', 'none');
        xlabel('Time (s)'); ylabel('Normalized Amplitude'); grid on; xlim([zoom_start, zoom_end]); hold off;

        subplot(4, 1, 3);
        shade_gaps(gap_blocks_t); hold on;
        plot(t_rel(idx_zoom), energy_full(idx_zoom), 'Color', [0.2 0.6 0.3], 'LineWidth', 1.0);
        yline(threshold, 'r--', sprintf('Threshold (%.2e)', threshold), 'LineWidth', 2);
        xregion(seizures.start_s(k), seizures.end_s(k), 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 2);
        xline(seizures.start_s(k), 'r--', 'LineWidth', 2);
        xline(seizures.end_s(k), 'r--', 'LineWidth', 2);
        title('Energy Metric (Stage 1)');
        xlabel('Time (s)'); ylabel('Energy (a.u.)'); grid on; xlim([zoom_start, zoom_end]); hold off;

        subplot(4, 1, 4);
        shade_gaps(gap_blocks_t); hold on;
        plot(t_rel(idx_zoom), ll_norm(idx_zoom), 'Color', [0.6 0.3 0.6], 'LineWidth', 1.0);
        yline(cfg.seizure_robust.ll_threshold, 'r--', sprintf('ll\\_threshold (%.2f)', cfg.seizure_robust.ll_threshold), 'LineWidth', 2);
        xregion(seizures.start_s(k), seizures.end_s(k), 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 2);
        xline(seizures.start_s(k), 'r--', 'LineWidth', 2);
        xline(seizures.end_s(k), 'r--', 'LineWidth', 2);
        title(sprintf('Line-Length Ratio: ll=%.2f | peak\\_energy=%.1f | hf\\_db=%.1f | env\\_cv=%.2f', ...
            seizures.ll_ratio(k), seizures.peak_energy_ratio(k), seizures.hf_ratio_db(k), seizures.envelope_cv(k)), ...
            'Interpreter', 'tex');
        xlabel('Time (s)'); ylabel('ll / ll\_median\_global'); grid on; xlim([zoom_start, zoom_end]); hold off;

        zoom_fig_file = fullfile(out_dir, sprintf('%s_seizure%d.fig', base, k));
        zoom_png_file = fullfile(out_dir, sprintf('%s_seizure%d.png', base, k));
        savefig(fig, zoom_fig_file);
        saveas(fig, zoom_png_file, 'png');
        close(fig);
    end
end
