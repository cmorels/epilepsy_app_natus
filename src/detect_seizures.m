function seizure_results = detect_seizures(data, cfg)
% DETECT_SEIZURES  Ports the seizure-detection method from
% complete_pipeline_seizures.m (STEP 3) verbatim: global min-max
% normalization, bandpass [5 75] Hz, 1 s edge trim, Hilbert envelope ^4,
% 2 s movmean, threshold = median(energy) x 10, minimum duration 15 s.
%
% Gap adaptation (see README for the full design rationale):
%   a) min-max normalization uses valid samples of the WHOLE recording.
%   b) bandpass runs per contiguous valid block (never across a NaN gap),
%      trimmed edge_trim_s at EACH block's own edges; blocks shorter than
%      2*edge_trim_s + min_seizure_duration are dropped and logged in QC.
%   c) Hilbert envelope / power / movmean computed per (trimmed) block.
%   d) threshold = median(energy) x median_factor over the CONCATENATION
%      of every surviving block's energy (one global statistic).
%   e) segmentation + duration filter per block; a seizure can never span
%      two blocks because it is found within one block's own energy trace.
%
%   seizure_results = detect_seizures(data, cfg)
%
%   data : struct from load_lfp_txt.m, loaded from a *_clean.txt file
%   cfg  : struct from pipeline_config.m
%
% OUTPUT (struct seizure_results):
%   .seizures : table (id, start_s, end_s, duration_s, start_abs, end_abs, block_id)
%   .metrics  : median_energy, threshold, pct_above, n_segments, n_rejected
%   .qc       : n_blocks_total, n_blocks_rejected_short, rejected_blocks (table)
%   .mat_file, .fig_file, .png_file : '' if nothing was saved (no seizures)

    signal = data.signal(:);
    fs = data.fs;
    t_rel = data.t_rel(:);
    valid_mask = data.valid_mask(:);
    session_start = data.session_start;

    valid_signal = signal(valid_mask);
    if numel(valid_signal) < 2
        error('detect_seizures:NotEnoughValidSamples', '%s has fewer than 2 valid (non-gap) samples.', data.file);
    end

    %% a) global min-max normalization
    min_val = min(valid_signal);
    max_val = max(valid_signal);
    lfp_norm = (signal - min_val) / (max_val - min_val);

    %% b) bandpass per valid block, then per-block edge trim
    raw_blocks = mask_to_segments(valid_mask);
    lfp_bp_full = filter_by_blocks(lfp_norm, raw_blocks, ...
        @(v) bandpass(v, cfg.seizure.bandpass_band, fs, 'Steepness', cfg.seizure.bandpass_steepness));

    trim_samples = round(cfg.seizure.edge_trim_s * fs);
    window_samples = round(cfg.seizure.window_sec * fs);
    min_len_samples = 2 * trim_samples + ceil(cfg.seizure.min_seizure_duration * fs);

    blocks = struct('block_id', {}, 'trimmed_start', {}, 'trimmed_end', {}, 'energy', {});
    rejected_rows = cell(0, 5);

    for i = 1:size(raw_blocks, 1)
        s = raw_blocks(i, 1);
        e = raw_blocks(i, 2);
        if (e - s + 1) < min_len_samples
            rejected_rows(end+1, :) = {i, t_rel(s), t_rel(e), (e - s + 1) / fs, 'too_short_after_trim'}; %#ok<AGROW>
            continue;
        end

        ts = s + trim_samples;
        te = e - trim_samples;

        %% c) Hilbert envelope, power, smoothing -- per block
        bp_trimmed = lfp_bp_full(ts:te);
        envelope = abs(hilbert(bp_trimmed));
        power = envelope .^ cfg.seizure.power_exponent;
        energy = movmean(power, window_samples);

        idx = numel(blocks) + 1;
        blocks(idx).block_id = i;
        blocks(idx).trimmed_start = ts;
        blocks(idx).trimmed_end = te;
        blocks(idx).energy = energy(:);
    end

    rejected_blocks = cell2table(rejected_rows, ...
        'VariableNames', {'block_id', 'start_s', 'end_s', 'duration_s', 'reason'});

    %% d) global threshold over the concatenation of surviving blocks
    if isempty(blocks)
        energy_all = zeros(0, 1);
    else
        energy_all = vertcat(blocks.energy);
    end

    if isempty(energy_all)
        median_energy = NaN;
        threshold = NaN;
        pct_above = 0;
    else
        median_energy = median(energy_all);
        threshold = median_energy * cfg.seizure.median_factor;
        pct_above = 100 * sum(energy_all > threshold) / numel(energy_all);
    end

    %% e) per-block segmentation + duration filter
    seizure_rows = cell(0, 7);
    n_segments = 0;
    seiz_id = 0;

    energy_full_plot = NaN(size(signal));
    bp_full_plot = NaN(size(signal));

    for b = 1:numel(blocks)
        blk = blocks(b);
        energy_full_plot(blk.trimmed_start:blk.trimmed_end) = blk.energy;
        bp_full_plot(blk.trimmed_start:blk.trimmed_end) = lfp_bp_full(blk.trimmed_start:blk.trimmed_end);

        above = blk.energy > threshold;
        local_segments = mask_to_segments(above);
        n_segments = n_segments + size(local_segments, 1);

        for k = 1:size(local_segments, 1)
            ls = local_segments(k, 1);
            le = local_segments(k, 2);
            duration_s = (le - ls + 1) / fs;
            if duration_s < cfg.seizure.min_seizure_duration
                continue;
            end
            gs = blk.trimmed_start + ls - 1;
            ge = blk.trimmed_start + le - 1;
            seiz_id = seiz_id + 1;
            start_s = t_rel(gs);
            end_s = t_rel(ge);
            seizure_rows(end+1, :) = {seiz_id, start_s, end_s, duration_s, ...
                abs_time_or_nat(start_s, session_start), abs_time_or_nat(end_s, session_start), blk.block_id}; %#ok<AGROW>
        end
    end

    n_rejected = n_segments - seiz_id;

    if isempty(seizure_rows)
        seizures = table('Size', [0 7], ...
            'VariableTypes', {'double', 'double', 'double', 'double', 'datetime', 'datetime', 'double'}, ...
            'VariableNames', {'id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', 'block_id'});
        seizures.start_abs.TimeZone = session_start.TimeZone;
        seizures.end_abs.TimeZone = session_start.TimeZone;
    else
        seizures = cell2table(seizure_rows, ...
            'VariableNames', {'id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', 'block_id'});
    end

    n_seizures = height(seizures);
    fprintf('detect_seizures: %s -> %d seizure(s) (median_energy=%.3e, threshold=%.3e, %d/%d blocks usable)\n', ...
        data.file, n_seizures, median_energy, threshold, numel(blocks), size(raw_blocks, 1));

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

    mat_file = '';
    fig_file = '';
    png_file = '';

    gap_blocks_t = gap_segments_seconds(valid_mask, t_rel);

    seizure_results = struct();
    seizure_results.seizures = seizures;
    seizure_results.metrics = struct( ...
        'median_energy', median_energy, 'threshold', threshold, 'pct_above', pct_above, ...
        'n_segments', n_segments, 'n_rejected', n_rejected);
    seizure_results.qc = struct( ...
        'n_blocks_total', size(raw_blocks, 1), 'n_blocks_rejected_short', height(rejected_blocks), ...
        'rejected_blocks', rejected_blocks);
    seizure_results.fs = fs;
    seizure_results.file = data.file;

    if n_seizures > 0
        if cfg.general.overwrite || exist(fullfile(out_dir, [base '_seizures.mat']), 'file') ~= 2
            mat_file = fullfile(out_dir, [base '_seizures.mat']);
            save(mat_file, 'seizure_results');
        end

        [fig_file, png_file] = save_panorama_figure(out_dir, base, cfg, data.meta, ...
            t_rel, signal, bp_full_plot, energy_full_plot, threshold, seizures, gap_blocks_t, n_seizures);

        save_zoom_figures(out_dir, base, cfg, t_rel, signal, bp_full_plot, energy_full_plot, ...
            threshold, seizures, gap_blocks_t);
    end

    seizure_results.mat_file = mat_file;
    seizure_results.fig_file = fig_file;
    seizure_results.png_file = png_file;
end

%% ======================================================================
function abs_t = abs_time_or_nat(t_rel_s, session_start)
    if isnat(session_start)
        abs_t = NaT;
        abs_t.TimeZone = session_start.TimeZone;  % keep zoned/unzoned consistent with session_start for vertcat downstream
    else
        abs_t = rel_to_abs_time(t_rel_s, session_start);
    end
end

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

function shade_seizures(seizures)
    for k = 1:height(seizures)
        xregion(seizures.start_s(k), seizures.end_s(k), 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
            'EdgeColor', 'r', 'LineWidth', 2);
    end
end

function [fig_file, png_file] = save_panorama_figure(out_dir, base, cfg, meta, ...
        t_rel, signal, bp_full_plot, energy_full_plot, threshold, seizures, gap_blocks_t, n_seizures)

    fig = figure('Position', [50, 50, 1400, 900], 'Visible', 'off');

    subplot(3, 1, 1);
    shade_gaps(gap_blocks_t); hold on;
    plot(t_rel, signal, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5);
    shade_seizures(seizures);
    title(sprintf('%s - Cleaned LFP (%d seizures)', region_label(meta, base), n_seizures), 'Interpreter', 'none');
    xlabel('Time (s)'); ylabel('uV'); grid on; hold off;

    subplot(3, 1, 2);
    shade_gaps(gap_blocks_t); hold on;
    plot(t_rel, bp_full_plot, 'Color', [0.2 0.4 0.7], 'LineWidth', 0.5);
    shade_seizures(seizures);
    title(sprintf('Band-passed [%g-%g Hz]', cfg.seizure.bandpass_band(1), cfg.seizure.bandpass_band(2)));
    xlabel('Time (s)'); ylabel('Normalized'); grid on; hold off;

    subplot(3, 1, 3);
    shade_gaps(gap_blocks_t); hold on;
    plot(t_rel, energy_full_plot, 'Color', [0.2 0.6 0.3], 'LineWidth', 0.5);
    yline(threshold, 'r--', sprintf('Threshold (median x %.1f)', cfg.seizure.median_factor), 'LineWidth', 2);
    above_mask = energy_full_plot > threshold;
    plot(t_rel(above_mask), energy_full_plot(above_mask), 'r.', 'MarkerSize', 3);
    shade_seizures(seizures);
    title('Energy Metric'); xlabel('Time (s)'); ylabel('Energy (a.u.)'); grid on; hold off;

    fig_file = fullfile(out_dir, [base '_seizures.fig']);
    png_file = fullfile(out_dir, [base '_seizures.png']);
    savefig(fig, fig_file);
    saveas(fig, png_file, 'png');
    close(fig);
end

function save_zoom_figures(out_dir, base, cfg, t_rel, signal, bp_full_plot, energy_full_plot, ...
        threshold, seizures, gap_blocks_t)

    for k = 1:height(seizures)
        zoom_start = max(seizures.start_s(k) - cfg.seizure.zoom_margin_s, t_rel(1));
        zoom_end = min(seizures.end_s(k) + cfg.seizure.zoom_margin_s, t_rel(end));
        idx_zoom = (t_rel >= zoom_start) & (t_rel <= zoom_end);

        fig = figure('Position', [50, 50, 1400, 900], 'Visible', 'off');

        subplot(3, 1, 1);
        shade_gaps(gap_blocks_t); hold on;
        plot(t_rel(idx_zoom), signal(idx_zoom), 'k', 'LineWidth', 0.8);
        xregion(seizures.start_s(k), seizures.end_s(k), 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 2);
        xline(seizures.start_s(k), 'r--', 'Start', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
        xline(seizures.end_s(k), 'r--', 'End', 'LineWidth', 2, 'LabelHorizontalAlignment', 'right');
        title(sprintf('Seizure #%d (full LFP)', k), 'Interpreter', 'none');
        xlabel('Time (s)'); ylabel('Voltage (uV)'); grid on; xlim([zoom_start, zoom_end]); hold off;

        subplot(3, 1, 2);
        shade_gaps(gap_blocks_t); hold on;
        plot(t_rel(idx_zoom), bp_full_plot(idx_zoom), 'k', 'LineWidth', 0.8);
        xregion(seizures.start_s(k), seizures.end_s(k), 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 2);
        xline(seizures.start_s(k), 'r--', 'LineWidth', 2);
        xline(seizures.end_s(k), 'r--', 'LineWidth', 2);
        title(sprintf('Seizure %d (band-passed [%g-%g Hz])', k, cfg.seizure.bandpass_band(1), cfg.seizure.bandpass_band(2)), 'Interpreter', 'none');
        xlabel('Time (s)'); ylabel('Normalized Amplitude'); grid on; xlim([zoom_start, zoom_end]); hold off;

        subplot(3, 1, 3);
        shade_gaps(gap_blocks_t); hold on;
        plot(t_rel(idx_zoom), energy_full_plot(idx_zoom), 'Color', [0.2 0.6 0.3], 'LineWidth', 1.0);
        yline(threshold, 'r--', sprintf('Threshold (%.2e)', threshold), 'LineWidth', 2);
        above_idx = idx_zoom & (energy_full_plot > threshold);
        if any(above_idx)
            plot(t_rel(above_idx), energy_full_plot(above_idx), 'r.', 'MarkerSize', 6);
        end
        xregion(seizures.start_s(k), seizures.end_s(k), 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 2);
        xline(seizures.start_s(k), 'r--', 'LineWidth', 2);
        xline(seizures.end_s(k), 'r--', 'LineWidth', 2);
        title(sprintf('Hilbert envelope(^%d) over a %ds smoothed window', cfg.seizure.power_exponent, cfg.seizure.window_sec));
        xlabel('Time (s)'); ylabel('Energy (a.u.)'); grid on; xlim([zoom_start, zoom_end]); hold off;

        zoom_fig_file = fullfile(out_dir, sprintf('%s_seizure%d.fig', base, k));
        zoom_png_file = fullfile(out_dir, sprintf('%s_seizure%d.png', base, k));
        savefig(fig, zoom_fig_file);
        saveas(fig, zoom_png_file, 'png');
        close(fig);
    end
end

function lbl = region_label(meta, fallback)
    if isfield(meta, 'region') && ~isempty(meta.region)
        lbl = meta.region;
    else
        lbl = fallback;
    end
end
