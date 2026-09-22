function iid_results = detect_iid(data, seizures, cfg)
% DETECT_IID  Ports the interictal-spike/polyspike/burst detector from
% IID_detection_FINAL.m verbatim: quantized baseline search, lower/upper
% amplitude thresholds, bandpass [15 70] Hz, findpeaks (distance/width/
% prominence), <=150ms ISI grouping into complexes, polyspike if
% N_spikes >= 2, burst grouping on <=5s inter-complex gaps with
% N>3 & 4<duration<40s kept.
%
% Differences from the original script (both spec-mandated):
%   - Exclusion zones are no longer a hand-typed list: they are built from
%     the detected seizures (+ exclusion_buffer_s) and the recording's own
%     gaps, unioned with any extra manual zones in cfg.iid.exclusion_zones_manual.
%   - Baseline, bandpass/findpeaks, and every rate are computed over
%     INCLUDED samples only (inclusion_mask = ~excluded & ~gap), and rates
%     always use analyzed_duration_min = sum(inclusion_mask)/fs/60.
%
% Gap adaptation: bandpass and findpeaks run per contiguous valid (gap-free)
% block, never across a NaN. Spike-to-complex and complex-to-burst grouping
% additionally break whenever consecutive spikes/complexes fall in
% different blocks, so no complex or burst can span a real gap even if its
% duration would otherwise pass the ISI/gap thresholds.
%
%   iid_results = detect_iid(data, seizures, cfg)
%
%   data     : struct from load_lfp_txt.m, loaded from a *_clean.txt file
%   seizures : seizures table from detect_seizures.m's output (or [] / omitted
%              to run with only gap+manual exclusion zones)
%   cfg      : struct from pipeline_config.m
%
% OUTPUT (struct iid_results): spike_complex_table, single_spike_table,
% polyspike_table, burst_table, individual_peaks_table, exclusion_table,
% summary (struct), mat_file/fig_file/png_file ('' if nothing saved).

    if nargin < 2 || isempty(seizures)
        seizures = empty_seizures_table();
    end
    ensure_findpeaks_signal_toolbox();

    signal = data.signal(:);
    fs = data.fs;
    t_rel = data.t_rel(:);
    valid_mask = data.valid_mask(:);
    session_start = data.session_start;

    %% ---- exclusion zones (seizures+buffer, gaps, manual) -----------------
    [exclusion_table, inclusion_mask] = build_exclusion_zones(seizures, valid_mask, t_rel, cfg);
    analyzed_duration_min = sum(inclusion_mask) / fs / 60;
    total_duration_min = numel(t_rel) / fs / 60;

    %% ---- baseline: quantized search over INCLUDED samples only -----------
    Voltage_mV = signal / 1000;
    included_abs_mV = abs(Voltage_mV(inclusion_mask));
    [baseline_uV, lower_uV, upper_uV] = compute_baseline(included_abs_mV, cfg);
    lower_mV = lower_uV / 1000;
    upper_mV = upper_uV / 1000;

    %% ---- bandpass [15 70] Hz per gap-free block ---------------------------
    raw_blocks = mask_to_segments(valid_mask);
    x_mV = (signal - median(signal(valid_mask))) / 1000;
    x_bp_full = filter_by_blocks(x_mV, raw_blocks, ...
        @(v) bandpass(v, cfg.iid.bandpass_band, fs, 'Steepness', cfg.iid.bandpass_steepness));
    abs_mV_full = abs(x_bp_full);

    %% ---- peak detection per block ----------------------------------------
    [spikes, pos, width, prom, spike_block_id] = detect_peaks_per_block(abs_mV_full, raw_blocks, fs, lower_mV, cfg);

    keep = spikes <= upper_mV;
    spikes = spikes(keep); pos = pos(keep); width = width(keep); prom = prom(keep); spike_block_id = spike_block_id(keep);

    keep2 = inclusion_mask(pos);
    spikes = spikes(keep2); pos = pos(keep2); width = width(keep2); prom = prom(keep2); spike_block_id = spike_block_id(keep2);

    pos_s = (pos - 1) / fs;

    %% ---- group spikes into complexes (ISI <= 150ms, never across a block)
    [complex_mat, complex_block_id, spike_complex_ids] = group_into_complexes(pos, pos_s, spikes, spike_block_id, fs, cfg);

    %% ---- group complexes into bursts (<=5s gap, never across a block) ----
    burst_mat = group_into_bursts(complex_mat, complex_block_id, cfg);

    %% ---- assemble tables ---------------------------------------------
    Individual_Peaks_Table = table(pos_s(:), spikes(:), width(:), prom(:), spike_block_id(:), ...
        'VariableNames', {'Time_s', 'Amplitude_mV', 'Width_samples', 'Prominence_mV', 'block_id'});

    if ~isempty(complex_mat)
        Spike_Complex_Table = array2table(complex_mat, 'VariableNames', ...
            {'Complex_ID', 'Start_s', 'End_s', 'Duration_ms', 'N_spikes', 'Max_amplitude_mV', 'Mean_amplitude_mV', 'Is_polyspike'});
        Spike_Complex_Table.Is_polyspike = logical(Spike_Complex_Table.Is_polyspike);
        Spike_Complex_Table.block_id = complex_block_id(:);
    else
        Spike_Complex_Table = empty_complex_table();
    end
    Spike_Complex_Table = add_abs_columns(Spike_Complex_Table, 'Start_s', 'End_s', session_start);
    Single_Spike_Table = Spike_Complex_Table(~Spike_Complex_Table.Is_polyspike, :);
    Polyspike_Table = Spike_Complex_Table(Spike_Complex_Table.Is_polyspike, :);

    if ~isempty(burst_mat)
        Burst_Table = array2table(burst_mat, 'VariableNames', {'N_complexes', 'Duration_s', 'Start_s', 'End_s', 'block_id'});
        Burst_Table = Burst_Table(:, {'Start_s', 'End_s', 'N_complexes', 'Duration_s', 'block_id'});
    else
        Burst_Table = empty_burst_table();
    end
    Burst_Table = add_abs_columns(Burst_Table, 'Start_s', 'End_s', session_start);

    %% ---- rates (always over analyzed_duration_min) ------------------------
    num_spike_complexes = height(Spike_Complex_Table);
    num_polyspikes = height(Polyspike_Table);
    num_single = num_spike_complexes - num_polyspikes;
    burst_count = height(Burst_Table);

    summary = struct( ...
        'file', data.file, ...
        'baseline_uV', baseline_uV, 'lower_threshold_uV', lower_uV, 'upper_threshold_uV', upper_uV, ...
        'total_duration_min', total_duration_min, 'analyzed_duration_min', analyzed_duration_min, ...
        'excluded_duration_min', total_duration_min - analyzed_duration_min, ...
        'n_exclusion_zones', height(exclusion_table), ...
        'total_peaks', numel(pos), 'total_complexes', num_spike_complexes, ...
        'n_single', num_single, 'n_polyspike', num_polyspikes, ...
        'pct_polyspike', 100 * num_polyspikes / max(1, num_spike_complexes), ...
        'complexes_per_min', num_spike_complexes / analyzed_duration_min, ...
        'single_per_min', num_single / analyzed_duration_min, ...
        'polyspikes_per_min', num_polyspikes / analyzed_duration_min, ...
        'n_bursts', burst_count, 'bursts_per_hour', burst_count / (analyzed_duration_min / 60), ...
        'mean_spikes_per_polyspike', mean_or_nan(Polyspike_Table.N_spikes));

    fprintf('detect_iid: %s -> %d complexes (%d single, %d polyspike), %d bursts (analyzed=%.1f min)\n', ...
        data.file, num_spike_complexes, num_single, num_polyspikes, burst_count, analyzed_duration_min);

    %% ---- output paths + figures -------------------------------------------
    out_dir = cfg.iid.output_dir;
    if isempty(out_dir)
        out_dir = pwd;
    end
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    [~, fname, ~] = fileparts(data.file);
    base = strrep(fname, '_clean', '');

    mat_file = ''; fig_file = ''; png_file = '';

    iid_results = struct();
    iid_results.spike_complex_table = Spike_Complex_Table;
    iid_results.single_spike_table = Single_Spike_Table;
    iid_results.polyspike_table = Polyspike_Table;
    iid_results.burst_table = Burst_Table;
    iid_results.individual_peaks_table = Individual_Peaks_Table;
    iid_results.exclusion_table = exclusion_table;
    iid_results.summary = summary;
    iid_results.fs = fs;
    iid_results.file = data.file;

    if num_spike_complexes > 0
        if cfg.general.overwrite || exist(fullfile(out_dir, [base '_IID_results.mat']), 'file') ~= 2
            mat_file = fullfile(out_dir, [base '_IID_results.mat']);
            save(mat_file, 'iid_results');
        end

        try
            [fig_file, png_file] = save_panorama_figure(out_dir, base, data.meta, t_rel, x_bp_full, ...
                lower_mV, exclusion_table, Single_Spike_Table, Polyspike_Table, Burst_Table, pos_s, spike_complex_ids, Spike_Complex_Table);

            save_polyspike_examples(out_dir, base, t_rel, x_bp_full, lower_mV, Polyspike_Table, pos_s, spike_complex_ids);
        catch ME
            warning('detect_iid:FigureSaveFailed', ...
                'Could not generate/save IID figures for %s: %s. IID detection results are unaffected.', ...
                data.file, ME.message);
        end
    end

    iid_results.mat_file = mat_file;
    iid_results.fig_file = fig_file;
    iid_results.png_file = png_file;
end

%% ======================================================================
function [baseline_uV, lower_uV, upper_uV] = compute_baseline(abs_mV_included, cfg)
    if isempty(abs_mV_included)
        error('detect_iid:NoIncludedSamples', 'No included (non-excluded, non-gap) samples remain for baseline calculation.');
    end

    sv = cfg.iid.baseline.start_uV;
    Percent = 0;
    while Percent < cfg.iid.baseline.target_fraction && sv <= cfg.iid.baseline.max_uV
        thr_mV = sv / 1000;
        Percent = mean(abs_mV_included < thr_mV);
        if Percent < cfg.iid.baseline.target_fraction
            sv = sv + cfg.iid.baseline.step_uV;
        end
    end

    if sv > cfg.iid.baseline.max_uV
        baseline_uV = cfg.iid.baseline.max_uV - cfg.iid.baseline.step_uV;
    else
        baseline_uV = sv - cfg.iid.baseline.step_uV;
    end

    lower_uV = cfg.iid.lower_threshold_factor * baseline_uV;
    upper_uV = cfg.iid.upper_threshold_uV;
end

function [exclusion_table, inclusion_mask] = build_exclusion_zones(seizures, valid_mask, t_rel, cfg)
    rows = cell(0, 4);  % {reason, source_id, start_s, end_s}

    for i = 1:height(seizures)
        rows(end+1, :) = {'seizure', seizures.id(i), ...
            seizures.start_s(i) - cfg.iid.exclusion_buffer_s, seizures.end_s(i) + cfg.iid.exclusion_buffer_s}; %#ok<AGROW>
    end

    gap_idx = mask_to_segments(~valid_mask);
    for i = 1:size(gap_idx, 1)
        rows(end+1, :) = {'gap', i, t_rel(gap_idx(i, 1)), t_rel(gap_idx(i, 2))}; %#ok<AGROW>
    end

    manual = cfg.iid.exclusion_zones_manual;
    for i = 1:size(manual, 1)
        rows(end+1, :) = {'manual', i, manual(i, 1), manual(i, 2)}; %#ok<AGROW>
    end

    inclusion_mask = valid_mask;

    if isempty(rows)
        exclusion_table = table('Size', [0 8], ...
            'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double', 'datetime', 'datetime'}, ...
            'VariableNames', {'zone_id', 'reason', 'source_id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs'});
        return;
    end

    t_min = t_rel(1);
    t_max = t_rel(end);
    reason = string(rows(:, 1));
    source_id = cell2mat(rows(:, 2));
    start_s = max(cell2mat(rows(:, 3)), t_min);
    end_s = min(cell2mat(rows(:, 4)), t_max);

    keep = end_s > start_s;
    reason = reason(keep); source_id = source_id(keep); start_s = start_s(keep); end_s = end_s(keep);

    for i = 1:numel(start_s)
        inclusion_mask(t_rel >= start_s(i) & t_rel <= end_s(i)) = false;
    end

    zone_id = (1:numel(start_s))';
    exclusion_table = table(zone_id, reason, source_id, start_s, end_s, end_s - start_s, ...
        'VariableNames', {'zone_id', 'reason', 'source_id', 'start_s', 'end_s', 'duration_s'});
end

function [spikes, pos, width, prom, block_id] = detect_peaks_per_block(abs_mV_full, raw_blocks, fs, lower_mV, cfg)
    MinPeakDistance = max(1, round(cfg.iid.min_peak_distance_s * fs));
    MaxPeakWidth = max(1, round(cfg.iid.max_peak_width_s * fs));

    spikes = zeros(0, 1); pos = zeros(0, 1); width = zeros(0, 1); prom = zeros(0, 1); block_id = zeros(0, 1);
    for i = 1:size(raw_blocks, 1)
        s = raw_blocks(i, 1);
        e = raw_blocks(i, 2);
        blk = abs_mV_full(s:e);
        [sp, po, wi, pr] = findpeaks(blk, ...
            'MinPeakDistance', MinPeakDistance, ...
            'MinPeakProminence', cfg.iid.min_peak_prominence_mV, ...
            'MinPeakHeight', lower_mV, ...
            'MaxPeakWidth', MaxPeakWidth);
        if isempty(po)
            continue;
        end
        spikes = [spikes; sp(:)]; pos = [pos; s - 1 + po(:)]; width = [width; wi(:)]; prom = [prom; pr(:)]; %#ok<AGROW>
        block_id = [block_id; repmat(i, numel(po), 1)]; %#ok<AGROW>
    end
end

function [complex_mat, complex_block_id, spike_complex_ids] = group_into_complexes(pos, pos_s, spikes, spike_block_id, fs, cfg)
    if isempty(pos)
        complex_mat = zeros(0, 8);
        complex_block_id = zeros(0, 1);
        spike_complex_ids = zeros(0, 1);
        return;
    end

    max_interspike_samples = round((cfg.iid.max_interspike_ms / 1000) * fs);
    ISI_samples = diff(pos);
    same_block = diff(spike_block_id) == 0;

    complex_id = 1;
    spike_complex_ids = ones(size(pos));
    for i = 1:numel(ISI_samples)
        if ISI_samples(i) > max_interspike_samples || ~same_block(i)
            complex_id = complex_id + 1;
        end
        spike_complex_ids(i+1) = complex_id;
    end

    num_complexes = max(spike_complex_ids);
    complex_mat = zeros(num_complexes, 8);
    complex_block_id = zeros(num_complexes, 1);

    for c = 1:num_complexes
        idx = find(spike_complex_ids == c);
        n_spikes = numel(idx);
        first_i = idx(1);
        last_i = idx(end);
        start_s = pos_s(first_i);
        end_s = pos_s(last_i);
        duration_ms = (end_s - start_s) * 1000;
        max_amp = max(spikes(idx));
        mean_amp = mean(spikes(idx));
        is_poly = n_spikes >= cfg.iid.polyspike_min_n;

        complex_mat(c, :) = [c, start_s, end_s, duration_ms, n_spikes, max_amp, mean_amp, is_poly];
        complex_block_id(c) = spike_block_id(first_i);
    end
end

function burst_mat = group_into_bursts(complex_mat, complex_block_id, cfg)
    if isempty(complex_mat)
        burst_mat = zeros(0, 5);
        return;
    end

    complex_times_s = complex_mat(:, 2);  % Start_s
    Interval_dt = diff(complex_times_s);
    same_block = diff(complex_block_id) == 0;
    Status = (Interval_dt <= cfg.iid.burst_max_gap_s) & same_block;

    D = zeros(0, 5);
    complex_tally = 1;
    dur = 0;
    t_ini = complex_times_s(1);
    blk = complex_block_id(1);

    for k = 1:numel(Interval_dt)
        if Status(k)
            complex_tally = complex_tally + 1;
            dur = dur + Interval_dt(k);
        else
            D = [D; complex_tally, dur, t_ini, t_ini + dur, blk]; %#ok<AGROW>
            complex_tally = 1;
            dur = 0;
            t_ini = complex_times_s(k+1);
            blk = complex_block_id(k+1);
        end
    end
    D = [D; complex_tally, dur, t_ini, t_ini + dur, blk];
    burst_mat = D;

    if ~isempty(burst_mat)
        mask = burst_mat(:, 1) > cfg.iid.burst_min_complexes & ...
               burst_mat(:, 2) > cfg.iid.burst_min_duration_s & ...
               burst_mat(:, 2) < cfg.iid.burst_max_duration_s;
        burst_mat = burst_mat(mask, :);
    end
end

function T = add_abs_columns(T, start_col, end_col, session_start)
    n = height(T);
    if n == 0 || isnat(session_start)
        T.start_abs = NaT(n, 1);
        T.end_abs = NaT(n, 1);
        T.start_abs.TimeZone = session_start.TimeZone;  % keep zoned/unzoned consistent with session_start for vertcat downstream
        T.end_abs.TimeZone = session_start.TimeZone;
    else
        T.start_abs = rel_to_abs_time(T.(start_col), session_start);
        T.end_abs = rel_to_abs_time(T.(end_col), session_start);
    end
end

function T = empty_complex_table()
    T = table('Size', [0 9], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'double', 'double', 'logical', 'double'}, ...
        'VariableNames', {'Complex_ID', 'Start_s', 'End_s', 'Duration_ms', 'N_spikes', 'Max_amplitude_mV', 'Mean_amplitude_mV', 'Is_polyspike', 'block_id'});
end

function T = empty_burst_table()
    T = table('Size', [0 5], 'VariableTypes', {'double', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'Start_s', 'End_s', 'N_complexes', 'Duration_s', 'block_id'});
end

function T = empty_seizures_table()
    T = table('Size', [0 7], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'datetime', 'datetime', 'double'}, ...
        'VariableNames', {'id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', 'block_id'});
end

function m = mean_or_nan(x)
    if isempty(x)
        m = NaN;
    else
        m = mean(x);
    end
end

%% ======================================================================
function shade_zones(exclusion_table)
    colors = struct('seizure', [1 0.85 0.6], 'gap', [0.85 0.85 0.85], 'manual', [0.75 0.7 0.95]);
    labels = struct('seizure', 'Excluded (seizure)', 'gap', 'Excluded (gap)', 'manual', 'Excluded (manual)');
    reasons = {'gap', 'seizure', 'manual'};  % draw gaps first (background), then events
    for r = 1:numel(reasons)
        reason = reasons{r};
        rows = find(exclusion_table.reason == reason);
        for i = rows(:)'
            xregion(exclusion_table.start_s(i), exclusion_table.end_s(i), ...
                'FaceColor', colors.(reason), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        if ~isempty(rows)
            % single legend proxy per reason actually present, instead of one per xregion call
            plot(NaN, NaN, 's', 'Color', colors.(reason), 'MarkerFaceColor', colors.(reason), ...
                'MarkerSize', 10, 'DisplayName', labels.(reason));
        end
    end
end

function [fig_file, png_file] = save_panorama_figure(out_dir, base, meta, t_rel, x_bp_full, lower_mV, ...
        exclusion_table, Single_Spike_Table, Polyspike_Table, Burst_Table, pos_s, spike_complex_ids, Spike_Complex_Table)

    % Min-max decimated for display only (see decimate_minmax.m) -- a long
    % recording's full-resolution line is what previously made this figure
    % too large to save; spike markers below are plotted individually and
    % unaffected, since there are orders of magnitude fewer of those than
    % raw samples.
    [t_plot, y_plot] = decimate_minmax(t_rel, x_bp_full, 20000);

    fig = figure('Name', 'IID Detection with Exclusion Zones', 'Position', [100 100 1400 600], 'Visible', 'off');
    shade_zones(exclusion_table); hold on;
    plot(t_plot, y_plot, 'k', 'DisplayName', 'Signal');

    single_ids = Spike_Complex_Table.Complex_ID(~Spike_Complex_Table.Is_polyspike);
    poly_ids = Spike_Complex_Table.Complex_ID(Spike_Complex_Table.Is_polyspike);
    single_mask = ismember(spike_complex_ids, single_ids);
    poly_mask = ismember(spike_complex_ids, poly_ids);

    if any(single_mask)
        plot(pos_s(single_mask), x_bp_full_at(x_bp_full, pos_s(single_mask), t_rel), 'bo', ...
            'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'Single spikes');
    end
    if any(poly_mask)
        plot(pos_s(poly_mask), x_bp_full_at(x_bp_full, pos_s(poly_mask), t_rel), 'ro', ...
            'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'Polyspikes');
    end

    first_box = true;
    for i = 1:height(Polyspike_Table)
        t_start = Polyspike_Table.Start_s(i);
        t_end = Polyspike_Table.End_s(i);
        idx = (t_rel >= t_start) & (t_rel <= t_end);
        y_range = [min(x_bp_full(idx), [], 'omitnan'), max(x_bp_full(idx), [], 'omitnan')];
        y_margin = 0.1 * diff(y_range);
        rectangle('Position', [t_start, y_range(1) - y_margin, t_end - t_start, diff(y_range) + 2 * y_margin], ...
            'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
        if first_box
            plot(NaN, NaN, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Polyspike boundaries');
            first_box = false;
        end
    end

    if height(Burst_Table) > 0
        plot([Burst_Table.Start_s(1) Burst_Table.End_s(1)], [1 1] * 1.2, 'm', 'LineWidth', 3, 'DisplayName', 'Bursts');
        for i = 2:height(Burst_Table)
            plot([Burst_Table.Start_s(i) Burst_Table.End_s(i)], [1 1] * 1.2, 'm', 'LineWidth', 3, 'HandleVisibility', 'off');
        end
    end

    yline(lower_mV, 'r-', 'LineWidth', 1, 'HandleVisibility', 'off');
    yline(-lower_mV, 'r-', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlabel('Time (s)'); ylabel('Voltage (mV)');
    title(sprintf('%s | Blue = single, Red = polyspikes, shaded = excluded (gap/seizure/manual)', region_label(meta, base)), 'Interpreter', 'none');
    legend('Location', 'best'); grid on; hold off;

    fig_file = fullfile(out_dir, [base '_IID.fig']);
    png_file = fullfile(out_dir, [base '_IID.png']);
    savefig(fig, fig_file);
    saveas(fig, png_file, 'png');
    close(fig);
end

function y = x_bp_full_at(x_bp_full, t_query_s, t_rel)
    fs_local = 1 / (t_rel(2) - t_rel(1));
    idx = round(t_query_s * fs_local) + 1;
    y = x_bp_full(idx);
end

function save_polyspike_examples(out_dir, base, t_rel, x_bp_full, lower_mV, Polyspike_Table, pos_s, spike_complex_ids)
    if isempty(Polyspike_Table) || height(Polyspike_Table) == 0
        return;
    end
    n_examples = min(4, height(Polyspike_Table));
    fig = figure('Name', 'Polyspike Complex Examples', 'Position', [200 200 1400 800], 'Visible', 'off');

    for i = 1:n_examples
        subplot(2, 2, i);
        t_center = (Polyspike_Table.Start_s(i) + Polyspike_Table.End_s(i)) / 2;
        window_s = 0.5;
        t_start_plot = max(0, t_center - window_s);
        t_end_plot = min(t_rel(end), t_center + window_s);
        idx_plot = (t_rel >= t_start_plot) & (t_rel <= t_end_plot);

        plot(t_rel(idx_plot), x_bp_full(idx_plot), 'k'); hold on;

        this_complex = Polyspike_Table.Complex_ID(i);
        spike_idx = find(spike_complex_ids == this_complex);
        for j = spike_idx(:)'
            plot(pos_s(j), x_bp_full_at(x_bp_full, pos_s(j), t_rel), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        end

        xline(Polyspike_Table.Start_s(i), '--r', 'LineWidth', 1.5);
        xline(Polyspike_Table.End_s(i), '--r', 'LineWidth', 1.5);
        yline(lower_mV, ':r');
        yline(-lower_mV, ':r');

        title(sprintf('Polyspike #%d: %d spikes, %.1f ms', i, Polyspike_Table.N_spikes(i), Polyspike_Table.Duration_ms(i)));
        xlabel('Time (s)'); ylabel('Voltage (mV)'); grid on; hold off;
    end

    fig_file = fullfile(out_dir, [base '_IID_polyspike_examples.fig']);
    png_file = fullfile(out_dir, [base '_IID_polyspike_examples.png']);
    savefig(fig, fig_file);
    saveas(fig, png_file, 'png');
    close(fig);
end

function lbl = region_label(meta, fallback)
    if isfield(meta, 'region') && ~isempty(meta.region)
        lbl = meta.region;
    else
        lbl = fallback;
    end
end
