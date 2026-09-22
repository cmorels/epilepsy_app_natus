function result = run_pipeline_edf(input_folder, cfg)
% RUN_PIPELINE_EDF  Batch orchestrator: EDF -> txt -> clean -> seizures ->
% IID -> consolidated summaries, for every .edf in input_folder.
%
%   result = run_pipeline_edf(input_folder, cfg)
%
% Output layout under cfg.paths.output_root:
%   01_txt/        txt per channel + gaps CSV (edf_import.m)
%   02_clean/      *_clean.txt (clean_lfp.m)
%   03_seizures/   mat/fig/png per channel (detect_seizures.m)
%   04_iid/        mat/fig/png per channel (detect_iid.m)
%   05_summaries/  the 9 consolidated CSVs + pipeline_summary.xlsx
%   logs/          timestamped run log
%   config_used.mat / config_used.json  (the exact cfg this run used)
%
% Robustness: a failing file never aborts the batch. Every stage call is
% individually try/caught; failures are logged (file, region, stage,
% message, line) into qc_report.csv and the run log, and the loop moves on
% with whatever partial results are available for that channel.
%
% Resumability (scoped, see README.md): 02_clean is genuinely skipped (not
% just not-overwritten) when its output txt already exists and
% cfg.general.overwrite is false -- that's the expensive-enough stage
% where re-running is worth avoiding and its output filename is a
% trivial, low-drift-risk one-line rule. 01_txt/03_seizures/04_iid still
% run every time; they internally refuse to overwrite an existing file,
% but do not skip the underlying computation.
%
% To combine several separate runs (e.g. one per subject) into one set of
% summaries, see merge_pipeline_runs.m -- this function only consolidates
% files processed within THIS call; a second call with the same
% cfg.paths.output_root overwrites 05_summaries/, it does not append.

    if nargin < 2 || isempty(cfg)
        cfg = pipeline_config();
    end
    if ~isfolder(input_folder)
        error('run_pipeline_edf:BadInput', 'Not a folder: %s', input_folder);
    end

    dirs = make_output_dirs(cfg.paths.output_root);
    log_file = start_log(dirs.logs);
    save_config_used(cfg, cfg.paths.output_root);

    files = discover_edf_files(input_folder);
    log_line(log_file, sprintf('run_pipeline_edf: %d EDF file(s) found in %s', numel(files), input_folder));
    if isempty(files)
        error('run_pipeline_edf:NoFilesFound', 'No .edf files found in: %s', input_folder);
    end

    seizure_event_parts = {}; seizure_summary_parts = {};
    iid_event_parts = {}; iid_summary_parts = {}; iid_burst_parts = {};
    gap_parts = {}; qc_parts = {};
    file_subject_map = containers.Map('KeyType', 'char', 'ValueType', 'char');

    for f = 1:numel(files)
        file_path = fullfile(files(f).folder, files(f).name);
        log_line(log_file, sprintf('FILE %d/%d: %s', f, numel(files), files(f).name));

        file_cfg = cfg;
        file_cfg.edf.output_dir = dirs.txt;

        try
            [manifest, gaps_table, ~] = edf_import(file_path, file_cfg);
        catch ME
            log_error(log_file, files(f).name, '', 'edf_import', ME);
            qc_parts{end+1} = build_qc_row('', '', files(f).name, NaT, {}, {sprintf('edf_import: %s', ME.message)}, {}, NaN, NaN); %#ok<AGROW>
            continue;
        end
        gap_parts{end+1} = gaps_table; %#ok<AGROW>
        if height(manifest) > 0
            file_subject_map(files(f).name) = manifest.subject_id{1};
        end

        for c = 1:height(manifest)
            row = manifest(c, :);
            region = row.region{1};
            stages = {}; errors = {}; warnings_list = {};
            raw_data = []; clean_data = []; clean_stats = []; seizure_results = []; iid_results = [];

            try
                raw_data = load_lfp_txt(row.txt_file{1});
                stages{end+1} = 'import'; %#ok<AGROW>
            catch ME
                log_error(log_file, files(f).name, region, 'load_raw', ME);
                errors{end+1} = sprintf('load_raw: %s', ME.message); %#ok<AGROW>
            end

            if ~isempty(raw_data)
                try
                    [clean_data, clean_stats, was_skipped] = run_clean_stage(raw_data, file_cfg, dirs.clean);
                    stages{end+1} = 'clean'; %#ok<AGROW>
                    if was_skipped
                        log_line(log_file, sprintf('  [%s] clean_lfp SKIPPED (output exists)', region));
                    elseif clean_stats.auto_switch
                        warnings_list{end+1} = clean_stats.switched_reason; %#ok<AGROW>
                    end
                catch ME
                    log_error(log_file, files(f).name, region, 'clean_lfp', ME);
                    errors{end+1} = sprintf('clean_lfp: %s', ME.message); %#ok<AGROW>
                end
            end

            if ~isempty(clean_data)
                file_cfg.seizure.output_dir = dirs.seizures;
                try
                    seizure_results = detect_seizures(clean_data, file_cfg);
                    stages{end+1} = 'seizures'; %#ok<AGROW>
                    if seizure_results.qc.n_blocks_rejected_short > 0
                        warnings_list{end+1} = sprintf('%d block(s) rejected as too short for seizure detection', ...
                            seizure_results.qc.n_blocks_rejected_short); %#ok<AGROW>
                    end
                catch ME
                    log_error(log_file, files(f).name, region, 'detect_seizures', ME);
                    errors{end+1} = sprintf('detect_seizures: %s', ME.message); %#ok<AGROW>
                end

                file_cfg.iid.output_dir = dirs.iid;
                try
                    seizures_for_iid = table_or_empty(seizure_results, 'seizures');
                    iid_results = detect_iid(clean_data, seizures_for_iid, file_cfg);
                    stages{end+1} = 'iid'; %#ok<AGROW>
                catch ME
                    log_error(log_file, files(f).name, region, 'detect_iid', ME);
                    errors{end+1} = sprintf('detect_iid: %s', ME.message); %#ok<AGROW>
                end
            end

            session_start = pick_session_start(raw_data, clean_data, cfg.general.timezone);

            if ~isempty(seizure_results)
                seizure_event_parts{end+1} = build_seizure_event_rows( ...
                    seizure_results.seizures, region, row.subject_id{1}, session_start, row.source_file{1}, ...
                    clean_data.valid_mask, clean_data.t_rel, file_cfg.seizure.edge_trim_s); %#ok<AGROW>
                seizure_summary_parts{end+1} = build_seizure_summary_row(row, session_start, seizure_results, file_cfg); %#ok<AGROW>
            end
            if ~isempty(iid_results)
                iid_event_parts{end+1} = build_iid_event_rows( ...
                    iid_results.spike_complex_table, iid_results.burst_table, region, row.subject_id{1}, session_start, row.source_file{1}); %#ok<AGROW>
                iid_summary_parts{end+1} = build_iid_summary_row(row, session_start, iid_results); %#ok<AGROW>
                iid_burst_parts{end+1} = build_iid_burst_rows(iid_results.burst_table, region, row.subject_id{1}, row.source_file{1}, session_start.TimeZone); %#ok<AGROW>
            end

            outlier_pct = NaN;
            if ~isempty(clean_stats)
                outlier_pct = clean_stats.outlier_pct;
            end
            nan_pct = 100 * (1 - row.n_valid_samples / row.n_samples);
            n_blocks_rejected = NaN;
            if ~isempty(seizure_results)
                n_blocks_rejected = seizure_results.qc.n_blocks_rejected_short;
            end
            qc_parts{end+1} = build_qc_row(row.subject_id{1}, region, row.source_file{1}, session_start, ...
                stages, errors, warnings_list, outlier_pct, nan_pct, n_blocks_rejected); %#ok<AGROW>
        end
    end

    tz = cfg.general.timezone;
    seizures_events = vertcat_or_empty(seizure_event_parts, @() empty_seizure_events_table(tz));
    seizures_summary = vertcat_or_empty(seizure_summary_parts, @() empty_seizure_summary_table(tz));
    iid_events = vertcat_or_empty(iid_event_parts, @() empty_iid_events_table(tz));
    iid_summary = vertcat_or_empty(iid_summary_parts, @() empty_iid_summary_table(tz));
    iid_bursts = vertcat_or_empty(iid_burst_parts, @() empty_iid_bursts_table(tz));
    gaps_summary = vertcat_or_empty(gap_parts, @() empty_gaps_table(tz));
    qc_report = vertcat_or_empty(qc_parts, @() empty_qc_table(tz));

    natus_review_sheet = build_natus_review_sheet(seizures_events, iid_bursts, gaps_summary, file_subject_map, tz);

    paths = write_all_summaries(dirs.summaries, seizures_events, seizures_summary, ...
        iid_events, iid_summary, iid_bursts, gaps_summary, qc_report, natus_review_sheet);

    log_line(log_file, sprintf('DONE: %d seizure(s), %d IID complex(es), %d burst(s), %d error(s) across %d file(s)', ...
        height(seizures_events), height(iid_events), height(iid_bursts), sum(qc_report.n_errors), numel(files)));
    fclose(log_file);

    result = struct();
    result.dirs = dirs;
    result.seizures_events = seizures_events;
    result.seizures_summary = seizures_summary;
    result.iid_events = iid_events;
    result.iid_summary = iid_summary;
    result.iid_bursts = iid_bursts;
    result.gaps_summary = gaps_summary;
    result.qc_report = qc_report;
    result.natus_review_sheet = natus_review_sheet;
    result.paths = paths;
end

%% ======================================================================
function dirs = make_output_dirs(output_root)
    names = {'01_txt', '02_clean', '03_seizures', '04_iid', '05_summaries', 'logs'};
    for i = 1:numel(names)
        d = fullfile(output_root, names{i});
        if ~isfolder(d)
            mkdir(d);
        end
    end
    dirs = struct('root', output_root, 'txt', fullfile(output_root, '01_txt'), ...
        'clean', fullfile(output_root, '02_clean'), 'seizures', fullfile(output_root, '03_seizures'), ...
        'iid', fullfile(output_root, '04_iid'), 'summaries', fullfile(output_root, '05_summaries'), ...
        'logs', fullfile(output_root, 'logs'));
end

function fid = start_log(logs_dir)
    log_path = fullfile(logs_dir, sprintf('run_%s.log', datestr(now, 'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
    fid = fopen(log_path, 'w');
    log_line(fid, sprintf('run_pipeline_edf started %s', datestr(now))); %#ok<TNOW1,DATST>
end

function log_line(fid, msg)
    line = sprintf('[%s] %s', datestr(now, 'HH:MM:SS'), msg); %#ok<TNOW1,DATST>
    fprintf('%s\n', line);
    fprintf(fid, '%s\n', line);
end

function log_error(fid, source_file, region, stage, ME)
    if isempty(ME.stack)
        loc = 'n/a';
    else
        loc = sprintf('%s:%d', ME.stack(1).name, ME.stack(1).line);
    end
    log_line(fid, sprintf('  ERROR [%s | %s | %s] %s (at %s)', source_file, region, stage, ME.message, loc));
end

function save_config_used(cfg, output_root)
    save(fullfile(output_root, 'config_used.mat'), 'cfg');
    try
        json_txt = jsonencode(jsonify(cfg), 'PrettyPrint', true);
        fid = fopen(fullfile(output_root, 'config_used.json'), 'w');
        fprintf(fid, '%s', json_txt);
        fclose(fid);
    catch ME
        warning('run_pipeline_edf:ConfigJson', 'Could not write config_used.json: %s', ME.message);
    end
end

function v = jsonify(v)
% Recursively replace containers.Map (not supported by jsonencode) with a
% plain key/value struct array, so config_used.json can be written.
    if isa(v, 'containers.Map')
        k = keys(v); val = values(v);
        if isempty(k)
            v = struct('key', {}, 'value', {});
        else
            v = struct('key', k(:), 'value', val(:));
        end
    elseif isstruct(v)
        for fn = fieldnames(v)'
            for i = 1:numel(v)
                v(i).(fn{1}) = jsonify(v(i).(fn{1}));
            end
        end
    end
end

function files_out = discover_edf_files(folder)
    files = [dir(fullfile(folder, '*.edf')); dir(fullfile(folder, '*.EDF'))];
    if isempty(files)
        files_out = files;
        return;
    end
    full_paths = fullfile({files.folder}, {files.name});
    [~, ia] = unique(full_paths, 'stable');
    files_out = files(ia);
end

%% ======================================================================
function [clean_data, clean_stats, was_skipped] = run_clean_stage(raw_data, cfg, clean_dir)
    [~, name, ~] = fileparts(raw_data.file);
    expected_clean = fullfile(clean_dir, [name '_clean.txt']);

    if ~cfg.general.overwrite && exist(expected_clean, 'file') == 2
        clean_data = load_lfp_txt(expected_clean);
        clean_stats = struct('auto_switch', false, 'switched_reason', '', 'outlier_pct', NaN);
        if isfield(clean_data.meta, 'outlier_pct')
            clean_stats.outlier_pct = str2double(clean_data.meta.outlier_pct);
        end
        was_skipped = true;
        return;
    end

    cfg.clean.output_dir = clean_dir;
    clean_result = clean_lfp(raw_data, cfg);
    clean_data = load_lfp_txt(clean_result.txt_file);
    clean_stats = clean_result.stats;
    was_skipped = false;
end

function session_start = pick_session_start(raw_data, clean_data, tz)
% Always returns a TimeZone-bearing datetime, even when the value itself
% is NaT (unknown session start) -- otherwise vertcat-ing datetime columns
% across files where some have a real session_start and others don't
% fails ("cannot concatenate a zoned datetime with an unzoned one").
    if ~isempty(clean_data)
        session_start = clean_data.session_start;
    elseif ~isempty(raw_data)
        session_start = raw_data.session_start;
    else
        session_start = NaT;
    end
    if isempty(session_start.TimeZone)
        session_start.TimeZone = tz;
    end
end

function T = table_or_empty(seizure_results, field)
    if isempty(seizure_results)
        T = [];
    else
        T = seizure_results.(field);
    end
end

%% ======================================================================
function tf = compute_adjacent_to_gap(start_s, end_s, valid_mask, t_rel, tol_s)
    gap_idx = mask_to_segments(~valid_mask);
    n = numel(start_s);
    tf = false(n, 1);
    if isempty(gap_idx)
        return;
    end
    gap_starts_t = t_rel(gap_idx(:, 1));
    gap_ends_t = t_rel(gap_idx(:, 2));
    for i = 1:n
        near = any(abs(start_s(i) - gap_ends_t) <= tol_s) || any(abs(start_s(i) - gap_starts_t) <= tol_s) || ...
               any(abs(end_s(i) - gap_ends_t) <= tol_s) || any(abs(end_s(i) - gap_starts_t) <= tol_s);
        tf(i) = near;
    end
end

function T = build_seizure_event_rows(seizures, region, subject_id, session_start, source_file, valid_mask, t_rel, edge_trim_s)
    n = height(seizures);
    if n == 0
        T = empty_seizure_events_table(session_start.TimeZone);
        return;
    end
    T = seizures;
    T.Properties.VariableNames{'id'} = 'seizure_id';
    T.subject_id = repmat({subject_id}, n, 1);
    T.region = repmat({region}, n, 1);
    T.session_start = repmat(session_start, n, 1);
    T.source_file = repmat({source_file}, n, 1);
    T.adjacent_to_gap = compute_adjacent_to_gap(seizures.start_s, seizures.end_s, valid_mask, t_rel, edge_trim_s);
    T = T(:, {'subject_id', 'region', 'session_start', 'source_file', 'seizure_id', 'start_s', 'end_s', ...
        'duration_s', 'start_abs', 'end_abs', 'block_id', 'adjacent_to_gap'});
end

function T = build_seizure_summary_row(row, session_start, seizure_results, cfg)
    seizures = seizure_results.seizures;
    n_seizures = height(seizures);
    if n_seizures > 0
        total_seizure_time_s = sum(seizures.duration_s);
        mean_duration_s = mean(seizures.duration_s);
        min_duration_s = min(seizures.duration_s);
        max_duration_s = max(seizures.duration_s);
    else
        total_seizure_time_s = 0; mean_duration_s = NaN; min_duration_s = NaN; max_duration_s = NaN;
    end
    total_duration_min = row.total_duration_s / 60;
    valid_duration_min = row.valid_duration_s / 60;

    T = table(row.subject_id, row.region, session_start, row.source_file, ...
        total_duration_min, valid_duration_min, row.n_gaps, total_duration_min - valid_duration_min, ...
        n_seizures, total_seizure_time_s, 100 * total_seizure_time_s / (60 * total_duration_min), ...
        mean_duration_s, min_duration_s, max_duration_s, ...
        seizure_results.metrics.median_energy, seizure_results.metrics.threshold, seizure_results.metrics.pct_above, ...
        seizure_results.metrics.n_segments, seizure_results.metrics.n_rejected, ...
        cfg.seizure.bandpass_band(1), cfg.seizure.bandpass_band(2), cfg.seizure.power_exponent, ...
        cfg.seizure.window_sec, cfg.seizure.median_factor, cfg.seizure.min_seizure_duration, ...
        'VariableNames', {'subject_id', 'region', 'session_start', 'source_file', 'total_duration_min', 'valid_duration_min', ...
        'n_gaps', 'gap_duration_min', 'n_seizures', 'total_seizure_time_s', 'pct_time_in_seizure', ...
        'mean_duration_s', 'min_duration_s', 'max_duration_s', 'median_energy', 'threshold_value', 'pct_above_thr', ...
        'n_segments', 'n_rejected', 'bandpass_low', 'bandpass_high', 'power_exponent', 'window_s', 'median_factor', 'min_seizure_duration'});
end

function [in_burst, burst_id] = find_burst_membership(complex_start_s, burst_table)
    n = numel(complex_start_s);
    in_burst = false(n, 1);
    burst_id = nan(n, 1);
    for b = 1:height(burst_table)
        mask = complex_start_s >= burst_table.Start_s(b) & complex_start_s <= burst_table.End_s(b);
        in_burst(mask) = true;
        burst_id(mask) = b;
    end
end

function T = build_iid_event_rows(spike_complex_table, burst_table, region, subject_id, session_start, source_file)
    n = height(spike_complex_table);
    if n == 0
        T = empty_iid_events_table(session_start.TimeZone);
        return;
    end
    T = spike_complex_table;
    classification = repmat({'single'}, n, 1);
    classification(T.Is_polyspike) = {'polyspike'};
    T.classification = classification;

    [in_burst, burst_id] = find_burst_membership(T.Start_s, burst_table);
    T.in_burst = in_burst;
    T.burst_id = burst_id;

    T.subject_id = repmat({subject_id}, n, 1);
    T.region = repmat({region}, n, 1);
    T.session_start = repmat(session_start, n, 1);
    T.source_file = repmat({source_file}, n, 1);

    T.Properties.VariableNames{'Complex_ID'} = 'complex_id';
    T.Properties.VariableNames{'Start_s'} = 'start_s';
    T.Properties.VariableNames{'End_s'} = 'end_s';
    T.Properties.VariableNames{'Duration_ms'} = 'duration_ms';
    T.Properties.VariableNames{'N_spikes'} = 'n_spikes';
    T.Properties.VariableNames{'Max_amplitude_mV'} = 'max_amplitude_mV';
    T.Properties.VariableNames{'Mean_amplitude_mV'} = 'mean_amplitude_mV';

    T = T(:, {'subject_id', 'region', 'session_start', 'source_file', 'complex_id', 'start_s', 'end_s', ...
        'start_abs', 'end_abs', 'duration_ms', 'n_spikes', 'classification', 'max_amplitude_mV', 'mean_amplitude_mV', ...
        'in_burst', 'burst_id'});
end

function T = build_iid_summary_row(row, session_start, iid_results)
    s = iid_results.summary;
    T = table(row.subject_id, row.region, session_start, row.source_file, ...
        s.total_duration_min, s.analyzed_duration_min, s.excluded_duration_min, s.n_exclusion_zones, ...
        s.baseline_uV, s.lower_threshold_uV, s.upper_threshold_uV, ...
        s.total_peaks, s.total_complexes, s.n_single, s.n_polyspike, s.pct_polyspike, ...
        s.complexes_per_min, s.single_per_min, s.polyspikes_per_min, s.n_bursts, s.bursts_per_hour, s.mean_spikes_per_polyspike, ...
        'VariableNames', {'subject_id', 'region', 'session_start', 'source_file', 'total_duration_min', 'analyzed_duration_min', ...
        'excluded_duration_min', 'n_exclusion_zones', 'baseline_uV', 'lower_threshold_uV', 'upper_threshold_uV', ...
        'total_peaks', 'total_complexes', 'n_single', 'n_polyspike', 'pct_polyspike', 'complexes_per_min', 'single_per_min', ...
        'polyspikes_per_min', 'n_bursts', 'bursts_per_hour', 'mean_spikes_per_polyspike'});
end

function T = build_iid_burst_rows(burst_table, region, subject_id, source_file, tz)
    n = height(burst_table);
    if n == 0
        T = empty_iid_bursts_table(tz);
        return;
    end
    T = burst_table;
    T.subject_id = repmat({subject_id}, n, 1);
    T.region = repmat({region}, n, 1);
    T.Properties.VariableNames{'Start_s'} = 'start_s';
    T.Properties.VariableNames{'End_s'} = 'end_s';
    T.Properties.VariableNames{'N_complexes'} = 'n_complexes';
    T.Properties.VariableNames{'Duration_s'} = 'duration_s';
    T = T(:, {'subject_id', 'region', 'start_s', 'end_s', 'start_abs', 'end_abs', 'duration_s', 'n_complexes'});
end

function T = build_qc_row(subject_id, region, source_file, session_start, stages, errors, warnings_list, outlier_pct, nan_pct, n_blocks_rejected)
    if nargin < 10
        n_blocks_rejected = NaN;
    end
    T = table({subject_id}, {region}, {source_file}, session_start, {strjoin(stages, ',')}, ...
        numel(errors), {strjoin(errors, '; ')}, n_blocks_rejected, outlier_pct, nan_pct, ...
        numel(warnings_list), {strjoin(warnings_list, '; ')}, ...
        'VariableNames', {'subject_id', 'region', 'source_file', 'session_start', 'stages_completed', ...
        'n_errors', 'error_messages', 'n_blocks_rejected_short', 'outlier_pct', 'nan_pct', 'n_warnings', 'warning_messages'});
end
