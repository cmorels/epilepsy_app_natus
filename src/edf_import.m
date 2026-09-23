function [manifest, gaps_table, file_meta] = edf_import(edf_input, cfg)
% EDF_IMPORT  Convert a Natus EDF(+C/+D) recording into one continuous txt
% per selected channel, plus a gaps table, WITHOUT splitting the recording
% at discontinuities (see KNOWN_ISSUES.md / README.md for the rationale).
%
%   [manifest, gaps_table, file_meta] = edf_import(edf_path, cfg)
%   [manifest, gaps_table, file_meta] = edf_import(folder_path, cfg)
%
% A single sample index i (1-based) always corresponds to
% t_rel = (i-1)/fs seconds after the EDF's own start (session_start); real
% acquisition gaps are filled with NaN at that position instead of
% compressing the time axis. Samples are never lost or duplicated: every
% record's samples are copied exactly once (see place_samples below).
%
% INPUTS
%   edf_input : path to a single .edf file, or a folder containing .edf
%               files (processed one by one, outputs concatenated).
%   cfg       : struct from pipeline_config.m (default: pipeline_config()).
%
% OUTPUTS
%   manifest   : table, one row per exported channel (subject_id, region,
%                channel_label, fs, txt_file, n_samples, n_valid_samples,
%                n_gaps, total_duration_s, valid_duration_s, units,
%                source_file).
%   gaps_table : table, one row per detected gap across all processed
%                files (gap_id, start_s, end_s, duration_s, start_abs,
%                end_abs, prev_record_idx, next_record_idx, source_file).
%   file_meta  : struct array, one entry per processed EDF file (timing,
%                subject/session identity, QC counters).

    if nargin < 2 || isempty(cfg)
        cfg = pipeline_config();
    end
    validate_cfg(cfg);

    if ~(ischar(edf_input) || isstring(edf_input))
        error('edf_import:BadInput', 'edf_input must be a char/string path to a .edf file or a folder.');
    end
    edf_input = char(edf_input);

    if isfolder(edf_input)
        files = [dir(fullfile(edf_input, '*.edf')); dir(fullfile(edf_input, '*.EDF'))];
        files = unique_by_name(files);
        if isempty(files)
            error('edf_import:NoFilesFound', 'No .edf files found in: %s', edf_input);
        end

        manifest = table();
        gaps_table = table();
        file_meta = struct([]);
        for i = 1:numel(files)
            fpath = fullfile(files(i).folder, files(i).name);
            [m_i, g_i, meta_i] = edf_import_one(fpath, cfg);
            manifest = vertcat_tables(manifest, m_i);
            gaps_table = vertcat_tables(gaps_table, g_i);
            if isempty(file_meta)
                file_meta = meta_i;
            else
                file_meta(end+1) = meta_i; %#ok<AGROW>
            end
        end
        return;
    end

    if exist(edf_input, 'file') ~= 2
        error('edf_import:FileNotFound', 'File not found: %s', edf_input);
    end
    [manifest, gaps_table, file_meta] = edf_import_one(edf_input, cfg);
end

%% ======================================================================
function [manifest, gaps_table, file_meta] = edf_import_one(edf_path, cfg)

    [~, ~, ext] = fileparts(edf_path);
    if ~ismember(lower(ext), {'.edf'})
        error('edf_import:BadExtension', 'Expected a .edf file, got "%s" (%s).', ext, edf_path);
    end

    info = edfinfo(edf_path);
    [labels, regions] = resolve_channels(info, cfg);

    source_file_name = char(info.Filename);
    source_format = strtrim(char(info.Reserved));

    if isempty(cfg.edf.subject_id)
        subject_id = derive_subject_id(edf_path);
        warning('edf_import:SubjectIdDerived', ...
            'cfg.edf.subject_id is empty; derived subject_id "%s" from EDF filename. Set cfg.edf.subject_id explicitly for real runs.', ...
            subject_id);
    else
        subject_id = cfg.edf.subject_id;
    end

    base_datetime = datetime([char(info.StartDate) ' ' char(info.StartTime)], ...
        'InputFormat', cfg.general.session_start_input_format, ...
        'TimeZone', cfg.general.timezone);

    record_duration_s = seconds(info.DataRecordDuration(1));
    nrec = info.NumDataRecords;

    [hdr, ~] = edfread(edf_path, 'TimeOutputType', 'duration', 'SelectedSignals', labels);
    % edfread sanitizes signal labels into table variable names (e.g. "EEG A5C2"
    % -> "EEGA5C2", spaces stripped), so hdr.(label) with the ORIGINAL label can
    % fail even though that label matched info.SignalLabels correctly. SelectedSignals
    % preserves request order, so hdr_var_names{i} is always the right column for
    % labels{i} regardless of how edfread renamed it.
    hdr_var_names = hdr.Properties.VariableNames;
    record_times = hdr.("Record Time");
    rt0 = record_times(1);
    rt = seconds(record_times - rt0);
    session_start = base_datetime + rt0;

    dt = diff(rt);
    [gaps_table, is_real_gap, gap_duration_s, qc] = compute_gaps( ...
        rt, dt, record_duration_s, cfg.edf.gap_tol_s, cfg.edf.gap_min_s, session_start);
    gaps_table.source_file = repmat({source_file_name}, height(gaps_table), 1);

    out_dir = cfg.edf.output_dir;
    if isempty(out_dir)
        out_dir = pwd;
    end
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    date_str = string(session_start, 'yyyyMMdd');
    time_str = string(session_start, 'HHmmss');
    base_name = sprintf('%s_%s_%s', subject_id, date_str, time_str);

    gaps_file = fullfile(out_dir, sprintf('%s_gaps.csv', base_name));
    if cfg.general.overwrite || exist(gaps_file, 'file') ~= 2
        writetable(gaps_table, gaps_file);
    end

    ss_iso = session_start;
    ss_iso.Format = cfg.general.iso8601_format;
    iso_str = char(ss_iso);
    unix_val = posixtime(session_start);

    manifest_rows = cell(0, 12);
    for i = 1:numel(labels)
        label = labels{i};
        region = regions{i};

        ch_idx = find(strcmp(cellstr(info.SignalLabels), label), 1);
        spr = info.NumSamples(ch_idx);
        fs_ch = spr / record_duration_s;
        unit_raw = char(info.PhysicalDimensions(ch_idx));
        [scale, units_out, recognized] = resolve_unit_scale(unit_raw, cfg.edf.unit_aliases);
        if ~recognized
            warning('edf_import:UnrecognizedUnit', ...
                'Channel "%s" has PhysicalDimensions "%s", not a recognized voltage unit; exporting unscaled.', ...
                label, unit_raw);
        end

        [sig, n_samples, n_valid] = place_samples(hdr.(hdr_var_names{i}), spr, fs_ch, nrec, is_real_gap, gap_duration_s, scale);
        total_duration_s = n_samples / fs_ch;
        valid_duration_s = n_valid / fs_ch;
        n_gaps = height(gaps_table);

        digital_span = info.DigitalMax(ch_idx) - info.DigitalMin(ch_idx);
        physical_min_uV = info.PhysicalMin(ch_idx) * scale;
        physical_max_uV = info.PhysicalMax(ch_idx) * scale;
        if digital_span > 0
            % abs(): some EDFs declare an inverted calibration for a given
            % channel (PhysicalMin > PhysicalMax); the step SIZE is still
            % positive regardless of that polarity convention.
            quantization_step_uV = abs(physical_max_uV - physical_min_uV) / digital_span;
        else
            quantization_step_uV = NaN;
        end

        if recognized
            columns_label = 'amplitude_microvolts';
        else
            unit_tag = regexprep(lower(unit_raw), '[^a-z0-9]+', '_');
            unit_tag = regexprep(unit_tag, '(^_+|_+$)', '');
            if isempty(unit_tag)
                unit_tag = 'raw';
            end
            columns_label = ['amplitude_' unit_tag];
        end

        header_pairs = { ...
            'subject_id',             subject_id; ...
            'fs',                     fs_ch; ...
            'time_unit',              'seconds'; ...
            'region',                 region; ...
            'channel_label',          label; ...
            'source_file',            source_file_name; ...
            'source_format',          source_format; ...
            'session_start_datetime', iso_str; ...
            'session_start_unix',     unix_val; ...
            'timezone',               cfg.general.timezone; ...
            'n_samples',              n_samples; ...
            'n_valid_samples',        n_valid; ...
            'n_gaps',                 n_gaps; ...
            'total_duration_s',       total_duration_s; ...
            'valid_duration_s',       valid_duration_s; ...
            'units',                  units_out; ...
            'columns',                columns_label; ...
            'quantization_step_uV',   quantization_step_uV; ...
            'physical_min_uV',        physical_min_uV; ...
            'physical_max_uV',        physical_max_uV ...
        };

        txt_name = sprintf('%s_%s.txt', base_name, region);
        txt_path = fullfile(out_dir, txt_name);
        if cfg.general.overwrite || exist(txt_path, 'file') ~= 2
            write_lfp_txt(txt_path, header_pairs, sig);
        end

        manifest_rows(end+1, :) = {subject_id, region, label, fs_ch, txt_path, ...
            n_samples, n_valid, n_gaps, total_duration_s, valid_duration_s, units_out, source_file_name}; %#ok<AGROW>
    end

    manifest = cell2table(manifest_rows, 'VariableNames', ...
        {'subject_id', 'region', 'channel_label', 'fs', 'txt_file', 'n_samples', ...
         'n_valid_samples', 'n_gaps', 'total_duration_s', 'valid_duration_s', 'units', 'source_file'});

    file_meta = struct( ...
        'source_file', source_file_name, ...
        'source_path', edf_path, ...
        'source_format', source_format, ...
        'subject_id', subject_id, ...
        'session_start', session_start, ...
        'session_start_unix', unix_val, ...
        'timezone', cfg.general.timezone, ...
        'n_records', nrec, ...
        'record_duration_s', record_duration_s, ...
        'n_gaps', height(gaps_table), ...
        'n_jitter', qc.n_jitter, ...
        'n_anomaly', qc.n_anomaly, ...
        'output_dir', out_dir);
end

%% ======================================================================
function [labels, regions] = resolve_channels(info, cfg)
    labels_avail = cellstr(info.SignalLabels);

    switch cfg.edf.channels.mode
        case 'list'
            req = cfg.edf.channels.labels(:)';
            regions = req;
            labels = match_labels(req, labels_avail);
        case 'map'
            req = keys(cfg.edf.channels.map);
            regions = values(cfg.edf.channels.map, req);
            labels = match_labels(req, labels_avail);
        case 'all'
            labels = labels_avail(:)';
            regions = labels;
        otherwise
            error('edf_import:BadChannelMode', ...
                'cfg.edf.channels.mode must be ''list'', ''map'', or ''all'' (got ''%s'').', cfg.edf.channels.mode);
    end
end

function labels = match_labels(req, labels_avail)
    labels = cell(size(req));
    for i = 1:numel(req)
        idx = find(strcmp(labels_avail, req{i}), 1);
        if isempty(idx)
            idx = find(strcmpi(labels_avail, req{i}), 1);
            if ~isempty(idx)
                warning('edf_import:CaseMismatch', ...
                    'Channel label "%s" matched "%s" case-insensitively.', req{i}, labels_avail{idx});
            end
        end
        if isempty(idx)
            error('edf_import:ChannelNotFound', ...
                'Channel label "%s" not found in EDF. Available labels: %s', ...
                req{i}, strjoin(labels_avail, ', '));
        end
        labels{i} = labels_avail{idx};
    end
end

%% ======================================================================
function [gaps_table, is_real_gap, gap_duration_s, qc] = compute_gaps( ...
        rt, dt, record_duration_s, tol, gap_min_s, session_start)
% A candidate deviation follows edf_txt_conversion_Samara.m exactly:
% abs(dt - record_duration_s) > tol. Candidates are then split into real
% gaps (>= gap_min_s), jitter (logged only), and timing anomalies
% (negative deviation beyond gap_min_s: non-monotonic/overlapping record
% times -- logged only, never subtracted from the sample cursor).

    is_real_gap = false(size(dt));
    gap_duration_s = zeros(size(dt));
    n_jitter = 0;
    n_anomaly = 0;
    rows = cell(0, 8);
    gap_id = 0;

    deviation = dt - record_duration_s;
    is_candidate = abs(deviation) > tol;

    for k = 1:numel(dt)
        if ~is_candidate(k)
            continue;
        end
        d = deviation(k);
        if d >= gap_min_s
            is_real_gap(k) = true;
            gap_duration_s(k) = d;
            gap_id = gap_id + 1;
            start_s = rt(k) + record_duration_s;
            end_s = rt(k+1);
            rows(end+1, :) = {gap_id, start_s, end_s, end_s - start_s, ...
                rel_to_abs_time(start_s, session_start), rel_to_abs_time(end_s, session_start), ...
                k, k+1}; %#ok<AGROW>
        elseif d <= -gap_min_s
            n_anomaly = n_anomaly + 1;
            warning('edf_import:TimingAnomaly', ...
                'Record %d starts %.6f s before record %d ends (non-monotonic/overlapping record times); treated as contiguous.', ...
                k+1, -d, k);
        else
            n_jitter = n_jitter + 1;
        end
    end

    if isempty(rows)
        gaps_table = table('Size', [0 8], ...
            'VariableTypes', {'double', 'double', 'double', 'double', 'datetime', 'datetime', 'double', 'double'}, ...
            'VariableNames', {'gap_id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', 'prev_record_idx', 'next_record_idx'});
        % Stamp the same TimeZone real gap rows would carry, otherwise a
        % file with zero gaps produces an unzoned start_abs/end_abs that
        % fails to vertcat against another file's zoned gap rows.
        gaps_table.start_abs.TimeZone = session_start.TimeZone;
        gaps_table.end_abs.TimeZone = session_start.TimeZone;
    else
        gaps_table = cell2table(rows, 'VariableNames', ...
            {'gap_id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', 'prev_record_idx', 'next_record_idx'});
    end

    qc = struct('n_jitter', n_jitter, 'n_anomaly', n_anomaly);
end

%% ======================================================================
function [sig, n_samples, n_valid] = place_samples(col, spr, fs_ch, nrec, is_real_gap, gap_duration_s, scale)
% Places each record's samples contiguously, jumping the write cursor by
% round(gap_duration_s * fs_ch) NaN samples only at a real gap. This
% guarantees exact sample conservation (n_valid == nrec*spr, always) and
% exact time conservation (valid_duration_s + sum(gaps) == total_duration_s,
% always) regardless of sub-sample timestamp jitter.

    gap_samples = zeros(size(gap_duration_s));
    gap_samples(is_real_gap) = round(gap_duration_s(is_real_gap) * fs_ch);

    n_valid = nrec * spr;
    n_samples = n_valid + sum(gap_samples);

    sig = NaN(n_samples, 1);
    cursor = 0;
    for r = 1:nrec
        if r > 1 && is_real_gap(r-1)
            cursor = cursor + gap_samples(r-1);
        end
        sig(cursor+1:cursor+spr) = col{r} * scale;
        cursor = cursor + spr;
    end
end

%% ======================================================================
function [scale, units_out, recognized] = resolve_unit_scale(unit_raw, alias_map)
    key = lower(strtrim(unit_raw));
    key = strrep(key, char(181), 'u');  % µ (micro sign, U+00B5)
    key = strrep(key, char(956), 'u');  % μ (greek mu, U+03BC)

    if isKey(alias_map, key)
        scale = alias_map(key);
        units_out = 'microvolts';
        recognized = true;
    else
        scale = 1;
        units_out = unit_raw;
        recognized = false;
    end
end

%% ======================================================================
function sid = derive_subject_id(edf_path)
    [~, name, ~] = fileparts(edf_path);
    sid = regexprep(strtrim(name), '[^a-zA-Z0-9_-]', '_');
    if isempty(sid)
        sid = 'unknown_subject';
    end
end

function files_out = unique_by_name(files_in)
    if isempty(files_in)
        files_out = files_in;
        return;
    end
    full_paths = fullfile({files_in.folder}, {files_in.name});
    [~, ia] = unique(full_paths, 'stable');
    files_out = files_in(ia);
end

function out = vertcat_tables(a, b)
    if width(a) == 0
        out = b;
    elseif width(b) == 0
        out = a;
    else
        out = [a; b];
    end
end

function validate_cfg(cfg)
    assert(isfield(cfg, 'edf') && isfield(cfg.edf, 'channels') && isfield(cfg.edf.channels, 'mode'), ...
        'edf_import:BadConfig', 'cfg.edf.channels.mode is required.');
    assert(isfield(cfg, 'general') && ~isempty(cfg.general.timezone), ...
        'edf_import:BadConfig', 'cfg.general.timezone is required.');
    assert(isnumeric(cfg.edf.gap_tol_s) && cfg.edf.gap_tol_s >= 0, ...
        'edf_import:BadConfig', 'cfg.edf.gap_tol_s must be >= 0.');
    assert(isnumeric(cfg.edf.gap_min_s) && cfg.edf.gap_min_s >= 0, ...
        'edf_import:BadConfig', 'cfg.edf.gap_min_s must be >= 0.');
end
