function T = load_recording_log(log_file)
% LOAD_RECORDING_LOG  Read the per-recording channel log (EEG_recording_log.xlsx)
% used by edf_import.m's cfg.edf.channels.mode = 'log'.
%
%   T = load_recording_log(log_file)
%
% Expected columns (first sheet): Filename, Animal ID, Port, HPCr_channel,
% HPCl_channel (others are ignored). Returned T has char/cellstr columns
% filename, animal_id, port, hpcr, hpcl, whitespace-trimmed.
%
% The same Filename appears once per animal recorded simultaneously, so a
% row is identified by (animal_id, filename), never by filename alone.
%
% Cached per file (keyed by path + modification date) because edf_import
% calls this once per EDF and the log lives on a network share.

    persistent cache_key cache_T

    if exist(log_file, 'file') ~= 2
        error('load_recording_log:FileNotFound', 'Recording log not found: %s', log_file);
    end
    d = dir(log_file);
    key = sprintf('%s|%.10f', log_file, d.datenum);
    if ~isempty(cache_key) && strcmp(cache_key, key)
        T = cache_T;
        return;
    end

    opts = detectImportOptions(log_file, 'VariableNamingRule', 'preserve');
    required = {'Filename', 'Animal ID', 'Port', 'HPCr_channel', 'HPCl_channel'};
    missing = setdiff(required, opts.VariableNames);
    if ~isempty(missing)
        error('load_recording_log:BadColumns', 'Recording log %s is missing column(s): %s', ...
            log_file, strjoin(missing, ', '));
    end
    opts.SelectedVariableNames = required;
    opts = setvartype(opts, required, 'char');
    raw = readtable(log_file, opts);

    T = table(strtrim(raw.('Filename')), strtrim(raw.('Animal ID')), strtrim(raw.('Port')), ...
        strtrim(raw.('HPCr_channel')), strtrim(raw.('HPCl_channel')), ...
        'VariableNames', {'filename', 'animal_id', 'port', 'hpcr', 'hpcl'});
    T = T(~cellfun(@isempty, T.filename), :);

    cache_key = key;
    cache_T = T;
end
