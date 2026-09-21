function data = load_lfp_txt(file_path)
% LOAD_LFP_TXT  Generic reader for the pipeline's single-channel txt
% format (replaces functions/load_LFP_intan_txt.m). Works for both
% EDF-origin files (src/edf_import.m output) and legacy Intan-origin
% files, since it parses every "# key = value" header line dynamically
% instead of hard-coding a handful of field names.
%
% OUTPUT (struct data):
%   .signal      : column vector, amplitude (NaN preserved at gaps)
%   .fs          : sampling frequency (Hz), read from header
%   .t_rel       : column vector, seconds from 0 (t_rel(1) = 0)
%   .valid_mask  : column logical, ~isnan(signal)
%   .session_start : scalar datetime with TimeZone, or NaT if the header
%                    carries no session_start_unix / session_start_datetime
%   .meta        : struct, every header field (typed, dynamic fieldnames)
%   .header_pairs: Nx2 cell {key, raw_value_string}, original order, for
%                  lossless passthrough (e.g. clean_lfp.m rewriting a header)
%   .file        : file name
%   .folder      : folder

    if ~(ischar(file_path) || isstring(file_path))
        error('load_lfp_txt:BadInput', 'file_path must be a char/string with the full path to the .txt file.');
    end
    file_path = char(file_path);

    [folder, name, ext] = fileparts(file_path);
    if ~strcmpi(ext, '.txt')
        error('load_lfp_txt:InputError', 'Expected a .txt file, got "%s".', ext);
    end
    if exist(file_path, 'file') ~= 2
        error('load_lfp_txt:FileNotFound', 'File not found: %s', file_path);
    end

    fid = fopen(file_path, 'rt');
    if fid == -1
        error('load_lfp_txt:FileOpenError', 'Could not open file: %s', file_path);
    end

    hdr = parse_header(fid);
    C = textscan(fid, '%f');
    fclose(fid);

    signal = C{1};
    n = numel(signal);

    if ~isfield(hdr.fields, 'fs')
        error('load_lfp_txt:NoFsFound', 'No "fs = ..." found in header of %s. Refusing to guess.', file_path);
    end
    fs = str2double(hdr.fields.fs);
    if isnan(fs) || ~(fs > 0)
        error('load_lfp_txt:BadFs', 'Header "fs" in %s ("%s") is not a positive numeric value.', file_path, hdr.fields.fs);
    end

    t_rel = (0:n-1)' / fs;
    valid_mask = ~isnan(signal);

    session_start = resolve_session_start(hdr.fields);

    data = struct();
    data.signal = signal;
    data.fs = fs;
    data.t_rel = t_rel;
    data.valid_mask = valid_mask;
    data.session_start = session_start;
    data.meta = hdr.fields;
    data.header_pairs = hdr.pairs;
    data.file = [name ext];
    data.folder = folder;
end

function session_start = resolve_session_start(fields)
    iso_fmt = 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX';

    tz = '';
    if isfield(fields, 'timezone') && (ischar(fields.timezone) || isstring(fields.timezone)) ...
            && strlength(string(fields.timezone)) > 0
        tz = char(fields.timezone);
    end

    if isfield(fields, 'session_start_unix') && ~isempty(tz)
        unix_val = str2double(fields.session_start_unix);
        if ~isnan(unix_val)
            session_start = datetime(unix_val, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
            return;
        end
    end

    if isfield(fields, 'session_start_datetime') && (ischar(fields.session_start_datetime) || isstring(fields.session_start_datetime))
        raw = char(fields.session_start_datetime);
        try
            session_start = datetime(raw, 'InputFormat', iso_fmt, 'TimeZone', 'UTC');
            if ~isempty(tz)
                session_start.TimeZone = tz;
            end
            return;
        catch
            warning('load_lfp_txt:BadSessionStart', ...
                'Could not parse session_start_datetime "%s"; leaving session_start as NaT.', raw);
        end
    end

    session_start = NaT;
end
