function [T, notes] = spec_read_seizures_table(file, tz)
% SPEC_READ_SEIZURES_TABLE  seizures_events of the pipeline (05_summaries
% CSV, the 'seizures_events' sheet of pipeline_summary.xlsx) or a reviewed
% copy of it (+ user_check / comments), with forced types.
%
%   [T, notes] = spec_read_seizures_table(file, tz)
%
% T columns: row (data row, 1 = first row under the header), subject_id,
% region, source_file, seizure_id, start_s, end_s, duration_s, event_id,
% category, detection_status, ll_status, user_check (upper case, trimmed),
% comments, session_start, start_abs, end_abs (zoned datetime, NaT when
% unreadable).
%
% - subject_id and every identifier are read as TEXT ('097' stays '097').
% - A CSV is first read with read_pipeline_csv.m (the pipeline's own typed
%   reader). A copy re-saved by Excel usually breaks it (session_start /
%   start_abs become e.g. "00:48.6"); then the file is read as text and
%   those datetimes are NaT -- harmless here, because event times are
%   always rebuilt as txt-header session_start + start_s (see notes).
% - user_check / comments are found by header NAME (not by column letter:
%   reviewed files have them in AE/AF, AF/AG or an unnamed column): the
%   first 'user_check' column, then 'comments' / 'comment', else an
%   unnamed column right after user_check.

    if nargin < 2 || isempty(tz)
        tz = 'Europe/Paris';
    end
    if exist(file, 'file') ~= 2
        error('spec_read_seizures_table:FileNotFound', 'File not found: %s', file);
    end
    notes = {};
    [~, ~, ext] = fileparts(file);
    is_csv = strcmpi(ext, '.csv');

    if is_csv
        R = spec_read_text_csv(file);
    else
        sheets = sheetnames(file);
        sheet = sheets{1};
        if any(strcmp(sheets, 'seizures_events'))
            sheet = 'seizures_events';
        end
        opts = detectImportOptions(file, 'Sheet', sheet, 'VariableNamingRule', 'preserve');
        opts = setvartype(opts, opts.VariableNames, 'char');
        R = readtable(file, opts);
        notes{end+1} = sprintf('Excel sheet ''%s''', sheet);
    end
    names = R.Properties.VariableNames;
    lnames = lower(strtrim(names));

    req = {'subject_id', 'region', 'source_file', 'seizure_id', 'start_s', 'end_s'};
    miss = setdiff(req, lnames);
    if ~isempty(miss)
        error('spec_read_seizures_table:MissingColumn', 'Column(s) missing in %s: %s', file, strjoin(miss, ', '));
    end
    col = @(nm) R.(names{find(strcmp(lnames, nm), 1)});
    keep = ~cellfun(@(s) isempty(strtrim(s)), col('subject_id'));
    R = R(keep, :);
    col = @(nm) R.(names{find(strcmp(lnames, nm), 1)});
    n = height(R);

    T = table();
    T.row = find(keep);
    T.subject_id = strtrim(col('subject_id'));
    T.region = strtrim(col('region'));
    T.source_file = strtrim(col('source_file'));
    T.seizure_id = str2double(col('seizure_id'));
    T.start_s = str2double(col('start_s'));
    T.end_s = str2double(col('end_s'));
    T.duration_s = opt_num('duration_s', T.end_s - T.start_s);
    T.event_id = opt_num('event_id', NaN(n, 1));
    T.category = opt_text('category');
    T.detection_status = opt_text('detection_status');
    T.ll_status = opt_text('ll_status');

    % ---- labels by header name ---------------------------------------
    iu = find(strcmp(lnames, 'user_check'), 1);
    ic = find(strcmp(lnames, 'comments') | strcmp(lnames, 'comment'), 1);
    if isempty(ic) && ~isempty(iu) && iu < numel(names) && ~isempty(regexp(names{iu + 1}, '^(Var|ExtraVar)\d+$', 'once'))
        ic = iu + 1;
        notes{end+1} = sprintf('comments taken from the unnamed column %d (right after user_check)', ic);
    end
    if isempty(iu)
        T.user_check = repmat({''}, n, 1);
        notes{end+1} = 'no user_check column';
    else
        T.user_check = upper(strtrim(R.(names{iu})));
        notes{end+1} = sprintf('user_check = column %d (%s)', iu, excel_col(iu));
    end
    if isempty(ic)
        T.comments = repmat({''}, n, 1);
    else
        T.comments = strtrim(R.(names{ic}));
        notes{end+1} = sprintf('comments = column %d (%s)', ic, excel_col(ic));
    end

    % ---- datetimes: typed reader first, text fallback -----------------
    T.session_start = zoned_nat(n, tz);
    T.start_abs = zoned_nat(n, tz);
    T.end_abs = zoned_nat(n, tz);
    typed = false;
    if is_csv
        try
            P = read_pipeline_csv(file, 'seizures_events', tz);
            if height(P) == n
                T.session_start = P.session_start;
                T.start_abs = P.start_abs;
                T.end_abs = P.end_abs;
                typed = true;
            else
                notes{end+1} = sprintf('read_pipeline_csv returned %d rows, text reader %d; datetimes not taken from it', height(P), n);
            end
        catch ME
            notes{end+1} = sprintf(['read_pipeline_csv could not read it (%s) -- typical of a CSV re-saved by Excel; ' ...
                'session_start/start_abs ignored, times rebuilt from the txt header + start_s'], ME.message);
        end
    end
    if ~typed
        for nm = {'session_start', 'start_abs', 'end_abs'}
            if any(strcmp(lnames, nm{1}))
                T.(nm{1}) = parse_dt(col(nm{1}), tz);
            end
        end
    end

    %% nested helpers (share R / names / lnames / n)
    function v = opt_num(nm, default)
        i = find(strcmp(lnames, nm), 1);
        if isempty(i)
            v = default;
        else
            v = str2double(R.(names{i}));
        end
    end
    function v = opt_text(nm)
        i = find(strcmp(lnames, nm), 1);
        if isempty(i)
            v = repmat({''}, n, 1);
        else
            v = strtrim(R.(names{i}));
        end
    end
end

function d = parse_dt(c, tz)
    d = zoned_nat(numel(c), tz);
    fmts = {'dd-MMM-yyyy HH:mm:ss.SSS', 'dd-MMM-yyyy HH:mm:ss', 'dd/MM/yyyy HH:mm:ss.SSS', ...
        'dd/MM/yyyy HH:mm:ss', 'yyyy-MM-dd HH:mm:ss.SSS', 'yyyy-MM-dd HH:mm:ss'};
    for k = 1:numel(fmts)
        todo = isnat(d) & ~cellfun(@isempty, c);
        if ~any(todo)
            break;
        end
        try
            v = datetime(c(todo), 'InputFormat', fmts{k}, 'Locale', 'en_US', 'TimeZone', tz);
            d(todo) = v;
        catch
        end
    end
end

function d = zoned_nat(n, tz)
    d = NaT(n, 1);
    d.TimeZone = tz;
end

function s = excel_col(i)
    s = '';
    while i > 0
        r = mod(i - 1, 26);
        s = [char(65 + r) s]; %#ok<AGROW>
        i = floor((i - 1) / 26);
    end
end
