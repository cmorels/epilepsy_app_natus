function events = spec_read_event_list(src)
% SPEC_READ_EVENT_LIST  Events chosen by hand: the table in
% run_spectrograms.m (cell array) or an Excel / CSV file with the same
% columns.
%
%   events = spec_read_event_list(manual_cell)
%   events = spec_read_event_list('spectrogram_events.xlsx')
%
% Columns, in this order in the cell array (by NAME in a file):
%   id | tipo ('seizure'|'iid') | subject | sesion ('yyyyMMdd_HHmmss' or EDF
%   name) | region ('HPCr'|'HPCl'|'both') | inicio | fin | time_ref
%   ('natus'|'pipeline'|'t_rel')
% inicio / fin: 'dd/MM/yyyy HH:mm:ss.SSS' in the time_ref clock, or a
% number = t_rel in seconds.
%
% Optional per-event overrides of cfg.spec (see spec_run.m for the list):
%   cell array : a 9th column with a struct, e.g. struct('win_s', 4), or
%                text 'win_s=4; fmax=60'
%   file       : any extra column whose header is an override name
%                (win_s, fmax, baseline, spike_source, signal_source, ...);
%                empty cells are ignored.
% 'peaks' (manual spike times for spike_source = 'manual'): numbers (t_rel)
% or times in the event's time_ref clock, separated by ';'.

    if iscell(src)
        events = from_cell(src);
    elseif ischar(src) || isstring(src)
        events = from_file(char(src));
    else
        error('spec_read_event_list:BadInput', 'Give the manual cell table or a .xlsx/.csv path.');
    end
end

%% ======================================================================
function events = from_cell(C)
    events = struct([]);
    if isempty(C)
        return;
    end
    if size(C, 2) < 8
        error('spec_read_event_list:BadTable', ...
            'The manual table needs 8 columns: id tipo subject sesion region inicio fin time_ref (+ optional overrides).');
    end
    for i = 1:size(C, 1)
        ov = struct();
        if size(C, 2) >= 9 && ~isempty(C{i, 9})
            if isstruct(C{i, 9})
                ov = C{i, 9};
            else
                ov = parse_kv(C{i, 9});
            end
        end
        events = [events, make(C{i, 1}, C{i, 2}, C{i, 3}, C{i, 4}, C{i, 5}, C{i, 6}, C{i, 7}, C{i, 8}, ov, ...
            sprintf('manual row %d', i))]; %#ok<AGROW>
    end
end

function events = from_file(file)
    if exist(file, 'file') ~= 2
        error('spec_read_event_list:FileNotFound', 'File not found: %s', file);
    end
    [~, ~, ext] = fileparts(file);
    req = {'id', 'tipo', 'subject', 'sesion', 'region', 'inicio', 'fin', 'time_ref'};
    if strcmpi(ext, '.csv')
        R = spec_read_text_csv(file);
    else
        opts = detectImportOptions(file, 'VariableNamingRule', 'preserve');
        % everything as text except inicio/fin, which may be real Excel dates
        ln = lower(strtrim(opts.VariableNames));
        timecols = opts.VariableNames(ismember(ln, {'inicio', 'fin'}));
        opts = setvartype(opts, setdiff(opts.VariableNames, timecols, 'stable'), 'char');
        R = readtable(file, opts);
    end
    names = R.Properties.VariableNames;
    miss = setdiff(req, lower(strtrim(names)));
    if ~isempty(miss)
        error('spec_read_event_list:MissingColumn', 'Column(s) missing in %s: %s', file, strjoin(miss, ', '));
    end
    lnames = lower(strtrim(names));
    get = @(i, nm) R.(names{find(strcmp(lnames, nm), 1)})(i);

    events = struct([]);
    for i = 1:height(R)
        id = get(i, 'id');
        if isempty(strtrim(id{1}))
            continue;
        end
        ov = struct();
        for j = 1:numel(names)
            if ismember(lnames{j}, req)
                continue;
            end
            v = R.(names{j})(i);
            if iscell(v)
                v = v{1};
            end
            if (ischar(v) && isempty(strtrim(v))) || (isnumeric(v) && all(isnan(v)))
                continue;
            end
            ov.(lnames{j}) = parse_value(lnames{j}, v);
        end
        events = [events, make(id{1}, cell1(get(i, 'tipo')), cell1(get(i, 'subject')), cell1(get(i, 'sesion')), ...
            cell1(get(i, 'region')), cell1(get(i, 'inicio')), cell1(get(i, 'fin')), cell1(get(i, 'time_ref')), ov, ...
            sprintf('%s row %d', file, i + 1))]; %#ok<AGROW>
    end
end

function v = cell1(v)
    if iscell(v)
        v = v{1};
    end
end

function ev = make(id, tipo, subject, sesion, region, inicio, fin, time_ref, ov, origin)
    peaks = [];
    if isfield(ov, 'peaks')
        peaks = ov.peaks;
        ov = rmfield(ov, 'peaks');
    end
    ev = spec_event(id, tipo, subject, sesion, region, inicio, fin, time_ref, ov, origin);
    ev.peaks = peaks;
end

function ov = parse_kv(s)
    ov = struct();
    % split only before "key=", so a value may itself hold ';' (peaks=1.2;3.4)
    parts = regexp(char(s), ';\s*(?=\w+\s*=)', 'split');
    for k = 1:numel(parts)
        p = strtrim(regexprep(parts{k}, ';\s*$', ''));
        if isempty(p)
            continue;
        end
        kv = regexp(p, '^\s*(\w+)\s*=\s*(.*)$', 'tokens', 'once');
        if isempty(kv)
            error('spec_read_event_list:BadOverride', 'Override "%s" is not key=value.', p);
        end
        ov.(lower(kv{1})) = parse_value(lower(kv{1}), kv{2});
    end
end

function v = parse_value(name, v)
% Text cell -> number / vector / list of times, depending on the override.
    if ~ischar(v) && ~isstring(v)
        return;
    end
    v = strtrim(char(v));
    if any(strcmp(name, {'peaks', 'baseline'})) && contains(v, ';') && ~contains(v, '[')
        parts = strtrim(strsplit(v, ';'));
        parts = parts(~cellfun(@isempty, parts));
        nums = str2double(parts);
        if all(~isnan(nums))
            v = nums;
        else
            v = parts;
        end
        return;
    end
    if ~isempty(regexp(v, '^\[?[\s\d.eE+,-]+\]?$', 'once'))
        num = sscanf(regexprep(v, '[\[\],]', ' '), '%f')';
        if ~isempty(num)
            v = num;
        end
    end
end
