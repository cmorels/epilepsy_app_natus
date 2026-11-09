function [events, chosen] = spec_select_seizures(sel, cfg)
% SPEC_SELECT_SEIZURES  Seizure events from a seizures_events table (the
% pipeline's, or your reviewed copy), by ID or by filter.
%
%   [events, chosen] = spec_select_seizures(sel, cfg)
%
% sel.seizures_file : seizures_events .csv / .xlsx (spec_read_seizures_table.m)
% sel.region        : 'both' (default: one event per seizure, every region
%                     that passes the filters, one figure) | 'HPCr' | 'HPCl'
% sel.id_field      : 'seizure_id' (default) | 'event_id'
%
% B. by ID: sel.ids = {subject, sesion, [ids]; ...}   (one row per session)
%    sesion = 'yyyyMMdd_HHmmss' or the EDF name (source_file).
% C. by filter (every field optional, combined with AND):
%    sel.filter.subject     = {'005-s', '128-s'}
%    sel.filter.sesion      = {'20260706_130048', '<EDF name>'}
%    sel.filter.region      = 'HPCr' | {'HPCr','HPCl'}
%    sel.filter.user_check  = 'TP' | {'TP', ...}   (case-insensitive)
%    sel.filter.category    = {...}
%    sel.filter.detection_status = {...}
%    sel.filter.date_range  = {'dd/MM/yyyy HH:mm:ss', 'dd/MM/yyyy HH:mm:ss'}  (pipeline clock)
%    sel.filter.duration_s  = [min max]
%    sel.filter.max_n       = N per subject (earliest first)
%    A row passes when IT passes; with sel.region = 'both' the regions of a
%    seizure that pass are grouped (a seizure TP in HPCr only gives an
%    HPCr-only event).
% Give sel.ids OR sel.filter (ids win if both are set).
%
% Onset / offset are the rows' start_s / end_s (t_rel): the event window is
% [min start_s, max end_s] over its regions, time_ref = 't_rel', i.e. the
% pipeline clock (start_abs) -- never the Excel text columns. events(k).sz
% keeps each region's row (seizure_id, event_id, user_check, comments, ...).
%
% chosen: the selected rows of the table (+ session_key, onset_pipeline).

    if ~isfield(sel, 'seizures_file') || isempty(sel.seizures_file)
        error('spec_select_seizures:NoFile', 'sel.seizures_file is required.');
    end
    s = cfg.spec;
    region_mode = field_or(sel, 'region', 'both');
    id_field = field_or(sel, 'id_field', 'seizure_id');
    [T, notes] = spec_read_seizures_table(sel.seizures_file, s.tz);
    for k = 1:numel(notes)
        spec_log('seizures table: %s', notes{k});
    end

    % session key + pipeline onset of each row, via the txt headers (cheap)
    T.session_key = repmat({''}, height(T), 1);
    T.onset_pipeline = NaT(height(T), 1);
    T.onset_pipeline.TimeZone = s.tz;
    T.resolve_error = repmat({''}, height(T), 1);
    combos = unique(strcat(T.subject_id, '|', T.source_file, '|', T.region));
    for k = 1:numel(combos)
        parts = strsplit(combos{k}, '|');
        m = strcmp(T.subject_id, parts{1}) & strcmp(T.source_file, parts{2}) & strcmp(T.region, parts{3});
        try
            r = spec_resolve_files(s.pipeline_root, parts{1}, parts{2}, parts{3}, 'clean', s.tz);
            T.session_key(m) = {r.session_key};
            T.onset_pipeline(m) = rel_to_abs_time(T.start_s(m), r.session_start);
        catch ME
            T.resolve_error(m) = {ME.message};
        end
    end

    if ~strcmp(region_mode, 'both')
        T = T(strcmp(T.region, region_mode), :);
    end

    if isfield(sel, 'ids') && ~isempty(sel.ids)
        ids = sel.ids;
        if size(ids, 2) ~= 3
            error('spec_select_seizures:BadIds', 'sel.ids must be {subject, sesion, [ids]; ...}.');
        end
        keep = false(height(T), 1);
        for k = 1:size(ids, 1)
            keep = keep | (strcmp(T.subject_id, ids{k, 1}) & session_match(T, ids{k, 2}) & ...
                ismember(T.(id_field), ids{k, 3}));
            found = unique(T.(id_field)(strcmp(T.subject_id, ids{k, 1}) & session_match(T, ids{k, 2})));
            missing = setdiff(ids{k, 3}, found);
            if ~isempty(missing)
                spec_log('WARNING sel.ids: %s %s has no %s %s in %s', ids{k, 1}, ids{k, 2}, id_field, ...
                    mat2str(missing), sel.seizures_file);
            end
        end
        origin = 'seizures_ids';
    elseif isfield(sel, 'filter') && ~isempty(sel.filter)
        keep = apply_filter(T, sel.filter, s.tz);
        origin = 'seizures_filter';
    else
        error('spec_select_seizures:NoSelection', 'Give sel.ids or sel.filter.');
    end
    chosen = T(keep, :);

    % ---- group rows into events ----------------------------------------
    events = struct([]);
    if height(chosen) == 0
        spec_log('%s: no row of %s matches the selection', origin, sel.seizures_file);
        return;
    end
    idtxt = arrayfun(@(v) sprintf('%.0f', v), chosen.(id_field), 'UniformOutput', false);
    if strcmp(region_mode, 'both')
        gkey = strcat(chosen.subject_id, '|', chosen.source_file, '|', idtxt);
    else
        gkey = strcat(chosen.subject_id, '|', chosen.source_file, '|', idtxt, '|', chosen.region);
    end
    [ukeys, first] = unique(gkey, 'stable');
    order = sortrows(table(chosen.subject_id(first), chosen.source_file(first), chosen.start_s(first), (1:numel(first))', ...
        'VariableNames', {'s', 'f', 't', 'k'}), {'s', 'f', 't'});
    ukeys = ukeys(order.k);

    prefix = 'sz';
    if strcmp(id_field, 'event_id')
        prefix = 'ev';
    end
    for k = 1:numel(ukeys)
        rows = chosen(strcmp(gkey, ukeys{k}), :);
        rows = sortrows(rows, 'region', 'descend');   % HPCr before HPCl
        if height(rows) > 1
            region = 'both';
        else
            region = rows.region{1};
        end
        id = sprintf('%s%03d', prefix, rows.(id_field)(1));
        ev = spec_event(id, 'seizure', rows.subject_id{1}, rows.source_file{1}, region, ...
            min(rows.start_s), max(rows.end_s), 't_rel', struct(), origin);
        ev.sz = table2struct(rows(:, {'region', 'seizure_id', 'event_id', 'start_s', 'end_s', 'duration_s', ...
            'user_check', 'comments', 'category', 'detection_status', 'll_status', 'row', 'start_abs', 'source_file'}));
        events = [events, ev]; %#ok<AGROW>
    end

    if isfield(sel, 'filter') && isfield(sel.filter, 'max_n') && ~isempty(sel.filter.max_n) && ~(isfield(sel, 'ids') && ~isempty(sel.ids))
        subs = unique({events.subject}, 'stable');
        keep_ev = false(1, numel(events));
        for k = 1:numel(subs)
            idx = find(strcmp({events.subject}, subs{k}));
            keep_ev(idx(1:min(numel(idx), sel.filter.max_n))) = true;
        end
        events = events(keep_ev);
    end
    spec_log('%s: %d row(s) selected -> %d seizure event(s) (region mode %s)', origin, height(chosen), numel(events), region_mode);
end

%% ======================================================================
function m = session_match(T, sesion)
    sesion = strtrim(char(sesion));
    if ~isempty(regexp(sesion, '^\d{8}_\d{6}$', 'once'))
        m = strcmp(T.session_key, sesion);
    else
        [~, stem] = fileparts(sesion);
        stems = cellfun(@stem_of, T.source_file, 'UniformOutput', false);
        m = strcmpi(stems, stem);
    end
end

function s = stem_of(f)
    [~, s] = fileparts(f);
end

function keep = apply_filter(T, f, tz)
    keep = true(height(T), 1);
    if isfield(f, 'subject') && ~isempty(f.subject)
        keep = keep & ismember(T.subject_id, cellstr(f.subject));
    end
    if isfield(f, 'sesion') && ~isempty(f.sesion)
        ses = cellstr(f.sesion);
        m = false(height(T), 1);
        for k = 1:numel(ses)
            m = m | session_match(T, ses{k});
        end
        keep = keep & m;
    end
    if isfield(f, 'region') && ~isempty(f.region)
        keep = keep & ismember(T.region, cellstr(f.region));
    end
    if isfield(f, 'user_check') && ~isempty(f.user_check)
        keep = keep & ismember(T.user_check, upper(cellstr(f.user_check)));
    end
    if isfield(f, 'category') && ~isempty(f.category)
        keep = keep & ismember(T.category, cellstr(f.category));
    end
    if isfield(f, 'detection_status') && ~isempty(f.detection_status)
        keep = keep & ismember(T.detection_status, cellstr(f.detection_status));
    end
    if isfield(f, 'duration_s') && ~isempty(f.duration_s)
        d = T.end_s - T.start_s;
        keep = keep & d >= f.duration_s(1) & d <= f.duration_s(2);
    end
    if isfield(f, 'date_range') && ~isempty(f.date_range)
        a = to_zoned(f.date_range{1}, tz);
        b = to_zoned(f.date_range{2}, tz);
        keep = keep & ~isnat(T.onset_pipeline) & T.onset_pipeline >= a & T.onset_pipeline <= b;
    end
end

function d = to_zoned(v, tz)
    if isdatetime(v)
        d = v;
    else
        d = spec_resolve_times('parse', v);
    end
    if isempty(d.TimeZone)
        d.TimeZone = tz;
    end
end

function v = field_or(s, name, default)
    if isfield(s, name) && ~isempty(s.(name))
        v = s.(name);
    else
        v = default;
    end
end
