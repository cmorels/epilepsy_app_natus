function S = cs_build_seizures(cfg, ev, rv, edf_tbl, edf_lists)
% CS_BUILD_SEIZURES  Final seizure set (section 5.1 of the task), the
% crisette_focal set, the scope of reviewed mouse x campaign units, the
% FP event set for detector performance, plus exclusion_log and flags.
%
% Scope: a (campaign, mouse) unit is IN when every one of its detected
% events carries a TP/FP label (event level: an unlabeled row takes its
% sibling rows' label) or it has no detections at all; units with no or
% partial labels are OUT (logged). Events of excluded EDFs (rv.excluded,
% cs_excluded_edfs.m) count nowhere, not for the scope either; crisette_focal
% events need no label. Group: control if in cfg.controls,
% epileptic if in cfg.epileptic (or, when that is empty, any other IN mouse).

    tz = cfg.tz;
    S.flags = cell(0, 3);      % severity, topic, detail
    S.excl = cell(0, 8);       % type, id, mouse, campaign, start, end, rule, reason
    loc = @(p) char(datetime(p, 'ConvertFrom', 'posixtime', 'TimeZone', tz, 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));

    % ---------------- excluded EDFs (cs_excluded_edfs.m) ----------------------
    % Their events count nowhere, and not for the scope either (a unit of
    % cfg.edf_level_units stays in with its fully labeled EDFs).
    disc = rv.excluded;
    ev.discarded = false(height(ev), 1);
    for k = 1:height(disc)
        hit = strcmp(ev.mouse, disc.mouse{k}) & strcmp(ev.edf, disc.edf{k});
        ev.discarded = ev.discarded | hit;
        S.excl(end+1, :) = {['EDF excluido (' disc.kind{k} ')'], disc.id{k}, disc.mouse{k}, disc.campaign{k}, loc(disc.start_posix(k)), ...
            loc(disc.end_posix(k)), disc.id{k}, sprintf('EDF %s fuera de tiempo, crisis, IIDs y desempeno; %d eventos detectados dentro (%d TP). %s', ...
            disc.edf{k}, sum(hit), sum(hit & strcmp(ev.label, 'TP')), disc.reason{k})};
    end

    % ---------------- scope per (campaign, mouse) ---------------------------
    % crisette_focal events (cfg.crisette) need no AE (5_Observaciones classifies them)
    is_cr = false(height(ev), 1);
    for k = 1:size(cfg.crisette, 1)
        c = cfg.crisette(k, :);
        is_cr = is_cr | (strcmp(ev.campaign, c{3}) & strcmp(ev.mouse, c{2}) & strcmp(ev.edf, c{4}) & ev.seizure_id == c{5});
    end
    units = unique(edf_tbl(:, {'campaign', 'mouse'}), 'rows');
    n = height(units);
    units.n_edfs = zeros(n, 1);
    units.n_events = zeros(n, 1);
    units.n_events_unlabeled = zeros(n, 1);
    units.n_events_conflict = zeros(n, 1);
    units.status = repmat({''}, n, 1);
    for i = 1:n
        sel = strcmp(ev.campaign, units.campaign{i}) & strcmp(ev.mouse, units.mouse{i}) & ~ev.discarded;
        units.n_edfs(i) = sum(strcmp(edf_tbl.campaign, units.campaign{i}) & strcmp(edf_tbl.mouse, units.mouse{i}));
        units.n_events(i) = sum(sel);
        units.n_events_unlabeled(i) = sum(sel & cellfun(@isempty, ev.label) & ~is_cr);
        units.n_events_conflict(i) = sum(sel & strcmp(ev.label, 'CONFLICT'));
        if units.n_events(i) == 0
            units.status{i} = 'no_detections';
        elseif units.n_events_unlabeled(i) == 0
            units.status{i} = 'reviewed';
        elseif units.n_events_unlabeled(i) == units.n_events(i)
            units.status{i} = 'not_reviewed';
        else
            units.status{i} = 'partial';
        end
    end
    units.included = ismember(units.status, {'reviewed', 'no_detections'});
    units.group = repmat({''}, n, 1);
    for i = 1:n
        if ismember(units.mouse{i}, cfg.controls)
            units.group{i} = 'control';
        elseif isempty(cfg.epileptic) || ismember(units.mouse{i}, cfg.epileptic)
            units.group{i} = 'epileptic';
        else
            units.group{i} = 'unassigned';
            units.included(i) = false;
        end
        if ~units.included(i)
            S.excl(end+1, :) = {'unidad raton x campaign', [units.campaign{i} '|' units.mouse{i}], units.mouse{i}, ...
                units.campaign{i}, '', '', 'alcance', sprintf('%s: %d de %d eventos sin etiqueta AE; fuera de las estadisticas', ...
                units.status{i}, units.n_events_unlabeled(i), units.n_events(i))};
        end
    end
    for i = find(units.n_events_conflict > 0)'
        S.flags(end+1, :) = {'decidir', 'AE contradictorio entre hemisferios', sprintf('%s %s: %d eventos con TP y FP a la vez', ...
            units.campaign{i}, units.mouse{i}, units.n_events_conflict(i))};
    end
    S.units = units;
    in_unit = @(c, m) any(strcmp(units.campaign, c) & strcmp(units.mouse, m) & units.included);

    % ---------------- crisette_focal --------------------------------------
    ev.crisette = repmat({''}, height(ev), 1);
    cr = cell(0, 12);
    for k = 1:size(cfg.crisette, 1)
        c = cfg.crisette(k, :);
        j = find(strcmp(ev.campaign, c{3}) & strcmp(ev.mouse, c{2}) & strcmp(ev.edf, c{4}) & ev.seizure_id == c{5});
        if numel(j) ~= 1
            error('cs_build_seizures:Crisette', '%s: event %s %s %s #%d not found', c{1}, c{2}, c{3}, c{4}, c{5});
        end
        n_in_edf = sum(strcmp(ev.campaign, c{3}) & strcmp(ev.mouse, c{2}) & strcmp(ev.edf, c{4}));
        ev.crisette{j} = c{1};
        cr(end+1, :) = {c{1}, c{2}, c{3}, c{4}, c{5}, loc(ev.start_posix(j)), loc(ev.end_posix(j)), ev.duration_s(j), ...
            ev.label{j}, ev.laterality{j}, n_in_edf, ev.start_posix(j)}; %#ok<AGROW>
        S.excl(end+1, :) = {'evento', ev.event_id{j}, c{2}, c{3}, loc(ev.start_posix(j)), loc(ev.end_posix(j)), ...
            'crisette_focal', sprintf('%s: categoria aparte (AE=%s), no suma al total de crisis', c{1}, ev.label{j})};
    end
    S.crisette = cell2table(cr, 'VariableNames', {'obs_id', 'mouse', 'campaign', 'edf', 'seizure_id', 'start_local', ...
        'end_local', 'duration_s', 'AE_label', 'laterality', 'n_detections_in_edf', 'start_posix'});

    % ---------------- manual merges (O-02) ---------------------------------
    ev.merged_into = zeros(height(ev), 1);
    for k = 1:size(cfg.merge_events, 1)
        c = cfg.merge_events(k, :);
        j = find(strcmp(ev.campaign, c{3}) & strcmp(ev.mouse, c{2}) & strcmp(ev.edf, c{4}) & ismember(ev.seizure_id, c{5}));
        if numel(j) ~= numel(c{5})
            error('cs_build_seizures:Merge', '%s: events not found', c{1});
        end
        j0 = j(1);
        ev.start_posix(j0) = min(ev.start_posix(j));
        ev.end_posix(j0) = max(ev.end_posix(j));
        ev.duration_s(j0) = ev.end_posix(j0) - ev.start_posix(j0);
        ev.merged_into(j(2:end)) = j0;
        S.flags(end+1, :) = {'info', c{1}, sprintf('%s %s seizure_id %s: %d evento(s) bilateral(es), %d filas (fragmentos) -> 1 crisis %s - %s', ...
            c{2}, c{3}, mat2str(c{5}), numel(j), sum(ev.n_rows(j)), loc(ev.start_posix(j0)), loc(ev.end_posix(j0)))};
    end

    % ---------------- event sets -------------------------------------------
    ev.in_scope = arrayfun(@(i) in_unit(ev.campaign{i}, ev.mouse{i}), (1:height(ev))');
    ev.usable = ev.in_scope & ~ev.discarded & cellfun(@isempty, ev.crisette) & ev.merged_into == 0;
    tp = ev(ev.usable & ismember(ev.label, {'TP', 'CONFLICT'}), :);
    S.fp_events = ev(ev.usable & strcmp(ev.label, 'FP'), :);
    S.events = ev;
    for i = find(strcmp(tp.label, 'CONFLICT'))'
        S.flags(end+1, :) = {'decidir', 'R?', sprintf('%s: TP y FP en hemisferios distintos; contado como TP', tp.event_id{i})};
    end
    for i = find(~cellfun(@isempty, tp.comments))'
        S.flags(end+1, :) = {'revisar', 'R7 TP con texto en AF', sprintf('%s (%s): AF="%s"; se mantiene como crisis', ...
            tp.event_id{i}, loc(tp.start_posix(i)), tp.comments{i})};
    end
    for i = 1:height(S.fp_events)
        S.excl(end+1, :) = {'evento', S.fp_events.event_id{i}, S.fp_events.mouse{i}, S.fp_events.campaign{i}, ...
            loc(S.fp_events.start_posix(i)), loc(S.fp_events.end_posix(i)), 'FP (AE)', ...
            ['no es crisis; AF=' S.fp_events.comments{i}]};
    end
    for i = find(ev.discarded & ev.in_scope)'
        S.excl(end+1, :) = {'evento', ev.event_id{i}, ev.mouse{i}, ev.campaign{i}, loc(ev.start_posix(i)), ...
            loc(ev.end_posix(i)), 'EDF descartado', ['AE=' ev.label{i}]};
    end

    % ---------------- seizure table: detected TP ---------------------------
    sz = table(tp.mouse, tp.campaign, tp.edf, tp.start_posix, tp.end_posix, repmat({'detected_TP'}, height(tp), 1), ...
        tp.event_id, tp.laterality, 'VariableNames', {'mouse', 'campaign', 'edf', 'start_posix', 'end_posix', 'source', ...
        'origin_ids', 'laterality'});

    % cross-campaign / cross-EDF duplicates of the same mouse: merge overlaps
    sz = sortrows(sz, {'mouse', 'start_posix'});
    keep = true(height(sz), 1);
    for i = 2:height(sz)
        p = find(keep(1:i-1) & strcmp(sz.mouse(1:i-1), sz.mouse{i}), 1, 'last');
        if ~isempty(p) && sz.start_posix(i) < sz.end_posix(p)
            S.flags(end+1, :) = {'info', 'dedup', sprintf('%s solapa %s -> fusionadas', sz.origin_ids{i}, sz.origin_ids{p})};
            sz.end_posix(p) = max(sz.end_posix(p), sz.end_posix(i));
            sz.origin_ids{p} = [sz.origin_ids{p} ' + ' sz.origin_ids{i}];
            sz.campaign{p} = strjoin(unique([split(sz.campaign{p}, ','); sz.campaign(i)])', ',');
            keep(i) = false;
        end
    end
    sz = sz(keep, :);

    % ---------------- missed seizures (ND) -------------------------------
    nd = rv.nd;
    for k = 1:height(nd)
        id = char(nd.id(k));
        m = char(nd.mouse(k));
        js = find(strcmp(cfg.nd_source(:, 1), id), 1);
        if isempty(js)
            if contains(lower(nd.tipo(k)), 'video')
                error('cs_build_seizures:Video', '%s looks like a video seizure but is not in cfg.nd_source', id);
            end
            S.excl(end+1, :) = {'crisis no detectada', id, m, '', '', '', 'sin horarios / fuera de spec', char(nd.tipo(k))};
            continue;
        end
        src = cfg.nd_source{js, 2};
        if ~isfinite(nd.start_posix(k))
            error('cs_build_seizures:NoTimes', '%s has no times', id);
        end
        % campaign whose analyzed EDF of this mouse contains the onset
        je = find(strcmp(edf_tbl.mouse, m) & edf_tbl.start_posix <= nd.start_posix(k) & edf_tbl.end_posix >= nd.start_posix(k));
        camp = strjoin(unique(edf_tbl.campaign(je))', ',');
        edf = strjoin(unique(edf_tbl.edf(je))', ',');
        if ~any(arrayfun(@(q) in_unit(edf_tbl.campaign{q}, m), je))
            S.flags(end+1, :) = {'revisar', 'ND fuera de alcance', sprintf('%s (%s) cae en una unidad no revisada/no analizada', id, m)};
        end
        if ~strcmp(nd.incluir(k), "") && ~startsWith(lower(nd.incluir(k)), "s")
            S.flags(end+1, :) = {'info', 'ND pendiente en el Excel', sprintf('%s: "%s" en la hoja; se suma segun la spec 5.1.3', id, nd.incluir(k))};
        end
        ov = find(strcmp(ev.mouse, m) & ev.start_posix < nd.end_posix(k) & ev.end_posix > nd.start_posix(k));
        if ~isempty(ov)
            S.flags(end+1, :) = {'revisar', 'ND solapa una deteccion', sprintf('%s solapa %s (AE=%s)', id, ...
                strjoin(ev.event_id(ov)', ' / '), strjoin(ev.label(ov)', '/'))};
        end
        dup = strcmp(sz.mouse, m) & sz.start_posix < nd.end_posix(k) & sz.end_posix > nd.start_posix(k);
        if any(dup)
            S.flags(end+1, :) = {'revisar', 'ND duplicada', sprintf('%s ya esta como TP; no se suma de nuevo', id)};
            continue;
        end
        if strcmp(src, 'video')
            inv = any(strcmp(rv.dx.mouse, m) & rv.dx.video_reviewed & rv.dx.start_posix <= nd.start_posix(k) & rv.dx.end_posix >= nd.end_posix(k));
            if ~inv
                error('cs_build_seizures:VideoOutsideDx', '%s (video) is not inside a DX interval with Video revisado = Si', id);
            end
        end
        sz(end+1, :) = {m, camp, edf, nd.start_posix(k), nd.end_posix(k), src, id, 'n/a'}; %#ok<AGROW>
    end
    if sum(strcmp(sz.source, 'video')) ~= 2
        error('cs_build_seizures:VideoCount', 'Expected exactly 2 video seizures (ND-05, ND-06), found %d', sum(strcmp(sz.source, 'video')));
    end

    % ---------------- checks / annotations ---------------------------------
    sz = sortrows(sz, {'mouse', 'start_posix'});
    sz.duration_s = sz.end_posix - sz.start_posix;
    sz.group = cellfun(@(m) group_of(units, m), sz.mouse, 'UniformOutput', false);
    sz.in_dx = repmat({''}, height(sz), 1);
    for i = 1:height(sz)
        d = find(strcmp(rv.dx.mouse, sz.mouse{i}) & rv.dx.start_posix < sz.end_posix(i) & rv.dx.end_posix > sz.start_posix(i));
        sz.in_dx{i} = strjoin(cellstr(rv.dx.id(d))', ',');
        if ~isempty(d) && strcmp(sz.source{i}, 'detected_TP')
            S.flags(end+1, :) = {'revisar', 'TP dentro de desconexion', sprintf('%s (%s) solapa %s', sz.origin_ids{i}, ...
                loc(sz.start_posix(i)), sz.in_dx{i})};
        end
    end
    sz.start_local = arrayfun(loc, sz.start_posix, 'UniformOutput', false);
    sz.end_local = arrayfun(loc, sz.end_posix, 'UniformOutput', false);
    sz.seizure_uid = arrayfun(@(i) sprintf('%s-%03d', sz.mouse{i}, sum(strcmp(sz.mouse(1:i), sz.mouse{i}))), (1:height(sz))', 'UniformOutput', false);
    S.seizures = sz(:, {'seizure_uid', 'mouse', 'group', 'campaign', 'edf', 'start_local', 'end_local', 'duration_s', ...
        'source', 'origin_ids', 'laterality', 'in_dx', 'start_posix', 'end_posix'});

    % pairs of the same mouse closer than cfg.merge_review_gap_s (listed, never merged)
    cp = cell(0, 4);
    for i = 2:height(sz)
        if strcmp(sz.mouse{i}, sz.mouse{i-1})
            gap = sz.start_posix(i) - sz.end_posix(i-1);
            if gap < cfg.merge_review_gap_s
                cp(end+1, :) = {sz.mouse{i}, sz.seizure_uid{i-1}, sz.seizure_uid{i}, gap}; %#ok<AGROW>
                S.flags(end+1, :) = {'decidir', 'TP separados < 30 s', sprintf('%s: %s y %s separados %.1f s (no fusionados)', ...
                    sz.mouse{i}, sz.seizure_uid{i-1}, sz.seizure_uid{i}, gap)};
            end
        end
    end
    S.close_pairs = cell2table(cp, 'VariableNames', {'mouse', 'seizure_a', 'seizure_b', 'gap_s'});

    % ---------------- not analyzed EDFs of in-scope mice -------------------
    for k = 1:height(edf_lists)
        m = edf_lists.mouse{k};
        e = edf_lists.edf{k};
        if any(strcmp(edf_tbl.mouse, m) & strcmp(edf_tbl.edf, e)) || ~any(strcmp(units.mouse, m) & units.included)
            continue;
        end
        S.excl(end+1, :) = {'EDF no analizado', e, m, edf_lists.campaign{k}, '', '', 'no analizado', ...
            ['fuera de T_eeg: ' edf_lists.reason{k}]};
    end

    S.exclusion_log = cell2table(S.excl, 'VariableNames', {'tipo', 'id', 'mouse', 'campaign', 'inicio', 'fin', 'regla', 'motivo'});
end

function g = group_of(units, m)
    j = find(strcmp(units.mouse, m) & units.included, 1);
    if isempty(j)
        g = 'fuera_de_alcance';
    else
        g = units.group{j};
    end
end
