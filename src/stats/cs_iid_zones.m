function Z = cs_iid_zones(cfg, S, rv, ch_tbl, det_rows)
% CS_IID_ZONES  Exclusion zones of the IID re-detection (section 5.3), per
% mouse in absolute time (POSIX s). A zone applies to BOTH hemispheres of
% its mouse and to every channel of the mouse whose span it overlaps
% (crossed by mouse + absolute time, never by EDF name). Every zone gets
% +/- cfg.iid_excl_buffer_s; start/end here are WITHOUT the buffer.
%
%   rule      zone                                             pass_as
%   R1        final seizure (detected TP, ND-01..06)           seizures
%   R7        final seizure whose TP rows carry text in AF     seizures
%   crisette  crisette_focal (O-07, O-08, O-13)                seizures
%   R3 / R4   FP grooming/scratching / other                   manual
%   R5        FP disconnect; consecutive ones (no other event  manual
%             in between, <= cfg.disconnect_run_max_gap_s apart) form one
%             span from the first start to the last end
%   R6        FP with empty or unrecognized AF; every event of  manual
%             a unit without review
%   R8        interval of 2_Desconexiones                      manual
%   (R2)      FP IIDs: NOT a zone (Z.fp_iids, for the section-8 check)
%
% pass_as = 'seizures': detect_iid's seizures argument, which adds
% cfg.iid.exclusion_buffer_s itself; 'manual': cfg.iid.exclusion_zones_manual,
% buffer added by cs_iid_channel_zones.m (detect_iid adds none there).
% Discarded EDFs (DX 'EDF completo') are not run and contribute no zones.
%
% Z.zones    one row per zone
% Z.fp_iids  FP-IIDs windows of reviewed units (counted, not excluded)
% Z.channels channel plan: run or not, zones that apply, expected T_iid
% Z.rule_summary  zones and excluded hours per rule over the EDFs that run,
%            vs the campaign's own exclusion (det_rows +/- buffer)
% Z.flags    {severity, topic, detail}

    loc = @(p) char(datetime(p, 'ConvertFrom', 'posixtime', 'TimeZone', cfg.tz, 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
    rows = cell(0, 10);   % mouse, rule, origin, campaign, edf, start, end, pass_as, unit_in_scope, detail
    flags = cell(0, 3);
    ev = S.events;
    units = S.units;
    in_unit = @(c, m) any(strcmp(units.campaign, c) & strcmp(units.mouse, m) & units.included);

    % ---------------- R1 / R7: final seizures --------------------------------
    r7_ids = ev.event_id(ev.usable & ismember(ev.label, {'TP', 'CONFLICT'}) & ~cellfun(@isempty, ev.comments));
    for i = 1:height(S.seizures)
        s = S.seizures(i, :);
        ids = strtrim(split(s.origin_ids{1}, ' + '));
        rule = 'R1';
        detail = sprintf('crisis %s (%s): %s', s.seizure_uid{1}, s.source{1}, s.origin_ids{1});
        hit = ismember(ev.event_id, ids) & ismember(ev.event_id, r7_ids);
        if any(hit)
            rule = 'R7';
            detail = [detail '; AF=' strjoin(ev.comments(hit)', ',')];
        end
        rows(end+1, :) = {s.mouse{1}, rule, s.seizure_uid{1}, s.campaign{1}, s.edf{1}, s.start_posix, s.end_posix, ...
            'seizures', true, detail}; %#ok<AGROW>
    end

    % ---------------- crisette_focal --------------------------------------
    for i = 1:height(S.crisette)
        c = S.crisette(i, :);
        rows(end+1, :) = {c.mouse{1}, 'crisette', c.obs_id{1}, c.campaign{1}, c.edf{1}, c.start_posix, ...
            c.start_posix + c.duration_s, 'seizures', in_unit(c.campaign{1}, c.mouse{1}), ...
            sprintf('crisette_focal %s (#%d, AE=%s)', c.obs_id{1}, c.seizure_id, c.AE_label{1})}; %#ok<AGROW>
    end

    % ---------------- detector events: FP (reviewed) and unreviewed units ---
    % Events of discarded EDFs, crisette events (above) and fragments merged
    % into a seizure (covered by its R1 span) are skipped. TP of reviewed
    % units are the R1 zones; TP of units out of scope (partial review)
    % become R1 zones from their own span.
    cand = ev(~ev.discarded & cellfun(@isempty, ev.crisette) & ev.merged_into == 0, :);
    n_mixed = 0;
    fp_iids = cand([], :);
    for i = 1:height(cand)
        e = cand(i, :);
        lab = e.label{1};
        if e.in_scope && ~strcmp(lab, 'FP')
            continue;   % TP/CONFLICT of a reviewed unit: already an R1 zone
        end
        if strcmp(lab, 'TP') || strcmp(lab, 'CONFLICT')
            rule = 'R1';
            detail = sprintf('TP de unidad fuera de alcance (AE=%s)', lab);
        elseif isempty(lab)
            rule = 'R6';
            detail = 'sin etiqueta AE (unidad sin revision)';
        else
            [rule, detail, mixed] = fp_rule(e.comments{1});
            if mixed && e.in_scope
                n_mixed = n_mixed + 1;
                flags(end+1, :) = {'decidir', 'FP con AF mixto (IIDs + artefacto)', sprintf(['%s (%s): AF="%s" ' ...
                    '(un hemisferio IIDs, el otro artefacto); excluido como %s en ambos hemisferios'], ...
                    e.event_id{1}, loc(e.start_posix), e.comments{1}, rule)}; %#ok<AGROW>
            end
            if strcmp(rule, 'R2')
                if e.in_scope
                    fp_iids(end+1, :) = e; %#ok<AGROW>
                end
                continue;
            end
            if strcmp(rule, 'R5')
                continue;   % disconnect: spans built per EDF below
            end
        end
        rows(end+1, :) = {e.mouse{1}, rule, e.event_id{1}, e.campaign{1}, e.edf{1}, e.start_posix, e.end_posix, ...
            'manual', e.in_scope, detail}; %#ok<AGROW>
    end

    % ---------------- R5: disconnect runs per (campaign, mouse, EDF) --------
    base = ev(~ev.discarded, :);
    [g, key] = findgroups(base(:, {'campaign', 'mouse', 'edf'}));
    for k = 1:height(key)
        evf = base(g == k, :);
        runs = cs_disconnect_runs(evf, cfg.disconnect_run_max_gap_s);
        isd = strcmp(evf.label, 'FP') & contains(evf.comments, 'disconnect');
        for r = 1:size(runs, 1)
            members = evf.seizure_id(isd & evf.start_posix >= runs(r, 1) & evf.end_posix <= runs(r, 2));
            rows(end+1, :) = {key.mouse{k}, 'R5', sprintf('%s|%s|%s|%s', key.campaign{k}, key.mouse{k}, key.edf{k}, ...
                strjoin(arrayfun(@num2str, sort(members)', 'UniformOutput', false), '+')), key.campaign{k}, key.edf{k}, ...
                runs(r, 1), runs(r, 2), 'manual', in_unit(key.campaign{k}, key.mouse{k}), ...
                sprintf('FP disconnect: %d evento(s) unidos', numel(members))}; %#ok<AGROW>
        end
    end

    % ---------------- R8: 2_Desconexiones ------------------------------------
    for k = find(~rv.dx.whole_edf)'
        rows(end+1, :) = {char(rv.dx.mouse(k)), 'R8', char(rv.dx.id(k)), '', char(rv.dx.edf(k)), rv.dx.start_posix(k), ...
            rv.dx.end_posix(k), 'manual', true, sprintf('desconexion %s (limites: %s / %s)', rv.dx.id(k), ...
            rv.dx.limit_start(k), rv.dx.limit_end(k))}; %#ok<AGROW>
    end

    zones = cell2table(rows, 'VariableNames', {'mouse', 'rule', 'origin', 'campaign', 'edf', 'start_posix', 'end_posix', ...
        'pass_as', 'unit_in_scope', 'detail'});
    zones = sortrows(zones, {'mouse', 'start_posix', 'rule'});
    zones.zone_id = (1:height(zones))';
    zones.duration_s = zones.end_posix - zones.start_posix;
    zones.start_local = arrayfun(loc, zones.start_posix, 'UniformOutput', false);
    zones.end_local = arrayfun(loc, zones.end_posix, 'UniformOutput', false);
    zones = zones(:, {'zone_id', 'mouse', 'rule', 'origin', 'campaign', 'edf', 'start_local', 'end_local', 'duration_s', ...
        'pass_as', 'unit_in_scope', 'detail', 'start_posix', 'end_posix'});
    if any(zones.duration_s < 0)
        error('cs_iid_zones:NegativeZone', 'zone(s) with end < start: %s', strjoin(zones.origin(zones.duration_s < 0)', ', '));
    end

    % ---------------- channel plan -------------------------------------------
    disc = rv.excluded;   % cs_excluded_edfs.m
    n = height(ch_tbl);
    run = false(n, 1); reason = repmat({''}, n, 1); unit_ok = false(n, 1);
    zone_ids = cell(n, 1); n_sz = zeros(n, 1); n_man = zeros(n, 1);
    excl_h = zeros(n, 1); expected_T_iid_h = zeros(n, 1); foreign = zeros(n, 1);
    for i = 1:n
        ch = ch_tbl(i, :);
        unit_ok(i) = in_unit(ch.campaign{1}, ch.mouse{1});
        jd = find(strcmp(disc.mouse, ch.mouse{1}) & strcmp(disc.edf, ch.edf{1}), 1);
        if ~isempty(jd)
            reason{i} = sprintf('EDF excluido (%s, %s)', disc.id{jd}, disc.kind{jd});
        elseif ~unit_ok(i) && strcmp(cfg.iid_rerun_scope, 'in_scope')
            reason{i} = 'unidad fuera de alcance (sin revision)';
        else
            run(i) = true;
        end
        [sz, man, ids] = cs_iid_channel_zones(zones, ch, cfg.iid_excl_buffer_s);
        zone_ids{i} = ids;
        n_sz(i) = height(sz);
        n_man(i) = size(man, 1);
        dur = ch.duration_s;
        iv = [sz.start_s - cfg.iid_excl_buffer_s, sz.end_s + cfg.iid_excl_buffer_s; man];
        iv = cs_iv('intersect', iv, [0 dur]);
        g = ch.gaps{1};
        excl_h(i) = cs_iv('length', iv) / 3600;
        expected_T_iid_h(i) = (dur - cs_iv('length', cs_iv('union', [iv; g]))) / 3600;
        foreign(i) = sum(~zones.unit_in_scope(ismember(zones.zone_id, ids))) * unit_ok(i);
    end
    plan = ch_tbl(:, {'campaign', 'mouse', 'edf', 'region', 'file', 'start_posix', 'duration_s', 'n_gaps'});
    plan.unit_in_scope = unit_ok;
    plan.run = run;
    plan.not_run_reason = reason;
    plan.n_zones_seizures = n_sz;
    plan.n_zones_manual = n_man;
    plan.excluded_h = excl_h;
    plan.expected_T_iid_h = expected_T_iid_h;
    plan.zone_ids = zone_ids;
    for i = find(plan.run & plan.expected_T_iid_h < 1)'
        flags(end+1, :) = {'decidir', 'canal IID con < 1 h analizable', sprintf(['%s %s %s %s: %.2f h de %.1f h quedan tras las zonas; ' ...
            'baseline y tasas sobre muy poco registro'], plan.campaign{i}, plan.mouse{i}, plan.edf{i}, plan.region{i}, ...
            plan.expected_T_iid_h(i), plan.duration_s(i) / 3600)}; %#ok<AGROW>
    end
    if any(foreign > 0)
        flags(end+1, :) = {'revisar', 'zona de unidad no revisada sobre canal en alcance', sprintf('%d canal(es)', sum(foreign > 0))};
    end

    % ---------------- excluded time per rule ---------------------------------
    % Over the channels that run, one channel per (campaign, mouse, EDF):
    % both hemispheres get the same zones. Rules overlap (e.g. R5 inside R8);
    % TOTAL is their union. 'antes' = the campaign's own exclusion (every
    % reported detector row of that channel +/- buffer, as 04_iid did).
    buf = cfg.iid_excl_buffer_s;
    [~, first] = unique(strcat(plan.campaign, '|', plan.mouse, '|', plan.edf), 'stable');
    sel = first(plan.run(first));
    total_h = sum(plan.duration_s(sel)) / 3600;
    summ = cell(0, 4);
    for r = {'R1', 'R7', 'crisette', 'R3', 'R4', 'R5', 'R6', 'R8', 'TOTAL', 'antes'}
        nz = 0;
        h = 0;
        for i = sel'
            ch = ch_tbl(i, :);
            if strcmp(r{1}, 'antes')
                d = det_rows(strcmp(det_rows.campaign, ch.campaign{1}) & strcmp(det_rows.mouse, ch.mouse{1}) & ...
                    strcmp(det_rows.edf, ch.edf{1}) & strcmp(det_rows.region, ch.region{1}), :);
                iv = [d.start_s - buf, d.end_s + buf];
            else
                zz = zones(ismember(zones.zone_id, zone_ids{i}), :);
                if ~strcmp(r{1}, 'TOTAL')
                    zz = zz(strcmp(zz.rule, r{1}), :);
                end
                iv = [zz.start_posix - buf, zz.end_posix + buf] - ch.start_posix;
            end
            iv = reshape(iv, [], 2);
            nz = nz + sum(iv(:, 2) > 0 & iv(:, 1) < ch.duration_s);
            iv = cs_iv('intersect', iv, [0 ch.duration_s]);
            h = h + cs_iv('length', iv) / 3600;
        end
        summ(end+1, :) = {r{1}, nz, h, 100 * h / total_h}; %#ok<AGROW>
    end
    Z.rule_summary = cell2table(summ, 'VariableNames', {'regla', 'n_zonas', 'horas_excluidas', 'pct_del_registro'});
    Z.n_edf_run = numel(sel);
    Z.hours_run = total_h;

    Z.zones = zones;
    Z.fp_iids = fp_iids(:, {'event_id', 'mouse', 'campaign', 'edf', 'seizure_id', 'start_posix', 'end_posix', 'duration_s', 'comments'});
    Z.channels = plan;
    Z.n_mixed = n_mixed;
    Z.flags = flags;
end

function [rule, detail, mixed] = fp_rule(c)
% AF value set (';'-joined, normalized by cs_load_detections) -> rule.
% An event whose hemispheres disagree (e.g. 'IIDs;other') is excluded by
% its artifact label: the zone applies to both hemispheres.
    parts = split(string(c), ';');
    parts = parts(strlength(parts) > 0);
    mixed = any(parts == "IIDs") && any(parts ~= "IIDs");
    if isempty(parts) || any(startsWith(parts, "unrecognized:"))
        rule = 'R6';
        detail = sprintf('FP sin categoria reconocida (AF="%s")', c);
    elseif all(parts == "IIDs")
        rule = 'R2';
        detail = 'FP IIDs (no se excluye)';
    elseif any(parts == "disconnect")
        rule = 'R5';
        detail = 'FP disconnect';
    elseif any(parts == "grooming/scratching")
        rule = 'R3';
        detail = sprintf('FP grooming/scratching (AF="%s")', c);
    elseif any(parts == "other")
        rule = 'R4';
        detail = sprintf('FP other (AF="%s")', c);
    else
        rule = 'R6';
        detail = sprintf('FP sin categoria reconocida (AF="%s")', c);
    end
end
