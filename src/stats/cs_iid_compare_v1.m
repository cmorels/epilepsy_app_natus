function T = cs_iid_compare_v1(cfg, Z, labels_changed, v1_units)
% CS_IID_COMPARE_V1  Section A.5: which channels of this run's IID plan can
% reuse an earlier run's re-detection (v1 or v2; same zone intervals, same
% code; cs_iid_previous.m; accion = 'reutilizar <run>') and which must be
% re-detected, with the reason; plus v1 channels that no
% longer run. A reason counts as explained when the channel's mouse has
% label changes vs v1 (labels_changed: {campaign|mouse}), its unit was not
% in v1's scope (v1_units: {campaign|mouse} in scope in v1), or its EDF
% changed state (rv.excluded vs v1's plan).
%
%   T: campaign, mouse, edf, region, accion, motivo, explicado

    buf = cfg.iid_excl_buffer_s;
    is_cr = strcmp(Z.zones.rule, 'crisette');
    P = Z.channels;
    q = load(fullfile(cfg.v1_dir, 'phase2_results.mat'), 'Z');
    P1 = q.Z.channels;
    key = @(X) strcat(X.campaign, '|', X.mouse, '|', X.edf, '|', X.region);
    [in1, j1] = ismember(key(P), key(P1));
    rows = cell(0, 7);
    for i = find(P.run)'
        ch = P(i, :);
        unit = [ch.campaign{1} '|' ch.mouse{1}];
        for variant = {'principal', 'sin_crisette'}
            zones = Z.zones;
            if strcmp(variant{1}, 'sin_crisette')
                if ~any(ismember(ch.zone_ids{1}, Z.zones.zone_id(is_cr))), continue; end
                zones = zones(~is_cr, :);
            end
            [sz, man] = cs_iid_channel_zones(zones, ch, buf);
            [rr, src, why] = cs_iid_previous(cfg.iid_rerun_dirs.(variant{1})(2:end), ch, sz, man, cfg.provenance);
            if ~isempty(rr)
                prev = regexp(src, '_stats_revision_\d{8}_?(\w*)[\\/]', 'tokens', 'once');   % '' = v1's folder
                if isempty(prev) || isempty(prev{1}), prev = {'v1'}; end
                rows(end+1, :) = {ch.campaign{1}, ch.mouse{1}, ch.edf{1}, ch.region{1}, variant{1}, ['reutilizar ' prev{1}], 'zonas y codigo iguales'}; %#ok<AGROW>
                continue;
            end
            if ~ismember(unit, v1_units)
                why = 'canal nuevo: unidad fuera de alcance en v1';
            elseif in1(i) && ~P1.run(j1(i))
                why = sprintf('canal nuevo: en v1 no se corria (%s)', P1.not_run_reason{j1(i)});
            elseif ismember(unit, labels_changed)
                why = [why '; etiquetas de la unidad cambiaron'];
            end
            rows(end+1, :) = {ch.campaign{1}, ch.mouse{1}, ch.edf{1}, ch.region{1}, variant{1}, 're-detectar', why}; %#ok<AGROW>
        end
    end
    % v1 channels that no longer run
    [in2, j2] = ismember(key(P1), key(P));
    for i = find(P1.run)'
        if in2(i) && P.run(j2(i)), continue; end
        r = 'canal ausente en v2';
        if in2(i), r = P.not_run_reason{j2(i)}; end
        rows(end+1, :) = {P1.campaign{i}, P1.mouse{i}, P1.edf{i}, P1.region{i}, 'principal', 'ya no se corre', r}; %#ok<AGROW>
    end
    T = cell2table(rows, 'VariableNames', {'campaign', 'mouse', 'edf', 'region', 'variante', 'accion', 'motivo'});
    T.explicado = ~strcmp(T.accion, 're-detectar') | startsWith(T.motivo, 'canal nuevo') | ...
        contains(T.motivo, 'etiquetas de la unidad cambiaron') | contains(T.motivo, 'excluido');
end
