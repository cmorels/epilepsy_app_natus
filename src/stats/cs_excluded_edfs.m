function [X, flags] = cs_excluded_edfs(cfg, rv, edf_tbl, ev)
% CS_EXCLUDED_EDFS  EDFs left out of EVERYTHING (observation time, seizure
% set, detector performance, IID re-detection), one row per (mouse, EDF):
%   - 2_Desconexiones rows whose limits are 'EDF completo' (DX-20)
%   - cfg.edf_excluded: {ED id, mouse, EDF, kind}, each checked against
%     sheet 3_EDF_estado (ED-16 descartado, ED-17 pendiente de reprocesamiento)
%   - units in cfg.edf_level_units: every EDF holding an event without AE
%     (the unit stays in scope with its other EDFs); crisette_focal events
%     (cfg.crisette) do not need AE
% start/end = the EDF span (*_clean.txt header).

    flags = cell(0, 3);
    rows = cell(0, 6);   % id, mouse, campaign, edf, kind, reason
    d = rv.dx(rv.dx.whole_edf, :);
    for k = 1:height(d)
        rows(end+1, :) = {char(d.id(k)), char(d.mouse(k)), '', char(d.edf(k)), 'descartado', 'DX con limites "EDF completo"'}; %#ok<AGROW>
    end
    for k = 1:size(cfg.edf_excluded, 1)
        c = cfg.edf_excluded(k, :);
        j = find(rv.ed.id == c{1});
        if numel(j) ~= 1 || rv.ed.mouse(j) ~= c{2} || rv.ed.edf(j) ~= c{3}
            error('cs_excluded_edfs:Ed', '%s (%s %s) not found as such in 3_EDF_estado', c{1}, c{2}, c{3});
        end
        rows(end+1, :) = {c{1}, c{2}, '', c{3}, c{4}, char(rv.ed.estado(j))}; %#ok<AGROW>
    end
    % crisette_focal events are classified by 5_Observaciones whatever AE says
    is_cr = false(height(ev), 1);
    for k = 1:size(cfg.crisette, 1)
        c = cfg.crisette(k, :);
        is_cr = is_cr | (strcmp(ev.campaign, c{3}) & strcmp(ev.mouse, c{2}) & strcmp(ev.edf, c{4}) & ev.seizure_id == c{5});
    end
    for u = 1:size(cfg.edf_level_units, 1)
        camp = cfg.edf_level_units{u, 1};
        mouse = cfg.edf_level_units{u, 2};
        E = edf_tbl(strcmp(edf_tbl.campaign, camp) & strcmp(edf_tbl.mouse, mouse), :);
        for i = 1:height(E)
            s = strcmp(ev.campaign, camp) & strcmp(ev.mouse, mouse) & strcmp(ev.edf, E.edf{i}) & ~is_cr;
            nu = sum(s & cellfun(@isempty, ev.label));
            if nu > 0
                rows(end+1, :) = {'AE incompleto', mouse, camp, E.edf{i}, 'sin AE completo', ...
                    sprintf('%d de %d eventos sin AE en el snapshot', nu, sum(s))}; %#ok<AGROW>
                flags(end+1, :) = {'revisar', 'EDF excluido por AE incompleto', sprintf(['%s %s %s: %d de %d eventos sin AE ' ...
                    '(snapshot); solo este EDF queda fuera'], camp, mouse, E.edf{i}, nu, sum(s))}; %#ok<AGROW>
            end
        end
    end
    % one row per (mouse, EDF), ids joined
    key = strcat(rows(:, 2), '|', rows(:, 4));
    [u, first] = unique(key, 'stable');
    X = cell2table(rows(first, :), 'VariableNames', {'id', 'mouse', 'campaign', 'edf', 'kind', 'reason'});
    for i = 1:numel(u)
        X.id{i} = strjoin(unique(rows(strcmp(key, u{i}), 1), 'stable')', '/');
    end
    X.start_posix = NaN(height(X), 1);
    X.end_posix = NaN(height(X), 1);
    for i = 1:height(X)
        j = find(strcmp(edf_tbl.mouse, X.mouse{i}) & strcmp(edf_tbl.edf, X.edf{i}));
        if isempty(j)
            error('cs_excluded_edfs:NoEdf', '%s: EDF %s of %s is not among the analyzed EDFs', X.id{i}, X.edf{i}, X.mouse{i});
        end
        X.campaign{i} = strjoin(unique(edf_tbl.campaign(j))', ',');
        X.start_posix(i) = min(edf_tbl.start_posix(j));
        X.end_posix(i) = max(edf_tbl.end_posix(j));
    end
end
