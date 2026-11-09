function [BC, SM, SG, FL] = cs_iid_sens_stats(cfg, T, V, F, S, R, C, Z)
% CS_IID_SENS_STATS  Tables of phase C (v2 section 5), from cs_iid_sens_all.m.
%
%   BC  iid_threshold_by_channel: percentiles of |x| / |x_bp|, where 312.5 uV
%       falls among |x_bp|, the threshold of every variant and its complexes/min
%   SM  iid_threshold_sensitivity: per variant x mouse (both hemispheres,
%       C1+C2): complexes/min, singles/min, polyspikes/min, bursts/h
%   SG  iid_threshold_summary: per variant x metric, epileptic vs control:
%       medians, fold change with bootstrap CI, Mann-Whitney, Cliff's delta,
%       >= 10x, and Spearman between IID rate and seizure rate (epileptic)
%   FL  flatline_scan: every flat segment with both-hemisphere overlap and
%       how much of it an IID exclusion already covers
% The prominence of findpeaks (200 uV) is the same in every variant.

    grp_of = @(m) S.units.group{find(strcmp(S.units.mouse, m), 1)};
    tz = cfg.tz;
    loc = @(p) char(datetime(p, 'ConvertFrom', 'posixtime', 'TimeZone', tz, 'Format', 'yyyy-MM-dd HH:mm:ss'));
    vn = unique(V.variant, 'stable');

    % ---- by channel -------------------------------------------------------
    BC = T;
    BC.group = cellfun(grp_of, BC.mouse, 'UniformOutput', false);
    BC = movevars(BC, 'group', 'After', 'mouse');
    for k = 1:numel(vn)
        v = V(strcmp(V.variant, vn{k}), :);
        [~, j] = ismember(strcat(BC.campaign, '|', BC.mouse, '|', BC.edf, '|', BC.region), ...
            strcat(v.campaign, '|', v.mouse, '|', v.edf, '|', v.region));
        nm = matlab.lang.makeValidName(['umbral_' vn{k}]);
        BC.(nm) = v.thr_uV(j);
        BC.(matlab.lang.makeValidName(['cx_min_' vn{k}])) = v.n_complexes(j) ./ v.analyzed_min(j);
    end

    % ---- per variant x mouse ------------------------------------------------
    rows = {};
    mice = unique(V.mouse);
    for k = 1:numel(vn)
        v = V(strcmp(V.variant, vn{k}), :);
        for i = 1:numel(mice)
            x = v(strcmp(v.mouse, mice{i}), :);
            mn = sum(x.analyzed_min);
            rows(end+1, :) = {vn{k}, mice{i}, grp_of(mice{i}), median(x.thr_uV), sum(x.n_complexes), mn / 60, ...
                sum(x.n_complexes) / mn, sum(x.n_single) / mn, sum(x.n_polyspike) / mn, sum(x.n_bursts) / (mn / 60)}; %#ok<AGROW>
        end
    end
    SM = cell2table(rows, 'VariableNames', {'variante', 'mouse', 'group', 'umbral_mediano_uV', 'n_complejos', 'T_iid_h', ...
        'complexes_per_min', 'single_per_min', 'polyspikes_per_min', 'bursts_per_h'});

    % ---- per variant x group ------------------------------------------------
    sr = R.mouse(strcmp(R.mouse.scope, 'merged'), {'mouse', 'rate_per_24h'});
    metrics = {'complexes_per_min', 'single_per_min', 'polyspikes_per_min', 'bursts_per_h'};
    rows = {};
    for k = 1:numel(vn)
        x = SM(strcmp(SM.variante, vn{k}), :);
        for q = 1:numel(metrics)
            e = x.(metrics{q})(strcmp(x.group, 'epileptic'));
            c = x.(metrics{q})(strcmp(x.group, 'control'));
            g = cs_group_compare(e, c, cfg);
            xe = x(strcmp(x.group, 'epileptic'), :);
            [~, j] = ismember(xe.mouse, sr.mouse);
            [rho, prho] = corr(xe.(metrics{q}), sr.rate_per_24h(j), 'Type', 'Spearman');
            rows(end+1, :) = {vn{k}, metrics{q}, numel(e), numel(c), median(e), median(c), g.fold, g.ci(1), g.ci(2), g.p, ...
                g.cliff, g.fold >= 10, rho, prho}; %#ok<AGROW>
        end
    end
    SG = cell2table(rows, 'VariableNames', {'variante', 'metrica', 'n_epilepticos', 'n_controles', 'mediana_epi', 'mediana_ctrl', ...
        'fold_change_medianas', 'fc_ci95_lo', 'fc_ci95_hi', 'p_mann_whitney', 'cliffs_delta', 'cumple_10x', ...
        'spearman_iid_vs_crisis', 'p_spearman'});

    % ---- flat segments --------------------------------------------------------
    n = height(F);
    both = zeros(n, 1); cov = zeros(n, 1);
    Cm = C(strcmp(C.variant, 'principal'), :);
    for i = 1:n
        o = F(strcmp(F.campaign, F.campaign{i}) & strcmp(F.mouse, F.mouse{i}) & strcmp(F.edf, F.edf{i}) & ...
            ~strcmp(F.region, F.region{i}), :);
        seg = [F.start_posix(i) F.end_posix(i)];
        both(i) = cs_iv('length', cs_iv('intersect', [o.start_posix o.end_posix], seg)) / diff(seg);
        j = find(strcmp(Cm.campaign, F.campaign{i}) & strcmp(Cm.mouse, F.mouse{i}) & strcmp(Cm.edf, F.edf{i}) & ...
            strcmp(Cm.region, F.region{i}), 1);
        ex = Cm.excl{j} + Cm.start_posix(j);
        cov(i) = cs_iv('length', cs_iv('intersect', ex, seg)) / diff(seg);
    end
    FL = F(:, {'mouse', 'campaign', 'edf', 'region'});
    FL.inicio = arrayfun(loc, F.start_posix, 'UniformOutput', false);
    FL.fin = arrayfun(loc, F.end_posix, 'UniformOutput', false);
    FL.duracion_s = F.end_posix - F.start_posix;
    FL.frac_simultaneo_otro_hemisferio = both;
    FL.ambos_hemisferios = both > 0;
    FL.frac_cubierto_por_exclusion = cov;
    FL.cubierto = repmat({'no'}, n, 1);
    FL.cubierto(cov >= 0.99) = {'si'};
    FL.cubierto(cov > 0 & cov < 0.99) = {'parcial'};
    FL = sortrows(FL, 'duracion_s', 'descend');
end
