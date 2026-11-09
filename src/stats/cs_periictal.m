function [PM, PS, no_lead] = cs_periictal(cfg, S, O, K, lead_only)
% CS_PERIICTAL  Peri-ictal IID rate (v2 section 6.3), epileptic mice.
% Windows cfg.periictal_win_min: pre windows relative to the seizure onset,
% post windows relative to its end. A window is valid when it overlaps no
% other final seizure and its analyzed IID time (mean of the channels that
% cover it) is >= cfg.periictal_min_coverage of its length. Baseline =
% analyzed time >= cfg.baseline_min_dist_h from any final seizure (as v1).
%
% lead_only: only lead seizures -- no other final seizure in the previous
% cfg.lead_min_gap_h AND those hours observed (T_seizure_obs) >=
% cfg.lead_min_coverage; when that cannot be known the seizure is not a lead.
%
% PM: per mouse x window (n seizures used, windows excluded and why, rate,
%     baseline, ratio); PS: per window across mice (median and IQR of the
%     per-mouse ratio, signed rank of the ratios vs 1); no_lead: mice left
%     out for having no lead seizure (lead_only).

    W = cfg.periictal_win_min;
    wname = {'pre -60..-30 min', 'pre -30..0 min', 'post 0..30 min', 'post 30..60 min'};
    is_post = [false false true true];
    set_name = 'todas las crisis';
    if lead_only, set_name = 'solo lideres'; end
    om = O.mouse(strcmp({O.mouse.scope}, 'merged') & strcmp({O.mouse.group}, 'epileptic'));
    rows = {};
    no_lead = {};
    ratio = NaN(numel(om), numel(wname));
    for mi = 1:numel(om)
        mouse = om(mi).mouse;
        cs = find(strcmp(K.mouse, mouse));
        sz = S.seizures(strcmp(S.seizures.mouse, mouse), :);
        if isempty(cs) || isempty(sz), continue; end
        lead = true(height(sz), 1);
        gap = cfg.lead_min_gap_h * 3600;
        for q = 1:height(sz)
            prev = [sz.start_posix(q) - gap, sz.start_posix(q)];
            others = sz([1:q-1, q+1:end], :);
            hit = any(others.start_posix < prev(2) & others.end_posix > prev(1));
            cov = cs_iv('length', cs_iv('intersect', om(mi).T_sz, prev)) / gap;
            lead(q) = ~hit && cov >= cfg.lead_min_coverage;
        end
        use = true(height(sz), 1);
        if lead_only
            use = lead;
            if ~any(use)
                no_lead{end+1} = mouse; %#ok<AGROW>
                continue;
            end
        end
        far = [sz.start_posix - cfg.baseline_min_dist_h * 3600, sz.end_posix + cfg.baseline_min_dist_h * 3600];
        [nb, eb] = count_in(K, cs, far, true);
        base = nb / (eb / 60);
        for w = 1:numel(wname)
            nv = 0; nov = 0; nlow = 0; nw = 0; ew = 0;
            for q = find(use)'
                ref = sz.start_posix(q);
                if is_post(w), ref = sz.end_posix(q); end
                win = ref + W(w, :) * 60;
                others = sz([1:q-1, q+1:end], :);
                if any(others.start_posix < win(2) & others.end_posix > win(1))
                    nov = nov + 1;
                    continue;
                end
                ncov = sum(K.start_posix(cs) < win(2) & K.start_posix(cs) + K.duration_s(cs) > win(1));
                [k, e] = count_in(K, cs, win, false);
                if ncov == 0 || e < cfg.periictal_min_coverage * diff(win) * ncov
                    nlow = nlow + 1;
                    continue;
                end
                nv = nv + 1; nw = nw + k; ew = ew + e;
            end
            rt = NaN;
            if nv > 0, rt = nw / (ew / 60); end
            ratio(mi, w) = rt / base;
            rows(end+1, :) = {set_name, mouse, wname{w}, sum(use), sum(lead), nv, nov, nlow, ew / 60, nw, rt, base, rt / base}; %#ok<AGROW>
        end
    end
    PM = cell2table(rows, 'VariableNames', {'conjunto', 'mouse', 'ventana', 'n_crisis_usadas', 'n_crisis_lideres', ...
        'n_ventanas_validas', 'n_excl_solapa_crisis', 'n_excl_cobertura', 'exposicion_canal_min', 'n_complejos', ...
        'tasa_por_min', 'tasa_basal_por_min', 'razon_vs_basal'});
    s = {};
    for w = 1:numel(wname)
        x = ratio(:, w);
        x = x(isfinite(x));
        p = NaN;
        if numel(x) >= 2, p = signrank(x, 1); end
        pw = PM(strcmp(PM.ventana, wname{w}), :);
        s(end+1, :) = {set_name, wname{w}, sum(pw.n_ventanas_validas), numel(x), median(x), prctile(x, 25), prctile(x, 75), p}; %#ok<AGROW>
    end
    PS = cell2table(s, 'VariableNames', {'conjunto', 'ventana', 'n_crisis', 'n_ratones', 'mediana_razon', 'q25_razon', ...
        'q75_razon', 'p_signed_rank_vs_1'});
end

function [k, e] = count_in(K, cs, iv, complement)
% complexes and analyzed channel-seconds of channels cs inside iv (or,
% complement = true, outside the union of iv)
    k = 0; e = 0;
    for c = cs(:)'
        if complement
            x = cs_iv('subtract', K.inc{c}, iv);
        else
            x = cs_iv('intersect', K.inc{c}, iv);
        end
        e = e + cs_iv('length', x);
        k = k + sum(cs_in_iv(K.on_abs{c}, x));
    end
end
