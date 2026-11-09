function I = cs_iid_stats(cfg, S, R, Z, C, Corig)
% CS_IID_STATS  Section 6.4 from the IID re-detection (C, main variant)
% and the campaigns' original 04_iid (Corig), over the channels of the
% plan (reviewed units, discarded EDF excluded).
%
%   I.channels     per channel x EDF: new vs original summary and rates
%   I.before_after per channel: change in complexes decomposed into
%                  (a) original complexes inside disconnections (R5/R8 zones,
%                  +/- buffer; removed), (b) new complexes inside FP-IIDs
%                  windows (+/- buffer; added), (c) the rest
%                  (= delta + a - b: threshold / other mask changes);
%                  sorted by |change in complexes/min|
%   I.mouse        per mouse x scope (C1, C2, merged) x hemisphere (HPCr,
%                  HPCl, both): pooled rates (counts / analyzed time),
%                  asymmetry (R-L)/(R+L) of complexes/min, new and original
%   I.group        control vs epileptic per-mouse rates (merged, both
%                  hemispheres): medians, fold change of medians with
%                  stratified bootstrap CI, Mann-Whitney, Cliff's delta,
%                  >= 10x criterion; new and original data, with and
%                  without 048-s (O-29)
%   I.circ_hourly / I.circ_tests  IIDs per clock hour normalized by the
%                  analyzed channel-time in each hour; dark/light: pooled
%                  Poisson rate ratio and per-mouse ratios (signed rank)
%   I.periictal    rate in cfg.periictal_win_min windows (pre: rel. to
%                  onset; post: rel. to end) vs interictal baseline
%                  (>= cfg.baseline_min_dist_h from any seizure), per mouse
%                  and paired across mice (signed rank)
%   I.day          per mouse x local day: analyzed IID time and rate, joined
%                  with the seizure day table (section 6.5)
% Rates use analyzed channel-time: both hemispheres pooled = sum of counts
% over sum of the two channels' analyzed minutes.

    tz = cfg.tz;
    buf = cfg.iid_excl_buffer_s;
    grp_of = @(m) S.units.group{find(strcmp(S.units.mouse, m), 1)};
    Cm = C(strcmp(C.variant, 'principal'), :);
    key = @(T) strcat(T.campaign, '|', T.mouse, '|', T.edf, '|', T.region);
    [tf, j] = ismember(key(Cm), key(Corig));
    if ~all(tf)
        error('cs_iid_stats:Join', 'channels without original results');
    end
    Co = Corig(j, :);
    n = height(Cm);
    Cm.group = cellfun(grp_of, Cm.mouse, 'UniformOutput', false);

    % absolute included intervals and complex onsets per channel
    inc = cell(n, 1);
    on_abs = cell(n, 1);
    for i = 1:n
        inc{i} = cs_iv('subtract', [0 Cm.duration_s(i)], Cm.excl{i}) + Cm.start_posix(i);
        on_abs{i} = Cm.complexes{i}(:, 1) + Cm.start_posix(i);
    end

    %% ---- per channel ------------------------------------------------------
    rate = @(k, mins) k ./ mins;
    T = table(Cm.campaign, Cm.mouse, Cm.group, Cm.edf, Cm.region, Cm.duration_s / 3600, Cm.analyzed_min / 60, Co.analyzed_min / 60, ...
        Cm.baseline_uV, Co.baseline_uV, Cm.lower_uV, Co.lower_uV, Cm.p97_uV, Cm.frac_below_130uV, ...
        Cm.n_complexes, Co.n_complexes, Cm.n_single, Co.n_single, Cm.n_polyspike, Co.n_polyspike, Cm.n_bursts, Co.n_bursts, ...
        'VariableNames', {'campaign', 'mouse', 'group', 'edf', 'region', 'duration_h', 'T_iid_h', 'T_iid_h_orig', ...
        'baseline_uV', 'baseline_uV_orig', 'lower_uV', 'lower_uV_orig', 'p97_real_uV', 'frac_below_130uV', ...
        'n_complexes', 'n_complexes_orig', 'n_single', 'n_single_orig', 'n_polyspike', 'n_polyspike_orig', 'n_bursts', 'n_bursts_orig'});
    T.complexes_per_min = rate(T.n_complexes, Cm.analyzed_min);
    T.complexes_per_min_orig = rate(T.n_complexes_orig, Co.analyzed_min);
    T.single_per_min = rate(T.n_single, Cm.analyzed_min);
    T.polyspikes_per_min = rate(T.n_polyspike, Cm.analyzed_min);
    T.pct_polyspike = 100 * T.n_polyspike ./ T.n_complexes;
    T.bursts_per_h = rate(T.n_bursts, Cm.analyzed_min / 60);
    T.bursts_per_h_orig = rate(T.n_bursts_orig, Co.analyzed_min / 60);
    T.orig_source = Co.source;
    I.channels = T;

    %% ---- before / after decomposition ------------------------------------
    disc = Z.zones(ismember(Z.zones.rule, {'R5', 'R8'}), :);
    fpi = Z.fp_iids;
    a = NaN(n, 1); b = zeros(n, 1); disc_h = zeros(n, 1); fpi_n = zeros(n, 1);
    for i = 1:n
        t0 = Cm.start_posix(i);
        d = disc(strcmp(disc.mouse, Cm.mouse{i}), :);
        div = [d.start_posix - buf, d.end_posix + buf] - t0;
        div = cs_iv('intersect', reshape(div, [], 2), [0 Cm.duration_s(i)]);
        disc_h(i) = cs_iv('length', div) / 3600;
        if Co.complexes_known(i)
            a(i) = sum(in_iv(Co.complexes{i}(:, 1), div));
        end
        w = fpi(strcmp(fpi.mouse, Cm.mouse{i}), :);
        wiv = [w.start_posix - buf, w.end_posix + buf] - t0;
        wiv = cs_iv('intersect', reshape(wiv, [], 2), [0 Cm.duration_s(i)]);
        fpi_n(i) = size(wiv, 1);
        b(i) = sum(in_iv(Cm.complexes{i}(:, 1), wiv));
    end
    BA = T(:, {'campaign', 'mouse', 'group', 'edf', 'region', 'T_iid_h_orig', 'T_iid_h', 'baseline_uV_orig', 'baseline_uV', ...
        'lower_uV_orig', 'lower_uV', 'n_complexes_orig', 'n_complexes', 'complexes_per_min_orig', 'complexes_per_min'});
    BA.delta_complexes = T.n_complexes - T.n_complexes_orig;
    BA.delta_per_min = T.complexes_per_min - T.complexes_per_min_orig;
    BA.delta_pct = 100 * BA.delta_per_min ./ T.complexes_per_min_orig;
    BA.disconnected_h = disc_h;
    BA.a_orig_in_disconnections = a;
    BA.n_fp_iids_windows = fpi_n;
    BA.b_new_in_fp_iids = b;
    BA.c_rest = BA.delta_complexes + a - b;
    BA.orig_source = Co.source;
    [~, k] = sort(abs(BA.delta_per_min), 'descend', 'MissingPlacement', 'last');
    I.before_after = BA(k, :);

    %% ---- per mouse x scope x hemisphere ---------------------------------
    rows = {};
    mice = unique(Cm.mouse);
    for mi = 1:numel(mice)
        m = mice{mi};
        for sc = {'C1', 'C2', 'merged'}
            for hm = {'HPCr', 'HPCl', 'ambos'}
                s = strcmp(T.mouse, m);
                if ~strcmp(sc{1}, 'merged'), s = s & strcmp(T.campaign, sc{1}); end
                if ~strcmp(hm{1}, 'ambos'), s = s & strcmp(T.region, hm{1}); end
                if ~any(s), continue; end
                x = T(s, :);
                mn = sum(x.T_iid_h) * 60;
                mo = sum(x.T_iid_h_orig) * 60;
                rows(end+1, :) = {m, grp_of(m), sc{1}, hm{1}, numel(unique(x.edf)), sum(x.T_iid_h), ...
                    sum(x.n_complexes), sum(x.n_single), sum(x.n_polyspike), sum(x.n_bursts), ...
                    sum(x.n_complexes) / mn, sum(x.n_single) / mn, sum(x.n_polyspike) / mn, ...
                    100 * sum(x.n_polyspike) / sum(x.n_complexes), sum(x.n_bursts) / (mn / 60), NaN, ...
                    sum(x.T_iid_h_orig), sum(x.n_complexes_orig), sum(x.n_complexes_orig) / mo, sum(x.n_single_orig) / mo, ...
                    sum(x.n_polyspike_orig) / mo, sum(x.n_bursts_orig) / (mo / 60), NaN}; %#ok<AGROW>
            end
            % asymmetry on the 'ambos' row
            r = find(strcmp(rows(:, 1), m) & strcmp(rows(:, 3), sc{1}));
            if numel(r) == 3
                rr = rows{r(1), 11}; ll = rows{r(2), 11};
                rows{r(3), 16} = (rr - ll) / (rr + ll);
                rr = rows{r(1), 19}; ll = rows{r(2), 19};
                rows{r(3), 23} = (rr - ll) / (rr + ll);
            end
        end
    end
    M = cell2table(rows, 'VariableNames', {'mouse', 'group', 'scope', 'hemisphere', 'n_edf', 'T_iid_h', 'n_complexes', ...
        'n_single', 'n_polyspike', 'n_bursts', 'complexes_per_min', 'single_per_min', 'polyspikes_per_min', 'pct_polyspike', ...
        'bursts_per_h', 'asym_index_RL', 'T_iid_h_orig', 'n_complexes_orig', 'complexes_per_min_orig', 'single_per_min_orig', ...
        'polyspikes_per_min_orig', 'bursts_per_h_orig', 'asym_index_RL_orig'});
    I.mouse = M;

    %% ---- control vs epileptic ------------------------------------------
    mm = M(strcmp(M.scope, 'merged') & strcmp(M.hemisphere, 'ambos'), :);
    metrics = {'complexes_per_min', 'single_per_min', 'polyspikes_per_min', 'bursts_per_h'};
    rows = {};
    for data = {'re-deteccion', 'original'}
        for set = {'todos', 'sin 048-s'}
            x = mm;
            if strcmp(set{1}, 'sin 048-s'), x = x(~strcmp(x.mouse, '048-s'), :); end
            for k = 1:numel(metrics)
                v = metrics{k};
                if strcmp(data{1}, 'original'), v = [v '_orig']; end %#ok<AGROW>
                e = x.(v)(strcmp(x.group, 'epileptic'));
                c = x.(v)(strcmp(x.group, 'control'));
                g = group_compare(e, c, cfg);
                rows(end+1, :) = {metrics{k}, data{1}, set{1}, numel(e), numel(c), median(e), prctile(e, 25), prctile(e, 75), ...
                    median(c), prctile(c, 25), prctile(c, 75), g.fold, g.ci(1), g.ci(2), g.p, g.cliff, g.fold >= 10, g.ci(1) >= 10}; %#ok<AGROW>
            end
        end
    end
    I.group = cell2table(rows, 'VariableNames', {'metrica', 'datos', 'conjunto', 'n_epilepticos', 'n_controles', ...
        'mediana_epi', 'q25_epi', 'q75_epi', 'mediana_ctrl', 'q25_ctrl', 'q75_ctrl', 'fold_change_medianas', 'fc_ci95_lo', ...
        'fc_ci95_hi', 'p_mann_whitney', 'cliffs_delta', 'criterio_10x_fold', 'criterio_10x_ci_lo'});

    %% ---- circadian -------------------------------------------------------
    on = sscanf(cfg.lights_on, '%d:%d'); on = on(1) * 60 + on(2);
    off = sscanf(cfg.lights_off, '%d:%d'); off = off(1) * 60 + off(2);
    hr_s = zeros(n, 24); hr_n = zeros(n, 24); ph_s = zeros(n, 2); ph_n = zeros(n, 2);   % phases: dark, light
    for i = 1:n
        hr_s(i, :) = cs_clock_split(inc{i}, tz, 0:60:1440);
        p = cs_clock_split(inc{i}, tz, [0 on off 1440]);
        ph_s(i, :) = [p(1) + p(3), p(2)];
        if ~isempty(on_abs{i})
            dt = datetime(on_abs{i}, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
            mod_ = hour(dt) * 60 + minute(dt) + second(dt) / 60;
            hr_n(i, :) = accumarray(hour(dt) + 1, 1, [24 1])';
            dark = mod_ < on | mod_ >= off;
            ph_n(i, :) = [sum(dark), sum(~dark)];
        end
    end
    rows = {};
    tests = {};
    for gname = {'epileptic', 'control'}
        s = strcmp(Cm.group, gname{1});
        for h = 1:24
            e = sum(hr_s(s, h)) / 3600;
            rows(end+1, :) = {gname{1}, h - 1, sum(hr_n(s, h)), e, sum(hr_n(s, h)) / e, sum(hr_n(s, h)) / e / 60}; %#ok<AGROW>
        end
        rr = cs_rate_ratio(sum(ph_n(s, 1)), sum(ph_s(s, 1)) / 3600, sum(ph_n(s, 2)), sum(ph_s(s, 2)) / 3600);
        % per mouse dark/light ratio
        ms = unique(Cm.mouse(s));
        ratio = NaN(numel(ms), 1);
        for k = 1:numel(ms)
            sk = strcmp(Cm.mouse, ms{k});
            ratio(k) = (sum(ph_n(sk, 1)) / sum(ph_s(sk, 1))) / (sum(ph_n(sk, 2)) / sum(ph_s(sk, 2)));
        end
        ok = isfinite(ratio) & ratio > 0;
        p_sr = NaN;
        if sum(ok) >= 2, p_sr = signrank(log(ratio(ok))); end
        tests(end+1, :) = {gname{1}, 'oscuridad/luz, agrupado (Poisson; ignora sobredispersion)', rr.n1, rr.e1, rr.n0, rr.e0, ...
            rr.rr, rr.ci_lo, rr.ci_hi, rr.p, NaN, NaN, NaN, NaN}; %#ok<AGROW>
        tests(end+1, :) = {gname{1}, 'oscuridad/luz por raton (signed rank de log razon)', NaN, NaN, NaN, NaN, ...
            NaN, NaN, NaN, p_sr, sum(ok), median(ratio(ok)), prctile(ratio(ok), 25), prctile(ratio(ok), 75)}; %#ok<AGROW>
    end
    I.circ_hourly = cell2table(rows, 'VariableNames', {'group', 'hour', 'n_complexes', 'channel_h', 'rate_per_h', 'rate_per_min'});
    I.circ_tests = cell2table(tests, 'VariableNames', {'group', 'test', 'n_dark', 'channel_h_dark', 'n_light', 'channel_h_light', ...
        'rate_ratio', 'ci_lo', 'ci_hi', 'p', 'n_mice', 'median_ratio_mouse', 'q25_ratio_mouse', 'q75_ratio_mouse'});

    %% ---- peri-ictal ------------------------------------------------------
    W = cfg.periictal_win_min;
    wname = {'pre -60..-30 min', 'pre -30..0 min', 'post 0..30 min', 'post 30..60 min'};
    is_post = [false false true true];
    rows = {};
    epi = unique(Cm.mouse(strcmp(Cm.group, 'epileptic')));
    per_mouse = NaN(numel(epi), numel(wname) + 1);
    for mi = 1:numel(epi)
        m = epi{mi};
        cs = find(strcmp(Cm.mouse, m));
        sz = S.seizures(strcmp(S.seizures.mouse, m), :);
        % interictal baseline
        far = [sz.start_posix - cfg.baseline_min_dist_h * 3600, sz.end_posix + cfg.baseline_min_dist_h * 3600];
        [nb, eb] = count_in(cs, inc, on_abs, far, true);
        base = nb / (eb / 60);
        per_mouse(mi, end) = base;
        for w = 1:numel(wname)
            nv = 0; nlow = 0; nov = 0; nw = 0; ew = 0;
            for q = 1:height(sz)
                ref = sz.start_posix(q);
                if is_post(w), ref = sz.end_posix(q); end
                win = ref + W(w, :) * 60;
                others = sz([1:q-1, q+1:end], :);
                if any(others.start_posix < win(2) & others.end_posix > win(1))
                    nov = nov + 1;
                    continue;
                end
                ncov = sum(arrayfun(@(c) Cm.start_posix(c) < win(2) && Cm.start_posix(c) + Cm.duration_s(c) > win(1), cs));
                [k, e] = count_in(cs, inc, on_abs, win, false);
                if ncov == 0 || e < 0.5 * diff(win) * ncov
                    nlow = nlow + 1;
                    continue;
                end
                nv = nv + 1; nw = nw + k; ew = ew + e;
            end
            r = nw / (ew / 60);
            if nv == 0, r = NaN; end
            per_mouse(mi, w) = r;
            rows(end+1, :) = {m, wname{w}, height(sz), nv, nov, nlow, ew / 60, nw, r, base, r / base}; %#ok<AGROW>
        end
    end
    PI = cell2table(rows, 'VariableNames', {'mouse', 'ventana', 'n_crisis', 'n_ventanas_validas', 'n_excluidas_solapan_crisis', ...
        'n_excluidas_poca_exposicion', 'exposicion_canal_min', 'n_complejos', 'tasa_por_min', 'tasa_basal_por_min', 'razon_vs_basal'});
    % paired across mice: window rate vs the same mouse's baseline
    PI.n_ratones = NaN(height(PI), 1);
    PI.q25_razon = NaN(height(PI), 1);
    PI.q75_razon = NaN(height(PI), 1);
    PI.p_signed_rank = NaN(height(PI), 1);
    for w = 1:numel(wname)
        ok = isfinite(per_mouse(:, w)) & isfinite(per_mouse(:, end));
        p = NaN;
        if sum(ok) >= 2, p = signrank(per_mouse(ok, w), per_mouse(ok, end)); end
        ratio = per_mouse(ok, w) ./ per_mouse(ok, end);
        r = PI(1, :);
        r.mouse = {'(epilepticos, mediana por raton)'};
        r.ventana = wname(w);
        r{1, {'n_crisis', 'n_ventanas_validas', 'n_excluidas_solapan_crisis', 'n_excluidas_poca_exposicion', ...
            'exposicion_canal_min', 'n_complejos'}} = [sum(PI.n_crisis(strcmp(PI.ventana, wname{w}))), ...
            sum(PI.n_ventanas_validas(strcmp(PI.ventana, wname{w}))), sum(PI.n_excluidas_solapan_crisis(strcmp(PI.ventana, wname{w}))), ...
            sum(PI.n_excluidas_poca_exposicion(strcmp(PI.ventana, wname{w}))), NaN, NaN];
        r{1, {'tasa_por_min', 'tasa_basal_por_min', 'razon_vs_basal', 'n_ratones', 'q25_razon', 'q75_razon', 'p_signed_rank'}} = ...
            [median(per_mouse(ok, w)), median(per_mouse(ok, end)), median(ratio), sum(ok), prctile(ratio, 25), prctile(ratio, 75), p];
        PI = [PI; r]; %#ok<AGROW>
    end
    I.periictal = PI;

    %% ---- mouse x day (section 6.5) -----------------------------------------
    rows = {};
    for mi = 1:numel(mice)
        m = mice{mi};
        cs = find(strcmp(Cm.mouse, m));
        dk = zeros(0, 1); dh = zeros(0, 1); dn = zeros(0, 1); dc = zeros(0, 1);
        for c = cs'
            [~, dt] = cs_clock_split(inc{c}, tz, [0 1440]);
            if isempty(dt), continue; end
            [y, mo, d] = ymd(dt.date);
            kk = y * 10000 + mo * 100 + d;
            cnt = zeros(numel(kk), 1);
            if ~isempty(on_abs{c})
                t = datetime(on_abs{c}, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
                [y2, mo2, d2] = ymd(t);
                [tfk, loc] = ismember(y2 * 10000 + mo2 * 100 + d2, kk);
                cnt = accumarray(loc(tfk), 1, [numel(kk) 1]);
            end
            dk = [dk; kk]; dh = [dh; dt.total_s / 3600]; dn = [dn; ones(numel(kk), 1)]; dc = [dc; cnt]; %#ok<AGROW>
        end
        [u, ~, g] = unique(dk);
        for k = 1:numel(u)
            ch_h = sum(dh(g == k));
            rows(end+1, :) = {m, grp_of(m), u(k), ch_h / sum(dn(g == k)), sum(dn(g == k)), sum(dc(g == k)), ...
                sum(dc(g == k)) / ch_h}; %#ok<AGROW>
        end
    end
    D = cell2table(rows, 'VariableNames', {'mouse', 'group', 'daykey', 'T_iid_h', 'n_channels', 'n_complexes', 'iid_per_h'});
    Rd = R.daily;
    [y, mo, d] = ymd(Rd.date);
    Rd.daykey = y * 10000 + mo * 100 + d;
    D = outerjoin(D, Rd(:, {'mouse', 'daykey', 'T_seizure_obs_h', 'valid_day', 'n_seizures'}), 'Keys', {'mouse', 'daykey'}, ...
        'MergeKeys', true, 'Type', 'left');
    D.date = datetime(floor(D.daykey / 10000), mod(floor(D.daykey / 100), 100), mod(D.daykey, 100), 'Format', 'yyyy-MM-dd');
    D.valid_iid_day = D.T_iid_h >= cfg.min_valid_h_day;
    D.valid_both = D.valid_iid_day & D.valid_day == 1;
    I.day = D;
end

function tf = in_iv(t, iv)
% t inside any [start end] row of iv (sorted, non-overlapping)
    tf = false(size(t));
    for k = 1:size(iv, 1)
        tf = tf | (t >= iv(k, 1) & t <= iv(k, 2));
    end
end

function [k, e] = count_in(cs, inc, on_abs, iv, complement)
% Complexes and analyzed channel-seconds of channels cs inside iv (or,
% with complement = true, outside the union of iv).
    k = 0; e = 0;
    for c = cs(:)'
        if complement
            x = cs_iv('subtract', inc{c}, iv);
        else
            x = cs_iv('intersect', inc{c}, iv);
        end
        e = e + cs_iv('length', x);
        k = k + sum(in_iv(on_abs{c}, x));
    end
end

function g = group_compare(e, c, cfg)
% fold change of medians (epileptic / control) with a stratified bootstrap
% percentile CI, Mann-Whitney p and Cliff's delta (P(e>c) - P(e<c)).
    g.fold = median(e) / median(c);
    rng(cfg.rng_seed);
    fb = zeros(cfg.n_boot, 1);
    ne = numel(e); nc = numel(c);
    for b = 1:cfg.n_boot
        fb(b) = median(e(randi(ne, ne, 1))) / median(c(randi(nc, nc, 1)));
    end
    fb = fb(~isnan(fb));
    g.ci = prctile(fb, [2.5 97.5]);
    g.p = ranksum(e, c);
    d = sign(e(:) - c(:)');
    g.cliff = mean(d(:));
end
