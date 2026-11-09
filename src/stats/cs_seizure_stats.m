function R = cs_seizure_stats(cfg, S, O)
% CS_SEIZURE_STATS  Section 6.1: per mouse x scope (C1, C2, merged) seizure
% counts by source, rate per 24 h over T_seizure_obs, seizure-free days,
% durations (EEG vs video), inter-seizure intervals, clusters, laterality;
% the daily table; circadian analysis normalized by observed time.
%
% Conventions:
%   - a seizure belongs to scope Cx when its campaign field contains Cx
%     (missed/video seizures take the campaign of the EDF that contains
%     them); merged = all.
%   - rate numerator = seizures whose onset lies in T_seizure_obs (TPs inside
%     a disconnection are counted in n_total, reported in n_outside_obs and
%     flagged by cs_build_seizures, but not in the rate).
%   - valid day = local calendar day with >= cfg.min_valid_h_day of
%     T_seizure_obs; seizure-free days are counted among valid days only.
%   - ISI = onset-to-onset of consecutive seizures, kept only when >= 95 %
%     of the interval is observed (T_seizure_obs); CV = std/mean.
%   - cluster: consecutive onsets <= gap h apart (6, 12, 24 h).

    tz = cfg.tz;
    sz = S.seizures;
    sz.minute = arrayfun(@(p) minute_of_day(p, tz), sz.start_posix);
    sz.date = arrayfun(@(p) local_date(p, tz), sz.start_posix);
    rows = {};
    daily = {};
    for k = 1:numel(O.mouse)
        o = O.mouse(k);
        sel = strcmp(sz.mouse, o.mouse);
        if ~strcmp(o.scope, 'merged')
            sel = sel & contains(sz.campaign, o.scope);
        end
        x = sortrows(sz(sel, :), 'start_posix');
        in_obs = cs_iv('overlap', [x.start_posix x.start_posix + 1e-6], o.T_sz) > 0;
        xo = x(in_obs, :);
        n_det = sum(strcmp(x.source, 'detected_TP'));
        n_mis = sum(strcmp(x.source, 'missed_eeg'));
        n_vid = sum(strcmp(x.source, 'video'));
        cr = S.crisette(strcmp(S.crisette.mouse, o.mouse), :);
        if ~strcmp(o.scope, 'merged')
            cr = cr(strcmp(cr.campaign, o.scope), :);
        end
        Teeg = cs_iv('length', o.T_eeg) / 3600;
        Tsz = cs_iv('length', o.T_sz) / 3600;
        % days
        [~, dsz] = cs_clock_split(o.T_sz, tz, [0 1440]);
        valid = dsz.total_s / 3600 >= cfg.min_valid_h_day;
        nday = arrayfun(@(d) sum(xo.date == d), dsz.date);
        if strcmp(o.scope, 'merged')
            for i = 1:height(dsz)
                xi = xo(xo.date == dsz.date(i), :);
                daily(end+1, :) = {o.mouse, o.group, dsz.date(i), dsz.total_s(i) / 3600, valid(i), height(xi), ...
                    sum(strcmp(xi.source, 'detected_TP')), sum(strcmp(xi.source, 'missed_eeg')), sum(strcmp(xi.source, 'video'))}; %#ok<AGROW>
            end
        end
        % durations
        de = x.duration_s(~strcmp(x.source, 'video'));
        dv = x.duration_s(strcmp(x.source, 'video'));
        % ISI
        isi = [];
        for i = 2:height(x)
            a = x.start_posix(i-1); b = x.start_posix(i);
            if cs_iv('length', cs_iv('intersect', o.T_sz, [a b])) >= 0.95 * (b - a)
                isi(end+1) = (b - a) / 3600; %#ok<AGROW>
            end
        end
        % clusters
        nc = zeros(1, 3); maxc = NaN;
        gaps_h = diff(x.start_posix) / 3600;
        for gi = 1:3
            g = cfg.cluster_gaps_h(gi);
            if height(x) > 0
                lab = cumsum([1; gaps_h > g]);
                nc(gi) = max(lab);
                if g == cfg.cluster_gap_h
                    maxc = max(accumarray(lab, 1));
                end
            end
        end
        lat = x.laterality(strcmp(x.source, 'detected_TP'));
        nl = sum(ismember(lat, {'bilateral', 'HPCr_only', 'HPCl_only'}));
        rows(end+1, :) = {o.mouse, o.group, o.scope, n_det, n_mis, n_vid, height(x), height(x) - height(xo), height(cr), ...
            Teeg, Tsz, height(xo) / Tsz * 24, sum(valid), sum(valid & nday == 0), sum(valid & nday == 0) / max(1, sum(valid)), ...
            q(de, 50), q(de, 25), q(de, 75), mn(de), mx(de), q(dv, 50), mn(dv), mx(dv), ...
            numel(isi), q(isi, 50), q(isi, 25), q(isi, 75), std(isi) / mean(isi), nc(1), nc(2), nc(3), maxc, ...
            nl, 100 * sum(strcmp(lat, 'bilateral')) / nl, 100 * sum(strcmp(lat, 'HPCr_only')) / nl, ...
            100 * sum(strcmp(lat, 'HPCl_only')) / nl}; %#ok<AGROW>
    end
    R.mouse = cell2table(rows, 'VariableNames', {'mouse', 'group', 'scope', 'n_detected_TP', 'n_missed_eeg', 'n_video', ...
        'n_total', 'n_outside_obs', 'n_crisette', 'T_eeg_h', 'T_seizure_obs_h', 'rate_per_24h', 'n_valid_days', ...
        'n_seizure_free_days', 'frac_seizure_free_days', 'dur_eeg_median_s', 'dur_eeg_q25_s', 'dur_eeg_q75_s', ...
        'dur_eeg_min_s', 'dur_eeg_max_s', 'dur_video_median_s', 'dur_video_min_s', 'dur_video_max_s', 'n_isi', ...
        'isi_median_h', 'isi_q25_h', 'isi_q75_h', 'isi_cv', 'n_clusters_6h', 'n_clusters_12h', 'n_clusters_24h', ...
        'max_cluster_size_24h', 'n_lateralized', 'pct_bilateral', 'pct_HPCr_only', 'pct_HPCl_only'});
    R.daily = cell2table(daily, 'VariableNames', {'mouse', 'group', 'date', 'T_seizure_obs_h', 'valid_day', ...
        'n_seizures', 'n_detected_TP', 'n_missed_eeg', 'n_video'});

    % group summaries (pooled counts/time + median of per-mouse rates)
    g = {};
    for grp = {'epileptic', 'control'}
        for sc = {'C1', 'C2', 'merged'}
            t = R.mouse(strcmp(R.mouse.group, grp{1}) & strcmp(R.mouse.scope, sc{1}), :);
            if isempty(t)
                continue;
            end
            g(end+1, :) = {['GRUPO:' grp{1}], grp{1}, sc{1}, sum(t.n_detected_TP), sum(t.n_missed_eeg), sum(t.n_video), ...
                sum(t.n_total), sum(t.n_outside_obs), sum(t.n_crisette), sum(t.T_eeg_h), sum(t.T_seizure_obs_h), ...
                (sum(t.n_total) - sum(t.n_outside_obs)) / sum(t.T_seizure_obs_h) * 24, sum(t.n_valid_days), ...
                sum(t.n_seizure_free_days), sum(t.n_seizure_free_days) / max(1, sum(t.n_valid_days)), ...
                median(t.rate_per_24h, 'omitnan'), height(t)}; %#ok<AGROW>
        end
    end
    R.group = cell2table(g, 'VariableNames', {'mouse', 'group', 'scope', 'n_detected_TP', 'n_missed_eeg', 'n_video', ...
        'n_total', 'n_outside_obs', 'n_crisette', 'T_eeg_h', 'T_seizure_obs_h', 'pooled_rate_per_24h', 'n_valid_days', ...
        'n_seizure_free_days', 'frac_seizure_free_days', 'median_mouse_rate_per_24h', 'n_mice'});

    % ---------------- circadian (merged scope) ------------------------------
    on = hm(cfg.lights_on); off = hm(cfg.lights_off); w = cfg.transition_window_h * 60;
    hourly = {};
    tests = {};
    units = {'ALL_epileptic'};
    om = O.mouse(strcmp({O.mouse.scope}, 'merged'));
    ep = om(strcmp({om.group}, 'epileptic'));
    units = [units, {ep.mouse}];
    for u = 1:numel(units)
        if u == 1
            mm = {ep.mouse};
            H = sum(vertcat(ep.hour_sz_s), 1);
            P = sum(vertcat(ep.phase_sz_s), 1);
        else
            mm = units(u);
            H = ep(u-1).hour_sz_s;
            P = ep(u-1).phase_sz_s;
        end
        x = sz(ismember(sz.mouse, mm), :);
        x = x(arrayfun(@(i) cs_iv('overlap', [x.start_posix(i) x.start_posix(i) + 1e-6], ...
            om(strcmp({om.mouse}, x.mouse{i})).T_sz) > 0, (1:height(x))'), :);
        if height(x) == 0
            continue;
        end
        hcount = accumarray(floor(x.minute / 60) + 1, 1, [24 1])';
        for h = 1:24
            hourly(end+1, :) = {units{u}, h - 1, hcount(h), H(h) / 3600, hcount(h) / (H(h) / 3600), 24 * hcount(h) / (H(h) / 3600)}; %#ok<AGROW>
        end
        pb = discretize(x.minute, cfg.phase_edges_min);  % bins of cfg.phase_edges_min
        pc = accumarray(pb, 1, [numel(cfg.phase_edges_min) - 1 1])';
        E = P / 3600;
        dark = ismember(1:numel(pc), cfg.phase_dark_bins);
        r = cs_rate_ratio(sum(pc(dark)), sum(E(dark)), sum(pc(~dark)), sum(E(~dark)));
        tests(end+1, :) = test_row(units{u}, 'oscuridad / luz', r); %#ok<AGROW>
        bl = cfg.phase_bin_after_on; br = cfg.phase_bins_light_rest;
        r = cs_rate_ratio(pc(bl), E(bl), sum(pc(br)), sum(E(br)));
        tests(end+1, :) = test_row(units{u}, sprintf('%s-%s vs resto de luz', cfg.lights_on, clock_str(on + w)), r); %#ok<AGROW>
        bd = cfg.phase_bin_after_off; bdr = cfg.phase_bins_dark_rest;
        r = cs_rate_ratio(pc(bd), E(bd), sum(pc(bdr)), sum(E(bdr)));
        tests(end+1, :) = test_row(units{u}, sprintf('%s-%s vs resto de oscuridad', cfg.lights_off, clock_str(off + w)), r); %#ok<AGROW>
        ry = cs_rayleigh(x.minute, H, cfg.n_boot, cfg.rng_seed);
        tests(end+1, :) = {units{u}, 'Rayleigh (reloj 24 h)', ry.n, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ry.p, ...
            ry.rbar, ry.mean_clock_h, ry.Z, ry.p_exposure}; %#ok<AGROW>
    end
    R.circ_hourly = cell2table(hourly, 'VariableNames', {'unit', 'hour', 'n_seizures', 'observed_h', 'rate_per_h', 'rate_per_24h'});
    R.circ_tests = cell2table(tests, 'VariableNames', {'unit', 'test', 'n1', 'exposure1_h', 'n0', 'exposure0_h', ...
        'rate1_per_24h', 'rate0_per_24h', 'RR', 'RR_ci_lo', 'RR_ci_hi', 'p', 'rbar', 'mean_clock_h', 'Z', 'p_exposure_corrected'});
    R.seizures = sz;
end

function row = test_row(unit, name, r)
    row = {unit, name, r.n1, r.e1, r.n0, r.e0, 24 * r.rate1, 24 * r.rate0, r.rr, r.ci_lo, r.ci_hi, r.p, NaN, NaN, NaN, NaN};
end

function v = q(x, p)
    if isempty(x), v = NaN; else, v = prctile(x, p); end
end
function v = mn(x)
    if isempty(x), v = NaN; else, v = min(x); end
end
function v = mx(x)
    if isempty(x), v = NaN; else, v = max(x); end
end
function m = minute_of_day(p, tz)
    t = datetime(p, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
    m = hour(t) * 60 + minute(t) + second(t) / 60;
end
function d = local_date(p, tz)
    t = datetime(p, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
    d = datetime(year(t), month(t), day(t));
end
function m = hm(s)
    v = sscanf(s, '%d:%d');
    m = v(1) * 60 + v(2);
end
function s = clock_str(m)
    s = sprintf('%02d:%02d', floor(m / 60), mod(m, 60));
end
