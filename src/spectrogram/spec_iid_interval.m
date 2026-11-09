function res = spec_iid_interval(ev, r, ctx, cfg)
% SPEC_IID_INTERVAL  Time-frequency of interictal spikes in one interval
% [inicio, fin] of one channel (never the whole recording).
%
%   res = spec_iid_interval(ev, r, ctx, cfg)
%
%   ev  : spec_event (tipo 'iid'); r: spec_resolve_files struct + on, off
%   ctx : seizures (table or []), cache_dir;   cfg : cfg.spec effective
%
% Two modes, chosen by the interval length:
%   interval (<= iid.max_interval_s, default 600 s): continuous Morlet
%     scalogram of the whole interval + spike-aligned maps (below).
%   zone (longer, e.g. hours of a spiking period): no continuous
%     scalogram (it would need GBs of RAM); instead a min/max overview of
%     the 1-100 Hz trace and the spike rate per iid.rate_bin_s. The
%     spike-aligned maps are accumulated over ALL usable spikes of the zone,
%     computing the CWT in iid.zone_chunk_s blocks (identical values to a
%     CWT of the whole zone: each block carries 3 sigma_t of extra signal).
%
% 1. Loads [on + epoch(1) - margin, off + epoch(2) + margin], margin = 3
%    sigma_t of the lowest wavelet; skipped if NaN > max_nan_frac of [on, off].
% 2. Peaks from spec_get_spike_peaks (pipeline / detect_local / manual),
%    re-aligned to the extremum of the analysed signal within +-align_ms
%    with polarity 'auto' (sign of the median extremum) | 'neg' | 'pos',
%    sub-sample (parabola). A peak whose extremum falls on the edge of the
%    search window is counted (n_align_at_limit) and warned, not dropped.
% 3. An epoch is excluded when it (or its wavelet support) leaves the
%    loaded window ('edge'), touches a NaN ('nan'), overlaps an exclusion
%    zone / seizure ('zone'), or repeats a used re-aligned sample ('duplicate').
% 4. Spike-aligned maps:
%      induced = mean over spikes of |W|^2 (dB re baseline)
%      evoked  = |mean over spikes of W|^2 (= CWT of the average spike)
%      ITC     = |mean of W / |W||
%    Baseline: [-0.5 -0.3] s re each peak (induced and evoked), or an
%    absolute clean interval (induced and the continuous map; evoked then
%    in absolute dB). The continuous scalogram without an absolute baseline
%    is shown in dB re the mean power of the interval.
% A sharp transient alone gives broadband power without any oscillation:
% the figure keeps the average trace aligned above the TF maps.
%
% Features: n detected / used / excluded (by reason), rate per minute
% (detected / interval length), mean peak-to-peak (uV) and FWHM (ms) of
% the 1-100 Hz trace within +-p2p_window_ms, polarity, and the frequency
% of maximum induced power (dB re baseline) within +-maxfreq_window_ms
% (notch band excluded).

    s = cfg.spec;
    p = s.iid;
    fs = r.fs;
    a = r.on;
    b = r.off;
    zone_mode = (b - a) > p.max_interval_s;
    warns = r.notes;
    sig1 = p.morlet.cycles(1) / (2 * pi * p.morlet.fmin);
    margin = 3 * sig1 + 2 / fs;
    ep = p.epoch_s;
    seg = spec_load_segment(r, a + ep(1) - margin, b + ep(2) + margin, struct('spec', s), ctx.cache_dir);
    core = seg.t_rel >= a & seg.t_rel <= b;
    nan_frac = mean(~seg.valid(core));
    if nan_frac > s.max_nan_frac
        error('spec_iid_interval:TooManyNaN', '%s %s: %.1f%% of the interval is NaN (> %.1f%%).', ...
            ev.id, r.region, 100 * nan_frac, 100 * s.max_nan_frac);
    end
    if ~isempty(seg.outlier_note)
        warns{end+1} = seg.outlier_note;
    end
    [x_tf, x_disp] = spec_prepare_signal(seg.x, fs, struct('spec', s));
    n = numel(seg.x);
    [~, mo] = spec_morlet(zeros(0, 1), fs, p.morlet);   % frequencies / supports only
    hmax = max(mo.half_samples);

    % ---- peaks + re-alignment --------------------------------------------
    pk = spec_get_spike_peaks(ev, r, seg, a, b, ctx.seizures, struct('spec', s));
    if ~isempty(pk.note)
        warns{end+1} = pk.note;
    end
    if strcmp(p.align_on, 'display')
        y_al = x_disp;
    else
        y_al = x_tf;
    end
    idx0 = round(pk.t * fs) + 1 - seg.i0 + 1;     % position in the loaded window
    w = round(p.align_ms / 1000 * fs);
    pol = p.polarity;
    if strcmp(pol, 'auto')
        v = NaN(numel(idx0), 1);
        for k = 1:numel(idx0)
            seg_k = win(y_al, idx0(k), w);
            if ~isempty(seg_k) && ~all(isnan(seg_k))
                [~, j] = max(abs(seg_k));
                v(k) = seg_k(j);
            end
        end
        if median(v, 'omitnan') < 0
            pol = 'neg';
        else
            pol = 'pos';
        end
    end
    idx = idx0;
    at_limit = false(size(idx0));     % extremum on the edge of the +-align_ms search window
    for k = 1:numel(idx0)
        [seg_k, base] = win(y_al, idx0(k), w);
        if isempty(seg_k) || all(isnan(seg_k))
            continue;
        end
        if strcmp(pol, 'neg')
            [~, j] = min(seg_k);
        else
            [~, j] = max(seg_k);
        end
        idx(k) = base + j - 1;
        at_limit(k) = j == 1 || j == numel(seg_k);
    end
    if any(at_limit)
        warns{end+1} = sprintf(['%d of %d peak(s) re-aligned to the edge of the +-%g ms search window ' ...
            '(no %s extremum inside it): check align_ms / polarity'], sum(at_limit), numel(at_limit), p.align_ms, pol);
    end
    % sub-sample peak time: parabola through the extremum and its neighbours
    delta = zeros(size(idx));
    for k = 1:numel(idx)
        i = idx(k);
        if i > 1 && i < n && all(~isnan(y_al(i - 1:i + 1)))
            den = y_al(i - 1) - 2 * y_al(i) + y_al(i + 1);
            if den ~= 0
                delta(k) = max(-0.5, min(0.5, 0.5 * (y_al(i - 1) - y_al(i + 1)) / den));
            end
        end
    end
    t_al = seg.t_rel(1) + (idx - 1 + delta) / fs;
    shift_ms = (idx + delta - idx0) / fs * 1000;

    % ---- exclusions -------------------------------------------------------
    e0 = round(ep(1) * fs);
    e1 = round(ep(2) * fs);
    te = (e0:e1) / fs;
    nanm = isnan(x_tf);
    reason = repmat({''}, numel(idx), 1);
    seen = false(n, 1);
    for k = 1:numel(idx)
        lo = idx(k) + e0;
        hi = idx(k) + e1;
        if lo - hmax < 1 || hi + hmax > n
            reason{k} = 'edge';
        elseif any(nanm(lo - hmax:hi + hmax))     % NaN in the epoch or in its wavelet support
            reason{k} = 'nan';
        elseif any(pk.zones.end_s >= t_al(k) + ep(1) & pk.zones.start_s <= t_al(k) + ep(2))
            reason{k} = 'zone';
        elseif seen(idx(k))
            reason{k} = 'duplicate';
        else
            seen(idx(k)) = true;
        end
    end
    clear seen
    used = cellfun(@isempty, reason);
    nu = sum(used);
    if nu == 0
        warns{end+1} = 'no usable spike in this interval';
    end
    ui = find(used);

    % ---- spike-aligned CWT (whole interval, or chunks in zone mode) --------
    nf = numel(mo.f);
    ne = numel(te);
    acc = struct('P', zeros(nf, ne), 'C', complex(zeros(nf, ne)), 'U', complex(zeros(nf, ne)));
    if ~zone_mode
        W = spec_morlet(x_tf, fs, p.morlet);
        acc = acc_epochs(acc, W, idx(ui), e0, e1);
    else
        L = round(p.zone_chunk_s * fs);
        ii = idx(ui);
        c_start = 1;
        while c_start <= n
            c_end = min(n, c_start + L - 1);
            in_c = ii >= c_start & ii <= c_end;
            if any(in_c)
                r0 = max(1, c_start + e0 - hmax);
                r1 = min(n, c_end + e1 + hmax);
                Wc = spec_morlet(x_tf(r0:r1), fs, p.morlet);
                acc = acc_epochs(acc, Wc, ii(in_c) - r0 + 1, e0, e1);
            end
            c_start = c_end + 1;
        end
        W = [];
    end
    induced = acc.P / max(nu, 1);
    evoked = abs(acc.C / max(nu, 1)) .^ 2;
    itc = abs(acc.U / max(nu, 1));
    if nu == 0
        induced(:) = NaN; evoked(:) = NaN; itc(:) = NaN;
    end

    % ---- waveforms (1-100 Hz trace) + per-spike features --------------------
    traces = NaN(nu, ne);
    for k = 1:nu
        traces(k, :) = x_disp(idx(ui(k)) + (e0:e1))';
    end

    % ---- baselines + normalization ----------------------------------------
    [bi, bkind] = spec_baseline_interval(p.baseline, p.baseline_is_t_rel, p.baseline_time_ref, r);
    Pc = []; tc = []; Zc = []; cont_label = '';
    if ~zone_mode
        [Pc, tc] = block_mean_time(abs(W) .^ 2, seg.t_rel, round(p.morlet.t_step_s * fs));
    end
    if strcmp(bkind, 'relative')
        blc = te >= bi(1) & te <= bi(2);
        [Zind, label] = spec_normalize(induced, s.norm, induced(:, blc));
        Zev = spec_normalize(evoked, s.norm, evoked(:, blc));
        if ~zone_mode
            zones_mask = in_zones(tc, pk.zones);
            [Zc, ~, used_c] = spec_normalize(Pc, s.norm, Pc(:, tc >= a & tc <= b & ~zones_mask));
            cont_label = 'dB re interval mean';
            if strcmp(used_c, 'db')
                cont_label = 'power (dB re 1 uV^2)';
            end
        end
    else
        if bi(1) - margin >= seg.t_rel(1) && bi(2) + margin <= seg.t_rel(end)
            sel = seg.t_rel >= bi(1) - margin & seg.t_rel <= bi(2) + margin;
            Wb = spec_morlet(x_tf(sel), fs, p.morlet);
            tb = seg.t_rel(sel);
        else
            bseg = spec_load_segment(r, bi(1) - margin, bi(2) + margin, struct('spec', s), ctx.cache_dir);
            [bx, ~] = spec_prepare_signal(bseg.x, fs, struct('spec', s));
            Wb = spec_morlet(bx, fs, p.morlet);
            tb = bseg.t_rel;
        end
        Pbl = abs(Wb(:, tb >= bi(1) & tb <= bi(2))) .^ 2;
        [Zind, label] = spec_normalize(induced, s.norm, Pbl);
        Zev = spec_normalize(evoked, 'db');
        if ~zone_mode
            [Zc, ~] = spec_normalize(Pc, s.norm, Pbl);
            cont_label = label;
        end
        warns{end+1} = 'absolute baseline: evoked map in absolute dB';
    end
    clear W

    % ---- features ----------------------------------------------------------
    pw = round(p.p2p_window_ms / 1000 * fs);
    p2p = NaN(nu, 1);
    fwhm = NaN(nu, 1);
    c0 = -e0 + 1;    % column of the peak inside an epoch
    if strcmp(bkind, 'relative')
        bcols = find(te >= bi(1) & te <= bi(2));
    else
        bcols = find(te <= -p.p2p_window_ms / 1000);
    end
    wnd = max(1, c0 - pw):min(ne, c0 + pw);
    for k = 1:nu
        tr = traces(k, :);
        p2p(k) = max(tr(wnd)) - min(tr(wnd));
        fwhm(k) = fwhm_ms(tr, c0, median(tr(bcols), 'omitnan'), wnd, fs);
    end
    mw = abs(te) <= p.maxfreq_window_ms / 1000;
    [fex, fex_txt, notch_band] = spec_notch_band(r, seg.notch_applied, 0);
    prof = mean(Zind(:, mw), 2, 'omitnan');
    prof(fex(mo.f)) = NaN;
    [~, im] = max(prof);
    f_ind = NaN;
    if nu > 0 && ~isempty(im) && ~all(isnan(prof))
        f_ind = mo.f(im);
    end

    % ---- assemble --------------------------------------------------------------
    res = struct();
    res.kind = 'iid';
    res.mode = ternary(zone_mode, 'zone', 'interval');
    res.id = ev.id;
    res.subject = r.subject;
    res.session_key = r.session_key;
    res.source_file = r.source_file;
    res.region = r.region;
    res.region_tag = r.region;
    res.origin = ev.origin;
    res.on = a;
    res.off = b;
    [res.on_natus, res.on_pipeline] = clock_text(a, r);
    [res.off_natus, res.off_pipeline] = clock_text(b, r);
    if zone_mode
        % overview: min/max per bin (decimate_minmax.m) instead of hours of samples
        [res.t_trace, res.trace] = decimate_minmax(seg.t_rel(core), x_disp(core), p.overview_points);
        edges = a:p.rate_bin_s:b;
        if edges(end) < b
            edges(end+1) = b;
        end
        res.rate = struct('t', edges(1:end - 1) + diff(edges) / 2, 'bin_s', diff(edges), ...
            'detected', histcounts(pk.t, edges) ./ (diff(edges) / 60), ...
            'used', histcounts(t_al(used), edges) ./ (diff(edges) / 60));
    else
        res.t_trace = seg.t_rel(core);
        res.trace = x_disp(core);
        res.rate = [];
    end
    res.signal = seg.x;
    res.signal_t_rel = seg.t_rel;
    res.peaks = struct('t_source', pk.t, 't_aligned', t_al, 'shift_ms', shift_ms, 'used', used, 'at_limit', at_limit, ...
        'reason', {reason}, 'source', pk.source, 'detect', pk.detect, 'zones', pk.zones, ...
        'polarity_mode', p.polarity, 'polarity', pol);
    if zone_mode
        res.cont = [];
    else
        keep_c = tc >= a & tc <= b;
        res.cont = struct('t', tc(keep_c), 'f', mo.f, 'P', Pc(:, keep_c), 'Z', Zc(:, keep_c), 'label', cont_label);
    end
    keep_tr = round(linspace(1, nu, min(nu, p.max_traces_saved)));
    res.epochs = struct('t', te, 'traces', traces(keep_tr, :), 'trace_index', ui(keep_tr), ...
        'mean', mean(traces, 1, 'omitnan'), 'sem', std(traces, 0, 1, 'omitnan') / sqrt(max(nu, 1)), 'n', nu);
    res.f = mo.f;
    res.coi_s = mo.coi_s;
    res.induced = struct('P', induced, 'Z', Zind);
    res.evoked = struct('P', evoked, 'Z', Zev);
    res.itc = itc;
    res.norm_label = label;
    res.baseline = bi;
    res.baseline_kind = bkind;
    res.nan_frac = nan_frac;
    res.cfg = s;
    res.src = struct('txt', r.txt, 'clean', r.clean, 'raw', r.raw, 'iid_mat', r.iid_mat, ...
        'signal_source', seg.source, 'gain_applied', seg.gain_applied, 'notch_applied', seg.notch_applied, ...
        'case_applied', seg.case_applied, 'n_outliers_window', seg.n_outliers_segment, ...
        'excel_attenuation', r.excel_attenuation, 'header', r.hdr.meta, ...
        'cache_build_s', seg.load_info.build_s, 'cache_reused', seg.load_info.reused, 'notch_band', notch_band);
    res.clim_tf = pick_clim({Zind, Zev}, s);
    if zone_mode
        res.clim_cont = [];
    else
        res.clim_cont = pick_clim({res.cont.Z}, s);
    end
    res.warnings = warns;

    count = @(rs) sum(strcmp(reason, rs));
    ss = r.session_start;
    ss.Format = 'dd-MMM-yyyy HH:mm:ss.SSS';
    thr = NaN;
    if ~isempty(pk.detect)
        thr = pk.detect.threshold_uV;
    end
    res.features = table({ev.id}, {res.mode}, {r.subject}, {r.region}, {r.session_key}, {r.source_file}, {char(ss)}, ...
        a, b, b - a, {res.on_natus}, {res.on_pipeline}, {res.off_natus}, {res.off_pipeline}, ...
        {pk.source}, thr, numel(pk.t), nu, numel(pk.t) - nu, count('nan'), count('zone'), count('edge'), count('duplicate'), ...
        numel(pk.t) / ((b - a) / 60), mean(p2p, 'omitnan'), mean(fwhm, 'omitnan'), {pol}, {p.polarity}, ...
        median(abs(shift_ms), 'omitnan'), sum(at_limit), f_ind, {fex_txt}, {s.signal_source}, {label}, {seg.case_applied}, seg.gain_applied, ...
        seg.notch_applied, seg.n_outliers_segment, nan_frac, {r.excel_attenuation}, {strjoin(warns, ' | ')}, ...
        'VariableNames', {'id', 'mode', 'subject_id', 'region', 'session', 'source_file', 'session_start', ...
        'start_t_rel', 'end_t_rel', 'duration_s', 'start_natus', 'start_pipeline', 'end_natus', 'end_pipeline', ...
        'spike_source', 'threshold_uV', 'n_detected', 'n_used', 'n_excluded', 'n_excl_nan', 'n_excl_zone', ...
        'n_excl_edge', 'n_excl_duplicate', 'rate_per_min', 'p2p_mean_uV', 'fwhm_mean_ms', 'polarity', ...
        'polarity_mode', 'median_abs_shift_ms', 'n_align_at_limit', 'induced_peak_freq_hz', 'freq_excluded', 'signal_source', 'norm', 'case_applied', ...
        'gain_applied', 'notch_applied', 'n_outliers_window', 'nan_frac', 'excel_attenuation', 'warnings'});
end

%% ======================================================================
function acc = acc_epochs(acc, W, pos, e0, e1)
% Add the epochs of W centred on columns pos to the running sums.
    for k = 1:numel(pos)
        C = W(:, pos(k) + (e0:e1));
        A = abs(C);
        acc.P = acc.P + A .^ 2;
        acc.C = acc.C + C;
        acc.U = acc.U + C ./ A;
    end
end

function [v, base] = win(y, i, w)
    base = max(1, i - w);
    last = min(numel(y), i + w);
    if i < 1 || i > numel(y) || last < base
        v = [];
        return;
    end
    v = y(base:last);
end

function m = in_zones(t, Z)
    m = false(size(t));
    for k = 1:height(Z)
        m = m | (t >= Z.start_s(k) & t <= Z.end_s(k));
    end
end

function v = fwhm_ms(tr, c0, base, wnd, fs)
% Full width at half the peak deviation from base, linear interpolation at
% both crossings; NaN if a crossing is not found within wnd.
    v = NaN;
    A = tr(c0) - base;
    if isnan(A) || A == 0
        return;
    end
    half = base + A / 2;
    g = sign(A) * (tr - half);        % > 0 above half (towards the peak)
    L = c0;
    while L > wnd(1) && g(L - 1) > 0
        L = L - 1;
    end
    R = c0;
    while R < wnd(end) && g(R + 1) > 0
        R = R + 1;
    end
    if L == wnd(1) || R == wnd(end) || isnan(g(L - 1)) || isnan(g(R + 1))
        return;
    end
    xl = (L - 1) + g(L - 1) / (g(L - 1) - g(L));
    xr = R + g(R) / (g(R) - g(R + 1));
    v = (xr - xl) / fs * 1000;
end

function [P2, t2] = block_mean_time(P, t_rel, k)
    k = max(1, k);
    n = floor(size(P, 2) / k) * k;
    P2 = reshape(mean(reshape(P(:, 1:n), size(P, 1), k, []), 2), size(P, 1), []);
    t2 = mean(reshape(t_rel(1:n), k, []), 1);
end

function [n, p] = clock_text(t, r)
    [a, b] = spec_resolve_times('clocks', t, r);
    n = spec_resolve_times('text', a);
    p = spec_resolve_times('text', b);
end

function v = ternary(c, a, b)
    if c
        v = a;
    else
        v = b;
    end
end

function cl = pick_clim(Zs, s)
    if isnumeric(s.clim) && numel(s.clim) == 2
        cl = s.clim(:)';
        return;
    end
    v = cellfun(@(Z) Z(isfinite(Z)), Zs, 'UniformOutput', false);
    v = vertcat(v{:});
    if isempty(v)
        cl = [0 1];
        return;
    end
    cl = prctile(v, s.clim_prctile);
    if cl(2) <= cl(1)
        cl(2) = cl(1) + 1;
    end
end
