function res = spec_seizure_event(ev, R, ctx, cfg)
% SPEC_SEIZURE_EVENT  Spectral analysis of one seizure (one or both
% hippocampal channels), on the window [onset - pad_pre, offset + pad_post]
% only -- never the whole recording.
%
%   res = spec_seizure_event(ev, R, ctx, cfg)
%
%   ev  : spec_event struct
%   R   : 1 x n_regions struct (spec_resolve_files output + on, off = this
%         region's onset / offset in t_rel)
%   ctx : seizures (table from spec_read_seizures_table or []), cache_dir
%   cfg : cfg.spec effective for this event
%
% Per region:
%   - window loaded with an extra margin = max(win_s / 2, 3 sigma_t of the
%     lowest-frequency wavelet), dropped before plotting / features;
%     skipped (error) when NaN > max_nan_frac of [on - pad_pre, off + pad_post]
%   - signal for TF: mean / trend removed per block, no band-pass;
%     trace: 1-100 Hz zero-phase
%   - multitaper (spec_multitaper) and Morlet (spec_morlet) power, both
%     normalized (spec_normalize) with the baseline: [-30 -5] s re onset by
%     default, or an absolute interval (loaded separately if outside)
%   - baseline warnings: NaN in it, another seizure of seizures_events in
%     it, an IID exclusion zone (04_iid) in it
%   - PSD (mean multitaper spectrum of the windows fully inside):
%     baseline, ictal [on, off], post-ictal [off, off + postictal_s], dB
%   - features: peak frequency of the ictal PSD (absolute and re baseline),
%     band power ictal / baseline (dB) for delta, theta, beta, gamma, and
%     the dominant-frequency curve (argmax over f of the multitaper map in
%     dB re baseline, each window inside [on, off]; its mean and linear slope)
% t = 0 at the event onset (earliest onset over its regions).

    s = cfg.spec;
    z = s.seizure;
    on = min([R.on]);
    off = max([R.off]);
    t_lo = on - z.pad_pre_s;
    t_hi = off + z.pad_post_s;
    sig1 = z.morlet.cycles(1) / (2 * pi * z.morlet.fmin);
    margin = max(z.mt.win_s / 2, 3 * sig1) + 2 / R(1).fs;

    res = struct('kind', 'seizure', 'id', ev.id, 'subject', R(1).subject, 'session_key', R(1).session_key, ...
        'source_file', R(1).source_file, 'region_tag', ev.region, 'origin', ev.origin, 'on', on, 'off', off, ...
        'window', [t_lo t_hi] - on, 'margin_s', margin, 'cfg', s, 'ch', []);
    [res.on_natus, res.on_pipeline] = clock_text(on, R(1));
    [res.off_natus, res.off_pipeline] = clock_text(off, R(1));

    ch = struct([]);
    for k = 1:numel(R)
        ch = [ch, one_region(R(k), ev, ctx, s, on, t_lo, t_hi, margin)]; %#ok<AGROW>
    end
    res.ch = ch;

    % common colour limits per map type (same scale for both regions)
    res.clim_mt = pick_clim(arrayfun(@(c) c.mt.Z, ch, 'UniformOutput', false), s);
    res.clim_morlet = pick_clim(arrayfun(@(c) c.morlet.Z, ch, 'UniformOutput', false), s);
    res.norm_label = ch(1).norm_label;
    res.features = vertcat(ch.features);
end

%% ======================================================================
function c = one_region(r, ev, ctx, s, on, t_lo, t_hi, margin)
    z = s.seizure;
    fs = r.fs;
    warns = {};
    seg = spec_load_segment(r, t_lo - margin, t_hi + margin, struct('spec', s), ctx.cache_dir);
    core = seg.t_rel >= t_lo & seg.t_rel <= t_hi;
    nan_frac = mean(~seg.valid(core));
    if nan_frac > s.max_nan_frac
        error('spec_seizure_event:TooManyNaN', '%s %s: %.1f%% of the window is NaN (> %.1f%%).', ...
            ev.id, r.region, 100 * nan_frac, 100 * s.max_nan_frac);
    end
    if nan_frac > 0
        warns{end+1} = sprintf('%.2f%% of the window is NaN (gap): those TF windows are blank', 100 * nan_frac);
    end
    if ~isempty(seg.outlier_note)
        warns{end+1} = seg.outlier_note;
    end
    for k = 1:numel(r.notes)
        warns{end+1} = r.notes{k}; %#ok<AGROW>
    end

    [x_tf, x_disp] = spec_prepare_signal(seg.x, fs, struct('spec', s));

    % ---- multitaper ---------------------------------------------------
    [S, tm, f_mt, mt_info] = spec_multitaper(x_tf, fs, z.mt);
    tm = seg.t_rel(1) + tm;                       % t_rel of window centres
    % ---- Morlet --------------------------------------------------------
    [W, mo_info] = spec_morlet(x_tf, fs, z.morlet);
    [P, tp] = block_mean_time(abs(W) .^ 2, seg.t_rel, round(z.morlet.t_step_s * fs));
    clear W

    % ---- baseline -------------------------------------------------------
    [bi, bkind] = spec_baseline_interval(z.baseline, z.baseline_is_t_rel, z.baseline_time_ref, r);
    if strcmp(bkind, 'relative')
        bl = on + bi;
    else
        bl = bi;
    end
    half_w = z.mt.win_s / 2;
    if bl(1) - margin >= seg.t_rel(1) && bl(2) + margin <= seg.t_rel(end)
        Sbl = S(:, tm - half_w >= bl(1) & tm + half_w <= bl(2));
        Pbl = P(:, tp >= bl(1) & tp <= bl(2));
        bl_nan = mean(~seg.valid(seg.t_rel >= bl(1) & seg.t_rel <= bl(2)));
    else
        bseg = spec_load_segment(r, bl(1) - margin, bl(2) + margin, struct('spec', s), ctx.cache_dir);
        [bx, ~] = spec_prepare_signal(bseg.x, fs, struct('spec', s));
        [Sb, tb] = spec_multitaper(bx, fs, z.mt);
        tb = bseg.t_rel(1) + tb;
        Sbl = Sb(:, tb - half_w >= bl(1) & tb + half_w <= bl(2));
        [Wb, ~] = spec_morlet(bx, fs, z.morlet);
        Pb = abs(Wb) .^ 2;
        Pbl = Pb(:, bseg.t_rel >= bl(1) & bseg.t_rel <= bl(2));
        bl_nan = mean(~bseg.valid(bseg.t_rel >= bl(1) & bseg.t_rel <= bl(2)));
    end
    if bl_nan > 0
        warns{end+1} = sprintf('baseline [%.1f, %.1f] s t_rel: %.1f%% NaN (gap or outside the recording)', bl(1), bl(2), 100 * bl_nan);
    end
    if isempty(Sbl) || all(isnan(Sbl(:)))
        warns{end+1} = 'no complete multitaper window inside the baseline: maps shown in absolute dB';
    end
    warns = [warns, baseline_conflicts(bl, r, ev, ctx)];

    % ---- normalized maps, cropped to the display window -----------------
    [Zmt, label, used] = spec_normalize(S, s.norm, Sbl);
    [Zmo, ~, ~] = spec_normalize(P, s.norm, Pbl);
    if ~strcmp(used, s.norm)
        warns{end+1} = sprintf('norm ''%s'' fell back to ''%s'' (no usable baseline)', s.norm, used);
    end
    keep_mt = tm >= t_lo & tm <= t_hi;
    keep_mo = tp >= t_lo & tp <= t_hi;
    keep_x = seg.t_rel >= t_lo & seg.t_rel <= t_hi;

    % ---- PSDs (dB re 1 uV^2/Hz) -----------------------------------------
    ron = r.on; roff = r.off;
    inside = @(a, b) tm - half_w >= a & tm + half_w <= b;
    psd_bl = 10 * log10(mean(Sbl, 2, 'omitnan'));
    psd_ic = 10 * log10(mean(S(:, inside(ron, roff)), 2, 'omitnan'));
    psd_po = NaN(size(psd_ic));
    if z.postictal_s > 0
        psd_po = 10 * log10(mean(S(:, inside(roff, roff + z.postictal_s)), 2, 'omitnan'));
    end
    if ~any(inside(ron, roff))
        warns{end+1} = sprintf('seizure shorter than one multitaper window (%.1f s): ictal PSD empty', z.mt.win_s);
    end

    % ---- features --------------------------------------------------------
    % the notch dip (+- its half-width + half the multitaper resolution) is no frequency of the seizure:
    % excluded from peak / dominant-frequency searches (band powers keep it, it is in both terms)
    [fex, fex_txt] = spec_notch_band(r, seg.notch_applied, mt_info.resolution_hz / 2);
    fsel = f_mt >= z.feature_fmin & f_mt <= z.mt.fmax & ~fex(f_mt);
    feat = struct();
    [~, i1] = max(psd_ic(fsel));
    fv = f_mt(fsel);
    feat.peak_freq_hz = pick(fv, i1, psd_ic(fsel));
    rel = psd_ic - psd_bl;
    [~, i2] = max(rel(fsel));
    feat.peak_freq_rel_hz = pick(fv, i2, rel(fsel));
    res_ic = aperiodic_residual(psd_ic(fsel) / 10, fv);   % log10 power minus its 1/f fit
    [~, i3] = max(res_ic);
    feat.peak_freq_osc_hz = pick(fv, i3, res_ic);
    lin_bl = 10 .^ (psd_bl / 10);
    lin_ic = 10 .^ (psd_ic / 10);
    for b = 1:numel(z.bands)
        m = f_mt >= z.bands(b).range(1) & f_mt <= z.bands(b).range(2);
        feat.(['rel_power_' z.bands(b).name '_db']) = 10 * log10(sum(lin_ic(m), 'omitnan') / sum(lin_bl(m), 'omitnan'));
    end
    ictal_cols = inside(ron, roff) & ~all(isnan(S), 1);
    has_bl = ~isempty(Sbl) && ~all(isnan(Sbl(:)));
    switch z.domfreq_method
        case 'aperiodic'       % peak above the 1/f background of each window
            D = aperiodic_residual(log10(S(fsel, :)), fv);
        case 'baseline_ratio'  % max of the power re baseline
            if has_bl
                D = 10 * log10(S(fsel, :) ./ mean(Sbl(fsel, :), 2, 'omitnan'));
            else
                D = 10 * log10(S(fsel, :));
            end
        case 'absolute_band'   % max of the absolute power within domfreq_band
            D = 10 * log10(S(fsel, :));
            D(fv < z.domfreq_band(1) | fv > z.domfreq_band(2), :) = -Inf;
        otherwise
            error('spec_seizure_event:BadDomfreq', ...
                'seizure.domfreq_method must be aperiodic, baseline_ratio or absolute_band (got ''%s'').', z.domfreq_method);
    end
    [~, im] = max(D, [], 1);
    dom = fv(im)';
    dom(all(isnan(D) | isinf(D), 1)) = NaN;
    dom_t = tm(ictal_cols) - on;
    dom_f = dom(ictal_cols);
    feat.domfreq_mean_hz = mean(dom_f, 'omitnan');
    ok = ~isnan(dom_f);
    if sum(ok) >= 2
        pf = polyfit(dom_t(ok), dom_f(ok), 1);
        feat.domfreq_slope_hz_per_s = pf(1);
    else
        feat.domfreq_slope_hz_per_s = NaN;
    end

    % ---- assemble -------------------------------------------------------------
    c = struct();
    c.region = r.region;
    c.on = ron;
    c.off = roff;
    [c.on_natus, c.on_pipeline] = clock_text(ron, r);
    [c.off_natus, c.off_pipeline] = clock_text(roff, r);
    c.baseline = bl;
    c.baseline_kind = bkind;
    c.t_trace = seg.t_rel(keep_x) - on;
    c.trace = x_disp(keep_x);
    c.signal = seg.x;                  % the loaded window as read (uV), margins included
    c.signal_t_rel = seg.t_rel;
    c.mt = struct('t', tm(keep_mt) - on, 'f', f_mt, 'S', S(:, keep_mt), 'Z', Zmt(:, keep_mt), 'info', mt_info);
    c.morlet = struct('t', tp(keep_mo) - on, 'f', mo_info.f, 'P', P(:, keep_mo), 'Z', Zmo(:, keep_mo), ...
        'coi_s', mo_info.coi_s, 'cycles', mo_info.cycles);
    c.psd = struct('f', f_mt, 'baseline_db', psd_bl, 'ictal_db', psd_ic, 'postictal_db', psd_po, ...
        'excluded', fex(f_mt), 'excluded_txt', fex_txt);
    c.domfreq = struct('t', dom_t, 'f', dom_f);
    c.norm_label = label;
    c.norm_used = used;
    c.nan_frac = nan_frac;
    c.src = struct('txt', r.txt, 'clean', r.clean, 'raw', r.raw, 'iid_mat', r.iid_mat, ...
        'signal_source', seg.source, 'gain_applied', seg.gain_applied, 'notch_applied', seg.notch_applied, ...
        'case_applied', seg.case_applied, 'quality_class', meta_text(r.clean_hdr, 'quality_class'), ...
        'n_outliers_session', str2double(meta_text(r.clean_hdr, 'n_outliers')), ...
        'n_outliers_window', seg.n_outliers_segment, 'excel_attenuation', r.excel_attenuation, ...
        'header', r.hdr.meta, 'cache_build_s', seg.load_info.build_s, 'cache_reused', seg.load_info.reused);
    c.sz = region_row(ev, r.region);
    c.warnings = warns;
    c.features = feature_row(c, r, ev, feat, s);
end

%% ======================================================================
function R = aperiodic_residual(L, f)
% L: log10 power (n_f x n_windows). Fits L = b0 - chi * log10(f) per column
% (least squares over every non-NaN frequency) and returns L minus the fit:
% > 0 where the spectrum rises above its 1/f background (an oscillation).
    X = [ones(numel(f), 1), log10(f(:))];
    R = NaN(size(L));
    ok = all(isfinite(L), 1);
    if any(ok)
        B = X \ L(:, ok);
        R(:, ok) = L(:, ok) - X * B;
    end
end

function [P2, t2] = block_mean_time(P, t_rel, k)
% Average power over blocks of k samples (no aliasing of the power envelope).
    k = max(1, k);
    n = floor(size(P, 2) / k) * k;
    P2 = squeeze(mean(reshape(P(:, 1:n), size(P, 1), k, []), 2));
    if size(P, 1) == 1
        P2 = P2(:)';
    end
    t2 = mean(reshape(t_rel(1:n), k, []), 1);
end

function w = baseline_conflicts(bl, r, ev, ctx)
    w = {};
    if ~isempty(ctx.seizures)
        Sz = ctx.seizures;
        [~, stem] = fileparts(r.source_file);
        stems = cellfun(@stem_of, Sz.source_file, 'UniformOutput', false);
        m = strcmp(Sz.subject_id, r.subject) & strcmpi(stems, stem) & Sz.end_s >= bl(1) & Sz.start_s <= bl(2);
        own = [];
        if ~isempty(ev.sz)
            own = [ev.sz.seizure_id];
        end
        m = m & ~ismember(Sz.seizure_id, own);
        if any(m)
            ids = unique(Sz.seizure_id(m));
            w{end+1} = sprintf('baseline overlaps other seizure(s) of seizures_events: seizure_id %s', mat2str(ids'));
        end
    end
    if ~isempty(r.iid_mat)
        try
            S = load(r.iid_mat, 'iid_results');
            ex = S.iid_results.exclusion_table;
            m = ex.end_s >= bl(1) & ex.start_s <= bl(2) & string(ex.reason) ~= "seizure";
            if any(m)
                w{end+1} = sprintf('baseline overlaps %d IID exclusion zone(s) (%s)', sum(m), strjoin(unique(cellstr(string(ex.reason(m))))', ', '));
            end
        catch ME
            w{end+1} = sprintf('could not read the IID exclusion zones: %s', ME.message);
        end
    end
end

function s = stem_of(f)
    [~, s] = fileparts(f);
end

function v = pick(fv, i, vals)
    if isempty(i) || all(isnan(vals))
        v = NaN;
    else
        v = fv(i);
    end
end

function [n, p] = clock_text(t, r)
    [a, b] = spec_resolve_times('clocks', t, r);
    n = spec_resolve_times('text', a);
    p = spec_resolve_times('text', b);
end

function s = meta_text(h, name)
    s = '';
    if ~isempty(h) && isfield(h.meta, name)
        s = strtrim(char(h.meta.(name)));
    end
end

function row = region_row(ev, region)
    row = struct('seizure_id', NaN, 'event_id', NaN, 'user_check', '', 'comments', '', 'category', '', ...
        'detection_status', '', 'll_status', '', 'row', NaN);
    if isempty(ev.sz)
        return;
    end
    k = find(strcmp({ev.sz.region}, region), 1);
    if isempty(k)
        return;
    end
    for f = fieldnames(row)'
        if isfield(ev.sz, f{1})
            row.(f{1}) = ev.sz(k).(f{1});
        end
    end
end

function T = feature_row(c, r, ev, feat, s)
    ss = r.session_start;
    ss.Format = 'dd-MMM-yyyy HH:mm:ss.SSS';
    T = table({ev.id}, {r.subject}, {r.region}, {r.session_key}, {r.source_file}, {char(ss)}, ...
        c.sz.seizure_id, c.sz.event_id, {c.sz.user_check}, {c.sz.comments}, {c.sz.category}, {c.sz.detection_status}, ...
        c.on, c.off, c.off - c.on, {c.on_natus}, {c.on_pipeline}, {c.off_natus}, {c.off_pipeline}, ...
        c.baseline(1), c.baseline(2), {s.signal_source}, {c.norm_used}, ...
        {c.src.case_applied}, c.src.gain_applied, c.src.notch_applied, c.src.n_outliers_session, c.src.n_outliers_window, ...
        c.nan_frac, feat.peak_freq_hz, feat.peak_freq_rel_hz, ...
        'VariableNames', {'id', 'subject_id', 'region', 'session', 'source_file', 'session_start', ...
        'seizure_id', 'event_id', 'user_check', 'comments', 'category', 'detection_status', ...
        'onset_t_rel', 'offset_t_rel', 'duration_s', 'onset_natus', 'onset_pipeline', 'offset_natus', 'offset_pipeline', ...
        'baseline_start_t_rel', 'baseline_end_t_rel', 'signal_source', 'norm', ...
        'case_applied', 'gain_applied', 'notch_applied', 'n_outliers_session', 'n_outliers_window', ...
        'nan_frac', 'peak_freq_hz', 'peak_freq_rel_hz'});
    for b = 1:numel(s.seizure.bands)
        nm = ['rel_power_' s.seizure.bands(b).name '_db'];
        T.(nm) = feat.(nm);
    end
    T.peak_freq_osc_hz = feat.peak_freq_osc_hz;
    T.domfreq_method = {s.seizure.domfreq_method};
    T.domfreq_mean_hz = feat.domfreq_mean_hz;
    T.domfreq_slope_hz_per_s = feat.domfreq_slope_hz_per_s;
    T.freq_excluded = {c.psd.excluded_txt};
    T.excel_attenuation = {c.src.excel_attenuation};
    T.warnings = {strjoin(c.warnings, ' | ')};
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
