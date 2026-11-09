function r = cs_iid_sens_channel(cfg, ch, sz, man)
% CS_IID_SENS_CHANNEL  Phase C of v2 (section 5): ONE read and ONE filtering
% of a channel's *_clean.txt, then, in memory:
%
% 5.1 flat segments: |x| < cfg.flat_thr_uV continuously for more than
%     cfg.flat_min_s, over valid (non-NaN) samples, excluded zones included.
%     Nothing is excluded because of them (list only).
% 5.2 IID threshold sensitivity (descriptive): detect_iid.m's processing
%     reproduced step by step with the verbatim copies in src/stats/iid_sens/
%     -- the same inclusion mask as the main run (sz, man: the zones
%     cs_iid_channel_zones.m gives), median subtraction, 15-70 Hz bandpass
%     per gap-free block, rectification, findpeaks (prominence 200 uV,
%     distance 30 ms, width 60 ms), upper 2000 uV, complexes (ISI <= 150 ms)
%     and bursts -- once per lower-threshold variant: 'actual' (2.5 x the
%     channel's baseline), cfg.sens_fixed_uV, and cfg.sens_pct_bp percentiles
%     of |x_bp| over the included samples. findpeaks is re-run for each
%     distinct threshold (filtering the lowest variant's peaks by height is
%     not equivalent because of MinPeakDistance); identical thresholds share
%     the run. The prominence (200 uV) is the same in every variant.
%
% r: channel keys, analyzed_min, baseline/actual threshold, P97 of |x|
% (unfiltered, uV), P99/P99.5/P99.9 of |x_bp| (uV), percentile of |x_bp| at
% 312.5 uV, flat segments (t_rel s), and r.variants (name, thr_uV, n_complexes,
% n_single, n_polyspike, n_bursts).

    t0 = tic;
    pc = pipeline_config();
    pc.iid.exclusion_zones_manual = reshape(man, [], 2);
    ensure_findpeaks_signal_toolbox();
    data = load_lfp_txt(ch.file{1});
    signal = data.signal(:);
    fs = data.fs;
    t_rel = data.t_rel(:);
    valid = data.valid_mask(:);
    t_load = toc(t0);

    % ---- 5.1 flat segments -------------------------------------------------
    seg = mask_to_segments(valid & abs(signal) < cfg.flat_thr_uV);
    dur = (seg(:, 2) - seg(:, 1) + 1) / fs;
    seg = seg(dur > cfg.flat_min_s, :);
    flat = [t_rel(seg(:, 1)), t_rel(seg(:, 2)) + 1 / fs];

    % ---- 5.2 detect_iid's mask, baseline and band-passed signal --------------
    [~, incl] = iids_build_exclusion_zones(sz, valid, t_rel, pc);
    Voltage_mV = signal / 1000;
    [baseline_uV, lower_uV, upper_uV] = iids_compute_baseline(abs(Voltage_mV(incl)), pc);
    upper_mV = upper_uV / 1000;
    raw_blocks = mask_to_segments(valid);
    x_mV = (signal - median(signal(valid))) / 1000;
    x_bp_full = filter_by_blocks(x_mV, raw_blocks, ...
        @(v) bandpass(v, pc.iid.bandpass_band, fs, 'Steepness', pc.iid.bandpass_steepness));
    abs_mV_full = abs(x_bp_full);
    clear x_mV x_bp_full Voltage_mV
    a_in = abs_mV_full(incl);
    pbp = prctile(a_in, cfg.sens_pct_bp) * 1000;
    r.p97_x_uV = prctile(abs(signal(incl)), 97);
    r.pct_bp_uV = pbp(:)';
    r.pctile_of_312_5 = 100 * mean(a_in < 0.3125);
    clear a_in signal

    names = [{'actual'}, arrayfun(@(v) sprintf('fijo %g uV', v), cfg.sens_fixed_uV, 'UniformOutput', false), ...
        arrayfun(@(p) sprintf('P%g |x_bp|', p), cfg.sens_pct_bp, 'UniformOutput', false)];
    thr_uV = [lower_uV, cfg.sens_fixed_uV, pbp(:)'];
    nv = numel(thr_uV);
    res = NaN(nv, 4);   % n_complexes, n_single, n_polyspike, n_bursts
    ws = warning('off', 'signal:findpeaks:largeMinPeakHeight');   % blocks with nothing above a high threshold: 0 peaks, expected
    restore = onCleanup(@() warning(ws));
    t_var = zeros(nv, 1);
    for v = 1:nv
        same = find(thr_uV(1:v-1) == thr_uV(v), 1);
        if ~isempty(same)
            res(v, :) = res(same, :);
            continue;
        end
        tv = tic;
        lower_mV = thr_uV(v) / 1000;
        if v == 1
            lower_mV = lower_uV / 1000;   % exactly detect_iid's value
        end
        [spikes, pos, ~, ~, blk] = iids_detect_peaks_per_block(abs_mV_full, raw_blocks, fs, lower_mV, pc);
        keep = spikes <= upper_mV;
        spikes = spikes(keep); pos = pos(keep); blk = blk(keep);
        keep2 = incl(pos);
        spikes = spikes(keep2); pos = pos(keep2); blk = blk(keep2);
        pos_s = (pos - 1) / fs;
        [cm, cb] = iids_group_into_complexes(pos, pos_s, spikes, blk, fs, pc);
        bm = iids_group_into_bursts(cm, cb, pc);
        n_poly = 0;
        if ~isempty(cm), n_poly = sum(cm(:, 8) ~= 0); end
        res(v, :) = [size(cm, 1), size(cm, 1) - n_poly, n_poly, size(bm, 1)];
        t_var(v) = toc(tv);
    end

    r.campaign = ch.campaign{1};
    r.mouse = ch.mouse{1};
    r.edf = ch.edf{1};
    r.region = ch.region{1};
    r.analyzed_min = sum(incl) / fs / 60;
    r.baseline_uV = baseline_uV;
    r.lower_uV = lower_uV;
    r.flat = flat;
    r.variants = table(names(:), thr_uV(:), res(:, 1), res(:, 2), res(:, 3), res(:, 4), t_var, 'VariableNames', ...
        {'variant', 'thr_uV', 'n_complexes', 'n_single', 'n_polyspike', 'n_bursts', 'findpeaks_s'});
    r.sz = sz;
    r.man = reshape(man, [], 2);
    r.load_s = t_load;
    r.total_s = toc(t0);
end
