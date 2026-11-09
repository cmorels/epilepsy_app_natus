function tests = test_spec_synthetic
% TEST_SPEC_SYNTHETIC  Known signals through the spectral core
% (run: runtests('tests/spectrogram')).
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(repo, 'src'), fullfile(repo, 'src', 'utils'), fullfile(repo, 'src', 'spectrogram'), ...
        fileparts(mfilename('fullpath')));
    tc.TestData.fs = 512;
    tc.TestData.cfg = spectrogram_config();
    rng(1);
end

%% 4 Hz sine + 30 Hz burst --------------------------------------------------
function testSineAndBurst(tc)
    fs = tc.TestData.fs;
    t = (0:60 * fs - 1)' / fs;
    x = 100 * sin(2 * pi * 4 * t) + 2 * randn(size(t));
    burst = t >= 30 & t < 35;
    x(burst) = x(burst) + 100 * sin(2 * pi * 30 * t(burst));
    z = tc.TestData.cfg.spec.seizure;

    [S, tm, f] = spec_multitaper(x, fs, z.mt);
    out = tm < 25 | tm > 40;
    [~, i] = max(mean(S(:, out), 2));
    tc.verifyEqual(f(i), 4, 'AbsTol', 0.5, 'multitaper: 4 Hz peak outside the burst');
    in = tm > 31 & tm < 34;
    hi = f > 10;
    fh = f(hi);
    [~, i] = max(mean(S(hi, in), 2));
    tc.verifyEqual(fh(i), 30, 'AbsTol', 0.5, 'multitaper: 30 Hz peak in the burst');
    [on_mt, off_mt] = burst_edges(tm, S(abs(f - 30) < 0.3, :));
    tc.verifyEqual([on_mt off_mt], [30 35], 'AbsTol', 1.1, 'multitaper: burst timing (2 s windows)');

    [W, info] = spec_morlet(x, fs, z.morlet);
    P = abs(W) .^ 2;
    [~, i4] = min(abs(info.f - 4));
    tc.verifyEqual(sqrt(median(P(i4, t > 5 & t < 25), 'omitnan')), 100, 'RelTol', 0.15, ...
        'Morlet: amplitude normalisation (|W| = sinusoid amplitude)');
    [~, i] = max(mean(P(info.f > 10, burst & t > 31 & t < 34), 2));
    ff = info.f(info.f > 10);
    tc.verifyEqual(ff(i), 30, 'AbsTol', 2, 'Morlet: 30 Hz peak in the burst');
    [~, i30] = min(abs(info.f - 30));
    [on_mo, off_mo] = burst_edges(t', P(i30, :));
    tc.verifyEqual([on_mo off_mo], [30 35], 'AbsTol', 0.15, 'Morlet: burst timing');
end

function [a, b] = burst_edges(t, p)
    p = p(:)';
    m = p > 0.5 * max(p);
    k = find(m);
    a = t(k(1));
    b = t(k(end));
end

%% chirp 2 -> 20 Hz followed by the dominant frequency -----------------------
function testChirpDominantFrequency(tc)
    fs = tc.TestData.fs;
    T = 20;
    t = (0:T * fs - 1)' / fs;
    finst = 2 + (20 - 2) * t / T;
    x = 100 * chirp(t, 2, T, 20) + 2 * randn(size(t));
    z = tc.TestData.cfg.spec.seizure;
    [S, tm, f] = spec_multitaper(x, fs, z.mt);
    [~, i] = max(S(f >= 1, :), [], 1);
    fsel = f(f >= 1);
    dom = fsel(i);
    truth = 2 + 18 * tm(:) / T;
    inner = tm > 1.5 & tm < T - 1.5;
    err = abs(dom(inner) - truth(inner));
    tc.verifyLessThan(median(err), 0.75, 'multitaper dominant frequency follows the chirp');
    tc.verifyLessThan(max(err), 2, 'multitaper dominant frequency: worst window');

    [W, info] = spec_morlet(x, fs, z.morlet);
    P = abs(W) .^ 2;
    ok = all(~isnan(P), 1) & t' > 1.5 & t' < T - 1.5;
    [~, i] = max(P(:, ok), [], 1);
    errm = abs(info.f(i) - finst(ok));
    tc.verifyLessThan(median(errm ./ finst(ok)), 0.08, 'Morlet ridge follows the chirp (relative error)');
end

%% Gaussian spikes + pink noise: re-alignment within +-2 ms --------------------
function testSpikeRealignment(tc)
    % default (align_on = 'analysis', broadband signal) with pink noise of
    % 5 uV std, and align_on = 'display' (1-100 Hz) with 15 uV
    for c = {{'analysis', 5}, {'display', 15}}
        [err_ms, res] = realign_case(tc, c{1}{1}, c{1}{2});
        tc.verifyEqual(res.peaks.polarity, 'neg');
        tc.verifyGreaterThan(sum(res.peaks.used), 45);
        tc.verifyLessThanOrEqual(max(err_ms), 2, sprintf('align_on=%s, pink noise %g uV: max error %.2f ms', ...
            c{1}{1}, c{1}{2}, max(err_ms)));
        tc.verifyEqual(res.features.fwhm_mean_ms, 30, 'AbsTol', 6, 'FWHM of 30 ms spikes (1-100 Hz trace)');
    end
end

function [err_ms, res] = realign_case(tc, align_on, noise_uV)
    fs = tc.TestData.fs;
    root = tempname;
    cleanup = onCleanup(@() rmdir(root, 's'));
    rng(3);
    n = 120 * fs;
    t = (0:n - 1)' / fs;
    x = noise_uV * pink(n);
    peaks = (10:2:110)' + 0.3 * rand(51, 1);          % true peak times (off the sample grid)
    sig = 0.030 / (2 * sqrt(2 * log(2)));              % FWHM 30 ms
    for k = 1:numel(peaks)
        x = x - 400 * exp(-(t - peaks(k)) .^ 2 / (2 * sig ^ 2));
    end
    ss = datetime(2026, 7, 1, 10, 0, 0, 'TimeZone', 'Europe/Paris');
    make_spec_session(fullfile(root, 'S1', 'EDFA'), 'S1', ss, 'HPCr', x, fs);
    spec_resolve_files('reset');
    spec_signal_cache('reset');
    r = spec_resolve_files(root, 'S1', '20260701_100000', 'HPCr', 'clean', 'Europe/Paris');
    r.on = 5; r.off = 115; r.excel_attenuation = 'NA';
    jitter = (rand(size(peaks)) - 0.5) * 0.030;         % given peaks off by up to +-15 ms
    ev = spec_event('iid_t', 'iid', 'S1', '20260701_100000', 'HPCr', 5, 115, 't_rel', struct(), 'test');
    ev.peaks = peaks + jitter;
    cfg = tc.TestData.cfg;
    cfg.spec.iid.spike_source = 'manual';
    cfg.spec.iid.align_on = align_on;
    cfg.spec.pipeline_root = root;
    ctx = struct('seizures', [], 'cache_dir', fullfile(root, 'cache'));
    res = spec_iid_interval(ev, r, ctx, cfg);
    used = res.peaks.used;
    err_ms = abs(res.peaks.t_aligned(used) - peaks(used)) * 1000;
    spec_signal_cache('reset');
end

%% zone mode (CWT in blocks) == interval mode (CWT of the whole interval) --------
function testZoneModeEqualsIntervalMode(tc)
    fs = tc.TestData.fs;
    root = tempname;
    cleanup = onCleanup(@() rmdir(root, 's'));
    rng(4);
    n = 200 * fs;
    t = (0:n - 1)' / fs;
    x = 20 * pink(n) + 30 * sin(2 * pi * 9 * t);
    peaks = (12:1.7:188)';
    for k = 1:numel(peaks)
        x = x - 350 * exp(-(t - peaks(k)) .^ 2 / (2 * 0.010 ^ 2));
    end
    x(round(100.2 * fs):round(100.6 * fs)) = NaN;       % a gap: same exclusions in both modes
    ss = datetime(2026, 7, 1, 10, 0, 0, 'TimeZone', 'Europe/Paris');
    make_spec_session(fullfile(root, 'S1', 'EDFA'), 'S1', ss, 'HPCr', x, fs);
    spec_resolve_files('reset');
    spec_signal_cache('reset');
    r = spec_resolve_files(root, 'S1', '20260701_100000', 'HPCr', 'clean', 'Europe/Paris');
    r.on = 10; r.off = 190; r.excel_attenuation = 'NA';
    ev = spec_event('z', 'iid', 'S1', '20260701_100000', 'HPCr', 10, 190, 't_rel', struct(), 'test');
    ev.peaks = peaks;
    cfg = tc.TestData.cfg;
    cfg.spec.iid.spike_source = 'manual';
    ctx = struct('seizures', [], 'cache_dir', fullfile(root, 'cache'));
    a = spec_iid_interval(ev, r, ctx, cfg);
    cfg.spec.iid.max_interval_s = 30;
    cfg.spec.iid.zone_chunk_s = 7;                      % many blocks, spikes on block edges
    z = spec_iid_interval(ev, r, ctx, cfg);
    tc.verifyEqual(a.mode, 'interval');
    tc.verifyEqual(z.mode, 'zone');
    tc.verifyEqual(z.peaks.used, a.peaks.used);
    tc.verifyGreaterThan(sum(z.peaks.used), 90);
    tc.verifyGreaterThan(sum(strcmp(z.peaks.reason, 'nan')), 0, 'the gap excludes some epochs');
    tc.verifyEqual(z.induced.P, a.induced.P, 'RelTol', 1e-9);
    tc.verifyEqual(z.evoked.P, a.evoked.P, 'RelTol', 1e-9);
    tc.verifyEqual(z.itc, a.itc, 'AbsTol', 1e-9);
    tc.verifyEqual(z.epochs.mean, a.epochs.mean, 'AbsTol', 1e-9);
    tc.verifyEqual(z.features.fwhm_mean_ms, a.features.fwhm_mean_ms, 'AbsTol', 1e-9);
    tc.verifyEqual(sum(z.rate.detected .* z.rate.bin_s / 60), numel(peaks), 'AbsTol', 1e-9);
    fig = spec_plot_iid(z);
    tc.verifyTrue(isgraphics(fig));
    close(fig);
    spec_signal_cache('reset');
end

function y = pink(n)
    X = fft(randn(n, 1));
    k = (0:n - 1)';
    k = min(k, n - k);
    k(1) = 1;
    y = real(ifft(X ./ sqrt(k)));
    y = y / std(y);
end

%% multitaper PSD vs pwelch ----------------------------------------------------
function testPsdMatchesPwelch(tc)
    fs = tc.TestData.fs;
    n = 600 * fs;
    [b, a] = butter(1, 30 / (fs / 2));
    x = 50 * filter(b, a, randn(n, 1)) + 5 * randn(n, 1);   % smooth coloured spectrum
    z = tc.TestData.cfg.spec.seizure.mt;
    z.step_s = 1;
    [S, ~, f] = spec_multitaper(x, fs, z);
    psd_mt = mean(S, 2);
    [pw, fw] = pwelch(x, hann(1024), 512, 1024, fs);
    pw_i = interp1(fw, pw, f);
    band = f >= 2 & f <= 100;
    d = abs(10 * log10(psd_mt(band) ./ pw_i(band)));
    tc.verifyLessThan(max(d), 1, sprintf('multitaper vs pwelch: max %.2f dB in 2-100 Hz', max(d)));
end

%% NaN gap is never crossed ------------------------------------------------------
function testGapNotCrossed(tc)
    fs = tc.TestData.fs;
    n = 40 * fs;
    t = (0:n - 1)' / fs;
    x = 100 * sin(2 * pi * 8 * t) + 10 * randn(n, 1);
    gap = t >= 20 & t < 22;
    x(gap) = NaN;
    x(~gap & t >= 22) = x(~gap & t >= 22) + 5000;     % step after the gap: must not leak backwards
    cfg = tc.TestData.cfg;
    [x_tf, x_disp] = spec_prepare_signal(x, fs, cfg);
    tc.verifyTrue(all(isnan(x_disp(gap))) && all(isnan(x_tf(gap))));
    b1 = find(~gap & t < 20);
    [~, d1] = spec_prepare_signal(x(b1), fs, cfg);
    tc.verifyEqual(x_disp(b1), d1, 'AbsTol', 1e-9, 'block before the gap filtered on its own');
    tc.verifyLessThan(max(abs(x_tf(b1))), 500, 'detrend per block: the step after the gap does not leak');

    [S, tm] = spec_multitaper(x, fs, cfg.spec.seizure.mt);
    touches = tm + 1 > 20 & tm - 1 < 22;
    tc.verifyTrue(all(all(isnan(S(:, touches)))) && ~any(any(isnan(S(:, ~touches)))));

    [W, info] = spec_morlet(x, fs, cfg.spec.seizure.morlet);
    k = 1;                                             % lowest frequency, widest support
    h = info.half_samples(k);
    i0 = find(gap, 1);
    tc.verifyTrue(all(isnan(W(k, i0 - h:i0 - 1))), 'Morlet: cone of influence before the gap is NaN');
    tc.verifyFalse(any(isnan(W(k, h + 1:i0 - h - 1))), 'Morlet: valid away from the gap');
end
