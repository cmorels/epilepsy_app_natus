function [S, t, f, info] = spec_multitaper(x, fs, p)
% SPEC_MULTITAPER  Multitaper spectrogram (DPSS tapers, Signal Processing
% Toolbox dpss) of a column signal that may contain NaN gaps.
%
%   [S, t, f, info] = spec_multitaper(x, fs, p)
%
%   p : win_s, step_s, TW, K, fmin, fmax  (spectrogram_config: cfg.spec.seizure.mt)
%   S : one-sided PSD, uV^2/Hz, numel(f) x numel(t); a window touching a
%       NaN is NaN (never interpolated across a gap)
%   t : window centres, s from x(1) (sample i is at (i-1)/fs)
%   f : Hz, within [fmin, fmax]
%
% Per window: linear detrend, K tapers of time-half-bandwidth TW, nfft =
% next power of 2, PSD = mean over tapers of |FFT|^2 / fs, doubled except
% at DC / Nyquist -- the same scaling as pwelch, so mean(S, 2) over a
% stationary stretch matches pwelch. Frequency resolution 2W = 2 TW / win_s
% (TW = 3, 2 s -> 3 Hz).

    x = x(:);
    n = numel(x);
    win_n = round(p.win_s * fs);
    step_n = max(1, round(p.step_s * fs));
    nfft = 2 ^ nextpow2(win_n);
    fall = (0:nfft / 2)' * fs / nfft;
    keep = fall >= p.fmin & fall <= p.fmax;
    f = fall(keep);

    info = struct('win_n', win_n, 'step_n', step_n, 'nfft', nfft, 'TW', p.TW, 'K', p.K, ...
        'resolution_hz', 2 * p.TW / p.win_s);

    if n < win_n
        S = zeros(numel(f), 0);
        t = zeros(1, 0);
        return;
    end
    starts = 1:step_n:(n - win_n + 1);
    t = (starts - 1 + (win_n - 1) / 2) / fs;
    tapers = get_tapers(win_n, p.TW, p.K);

    nw = numel(starts);
    S = NaN(numel(f), nw);
    chunk = 1000;
    scale = ones(nfft / 2 + 1, 1) * 2;
    scale(1) = 1;
    scale(end) = 1;
    for c0 = 1:chunk:nw
        cols = c0:min(nw, c0 + chunk - 1);
        idx = (0:win_n - 1)' + starts(cols);
        X = x(idx);
        bad = any(isnan(X), 1);
        X(:, bad) = 0;
        X = detrend(X);
        P = zeros(nfft / 2 + 1, numel(cols));
        for k = 1:size(tapers, 2)
            F = fft(X .* tapers(:, k), nfft);
            P = P + abs(F(1:nfft / 2 + 1, :)) .^ 2;
        end
        P = P / size(tapers, 2) / fs .* scale;
        P(:, bad) = NaN;
        S(:, cols) = P(keep, :);
    end
end

function E = get_tapers(n, TW, K)
    persistent cache
    if isempty(cache)
        cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end
    key = sprintf('%d_%.6g_%d', n, TW, K);
    if ~isKey(cache, key)
        cache(key) = dpss(n, TW, K);
    end
    E = cache(key);
end
