function [W, info] = spec_morlet(x, fs, p)
% SPEC_MORLET  Complex Morlet transform by FFT convolution (no Wavelet
% Toolbox), for a column signal that may contain NaN gaps.
%
%   [W, info] = spec_morlet(x, fs, p)
%
%   p : fmin, fmax, nf, cycles ([n_min n_max] or a scalar), cycles_scale
%       ('linear': n grows linearly with f | 'log': linearly with log f)
%   W : complex, nf x numel(x). abs(W).^2 is power in uV^2: a sinusoid of
%       amplitude a at the wavelet frequency gives abs(W) = a.
%   info : f (log-spaced, column), cycles, sigma_t (= n / (2 pi f), s),
%          coi_s (= 3 sigma_t, the cone of influence half-width)
%
% Each wavelet is exp(2 pi i f t) exp(-t^2 / (2 sigma_t^2)) truncated at
% +-3 sigma_t and scaled by 2 / sum(gaussian). An output sample is NaN
% when its +-3 sigma_t support touches a NaN or the ends of x (that is the
% cone of influence): callers load an extra margin of 3 sigma_t(fmin) on
% each side and drop it before plotting.

    x = x(:);
    n = numel(x);
    f = logspace(log10(p.fmin), log10(p.fmax), p.nf)';
    c = p.cycles;
    if isscalar(c)
        cyc = repmat(c, p.nf, 1);
    else
        switch p.cycles_scale
            case 'linear'
                cyc = c(1) + (c(2) - c(1)) * (f - f(1)) / max(f(end) - f(1), eps);
            case 'log'
                cyc = c(1) + (c(2) - c(1)) * (log(f) - log(f(1))) / max(log(f(end)) - log(f(1)), eps);
            otherwise
                error('spec_morlet:BadScale', 'cycles_scale must be ''linear'' or ''log''.');
        end
    end
    sigma_t = cyc ./ (2 * pi * f);
    half = ceil(3 * sigma_t * fs);

    info = struct('f', f, 'cycles', cyc, 'sigma_t', sigma_t, 'coi_s', 3 * sigma_t, 'half_samples', half);

    W = complex(NaN(p.nf, n), NaN(p.nf, n));
    if n == 0
        return;
    end
    nanm = isnan(x);
    x0 = x;
    x0(nanm) = 0;
    nfft = 2 ^ nextpow2(n + 2 * max(half));
    X = fft(x0, nfft);

    for k = 1:p.nf
        h = half(k);
        tw = (-h:h)' / fs;
        g = exp(-tw .^ 2 / (2 * sigma_t(k) ^ 2));
        w = exp(2i * pi * f(k) * tw) .* g * (2 / sum(g));
        y = ifft(X .* fft(w, nfft));
        y = y(h + 1:h + n);
        bad = movsum(double(nanm), [h h]) > 0;
        bad(1:min(h, n)) = true;
        bad(max(1, n - h + 1):n) = true;
        y(bad) = NaN;
        W(k, :) = y.';
    end
end
