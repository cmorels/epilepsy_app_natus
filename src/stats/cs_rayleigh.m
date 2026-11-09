function r = cs_rayleigh(minute_of_day, hour_exposure_s, n_boot, seed)
% CS_RAYLEIGH  Rayleigh test of 24-h clock uniformity of event onsets
% (no toolbox), plus an exposure-corrected version.
%
%   Standard: theta = 2*pi*minute/1440, Rbar = |mean(exp(i*theta))|,
%   Z = n*Rbar^2, p from Zar (1999) eq. 27.4. It assumes uniform
%   observation over the clock, which recordings do not guarantee.
%   Exposure-corrected: p_exposure = fraction of n_boot simulated samples
%   of n onsets, drawn from the observed hourly exposure (uniform within the
%   hour), whose Rbar >= the observed Rbar.
%
%   r: n, rbar, mean_clock_h (circular mean, hours), Z, p, p_exposure

    th = 2 * pi * minute_of_day(:) / 1440;
    n = numel(th);
    r = struct('n', n, 'rbar', NaN, 'mean_clock_h', NaN, 'Z', NaN, 'p', NaN, 'p_exposure', NaN);
    if n < 2
        return;
    end
    m = mean(exp(1i * th));
    r.rbar = abs(m);
    r.mean_clock_h = mod(angle(m), 2 * pi) / (2 * pi) * 24;
    r.Z = n * r.rbar ^ 2;
    r.p = min(1, exp(sqrt(1 + 4 * n + 4 * (n ^ 2 - (n * r.rbar) ^ 2)) - (1 + 2 * n)));

    w = hour_exposure_s(:)' / sum(hour_exposure_s);
    if any(~isfinite(w)) || n_boot < 1
        return;
    end
    rng(seed);
    hrs = find(w > 0) - 1;  % hours with exposure (0-based)
    edges = [0 cumsum(w(w > 0))];
    edges(end) = 1 + eps;
    hits = 0;
    for b = 1:n_boot
        h = hrs(discretize(rand(n, 1), edges));
        h = h(:);
        ths = 2 * pi * (h + rand(n, 1)) / 24;
        hits = hits + (abs(mean(exp(1i * ths))) >= r.rbar);
    end
    r.p_exposure = (hits + 1) / (n_boot + 1);
end
