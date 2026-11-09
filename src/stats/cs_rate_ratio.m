function r = cs_rate_ratio(n1, e1, n0, e0, alpha)
% CS_RATE_RATIO  Poisson rate ratio (n1/e1)/(n0/e0) with the exact
% conditional CI: given n = n1+n0, n1 ~ Binomial(n, pi) with
% pi = e1*RR/(e1*RR + e0); Clopper-Pearson bounds on pi are mapped back to
% RR. p = two-sided exact binomial test of pi0 = e1/(e1+e0) (RR = 1).
% Exposures in any common unit (never assumed equal: light 11.5 h vs dark
% 12.5 h enter through e1/e0).

    if nargin < 5
        alpha = 0.05;
    end
    r = struct('n1', n1, 'e1', e1, 'n0', n0, 'e0', e0, 'rate1', n1 / e1, 'rate0', n0 / e0, ...
        'rr', NaN, 'ci_lo', NaN, 'ci_hi', NaN, 'p', NaN);
    n = n1 + n0;
    if n == 0 || e1 <= 0 || e0 <= 0
        return;
    end
    r.rr = (n1 / e1) / (n0 / e0);
    if n1 == 0
        plo = 0;
    else
        plo = betaincinv(alpha / 2, n1, n - n1 + 1);
    end
    if n1 == n
        phi = 1;
    else
        phi = betaincinv(1 - alpha / 2, n1 + 1, n - n1);
    end
    r.ci_lo = plo / (1 - plo) * e0 / e1;
    r.ci_hi = phi / (1 - phi) * e0 / e1;
    pi0 = e1 / (e1 + e0);
    pk = binopdf(0:n, n, pi0);
    r.p = min(1, sum(pk(pk <= binopdf(n1, n, pi0) * (1 + 1e-7))));
end
