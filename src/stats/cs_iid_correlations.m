function [CR, flags] = cs_iid_correlations(cfg, I, R)
% CS_IID_CORRELATIONS  Section 6.5, epileptic mice only.
%
% Between animals (n = mice; merged scope, both hemispheres):
%   x = IID complexes/min, y = seizures per 24 h (R.mouse, T_seizure_obs).
%   Spearman rho, Kendall tau-b (explicit tie correction), Pearson r on
%   log10(v + eps_v), eps_v = half the smallest non-zero value of v.
%   95 % CI: bootstrap BCa (bootci, cfg.n_boot resamples of mice);
%   p: permutation of y (cfg.n_boot permutations), two-sided.
%   For the re-detected and the original IID rates, with and without 048-s.
%
% Within animals (n = mouse-days with >= cfg.min_valid_h_day of BOTH
% seizure observation and analyzed IID time):
%   rmcorr via ANCOVA fitlm(n_seizures ~ iid_per_h + mouse), mouse
%   categorical; r_rm = sign(b) * sqrt(SS_iid / (SS_iid + SS_error)),
%   p of the iid term, CI by Fisher z with se = 1/sqrt(df_error - 1).
%   Poisson GLMM n_seizures ~ iid10 + (1|mouse), offset log(T_seizure_obs_h),
%   DispersionFlag true (iid10 = IID/h divided by 10 -> IRR per +10 IID/h).
%   Lagged: IID of day d vs seizures of day d+1 (both days valid).

    flags = cell(0, 3);
    rows = {};
    mm = I.mouse(strcmp(I.mouse.scope, 'merged') & strcmp(I.mouse.hemisphere, 'ambos') & strcmp(I.mouse.group, 'epileptic'), :);
    sm = R.mouse(strcmp(R.mouse.scope, 'merged') & strcmp(R.mouse.group, 'epileptic'), :);
    [tf, j] = ismember(mm.mouse, sm.mouse);
    if ~all(tf)
        error('cs_iid_correlations:Join', 'epileptic mice without seizure rate');
    end
    sr = sm.rate_per_24h(j);
    for data = {'re-deteccion', 'original'}
        v = 'complexes_per_min';
        if strcmp(data{1}, 'original'), v = [v '_orig']; end %#ok<AGROW>
        for set = {'todos', 'sin 048-s'}
            keep = true(height(mm), 1);
            if strcmp(set{1}, 'sin 048-s'), keep = ~strcmp(mm.mouse, '048-s'); end
            x = mm.(v)(keep);
            y = sr(keep);
            ex = half_min_nonzero(x);
            ey = half_min_nonzero(y);
            meth = { ...
                'Spearman rho', @(a, b) corr(a, b, 'Type', 'Spearman'), x, y, ''; ...
                'Kendall tau-b', @(a, b) kendall_tau_b(a, b), x, y, ''; ...
                'Pearson r log10', @(a, b) corr(a, b), log10(x + ex), log10(y + ey), ...
                    sprintf('eps IID = %.4g /min, eps crisis = %.4g /24 h (mitad del minimo no-cero)', ex, ey)};
            for k = 1:size(meth, 1)
                f = meth{k, 2};
                a = meth{k, 3}; b = meth{k, 4};
                est = f(a, b);
                [lo, hi, ci_type] = boot_ci(f, a, b, cfg);
                p = perm_p(f, a, b, est, cfg);
                rows(end+1, :) = {'entre animales', meth{k, 1}, data{1}, set{1}, numel(a), NaN, est, lo, hi, ci_type, p, meth{k, 5}}; %#ok<AGROW>
            end
        end
    end

    % ---------------- within animals ------------------------------------
    D = I.day(strcmp(I.day.group, 'epileptic'), :);
    Dv = D(D.valid_both, :);
    rows = [rows; cs_within_assoc(Dv, 'mismo dia')];
    % lag: IID of day d (valid) vs seizures of day d+1 (valid)
    Rd = R.daily(strcmp(R.daily.group, 'epileptic') & R.daily.valid_day, :);
    [y, mo, d] = ymd(Rd.date);
    Rd.daykey = y * 10000 + mo * 100 + d;
    L = D(D.valid_iid_day, {'mouse', 'daykey', 'date', 'iid_per_h'});
    [y, mo, d] = ymd(L.date + caldays(1));
    L.next = y * 10000 + mo * 100 + d;
    [tf, j] = ismember(strcat(L.mouse, '|', string(L.next)), strcat(Rd.mouse, '|', string(Rd.daykey)));
    L = L(tf, :);
    L.n_seizures = Rd.n_seizures(j(tf));
    L.T_seizure_obs_h = Rd.T_seizure_obs_h(j(tf));
    rows = [rows; cs_within_assoc(L, 'IID dia d vs crisis dia d+1')];

    CR = cell2table(rows, 'VariableNames', {'analisis', 'metodo', 'datos', 'conjunto', 'n', 'n_ratones', 'estimacion', ...
        'ci95_lo', 'ci95_hi', 'tipo_ic', 'p', 'nota'});
    failed = startsWith(CR.nota, 'error:');
    if any(failed)
        flags(end+1, :) = {'revisar', 'modelo intra-animal no ajustado', strjoin(CR.nota(failed)', ' | ')};
    end
end

function e = half_min_nonzero(v)
    nz = v(v > 0);
    if isempty(nz)
        e = 1;
    else
        e = min(nz) / 2;
    end
end

function t = kendall_tau_b(x, y)
    n = numel(x);
    s = 0; tx = 0; ty = 0; n0 = n * (n - 1) / 2;
    for i = 1:n - 1
        dx = sign(x(i+1:end) - x(i));
        dy = sign(y(i+1:end) - y(i));
        s = s + sum(dx .* dy);
        tx = tx + sum(dx == 0);
        ty = ty + sum(dy == 0);
    end
    t = s / sqrt((n0 - tx) * (n0 - ty));
end

function [lo, hi, typ] = boot_ci(f, x, y, cfg)
    rng(cfg.rng_seed);
    typ = 'bootstrap BCa';
    try
        ci = bootci(cfg.n_boot, {f, x, y}, 'Type', 'bca');
    catch
        rng(cfg.rng_seed);
        ci = bootci(cfg.n_boot, {f, x, y}, 'Type', 'per');
        typ = 'bootstrap percentil (BCa no calculable)';
    end
    lo = ci(1);
    hi = ci(2);
end

function p = perm_p(f, x, y, est, cfg)
    rng(cfg.rng_seed);
    n = numel(y);
    k = 0;
    for b = 1:cfg.n_boot
        k = k + (abs(f(x, y(randperm(n)))) >= abs(est) - 1e-12);
    end
    p = (k + 1) / (cfg.n_boot + 1);
end
