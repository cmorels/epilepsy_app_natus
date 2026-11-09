function rows = cs_within_assoc(T, label)
% CS_WITHIN_ASSOC  Within-animal association of daily seizures with daily
% IID rate (n = mouse-days; table T: mouse, iid_per_h, n_seizures,
% T_seizure_obs_h):
%   rmcorr via ANCOVA fitlm(n_seizures ~ iid_per_h + mouse), mouse
%   categorical; r_rm = sign(b) * sqrt(SS_iid / (SS_iid + SS_error)), p of
%   the iid term, CI by Fisher z with se = 1/sqrt(df_error - 1).
%   Poisson GLMM n_seizures ~ iid10 + (1|mouse), offset log(T_seizure_obs_h),
%   DispersionFlag true; iid10 = IID/h / 10 -> IRR per +10 IID/h.
% rows: {analisis, metodo, label, conjunto, n, n_ratones, estimacion, ci_lo, ci_hi, tipo_ic, p, nota}

    rows = cell(0, 12);
    T.mouse = categorical(T.mouse);
    nm = numel(unique(T.mouse));
    if height(T) < nm + 3
        rows(end+1, :) = {'intra-animal', 'rmcorr', label, 'epilepticos', height(T), nm, NaN, NaN, NaN, '', NaN, 'muy pocos raton-dias'};
        return;
    end
    try
        mdl = fitlm(T, 'n_seizures ~ iid_per_h + mouse');
        a = anova(mdl, 'component');
        ss = a{'iid_per_h', 'SumSq'};
        se = a{'Error', 'SumSq'};
        b = mdl.Coefficients{'iid_per_h', 'Estimate'};
        r = sign(b) * sqrt(ss / (ss + se));
        df = mdl.DFE;
        z = atanh(r);
        h = 1.96 / sqrt(df - 1);
        rows(end+1, :) = {'intra-animal', 'rmcorr (ANCOVA)', label, 'epilepticos', height(T), nm, r, tanh(z - h), tanh(z + h), ...
            'Fisher z', a{'iid_per_h', 'pValue'}, sprintf('gl error = %d', df)};
    catch ME
        rows(end+1, :) = {'intra-animal', 'rmcorr (ANCOVA)', label, 'epilepticos', height(T), nm, NaN, NaN, NaN, '', NaN, ['error: ' ME.message]};
    end
    try
        T.iid10 = T.iid_per_h / 10;
        g = fitglme(T, 'n_seizures ~ iid10 + (1|mouse)', 'Distribution', 'Poisson', 'Link', 'log', ...
            'Offset', log(T.T_seizure_obs_h), 'DispersionFlag', true);
        k = strcmp(g.Coefficients.Name, 'iid10');
        ci = coefCI(g);
        rows(end+1, :) = {'intra-animal', 'GLMM Poisson, IRR por +10 IID/h', label, 'epilepticos', height(T), nm, ...
            exp(g.Coefficients.Estimate(k)), exp(ci(k, 1)), exp(ci(k, 2)), 'Wald', g.Coefficients.pValue(k), ...
            sprintf('dispersion = %.3g', g.Dispersion)};
    catch ME
        rows(end+1, :) = {'intra-animal', 'GLMM Poisson, IRR por +10 IID/h', label, 'epilepticos', height(T), nm, NaN, NaN, NaN, '', NaN, ...
            ['error: ' ME.message]};
    end
end
