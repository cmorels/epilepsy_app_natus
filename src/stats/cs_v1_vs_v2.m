function T = cs_v1_vs_v2(cfg, S, R, I, CR, DA, PSa, PSl, DW, rv)
% CS_V1_VS_V2  Key figures of v1 (cfg.v1_dir) and v2 side by side, with
% the cause of each change (sheet v1_vs_v2).

    q1 = load(fullfile(cfg.v1_dir, 'phase1_results.mat'), 'R', 'S');
    q3 = load(fullfile(cfg.v1_dir, 'phase3_results.mat'), 'I', 'CR');
    g1 = q1.R.group(strcmp(q1.R.group.group, 'epileptic') & strcmp(q1.R.group.scope, 'merged'), :);
    g2 = R.group(strcmp(R.group.group, 'epileptic') & strcmp(R.group.scope, 'merged'), :);
    fc = @(G, m) G.fold_change_medianas(strcmp(G.metrica, m) & strcmp(G.datos, 're-deteccion') & strcmp(G.conjunto, 'todos'));
    sp = @(C, m) C.estimacion(strcmp(C.analisis, 'entre animales') & strcmp(C.metodo, m) & strcmp(C.datos, 're-deteccion') & strcmp(C.conjunto, 'todos'));
    w1 = @(C, m) C(strcmp(C.analisis, 'intra-animal') & startsWith(C.metodo, m) & strcmp(C.datos, 'mismo dia'), :);
    w2 = @(m, v) DA(startsWith(DA.version, v) & startsWith(DA.metodo, m) & strcmp(DA.datos, 'mismo dia'), :);
    n_ep = @(S) numel(unique(S.units.mouse(S.units.included & strcmp(S.units.group, 'epileptic'))));
    pend = rv.excluded(strcmp(rv.excluded.kind, 'pendiente de reprocesamiento'), :);
    ae = rv.excluded(strcmp(rv.excluded.kind, 'sin AE completo'), :);
    c_scope = sprintf(['%s excluido (%.1f h, pendiente de reprocesamiento); 084-s C2 entra con sus EDF de AE completo ' ...
        '(%d EDF sin AE completo quedan fuera, regla A.2)'], strjoin(strcat(pend.mouse, {' '}, pend.edf)', ', '), ...
        sum(pend.end_posix - pend.start_posix) / 3600, height(ae));
    r1 = w1(q3.CR, 'rmcorr'); r2 = w2('rmcorr', '(i)');
    m1 = w1(q3.CR, 'GLMM'); m2 = w2('GLMM', '(i)');
    pv1 = q3.I.periictal(startsWith(q3.I.periictal.mouse, '('), :);
    rows = {
        'n crisis finales', height(q1.S.seizures), height(S.seizures), c_scope
        'n ratones epilepticos en alcance', n_ep(q1.S), n_ep(S), '084-s (C2) entra en alcance por la regla A.2'
        'T_seizure_obs grupo epileptico (h)', g1.T_seizure_obs_h, g2.T_seizure_obs_h, c_scope
        'tasa grupo epileptico (crisis / 24 h, agrupada)', g1.pooled_rate_per_24h, g2.pooled_rate_per_24h, c_scope
        'fold change complejos/min epi/control', fc(q3.I.group, 'complexes_per_min'), fc(I.group, 'complexes_per_min'), ...
            'medianas por raton con 084-s incluido y 128-s sin 20260729'
        'fold change polyspikes/min epi/control', fc(q3.I.group, 'polyspikes_per_min'), fc(I.group, 'polyspikes_per_min'), ...
            'idem'
        'Spearman IIDs vs crisis entre animales', sp(q3.CR, 'Spearman rho'), sp(CR, 'Spearman rho'), 'idem; n epilepticos cambia'
        'raton-dias en la asociacion diaria', r1.n(1), r2.n_raton_dias(1), sprintf(['v1: dias calendario con >= 12 h; v2: ventanas de ' ...
            '24 h desde la primera muestra de cada raton, >= %g h de T_iid y de T_seizure_obs'], cfg.min_valid_h_day)
        'rmcorr mismo dia (sin exclusion post-ictal)', r1.estimacion(1), r2.estimacion(1), 'definicion de dia (ver fila anterior)'
        'GLMM IRR por +10 IID/h, mismo dia (sin exclusion)', m1.estimacion(1), m2.estimacion(1), 'idem'
        };
    wn = {'pre -60..-30 min', 'pre -30..0 min', 'post 0..30 min', 'post 30..60 min'};
    for k = 1:numel(wn)
        a = pv1.razon_vs_basal(strcmp(pv1.ventana, wn{k}));
        b = PSa.mediana_razon(strcmp(PSa.ventana, wn{k}));
        c = PSl.mediana_razon(strcmp(PSl.ventana, wn{k}));
        rows(end+1, :) = {['peri-ictal ' wn{k} ': razon vs basal (mediana por raton)'], a, b, sprintf(['v2 exige >= %g %% de la ' ...
            'ventana observada (v1: 50 %%); solo crisis lideres: %.2f'], 100 * cfg.periictal_min_coverage, c)}; %#ok<AGROW>
    end
    rows(end+1, :) = {'raton-dias en el violin por dia', NaN, sum(DW.enters_violin & strcmp(DW.group, 'epileptic')), ...
        sprintf('nuevo en v2 (T_seizure_obs >= %g h)', cfg.violin_min_h)};
    T = cell2table(rows, 'VariableNames', {'metrica', 'v1', cfg.run_tag, 'causa_del_cambio'});   % 'v2' / 'final' (v2 methods)
end
