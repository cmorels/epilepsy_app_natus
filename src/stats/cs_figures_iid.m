function files = cs_figures_iid(cfg, I, CR, fig_dir)
% CS_FIGURES_IID  Section 7 IID figures (PNG, 300 dpi):
%   iid_rate_by_mouse.png          complexes/min per mouse (merged, both
%                                  hemispheres), control vs epileptic, log
%                                  axis, original (open) vs re-detection (filled)
%   iid_baseline_before_after.png  baseline / lower threshold per channel,
%                                  original vs re-detection, and the real P97
%                                  of the included samples vs the baseline grid
%   iid_vs_seizures.png            between animals: IID rate vs seizure rate
%                                  (epileptic, log-log, mouse labels)
%   iid_periictal.png              per-mouse rate in each peri-ictal window
%                                  relative to its interictal baseline
%   iid_circadian.png              rate per clock hour normalized by analyzed
%                                  channel-time, per group, dark phase shaded

    if ~isfolder(fig_dir), mkdir(fig_dir); end
    files = {};
    col = struct('epileptic', [0.80 0.25 0.20], 'control', [0.20 0.40 0.75]);

    % ---- rate by mouse ---------------------------------------------------
    mm = I.mouse(strcmp(I.mouse.scope, 'merged') & strcmp(I.mouse.hemisphere, 'ambos'), :);
    mm = sortrows(mm, {'group', 'complexes_per_min'}, {'ascend', 'ascend'});
    fig = figure('Visible', 'off', 'Position', [100 100 1000 420]);
    ax = axes(fig); hold(ax, 'on');
    for g = {'control', 'epileptic'}
        s = find(strcmp(mm.group, g{1}));
        plot(ax, s, mm.complexes_per_min_orig(s), 'o', 'Color', col.(g{1}), 'MarkerSize', 7, 'LineWidth', 1.2, ...
            'DisplayName', [g{1} ' original']);
        plot(ax, s, mm.complexes_per_min(s), 'o', 'Color', col.(g{1}), 'MarkerFaceColor', col.(g{1}), 'MarkerSize', 6, ...
            'DisplayName', [g{1} ' re-deteccion']);
        for k = s'
            plot(ax, [k k], [mm.complexes_per_min_orig(k) mm.complexes_per_min(k)], '-', 'Color', [col.(g{1}) 0.5], 'HandleVisibility', 'off');
        end
    end
    v = [mm.complexes_per_min; mm.complexes_per_min_orig];
    v = v(v > 0);
    set(ax, 'YScale', 'log', 'XTick', 1:height(mm), 'XTickLabel', mm.mouse, 'XTickLabelRotation', 60, ...
        'XLim', [0.5 height(mm) + 0.5], 'YLim', [min(v) / 1.6, max(v) * 1.6]);
    ylabel(ax, 'complejos / min (log)'); grid(ax, 'on'); legend(ax, 'Location', 'northwest');
    title(ax, 'Tasa de IIDs por raton (C1+C2, ambos hemisferios)');
    files{end+1} = save_fig(fig, fig_dir, 'iid_rate_by_mouse.png');

    % ---- baseline / threshold before vs after ------------------------------
    T = I.channels;
    fig = figure('Visible', 'off', 'Position', [100 100 1100 420]);
    ax = subplot(1, 2, 1, 'Parent', fig); hold(ax, 'on');
    rng(cfg.rng_seed);
    j = @(v) v + 0.8 * (rand(size(v)) - 0.5);
    plot(ax, j(T.baseline_uV_orig), j(T.baseline_uV), '.', 'MarkerSize', 10, 'Color', [0.3 0.3 0.3]);
    lim = [min([T.baseline_uV_orig; T.baseline_uV]) - 5, max([T.baseline_uV_orig; T.baseline_uV]) + 5];
    plot(ax, lim, lim, 'k--');
    xlim(ax, lim); ylim(ax, lim); axis(ax, 'square'); grid(ax, 'on');
    xlabel(ax, 'baseline original (\muV)'); ylabel(ax, 'baseline re-deteccion (\muV)');
    title(ax, sprintf('Baseline por canal (n = %d; umbral = 2.5 x baseline)', height(T)));
    ax = subplot(1, 2, 2, 'Parent', fig); hold(ax, 'on');
    histogram(ax, T.p97_real_uV, 30, 'FaceColor', [0.5 0.5 0.6]);
    xline(ax, 130, 'r-', 'tope de la grilla (130 \muV)', 'LineWidth', 1.2);
    xlabel(ax, 'P97 real de |senal| en muestras incluidas (\muV)'); ylabel(ax, 'canales'); grid(ax, 'on');
    title(ax, 'P97 real vs grilla 100:5:130 \muV');
    files{end+1} = save_fig(fig, fig_dir, 'iid_baseline_before_after.png');

    % ---- IIDs vs seizures, between animals --------------------------------
    e = mm(strcmp(mm.group, 'epileptic'), :);
    if ismember('seizure_rate_per_24h', e.Properties.VariableNames)
        fig = figure('Visible', 'off', 'Position', [100 100 600 520]);
        ax = axes(fig); hold(ax, 'on');
        plot(ax, e.complexes_per_min, e.seizure_rate_per_24h, 'o', 'MarkerFaceColor', col.epileptic, 'Color', col.epileptic);
        text(ax, e.complexes_per_min, e.seizure_rate_per_24h, strcat({'  '}, e.mouse), 'FontSize', 8);
        xr = e.complexes_per_min(e.complexes_per_min > 0);
        yr = e.seizure_rate_per_24h(e.seizure_rate_per_24h > 0);
        set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', [min(xr) / 1.5, max(xr) * 2], 'YLim', [min(yr) / 1.6, max(yr) * 1.6], ...
            'XTick', [0.5 1 2 5 10 20 50], 'YTick', [0.1 0.2 0.5 1 2 5 10 20]);
        grid(ax, 'on');
        xlabel(ax, 'IIDs: complejos / min (log)'); ylabel(ax, 'crisis / 24 h (log)');
        r = CR(strcmp(CR.analisis, 'entre animales') & strcmp(CR.datos, 're-deteccion') & strcmp(CR.conjunto, 'todos'), :);
        lbl = strjoin(arrayfun(@(k) sprintf('%s = %.2f [%.2f, %.2f], p = %.3g', r.metodo{k}, r.estimacion(k), r.ci95_lo(k), ...
            r.ci95_hi(k), r.p(k)), (1:height(r))', 'UniformOutput', false), newline);
        title(ax, sprintf('Entre animales (n = %d epilepticos)', height(e)));
        text(ax, 0.02, 0.98, lbl, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8);
        files{end+1} = save_fig(fig, fig_dir, 'iid_vs_seizures.png');
    end

    % ---- peri-ictal --------------------------------------------------------
    P = I.periictal(~startsWith(I.periictal.mouse, '('), :);
    wins = unique(P.ventana, 'stable');
    ms = unique(P.mouse);
    Y = NaN(numel(ms), numel(wins));
    for k = 1:numel(ms)
        for w = 1:numel(wins)
            r = P(strcmp(P.mouse, ms{k}) & strcmp(P.ventana, wins{w}), :);
            if ~isempty(r), Y(k, w) = r.razon_vs_basal(1); end
        end
    end
    fig = figure('Visible', 'off', 'Position', [100 100 700 450]);
    ax = axes(fig); hold(ax, 'on');
    Yp = Y; Yp(Yp <= 0) = NaN;
    plot(ax, 1:numel(wins), Yp', '-', 'Color', [0.7 0.7 0.7]);
    plot(ax, 1:numel(wins), median(Y, 1, 'omitnan'), 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
    yline(ax, 1, 'r--');
    yl = [min(Yp(:), [], 'omitnan') / 1.3, max(Yp(:), [], 'omitnan') * 1.3];
    set(ax, 'YScale', 'log', 'XTick', 1:numel(wins), 'XTickLabel', wins, 'XLim', [0.5 numel(wins) + 0.5], ...
        'YLim', yl, 'YTick', [0.125 0.25 0.5 1 2 4 8 16 32 64]);
    ylabel(ax, 'tasa en ventana / tasa basal interictal (log)'); grid(ax, 'on');
    title(ax, sprintf('Perfil peri-ictal de IIDs (gris: raton; negro: mediana; n = %d)', numel(ms)));
    files{end+1} = save_fig(fig, fig_dir, 'iid_periictal.png');

    % ---- circadian ---------------------------------------------------------
    on = hm(cfg.lights_on) / 60; off = hm(cfg.lights_off) / 60;
    H = I.circ_hourly;
    fig = figure('Visible', 'off', 'Position', [100 100 900 600]);
    gs = {'epileptic', 'control'};
    for gi = 1:2
        ax = subplot(2, 1, gi, 'Parent', fig); hold(ax, 'on');
        x = H(strcmp(H.group, gs{gi}), :);
        ymax = max(x.rate_per_min) * 1.15;
        if ~(ymax > 0), ymax = 1; end
        patch(ax, [off 24 24 off], [0 0 ymax ymax], [0.88 0.88 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
        patch(ax, [0 on on 0], [0 0 ymax ymax], [0.88 0.88 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
        bar(ax, x.hour + 0.5, x.rate_per_min, 1, 'FaceColor', col.(gs{gi}), 'EdgeColor', 'none');
        xlim(ax, [0 24]); ylim(ax, [0 ymax]); set(ax, 'XTick', 0:2:24, 'Layer', 'top');
        ylabel(ax, 'complejos / min'); title(ax, sprintf('%s: IIDs por hora normalizados por tiempo analizado', gs{gi}));
    end
    xlabel(ax, 'hora del dia (sombreado: oscuridad)');
    files{end+1} = save_fig(fig, fig_dir, 'iid_circadian.png');
end

function f = save_fig(fig, d, name)
    f = fullfile(d, name);
    exportgraphics(fig, f, 'Resolution', 300);
    close(fig);
end

function m = hm(s)
    v = sscanf(s, '%d:%d');
    m = v(1) * 60 + v(2);
end
