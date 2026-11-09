function f = cs_figure_sensitivity(cfg, SG, SM, fig_dir)
% CS_FIGURE_SENSITIVITY  iid_threshold_sensitivity.png (v2 section 8):
% (top) fold change epileptic/control of the medians with its bootstrap
% 95 % CI, per threshold variant and metric, log axis, line at 10x;
% (bottom) complexes/min per mouse across variants, control vs epileptic.

    if ~isfolder(fig_dir), mkdir(fig_dir); end
    vn = unique(SG.variante, 'stable');
    metrics = {'complexes_per_min', 'single_per_min', 'polyspikes_per_min', 'bursts_per_h'};
    mlab = {'complejos/min', 'singles/min', 'polyspikes/min', 'bursts/h'};
    mcol = lines(numel(metrics));
    fig = figure('Visible', 'off', 'Position', [100 100 1200 820], 'Color', 'w');
    ax = subplot(2, 1, 1, 'Parent', fig); hold(ax, 'on');
    off = linspace(-0.27, 0.27, numel(metrics));
    for q = 1:numel(metrics)
        x = SG(strcmp(SG.metrica, metrics{q}), :);
        [~, j] = ismember(vn, x.variante);
        xx = (1:numel(vn))' + off(q);
        errorbar(ax, xx, x.fold_change_medianas(j), x.fold_change_medianas(j) - x.fc_ci95_lo(j), ...
            x.fc_ci95_hi(j) - x.fold_change_medianas(j), 'o', 'Color', mcol(q, :), 'MarkerFaceColor', mcol(q, :), ...
            'LineWidth', 1.1, 'DisplayName', mlab{q});
    end
    yline(ax, 10, 'r--', '10x', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(ax, 1, 'k:', 'HandleVisibility', 'off');
    set(ax, 'YScale', 'log', 'XTick', 1:numel(vn), 'XTickLabel', vn, 'TickLabelInterpreter', 'none', 'XLim', [0.5 numel(vn) + 0.5], 'FontSize', 9);
    ylabel(ax, 'fold change epileptico/control (mediana, IC95 bootstrap)'); grid(ax, 'on');
    legend(ax, 'Location', 'eastoutside');
    title(ax, 'Sensibilidad del umbral IID (descriptiva; prominencia 200 \muV fija en todas las variantes)');
    subtitle(ax, ['Sin punto o sin limite superior: la mediana control (o la de un remuestreo) es 0 -> cociente infinito ' ...
        '(ver iid_threshold_summary)'], 'Interpreter', 'none', 'FontSize', 8);

    ax = subplot(2, 1, 2, 'Parent', fig); hold(ax, 'on');
    col = struct('epileptic', [0.80 0.25 0.20], 'control', [0.20 0.40 0.75]);
    mice = unique(SM.mouse);
    for i = 1:numel(mice)
        x = SM(strcmp(SM.mouse, mice{i}), :);
        [~, j] = ismember(vn, x.variante);
        y = x.complexes_per_min(j);
        y(y <= 0) = NaN;
        plot(ax, 1:numel(vn), y, '-o', 'Color', [col.(x.group{1}) 0.6], 'MarkerFaceColor', col.(x.group{1}), ...
            'MarkerSize', 4, 'HandleVisibility', 'off');
    end
    plot(ax, NaN, NaN, '-o', 'Color', col.epileptic, 'MarkerFaceColor', col.epileptic, 'DisplayName', 'epileptico');
    plot(ax, NaN, NaN, '-o', 'Color', col.control, 'MarkerFaceColor', col.control, 'DisplayName', 'control');
    set(ax, 'YScale', 'log', 'XTick', 1:numel(vn), 'XTickLabel', vn, 'TickLabelInterpreter', 'none', 'XLim', [0.5 numel(vn) + 0.5], 'FontSize', 9);
    ylabel(ax, 'complejos / min por raton (log)'); grid(ax, 'on');
    legend(ax, 'Location', 'eastoutside');
    f = fullfile(fig_dir, 'iid_threshold_sensitivity.png');
    exportgraphics(fig, f, 'Resolution', 300);
    close(fig);
end
