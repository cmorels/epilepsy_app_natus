function files = cs_iid_vs_seizures_v2(out_dir)
% CS_IID_VS_SEIZURES_V2  Reformatted version of the between-animals figure
% of cs_figures_iid.m (iid_vs_seizures.png), from a run's
% phase3_results.mat: same data (I.mouse, merged, both hemispheres,
% epileptic), text in English, no statistics text, fixed log-log limits
% and major grid only. Writes figures/iid_vs_seizures_v2.png (300 dpi) and
% .fig; the original figure and cs_figures_iid.m are not touched.
%
%   cs_iid_vs_seizures_v2('<out_dir>')

    q = load(fullfile(out_dir, 'phase3_results.mat'), 'I');
    mm = q.I.mouse(strcmp(q.I.mouse.scope, 'merged') & strcmp(q.I.mouse.hemisphere, 'ambos'), :);
    mm = sortrows(mm, {'group', 'complexes_per_min'}, {'ascend', 'ascend'});
    e = mm(strcmp(mm.group, 'epileptic'), :);
    x = e.complexes_per_min;
    y = e.seizure_rate_per_24h;
    ids = e.mouse;

    xl = [0.9 25];
    yl = [0.1 12];
    red = [0.80 0.22 0.17];
    fig = figure('Visible', 'off', 'Position', [100 100 600 520], 'Color', 'w');
    ax = axes(fig, 'Color', 'w', 'FontSize', 11); hold(ax, 'on');
    scatter(ax, x, y, 70, red, 'filled', 'MarkerEdgeColor', 'w');
    h = text(ax, x * 1.06, y, ids, 'FontSize', 8);
    % 072-s / 058-s labels overlap otherwise; 058-s goes to the left of its
    % point (below it, it would cover 001-s)
    set(h(strcmp(ids, '072-s')), 'VerticalAlignment', 'bottom');
    k058 = strcmp(ids, '058-s');
    set(h(k058), 'Position', [x(k058) / 1.06, y(k058), 0], 'HorizontalAlignment', 'right');
    set(ax, 'XScale', 'log', 'YScale', 'log', 'XLim', xl, 'YLim', yl, ...
        'XTick', [1 2 5 10 20], 'YTick', [0.1 0.2 0.5 1 2 5 10], 'YColor', 'k', ...
        'XMinorTick', 'off', 'YMinorTick', 'off');
    grid(ax, 'on');
    set(ax, 'XMinorGrid', 'off', 'YMinorGrid', 'off', 'GridLineStyle', '-', 'GridColor', 'k', 'GridAlpha', 0.10);
    box(ax, 'off');
    xlabel(ax, 'IIDs / min (log)'); ylabel(ax, 'Seizures / 24 h (log)');
    title(ax, 'IIDs vs seizures');
    subtitle(ax, sprintf('(n = %d)', sum(x > 0 & y > 0)));   % mice drawn (0 seizures cannot be shown on a log axis)

    fig_dir = fullfile(out_dir, 'figures');
    files = {fullfile(fig_dir, 'iid_vs_seizures_v2.png'), fullfile(fig_dir, 'iid_vs_seizures_v2.fig')};
    exportgraphics(fig, files{1}, 'Resolution', 300);
    set(fig, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');   % the .fig opens visible
    savefig(fig, files{2});
    close(fig);
end
