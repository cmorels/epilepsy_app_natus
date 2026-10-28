function plot_ll_band(ll_accept, ll_reject, ax)
% PLOT_LL_BAND  Draw the robust branch's review band on a line-length axis:
% two lines (ll_accept, ll_reject) with the in-between band shaded.
    if nargin < 3
        ax = gca;
    end
    if ll_accept > ll_reject
        yregion(ax, ll_reject, ll_accept, 'FaceColor', [1 0.85 0.55], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
    end
    yline(ax, ll_accept, 'r--', sprintf('ll\\_accept (%.2f)', ll_accept), 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    yline(ax, ll_reject, '--', sprintf('ll\\_reject (%.2f)', ll_reject), 'Color', [0.9 0.45 0], 'LineWidth', 1.5, ...
        'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');
end
