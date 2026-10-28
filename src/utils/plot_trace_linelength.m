function r = plot_trace_linelength(ax, t, ll_norm, win, cfg)
% PLOT_TRACE_LINELENGTH  Row 4 of every event figure: line length divided by
% ll_median_global, with the review band (ll_accept / ll_reject, shaded
% between) from resolve_ll_band.m. Same trace as detect_seizures_robust.m's
% zoom panel 4.
    hold(ax, 'on');
    [ll_accept, ll_reject] = resolve_ll_band(cfg);
    plot_ll_band(ll_accept, ll_reject, ax);
    [a, b] = plot_trace_segment(ax, t, ll_norm, win, [0.6 0.3 0.6]);
    r = [min(a, ll_reject), max(b, ll_accept)];
    ylabel(ax, 'll / ll\_median\_global'); xlabel(ax, 'Time (s)'); grid(ax, 'on');
end
