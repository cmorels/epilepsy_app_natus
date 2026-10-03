function [ymin, ymax] = plot_trace_segment(ax, t, y, win, color)
% PLOT_TRACE_SEGMENT  Plot y(t) over the time window win = [t0 t1] on ax,
% min-max decimated (decimate_minmax.m) above 20000 points. Returns the
% finite data range inside the window (for common y-limits across columns).
% Shared core of plot_trace_voltage / _bandpassed / _energy / _linelength.
    iz = t >= win(1) & t <= win(2);
    tz = t(iz);
    yz = y(iz);
    if numel(tz) > 20000
        [tz, yz] = decimate_minmax(tz, yz, 20000);
    end
    plot(ax, tz, yz, 'Color', color, 'LineWidth', 0.8);
    yf = yz(isfinite(yz));
    if isempty(yf)
        ymin = NaN; ymax = NaN;
    else
        ymin = min(yf); ymax = max(yf);
    end
end
