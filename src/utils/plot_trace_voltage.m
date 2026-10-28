function r = plot_trace_voltage(ax, t, signal, win)
% PLOT_TRACE_VOLTAGE  Row 1 of every event figure: cleaned LFP in uV.
% Same trace as detect_seizures_robust.m's zoom panel 1. Returns [ymin ymax].
    hold(ax, 'on');
    [a, b] = plot_trace_segment(ax, t, signal, win, [0 0 0]);
    r = [a b];
    ylabel(ax, 'Voltage (uV)'); grid(ax, 'on');
end
