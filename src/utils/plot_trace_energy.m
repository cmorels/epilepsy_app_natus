function r = plot_trace_energy(ax, t, energy, win, threshold)
% PLOT_TRACE_ENERGY  Row 3 of every event figure: Hilbert-envelope energy
% (seizure_energy_trace.m's energy_full) with its global threshold. Same
% trace as detect_seizures_robust.m's zoom panel 3.
    hold(ax, 'on');
    [a, b] = plot_trace_segment(ax, t, energy, win, [0.2 0.6 0.3]);
    r = [min(a, threshold), max(b, threshold)];
    yline(ax, threshold, 'r--', sprintf('Threshold (%.2e)', threshold), 'LineWidth', 1.5);
    ylabel(ax, 'Energy (a.u.)'); grid(ax, 'on');
end
