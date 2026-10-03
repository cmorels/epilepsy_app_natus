function r = plot_trace_bandpassed(ax, t, bp, win, cfg)
% PLOT_TRACE_BANDPASSED  Row 2 of every event figure: the min-max normalized,
% cfg.seizure.bandpass_band band-passed signal (seizure_energy_trace.m's
% bp_full). Same trace as detect_seizures_robust.m's zoom panel 2.
    hold(ax, 'on');
    [a, b] = plot_trace_segment(ax, t, bp, win, [0 0 0]);
    r = [a b];
    ylabel(ax, 'Normalized amplitude'); grid(ax, 'on');
    text(ax, 0.01, 0.95, sprintf('band-passed %g-%g Hz', cfg.seizure.bandpass_band(1), cfg.seizure.bandpass_band(2)), ...
        'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0.3 0.3 0.3]);
end
