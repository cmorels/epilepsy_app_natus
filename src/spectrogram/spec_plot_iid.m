function fig = spec_plot_iid(res)
% SPEC_PLOT_IID  Figure of one IID interval result (spec_iid_interval.m),
% drawn ONLY from res (S = load(mat); spec_plot_iid(S.res)).
%
% Title: animal, session, 'IID interval' (or 'IID zone'). Warnings and the
% notch band are not drawn (they stay in res / log / CSV).
%   1. 1-100 Hz trace of the interval (uV), spikes used (o)
%   2. continuous morlet scalogram of the interval, dB (zone mode: spike rate)
%   3. spike-aligned average +- SEM over the individual spikes (light)
%   4. induced TF (mean of each spike's power)      } x linked with 3,
%   5. evoked TF (TF of the average spike)          } t = 0 at the
%   6. inter-trial phase coherence (if cfg itc)      } re-aligned peak
% A sharp transient gives broadband power without any oscillation: read
% 4-5 against the trace in 3.

    s = res.cfg;
    vis = 'off';
    if s.visible
        vis = 'on';
    end
    fig = figure('Visible', vis, 'Color', 'w', 'Units', 'pixels', 'Position', [50 50 1300 1150]);
    tl = tiledlayout(fig, 4, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    t0 = res.on;
    pkk = res.peaks;
    % y axis of the four TF maps (continuous, induced, evoked, ITC): log or linear,
    % labels iid_freq_ticks (10 30 60 90 120) + the map's maximum (150)
    logy = ~isfield(s, 'iid_freq_scale') || ~strcmp(s.iid_freq_scale, 'linear');

    zone = isfield(res, 'mode') && strcmp(res.mode, 'zone');
    ax1 = nexttile(tl, 1, [1 2]);
    hold(ax1, 'on');
    if zone
        % hours: minutes on x, min/max overview of the trace, spike rate below
        tu = 60;
        xunit = 'min';
    else
        tu = 1;
        xunit = 's';
    end
    plot(ax1, (res.t_trace - t0) / tu, res.trace, 'k', 'LineWidth', 0.4);
    if ~zone
        yv = interp1(res.t_trace - t0, res.trace, pkk.t_aligned - t0, 'nearest', NaN);
        % only the spikes used in the averages / maps are marked
        plot(ax1, pkk.t_aligned(pkk.used) - t0, yv(pkk.used), 'o', 'Color', [0 0.45 0.74], 'MarkerSize', 4);
    end
    for k = 1:height(pkk.zones)
        xregion_safe(ax1, (pkk.zones.start_s(k) - t0) / tu, (pkk.zones.end_s(k) - t0) / tu);
    end
    xlim(ax1, [0 (res.off - t0) / tu]);
    ylabel(ax1, 'µV');

    ax2 = nexttile(tl, 3, [1 2]);
    if zone
        hold(ax2, 'on');
        stairs(ax2, (res.rate.t - res.rate.bin_s / 2 - t0) / tu, res.rate.detected, 'Color', [0.5 0.5 0.5], ...
            'DisplayName', 'detected');
        stairs(ax2, (res.rate.t - res.rate.bin_s / 2 - t0) / tu, res.rate.used, 'Color', [0 0.45 0.74], ...
            'LineWidth', 1.2, 'DisplayName', 'used');
        ylabel(ax2, 'spikes / min');
        legend(ax2, 'Location', 'northeast', 'Box', 'off');
        grid(ax2, 'on');
        title(ax2, 'spike rate', 'FontSize', 9, 'FontWeight', 'normal');
    else
        tf_image(ax2, res.cont.t - t0, res.f, res.cont.Z, res.clim_cont, s, logy);
        ylabel(ax2, 'Hz');
        cb = colorbar(ax2);
        cb.Label.String = 'dB';
        title(ax2, 'continuous morlet scalogram', 'FontSize', 9, 'FontWeight', 'normal');
    end
    xlim(ax2, [0 (res.off - t0) / tu]);
    xlabel(ax2, sprintf('time (%s)', xunit));
    linkaxes([ax1 ax2], 'x');

    te = res.epochs.t * 1000;
    nu = size(res.epochs.traces, 1);
    n_all = nu;
    if isfield(res.epochs, 'n')
        n_all = res.epochs.n;
    end
    ax3 = nexttile(tl, 5);
    hold(ax3, 'on');
    if nu > 0
        plot(ax3, te, res.epochs.traces', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.3);
        m = res.epochs.mean;
        e = res.epochs.sem;
        fill(ax3, [te fliplr(te)], [m + e fliplr(m - e)], [0.3 0.5 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.35);
        plot(ax3, te, m, 'b', 'LineWidth', 1.5);
    end
    xline(ax3, 0, 'k:');
    xlim(ax3, [te(1) te(end)]);
    ylabel(ax3, 'µV');
    title(ax3, sprintf('aligned average ± SEM (n = %d)', n_all), 'FontSize', 9, 'FontWeight', 'normal');

    ax4 = nexttile(tl, 6);
    tf_image(ax4, te, res.f, res.induced.Z, map_clim(res.induced.Z, s), s, logy);
    xline(ax4, 0, 'w:');
    cb = colorbar(ax4);
    cb.Label.String = 'dB';
    title(ax4, 'induced (mean of single-spike power)', 'FontSize', 9, 'FontWeight', 'normal');

    ax5 = nexttile(tl, 7);
    tf_image(ax5, te, res.f, res.evoked.Z, map_clim(res.evoked.Z, s), s, logy);
    xline(ax5, 0, 'w:');
    xlabel(ax5, 'time relative to peak (ms)');
    ylabel(ax5, 'Hz');
    cb = colorbar(ax5);
    cb.Label.String = 'dB';
    title(ax5, 'evoked', 'FontSize', 9, 'FontWeight', 'normal');

    ax6 = nexttile(tl, 8);
    if s.iid.itc
        tf_image(ax6, te, res.f, res.itc, [0 1], s, logy);
        xline(ax6, 0, 'w:');
        cb = colorbar(ax6);
        cb.Label.String = 'ITC';
        title(ax6, 'inter-trial phase coherence', 'FontSize', 9, 'FontWeight', 'normal');
    else
        axis(ax6, 'off');
    end
    xlabel(ax6, 'time relative to peak (ms)');
    linkaxes([ax3 ax4 ax5 ax6], 'x');

    % general title: animal, session, kind only. Warnings are not drawn: they stay in
    % res.warnings, the run log and the feature CSV.
    if zone
        kind_txt = 'IID zone';
    else
        kind_txt = 'IID interval';
    end
    title(tl, sprintf('%s   session %s   %s', res.subject, res.session_key, kind_txt), 'FontSize', 10, 'Interpreter', 'none');
end

%% ======================================================================
function tf_image(ax, t, f, Z, cl, s, logy)
    h = pcolor(ax, t, f, Z);
    h.EdgeColor = 'none';
    set(ax, 'Color', [0.8 0.8 0.8], 'Layer', 'top');
    colormap(ax, s.cmap);
    set(ax, 'CLim', cl);
    ylim(ax, [f(1) f(end)]);
    ticks = [10 30 60 90 120];                    % default (also for a .mat saved before iid_freq_ticks)
    if isfield(s, 'iid_freq_ticks') && ~isempty(s.iid_freq_ticks)
        ticks = s.iid_freq_ticks;
    end
    if logy
        spec_set_freq_ticks(ax, s, [f(1) f(end)], 'log', ticks);
    else
        spec_set_freq_ticks(ax, s, [f(1) f(end)], 'linear', ticks);
    end
end

function cl = map_clim(Z, s)
% Colour limits of ONE map: fixed cfg.spec.clim, or its own percentiles
% (induced and evoked differ by tens of dB: a shared scale flattens one).
    if isnumeric(s.clim) && numel(s.clim) == 2
        cl = s.clim(:)';
        return;
    end
    v = Z(isfinite(Z));
    if isempty(v)
        cl = [0 1];
        return;
    end
    cl = prctile(v, s.clim_prctile);
    if cl(2) <= cl(1)
        cl(2) = cl(1) + 1;
    end
end

function xregion_safe(ax, a, b)
    yl = ylim(ax);
    patch(ax, [a b b a], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

