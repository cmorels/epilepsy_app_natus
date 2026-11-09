function fig = spec_plot_iid(res)
% SPEC_PLOT_IID  Figure of one IID interval result (spec_iid_interval.m),
% drawn ONLY from res (S = load(mat); spec_plot_iid(S.res)).
%
%   1. 1-100 Hz trace of the interval: used peaks (o), excluded (x)
%   2. continuous Morlet scalogram of the interval (x linked with 1)
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
    logy = true;

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
        plot(ax1, pkk.t_aligned(pkk.used) - t0, yv(pkk.used), 'o', 'Color', [0 0.45 0.74], 'MarkerSize', 4);
        plot(ax1, pkk.t_aligned(~pkk.used) - t0, yv(~pkk.used), 'x', 'Color', [0.85 0.33 0.1], 'MarkerSize', 5);
        if isfield(pkk, 'at_limit') && any(pkk.at_limit)
            plot(ax1, pkk.t_aligned(pkk.at_limit) - t0, yv(pkk.at_limit), 's', 'Color', [0.8 0 0.8], 'MarkerSize', 7);
        end
    end
    for k = 1:height(pkk.zones)
        xregion_safe(ax1, (pkk.zones.start_s(k) - t0) / tu, (pkk.zones.end_s(k) - t0) / tu);
    end
    xlim(ax1, [0 (res.off - t0) / tu]);
    ylabel(ax1, 'uV (1-100 Hz)');
    if zone
        title(ax1, sprintf('trace overview (min/max per bin; shaded = exclusion zone)   peaks: %s', pkk.source), ...
            'FontSize', 8, 'FontWeight', 'normal', 'Interpreter', 'none');
    else
        title(ax1, sprintf('trace  (o used, x excluded, square = re-aligned at the +-%g ms limit; shaded = exclusion zone)   peaks: %s', ...
            s.iid.align_ms, pkk.source), 'FontSize', 8, 'FontWeight', 'normal', 'Interpreter', 'none');
    end

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
        title(ax2, sprintf('spike rate (%g s bins)', res.rate.bin_s(1)), 'FontSize', 8, 'FontWeight', 'normal');
    else
        tf_image(ax2, res.cont.t - t0, res.f, res.cont.Z, res.clim_cont, s.cmap, logy);
        ylabel(ax2, 'Hz');
        cb = colorbar(ax2);
        cb.Label.String = res.cont.label;
        title(ax2, 'continuous Morlet scalogram; grey = gap / cone of influence', 'FontSize', 8, 'FontWeight', 'normal');
    end
    xlim(ax2, [0 (res.off - t0) / tu]);
    xlabel(ax2, sprintf('%s from %s (Natus)', xunit, res.on_natus));
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
    ylabel(ax3, 'uV (1-100 Hz)');
    f = res.features;
    title(ax3, sprintf('aligned average +- SEM (n = %d of %d, %d drawn; polarity %s [%s])  p2p %.0f uV, FWHM %.1f ms', ...
        n_all, f.n_detected, nu, pkk.polarity, pkk.polarity_mode, f.p2p_mean_uV, f.fwhm_mean_ms), ...
        'FontSize', 8, 'FontWeight', 'normal');

    ax4 = nexttile(tl, 6);
    tf_image(ax4, te, res.f, res.induced.Z, res.clim_tf, s.cmap, logy);
    xline(ax4, 0, 'w:');
    cb = colorbar(ax4);
    cb.Label.String = res.norm_label;
    title(ax4, sprintf('induced (mean of single-spike power); max within +-%g ms at %.1f Hz', ...
        s.iid.maxfreq_window_ms, f.induced_peak_freq_hz), 'FontSize', 8, 'FontWeight', 'normal');

    ax5 = nexttile(tl, 7);
    tf_image(ax5, te, res.f, res.evoked.Z, res.clim_tf, s.cmap, logy);
    xline(ax5, 0, 'w:');
    xlabel(ax5, 'ms re peak');
    ylabel(ax5, 'Hz');
    cb = colorbar(ax5);
    cb.Label.String = res.norm_label;
    title(ax5, 'evoked (TF of the average spike)', 'FontSize', 8, 'FontWeight', 'normal');

    ax6 = nexttile(tl, 8);
    if s.iid.itc
        tf_image(ax6, te, res.f, res.itc, [0 1], s.cmap, logy);
        xline(ax6, 0, 'w:');
        cb = colorbar(ax6);
        cb.Label.String = 'ITC';
        title(ax6, 'inter-trial phase coherence', 'FontSize', 8, 'FontWeight', 'normal');
    else
        axis(ax6, 'off');
    end
    xlabel(ax6, 'ms re peak');
    linkaxes([ax3 ax4 ax5 ax6], 'x');

    src = res.src;
    if zone
        kind_txt = sprintf('IID zone, %.2f h', (res.off - res.on) / 3600);
    else
        kind_txt = sprintf('IID interval, %.1f s', res.off - res.on);
    end
    head = {sprintf('%s   %s   session %s   id %s   (%s)', res.subject, res.region, ...
        res.session_key, res.id, kind_txt), ...
        sprintf('%s -> %s (Natus) | %s (pipeline) | t_rel %.3f-%.3f s      signal %s, case %s, gain %g, notch %s', ...
        res.on_natus, res.off_natus, res.on_pipeline, res.on, res.off, src.signal_source, blank(src.case_applied), ...
        src.gain_applied, tf_text(src.notch_applied)), ...
        sprintf('spikes: %d detected, %d used (excluded: nan %d, zone %d, edge %d, duplicate %d), %.1f / min; %d re-aligned at the +-%g ms limit', ...
        f.n_detected, f.n_used, f.n_excl_nan, f.n_excl_zone, f.n_excl_edge, f.n_excl_duplicate, f.rate_per_min, ...
        f.n_align_at_limit, s.iid.align_ms)};
    title(tl, head, 'FontSize', 10, 'Interpreter', 'none');
    if ~isempty(res.warnings)
        tl.OuterPosition = [0 0.04 1 0.96];
        annotation(fig, 'textbox', [0.005 0 0.99 0.037], 'String', ['WARNING  ' strjoin(res.warnings, '  |  ')], ...
            'EdgeColor', 'none', 'FontSize', 7, 'Color', [0.7 0 0], 'Interpreter', 'none', ...
            'VerticalAlignment', 'bottom');
    end
end

%% ======================================================================
function tf_image(ax, t, f, Z, cl, cmap, logy)
    h = pcolor(ax, t, f, Z);
    h.EdgeColor = 'none';
    set(ax, 'Color', [0.8 0.8 0.8], 'Layer', 'top');
    colormap(ax, cmap);
    set(ax, 'CLim', cl);
    ylim(ax, [f(1) f(end)]);
    if logy
        set(ax, 'YScale', 'log');
    end
end

function xregion_safe(ax, a, b)
    yl = ylim(ax);
    patch(ax, [a b b a], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

function s = blank(v)
    if isempty(v)
        s = '-';
    else
        s = char(v);
    end
end

function s = tf_text(b)
    if b
        s = 'yes';
    else
        s = 'no';
    end
end
