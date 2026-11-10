function fig = spec_plot_seizure(res)
% SPEC_PLOT_SEIZURE  Figure of one seizure result (spec_seizure_event.m),
% drawn ONLY from res -- so  S = load(mat); spec_plot_seizure(S.res)
% regenerates it without reading any txt.
%
% Title: subject, session, id. One column per region (HPCr, HPCl with
% region 'both'), x axes linked:
%   1. region: 1-100 Hz trace (uV), onset / offset of that region, baseline shaded
%   2. 'spectrogram': multitaper (dB, norm of res.cfg.norm)
%   3. 'morlet': Morlet scalogram (log frequency axis, dB); grey = NaN: gap or
%      cone of influence
%   4. 'PSD': baseline vs ictal (vs post-ictal), dB re 1 uV^2/Hz
% t = 0 at the event onset. Warnings are not drawn (res.ch(k).warnings, run
% log and feature CSV); neither is the notch band (excluded from the
% frequency features all the same).

    s = res.cfg;
    nreg = numel(res.ch);
    vis = 'off';
    if s.visible
        vis = 'on';
    end
    fig = figure('Visible', vis, 'Color', 'w', 'Units', 'pixels', 'Position', [50 50 720 * nreg 1050]);
    tl = tiledlayout(fig, 4, nreg, 'TileSpacing', 'compact', 'Padding', 'compact');
    xl = res.window;
    axs = gobjects(0);
    for k = 1:nreg
        c = res.ch(k);
        on_k = c.on - res.on;
        off_k = c.off - res.on;
        bl = c.baseline - res.on;

        ax1 = nexttile(tl, k);
        hold(ax1, 'on');
        yl = robust_ylim(c.trace);
        patch(ax1, [bl(1) bl(2) bl(2) bl(1)], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.92 1], 'EdgeColor', 'none');
        plot(ax1, c.t_trace, c.trace, 'k', 'LineWidth', 0.5);
        xline(ax1, on_k, 'r', 'LineWidth', 1.2);
        xline(ax1, off_k, 'b', 'LineWidth', 1.2);
        ylim(ax1, yl);
        xlim(ax1, xl);
        ylabel(ax1, 'µV');
        title(ax1, c.region, 'FontSize', 9, 'FontWeight', 'normal', 'Interpreter', 'none');

        ax2 = nexttile(tl, nreg + k);
        tf_image(ax2, c.mt.t, c.mt.f, c.mt.Z, res.clim_mt, s, strcmp(s.freq_scale, 'log'));
        mark(ax2, on_k, off_k);
        xlim(ax2, xl);
        ylabel(ax2, 'Hz');
        title(ax2, 'spectrogram', 'FontSize', 9, 'FontWeight', 'normal');
        cb = colorbar(ax2);
        cb.Label.String = 'dB';

        ax3 = nexttile(tl, 2 * nreg + k);
        tf_image(ax3, c.morlet.t, c.morlet.f, c.morlet.Z, res.clim_morlet, s, true);
        mark(ax3, on_k, off_k);
        xlim(ax3, xl);
        ylabel(ax3, 'Hz');
        xlabel(ax3, 'Time (s)');
        title(ax3, 'morlet', 'FontSize', 9, 'FontWeight', 'normal');
        cb = colorbar(ax3);
        cb.Label.String = 'dB';

        ax4 = nexttile(tl, 3 * nreg + k);
        hold(ax4, 'on');
        plot(ax4, c.psd.f, c.psd.baseline_db, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, 'DisplayName', 'baseline');
        plot(ax4, c.psd.f, c.psd.ictal_db, 'r', 'LineWidth', 1.2, 'DisplayName', 'ictal');
        if any(~isnan(c.psd.postictal_db))
            plot(ax4, c.psd.f, c.psd.postictal_db, 'b', 'LineWidth', 1, 'DisplayName', 'post-ictal');
        end
        grid(ax4, 'on');
        xlim(ax4, [c.psd.f(1) c.psd.f(end)]);
        if strcmp(s.freq_scale, 'log')
            set(ax4, 'XScale', 'log');
        end
        xlabel(ax4, 'Hz');
        ylabel(ax4, 'dB re 1 uV^2/Hz');
        legend(ax4, 'Location', 'northeast', 'Box', 'off');
        title(ax4, 'PSD', 'FontSize', 9, 'FontWeight', 'normal');

        axs = [axs, ax1, ax2, ax3]; %#ok<AGROW>
    end
    linkaxes(axs, 'x');

    % warnings are not drawn: they stay in res.ch(k).warnings, the run log and the feature CSV
    title(tl, sprintf('%s   session %s   %s', res.subject, res.session_key, res.id), 'FontSize', 10, 'Interpreter', 'none');
end

%% ======================================================================
function tf_image(ax, t, f, Z, cl, s, logy)
    h = pcolor(ax, t, f, Z);
    h.EdgeColor = 'none';
    set(ax, 'Color', [0.8 0.8 0.8], 'Layer', 'top');
    colormap(ax, s.cmap);
    set(ax, 'CLim', cl);
    axis(ax, 'xy');
    ylim(ax, [f(1) f(end)]);
    if logy
        spec_set_freq_ticks(ax, s, [f(1) f(end)]);   % Hz labels, not 10^k
    end
end

function mark(ax, on_k, off_k)
    xline(ax, on_k, 'w--', 'LineWidth', 1);
    xline(ax, off_k, 'w--', 'LineWidth', 1);
end

function yl = robust_ylim(x)
    v = x(isfinite(x));
    if isempty(v)
        yl = [-1 1];
        return;
    end
    m = max(abs(v));
    yl = [-m m] * 1.05;
end
