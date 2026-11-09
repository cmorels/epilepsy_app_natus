function fig = spec_plot_seizure(res)
% SPEC_PLOT_SEIZURE  Figure of one seizure result (spec_seizure_event.m),
% drawn ONLY from res -- so  S = load(mat); spec_plot_seizure(S.res)
% regenerates it without reading any txt.
%
% One column per region (HPCr, HPCl with region 'both'), x axes linked:
%   1. 1-100 Hz trace, onset / offset of that region, baseline shaded
%   2. multitaper spectrogram (norm of res.cfg.norm)
%   3. Morlet scalogram (log frequency axis); grey = NaN: gap or cone of
%      influence
%   4. PSD baseline vs ictal (vs post-ictal), dB re 1 uV^2/Hz
% t = 0 at the event onset.

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
        ylabel(ax1, 'uV (1-100 Hz)');
        title(ax1, {sprintf('%s   onset %s (Natus) | %s (pipeline)', c.region, c.on_natus, c.on_pipeline), ...
            sprintf('offset %s (Natus)   duration %.1f s   user_check: %s', c.off_natus, c.off - c.on, blank(c.sz.user_check))}, ...
            'FontSize', 8, 'FontWeight', 'normal', 'Interpreter', 'none');

        ax2 = nexttile(tl, nreg + k);
        tf_image(ax2, c.mt.t, c.mt.f, c.mt.Z, res.clim_mt, s.cmap, strcmp(s.freq_scale, 'log'));
        mark(ax2, on_k, off_k);
        xlim(ax2, xl);
        ylabel(ax2, 'Hz');
        title(ax2, sprintf('multitaper  win %.1f s, step %.2f s, TW %g, K %d (%.1f Hz res.)', ...
            s.seizure.mt.win_s, s.seizure.mt.step_s, s.seizure.mt.TW, s.seizure.mt.K, c.mt.info.resolution_hz), ...
            'FontSize', 8, 'FontWeight', 'normal');
        cb = colorbar(ax2);
        cb.Label.String = res.norm_label;

        ax3 = nexttile(tl, 2 * nreg + k);
        tf_image(ax3, c.morlet.t, c.morlet.f, c.morlet.Z, res.clim_morlet, s.cmap, true);
        mark(ax3, on_k, off_k);
        xlim(ax3, xl);
        ylabel(ax3, 'Hz');
        xlabel(ax3, 's re onset');
        title(ax3, sprintf('Morlet  %g-%g Hz, %g-%g cycles (%s); grey = gap / cone of influence', ...
            s.seizure.morlet.fmin, s.seizure.morlet.fmax, s.seizure.morlet.cycles(1), s.seizure.morlet.cycles(end), ...
            s.seizure.morlet.cycles_scale), 'FontSize', 8, 'FontWeight', 'normal');
        cb = colorbar(ax3);
        cb.Label.String = res.norm_label;

        ax4 = nexttile(tl, 3 * nreg + k);
        hold(ax4, 'on');
        if isfield(c.psd, 'excluded') && any(c.psd.excluded)
            fx = c.psd.f(c.psd.excluded);
            yl4 = [min([c.psd.baseline_db; c.psd.ictal_db]) max([c.psd.baseline_db; c.psd.ictal_db])];
            if all(isfinite(yl4))
                patch(ax4, [fx(1) fx(end) fx(end) fx(1)], yl4([1 1 2 2]), [0.9 0.9 0.9], 'EdgeColor', 'none', ...
                    'DisplayName', 'notch (excluded from peak / dominant f)');
            end
        end
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
        f = c.features;
        title(ax4, sprintf('ictal peak above 1/f %.1f Hz (absolute %.1f, re baseline %.1f); dominant f (%s) %.1f Hz, %.2f Hz/s', ...
            f.peak_freq_osc_hz, f.peak_freq_hz, f.peak_freq_rel_hz, s.seizure.domfreq_method, f.domfreq_mean_hz, ...
            f.domfreq_slope_hz_per_s), ...
            'FontSize', 8, 'FontWeight', 'normal');

        axs = [axs, ax1, ax2, ax3]; %#ok<AGROW>
    end
    linkaxes(axs, 'x');

    c1 = res.ch(1);
    ids = sprintf('id %s', res.id);
    if ~isnan(c1.sz.seizure_id)
        ids = sprintf('%s   seizure_id %d', ids, c1.sz.seizure_id);
    end
    if ~isnan(c1.sz.event_id)
        ids = sprintf('%s   event_id %d', ids, c1.sz.event_id);
    end
    src = c1.src;
    head = {sprintf('%s   %s   session %s   %s', res.subject, res.region_tag, res.session_key, ids), ...
        sprintf('onset %s (Natus) | %s (pipeline) | t_rel %.3f s      signal %s, case %s, gain %g, notch %s', ...
        res.on_natus, res.on_pipeline, res.on, src.signal_source, blank(src.case_applied), src.gain_applied, ...
        tf_text(src.notch_applied))};
    title(tl, head, 'FontSize', 10, 'Interpreter', 'none');

    w = {};
    for k = 1:nreg
        w = [w, cellfun(@(x) [res.ch(k).region ': ' x], res.ch(k).warnings, 'UniformOutput', false)]; %#ok<AGROW>
    end
    if ~isempty(w)
        tl.OuterPosition = [0 0.045 1 0.955];
        annotation(fig, 'textbox', [0.005 0 0.99 0.042], 'String', ['WARNING  ' strjoin(w, '  |  ')], ...
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
    axis(ax, 'xy');
    ylim(ax, [f(1) f(end)]);
    if logy
        set(ax, 'YScale', 'log');
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

function s = blank(v)
    if isempty(v) || (isnumeric(v) && isnan(v))
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
