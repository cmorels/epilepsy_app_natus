function [files, chk] = cs_figures_v2(cfg, S, R, rv, DW, M, fig_dir)
% CS_FIGURES_V2  Section 7 figures of v2 (PNG, 300 dpi, no calendar dates):
%   violin_group.png       one violin, epileptic group: (a) seizures per mouse,
%                          (b) seizures per 24 h observed
%   violin_by_day.png      one violin per day (cs_day_windows.m days): seizures
%                          of each mouse that day (days with T_seizure_obs >=
%                          cfg.violin_min_h); violin_by_day_rate.png: per 24 h
%   raster_<mouse>.png     rows = days with any recording ("Dia k", empty days
%                          dropped), x = hours since the day start, dark phase
%                          at its real local time, layers as in the legend
% chk: {file, calendar-date strings found} for the section-9 check.

    if ~isfolder(fig_dir), mkdir(fig_dir); end
    files = {};
    chk = cell(0, 2);
    excl_mice = unique(rv.excluded.mouse);
    excl_note = strjoin(unique(strcat(rv.excluded.mouse, {' ('}, rv.excluded.kind, {')'}))', '; ');
    col = [0.80 0.30 0.25];

    % ---- 7.1 group violin -------------------------------------------------
    x = R.mouse(strcmp(R.mouse.scope, 'merged') & strcmp(R.mouse.group, 'epileptic'), :);
    n = height(x);
    vals = {x.n_total, x.n_total ./ x.T_seizure_obs_h * 24};
    ylab = {'crisis por raton (total)', 'crisis / 24 h observadas'};
    days_obs = x.T_seizure_obs_h / 24;
    fig = figure('Visible', 'off', 'Position', [100 100 1000 640], 'Color', 'w');
    for p = 1:2
        ax = axes(fig, 'Position', [0.08 + (p - 1) * 0.49, 0.27, 0.38, 0.65]);   % room below for the footnotes
        y = vals{p};
        o = struct('color', col, 'hollow', ismember(x.mouse, excl_mice), 'labels', {x.mouse}, ...
            'show_label', label_mask(y), 'width', 0.32);
        plot_violin(ax, 1, y, o);
        set(ax, 'XTick', 1, 'XTickLabel', {sprintf('Epilépticos (n = %d)', n)}, 'XLim', [0.4 1.8], 'FontSize', 10);
        yl = ylim(ax); ylim(ax, [0 yl(2)]);
        ylabel(ax, ylab{p}); grid(ax, 'on');
        title(ax, sprintf('(%s) %s', char('a' + p - 1), ylab{p}));
    end
    notes = {['Marcador hueco: raton con algun EDF excluido o pendiente: ' excl_note], ...
        sprintf('(a) Dias observados por raton: mediana %.1f (rango %.1f-%.1f); el total depende de la duracion del registro.', ...
        median(days_obs), min(days_obs), max(days_obs))};
    annotation(fig, 'textbox', [0.02 0.02 0.96 0.13], 'String', notes, 'EdgeColor', 'none', 'FontSize', 8, 'Interpreter', 'none', ...
        'VerticalAlignment', 'top');
    [files, chk] = save_fig(fig, fig_dir, 'violin_group.png', files, chk);

    % ---- 7.2 violin per day --------------------------------------------------
    d = DW(strcmp(DW.group, 'epileptic') & DW.enters_violin, :);
    for v = 1:2
        fig = figure('Visible', 'off', 'Position', [100 100 1200 520], 'Color', 'w');
        ax = axes(fig); hold(ax, 'on');
        ks = unique(d.day);
        lbl = cell(numel(ks), 1);
        for i = 1:numel(ks)
            s = d(d.day == ks(i), :);
            y = s.n_seizures;
            if v == 2, y = s.seizures_per_24h; end
            plot_violin(ax, i, y, struct('color', col, 'min_n', cfg.violin_min_n, 'width', 0.38, ...
                'hollow', ismember(s.mouse, excl_mice), 'marker_size', 18));
            lbl{i} = sprintf('Día %d\\newlinen = %d', ks(i), height(s));
        end
        set(ax, 'XTick', 1:numel(ks), 'XTickLabel', lbl, 'TickLabelInterpreter', 'tex', 'XLim', [0.4 numel(ks) + 0.6], 'FontSize', 9);
        yl = ylim(ax); ylim(ax, [0 yl(2)]); grid(ax, 'on');
        if v == 1
            ylabel(ax, 'crisis del raton en esas 24 h'); name = 'violin_by_day.png';
        else
            ylabel(ax, 'crisis / 24 h observadas'); name = 'violin_by_day_rate.png';
        end
        title(ax, sprintf(['Epilepticos: crisis por dia (dias de 24 h desde la primera muestra de cada raton; ' ...
            'solo raton-dias con >= %g h observadas)'], cfg.violin_min_h));
        [files, chk] = save_fig(fig, fig_dir, name, files, chk);
    end

    % ---- 7.3 raster per mouse -------------------------------------------------
    for m = 1:numel(M)
        files_m = raster(cfg, M(m), DW(strcmp(DW.mouse, M(m).mouse), :), fig_dir);
        [files, chk] = save_fig(files_m{1}, fig_dir, files_m{2}, files, chk);
    end
end

function show = label_mask(y)
% label every point unless two are closer than 2 % of the range; then only the 3 highest and 3 lowest
    y = y(:);
    show = true(size(y));
    s = sort(y);
    if numel(y) > 6 && any(diff(s) < 0.02 * max(range(y), eps))
        [~, o] = sort(y);
        show(:) = false;
        show(o([1:3, end-2:end])) = true;
    end
end

function out = raster(cfg, Mm, dw, fig_dir)
    tz = cfg.tz;
    t0 = Mm.t0;
    days = dw.day(dw.rec_h > 0);
    nr = numel(days);
    grp_es = struct('epileptic', 'epileptico', 'control', 'control');
    fig = figure('Visible', 'off', 'Position', [100 100 1150 max(320, 90 + 26 * nr)], 'Color', 'w');
    ax = axes(fig, 'Position', [0.08 0.10 0.68 0.80]); hold(ax, 'on');
    c = struct('eeg', [0.45 0.75 0.45], 'video', [0.20 0.55 0.20], 'dx', [0.88 0.25 0.25], 'excl', [0.55 0.55 0.55], ...
        'pend', [0.98 0.85 0.20], 'dark', [0.86 0.86 0.90]);
    vid = cs_iv('subtract', Mm.T_sz, Mm.T_eeg);
    dx_only = cs_iv('subtract', Mm.dx, vid);
    for i = 1:nr
        w = t0 + [(days(i) - 1) days(i)] * 86400;
        band(ax, dark_iv(w, tz), w, i, 0.5, c.dark, []);
        band(ax, Mm.excluded, w, i, 0.32, [0.88 0.88 0.88], c.excl);
        band(ax, Mm.pending, w, i, 0.32, c.pend, []);
        band(ax, Mm.T_eeg, w, i, 0.32, c.eeg, []);
        band(ax, vid, w, i, 0.32, [0.80 0.92 0.80], c.video);
        band(ax, dx_only, w, i, 0.32, c.dx, []);
    end
    sz = Mm.seizures;
    mk = {'detected_TP', 'o', [0.10 0.35 0.80], true; 'missed_eeg', '^', [0.95 0.50 0.05], true; 'video', 'p', [0.55 0.15 0.65], true};
    for q = 1:size(mk, 1)
        s = sz(strcmp(sz.source, mk{q, 1}), :);
        [xh, yi] = place(s.start_posix, t0, days);
        if mk{q, 4}
            plot(ax, xh, yi, mk{q, 2}, 'MarkerFaceColor', mk{q, 3}, 'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off');
        end
    end
    [xh, yi] = place(Mm.crisette.start_posix, t0, days);
    plot(ax, xh, yi, 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'MarkerSize', 7, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    % legend proxies
    L = {patch(ax, NaN, NaN, c.eeg, 'EdgeColor', 'none', 'DisplayName', 'observado EEG'), ...
        patch(ax, NaN, NaN, [0.80 0.92 0.80], 'EdgeColor', c.video, 'LineWidth', 1.5, 'DisplayName', 'observado solo por video (trama)'), ...
        patch(ax, NaN, NaN, c.dx, 'EdgeColor', 'none', 'DisplayName', 'desconexion'), ...
        patch(ax, NaN, NaN, [0.88 0.88 0.88], 'EdgeColor', c.excl, 'LineWidth', 1.5, 'DisplayName', 'EDF descartado / excluido (trama)'), ...
        patch(ax, NaN, NaN, c.pend, 'EdgeColor', 'none', 'DisplayName', 'pendiente de reprocesamiento'), ...
        patch(ax, NaN, NaN, c.dark, 'EdgeColor', 'none', 'DisplayName', 'oscuridad (19:30-08:00)'), ...
        plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', mk{1, 3}, 'MarkerEdgeColor', 'k', 'DisplayName', 'crisis TP detectada'), ...
        plot(ax, NaN, NaN, '^', 'MarkerFaceColor', mk{2, 3}, 'MarkerEdgeColor', 'k', 'DisplayName', 'no detectada validada'), ...
        plot(ax, NaN, NaN, 'p', 'MarkerFaceColor', mk{3, 3}, 'MarkerEdgeColor', 'k', 'DisplayName', 'crisis por video'), ...
        plot(ax, NaN, NaN, 'd', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'DisplayName', 'crisette')};
    legend(ax, [L{:}], 'Location', 'eastoutside', 'FontSize', 8);
    set(ax, 'YDir', 'reverse', 'YTick', 1:nr, 'YTickLabel', arrayfun(@(k) sprintf('Día %d', k), days, 'UniformOutput', false), ...
        'XLim', [0 24], 'XTick', 0:2:24, 'YLim', [0.4 nr + 0.6], 'FontSize', 9, 'Layer', 'top');
    xlabel(ax, 'Horas desde el inicio del día');
    title(ax, sprintf('%s (%s)', Mm.mouse, grp_es.(Mm.group)));
    subtitle(ax, sprintf('Cada día empieza a las %s (hora local de la primera muestra)', ...
        char(datetime(t0, 'ConvertFrom', 'posixtime', 'TimeZone', tz, 'Format', 'HH:mm'))));
    out = {fig, sprintf('raster_%s.png', Mm.mouse)};
end

function band(ax, X, w, row, hh, face, hatch)
% draw the part of interval set X inside window w on raster row `row`
    x = cs_iv('intersect', X, w);
    for k = 1:size(x, 1)
        a = (x(k, 1) - w(1)) / 3600;
        b = (x(k, 2) - w(1)) / 3600;
        patch(ax, [a b b a], row + [-hh -hh hh hh], face, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        if ~isempty(hatch)
            hatch_rect(ax, a, b, row - hh, row + hh, hatch, 0.35, 0.25);
        end
    end
end

function iv = dark_iv(w, tz)
% dark phase (cfg lights 19:30-08:00) at its real local time inside window w
    d0 = dateshift(datetime(w(1), 'ConvertFrom', 'posixtime', 'TimeZone', tz), 'start', 'day');
    iv = zeros(0, 2);
    for d = [d0 - caldays(1), d0, d0 + caldays(1)]
        [y, mo, dd] = ymd(d);
        iv(end+1, :) = posixtime(datetime(y, mo, dd, [0 8], 0, 0, 'TimeZone', tz)); %#ok<AGROW>
        iv(end+1, :) = [posixtime(datetime(y, mo, dd, 19, 30, 0, 'TimeZone', tz)), posixtime(datetime(y, mo, dd + 1, 0, 0, 0, 'TimeZone', tz))]; %#ok<AGROW>
    end
    iv = cs_iv('intersect', iv, w);
end

function [xh, yi] = place(t, t0, days)
% hours since the start of the day and raster row of event times t
    k = floor((t - t0) / 86400) + 1;
    [tf, yi] = ismember(k, days);
    xh = (t(tf) - t0 - (k(tf) - 1) * 86400) / 3600;
    yi = yi(tf);
end

function [files, chk] = save_fig(fig, d, name, files, chk)
    bad = cs_fig_dates(fig);
    f = fullfile(d, name);
    exportgraphics(fig, f, 'Resolution', 300);
    close(fig);
    files{end+1} = f;
    chk(end+1, :) = {name, strjoin(bad', ' | ')};
end
