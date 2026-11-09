function files = cs_figures_seizures(cfg, S, O, R, fig_dir)
% CS_FIGURES_SEIZURES  Section 7 seizure figures (PNG, 300 dpi):
%   raster_<mouse>.png   day x clock-hour raster: observed time (T_seizure_obs),
%                        disconnections, dark phase shading, seizures by source
%                        (not with cfg.raster_v2: cs_figures_v2.m draws them)
%   daily_counts.png     seizures per local day per epileptic mouse (invalid
%                        days, < cfg.min_valid_h_day observed, in grey)
%   circadian.png        pooled epileptic rate per hour normalized by observed
%                        hours, dark phase (lights_off-lights_on) shaded

    if ~isfolder(fig_dir), mkdir(fig_dir); end
    tz = cfg.tz;
    files = {};
    on = hm(cfg.lights_on) / 60; off = hm(cfg.lights_off) / 60;
    mk = struct('detected_TP', {{'o', [0.10 0.35 0.75]}}, 'missed_eeg', {{'s', [0.85 0.45 0.05]}}, 'video', {{'^', [0.55 0.10 0.60]}});
    om = O.mouse(strcmp({O.mouse.scope}, 'merged'));
    if isfield(cfg, 'raster_v2') && cfg.raster_v2
        om = om([]);   % v2: rasters are drawn by cs_figures_v2.m
    end
    for k = 1:numel(om)
        o = om(k);
        segs_obs = day_segments(o.T_sz, tz);
        segs_dx = day_segments(o.dx, tz);
        if isempty(segs_obs), continue; end
        d0 = min(segs_obs.date); d1 = max(segs_obs.date);
        nd = days(d1 - d0) + 1;
        fig = figure('Visible', 'off', 'Position', [100 100 900 max(300, 40 + 22 * nd)]);
        ax = axes(fig); hold(ax, 'on');
        for d = 0:nd-1
            patch(ax, [off 24 24 off], [d d d+1 d+1] - 0.5, [0.88 0.88 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
            patch(ax, [0 on on 0], [d d d+1 d+1] - 0.5, [0.88 0.88 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        draw(ax, segs_obs, d0, [0.55 0.80 0.55], 0.35);
        draw(ax, segs_dx, d0, [0.90 0.30 0.30], 0.55);
        x = S.seizures(strcmp(S.seizures.mouse, o.mouse), :);
        src = fieldnames(mk);
        for s = 1:numel(src)
            xs = x(strcmp(x.source, src{s}), :);
            if isempty(xs), continue; end
            t = datetime(xs.start_posix, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
            dd = days(datetime(year(t), month(t), day(t)) - d0);
            plot(ax, hour(t) + minute(t) / 60 + second(t) / 3600, dd, mk.(src{s}){1}, 'MarkerSize', 6, ...
                'MarkerFaceColor', mk.(src{s}){2}, 'MarkerEdgeColor', 'k', 'DisplayName', strrep(src{s}, '_', ' '));
        end
        cr = S.crisette(strcmp(S.crisette.mouse, o.mouse), :);
        if ~isempty(cr)
            t = datetime(cr.start_posix, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
            plot(ax, hour(t) + minute(t) / 60, days(datetime(year(t), month(t), day(t)) - d0), 'kx', 'MarkerSize', 9, ...
                'LineWidth', 1.5, 'DisplayName', 'crisette focal');
        end
        set(ax, 'YDir', 'reverse', 'XLim', [0 24], 'XTick', 0:2:24, 'YLim', [-0.5 nd - 0.5], 'YTick', 0:nd-1, ...
            'YTickLabel', cellstr(string(d0 + caldays(0:nd-1), 'dd/MM')), 'FontSize', 8, 'Layer', 'top');
        xlabel(ax, 'hora local'); ylabel(ax, 'dia');
        title(ax, sprintf('%s (%s) - verde: observado (T_{seizure\\_obs}); rojo: desconexion; gris: oscuridad %s-%s', ...
            o.mouse, o.group, cfg.lights_off, cfg.lights_on), 'FontSize', 9, 'Interpreter', 'tex');
        if ~isempty(x) || ~isempty(cr), legend(ax, 'Location', 'eastoutside'); end
        f = fullfile(fig_dir, sprintf('raster_%s.png', o.mouse));
        exportgraphics(fig, f, 'Resolution', 300); close(fig);
        files{end+1} = f; %#ok<AGROW>
    end

    % daily counts, epileptic mice
    D = R.daily(strcmp(R.daily.group, 'epileptic'), :);
    mice = unique(D.mouse);
    if ~isempty(mice)
        nc = 3; nr = ceil(numel(mice) / nc);
        fig = figure('Visible', 'off', 'Position', [100 100 1200 220 * nr]);
        tl = tiledlayout(fig, nr, nc, 'TileSpacing', 'compact');
        for i = 1:numel(mice)
            ax = nexttile(tl); hold(ax, 'on');
            t = D(strcmp(D.mouse, mice{i}), :);
            bar(ax, t.date(t.valid_day), t.n_seizures(t.valid_day), 'FaceColor', [0.10 0.35 0.75]);
            if any(~t.valid_day)
                bar(ax, t.date(~t.valid_day), t.n_seizures(~t.valid_day), 'FaceColor', [0.75 0.75 0.75]);
            end
            title(ax, mice{i}); ylabel(ax, 'crisis/dia'); ax.FontSize = 8;
        end
        title(tl, sprintf('Crisis por dia local (gris: dia con < %g h observadas)', cfg.min_valid_h_day));
        f = fullfile(fig_dir, 'daily_counts.png');
        exportgraphics(fig, f, 'Resolution', 300); close(fig);
        files{end+1} = f;
    end

    % circadian, pooled epileptic
    H = R.circ_hourly(strcmp(R.circ_hourly.unit, 'ALL_epileptic'), :);
    if ~isempty(H)
        fig = figure('Visible', 'off', 'Position', [100 100 800 420]);
        ax = axes(fig); hold(ax, 'on');
        ymax = max(H.rate_per_24h) * 1.15 + eps;
        patch(ax, [off 24 24 off], [0 0 ymax ymax], [0.88 0.88 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
        patch(ax, [0 on on 0], [0 0 ymax ymax], [0.88 0.88 0.92], 'EdgeColor', 'none', 'HandleVisibility', 'off');
        bar(ax, H.hour + 0.5, H.rate_per_24h, 1, 'FaceColor', [0.10 0.35 0.75], 'DisplayName', 'crisis / 24 h observadas');
        yyaxis(ax, 'right');
        plot(ax, H.hour + 0.5, H.observed_h, 'k.-', 'DisplayName', 'horas observadas');
        ylabel(ax, 'horas observadas en esa hora del reloj');
        yyaxis(ax, 'left');
        set(ax, 'XLim', [0 24], 'XTick', 0:2:24, 'YLim', [0 ymax]);
        xlabel(ax, 'hora local'); ylabel(ax, 'tasa (crisis por 24 h de observacion)');
        title(ax, sprintf('Epilepticos agrupados: tasa por hora normalizada por exposicion (n=%d crisis)', sum(H.n_seizures)));
        legend(ax, 'Location', 'northoutside', 'Orientation', 'horizontal');
        f = fullfile(fig_dir, 'circadian.png');
        exportgraphics(fig, f, 'Resolution', 300); close(fig);
        files{end+1} = f;
    end
end

function draw(ax, segs, d0, col, alpha)
    for i = 1:height(segs)
        d = days(segs.date(i) - d0);
        patch(ax, [segs.h0(i) segs.h1(i) segs.h1(i) segs.h0(i)], [d d d+1 d+1] - 0.4 + [0 0 -0.2 -0.2], col, ...
            'EdgeColor', 'none', 'FaceAlpha', alpha, 'HandleVisibility', 'off');
    end
end

function segs = day_segments(iv, tz)
% Split POSIX intervals into (local date, start hour, end hour) pieces.
    rows = {};
    for i = 1:size(iv, 1)
        t = datetime(iv(i, 1), 'ConvertFrom', 'posixtime', 'TimeZone', tz);
        tend = datetime(iv(i, 2), 'ConvertFrom', 'posixtime', 'TimeZone', tz);
        while t < tend
            d = datetime(year(t), month(t), day(t), 'TimeZone', tz);
            nxt = min(tend, datetime(year(t), month(t), day(t) + 1, 'TimeZone', tz));
            h0 = hours(t - d); h1 = h0 + hours(nxt - t);
            rows(end+1, :) = {datetime(year(t), month(t), day(t)), h0, min(24, h1)}; %#ok<AGROW>
            t = nxt;
        end
    end
    if isempty(rows)
        segs = table(NaT(0, 1), zeros(0, 1), zeros(0, 1), 'VariableNames', {'date', 'h0', 'h1'});
    else
        segs = cell2table(rows, 'VariableNames', {'date', 'h0', 'h1'});
    end
end

function m = hm(s)
    v = sscanf(s, '%d:%d');
    m = v(1) * 60 + v(2);
end
