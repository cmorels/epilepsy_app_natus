function h = plot_violin(ax, x0, y, opts)
% PLOT_VIOLIN  One violin at x = x0 for non-negative data y (MATLAB's
% violinplot is not used, to keep control of the support).
%
% Density: ksdensity of [y; -y] (reflection at 0) with the bandwidth of y,
% evaluated on [0, max + 1 bandwidth] and doubled, so there is no negative
% tail and the mass near 0 is kept. Half
% width scaled to opts.width. Inside: IQR as a thick bar and the median as
% a white dot; on top: the points with horizontal jitter.
%
%   opts.width      max half-width (default 0.35)
%   opts.color      face color
%   opts.min_n      draw the violin only when numel(y) >= min_n (default 3;
%                   below that only the points)
%   opts.hollow     logical per point -> open marker (default all filled)
%   opts.labels     cellstr per point (shown when opts.show_label(i))
%   opts.show_label logical per point
%   opts.jitter     max horizontal jitter (default 0.6 * width)
%   opts.marker_size
% h: handles (patch, points)

    if nargin < 4, opts = struct(); end
    d = struct('width', 0.35, 'color', [0.55 0.65 0.85], 'min_n', 3, 'hollow', false(numel(y), 1), ...
        'labels', {repmat({''}, numel(y), 1)}, 'show_label', false(numel(y), 1), 'jitter', [], 'marker_size', 28);
    for f = fieldnames(d)'
        if ~isfield(opts, f{1}), opts.(f{1}) = d.(f{1}); end
    end
    if isempty(opts.jitter), opts.jitter = 0.6 * opts.width; end
    y = y(:);
    hold(ax, 'on');
    h = struct('patch', [], 'points', []);
    if numel(y) >= opts.min_n && range(y) > 0
        [~, ~, bw] = ksdensity(y);                % bandwidth of the data itself (the mirrored set would inflate it)
        g = linspace(0, max(y) + bw, 256)';       % up to 1 bandwidth above the maximum: the shape closes
        f = ksdensity([y; -y], g, 'Bandwidth', bw) * 2;
        f = f / max(f) * opts.width;
        h.patch = patch(ax, [x0 - f; flipud(x0 + f)], [g; flipud(g)], opts.color, 'EdgeColor', opts.color * 0.7, ...
            'FaceAlpha', 0.45, 'HandleVisibility', 'off');
        q = prctile(y, [25 50 75]);
        plot(ax, [x0 x0], q([1 3]), '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 5, 'HandleVisibility', 'off');
        plot(ax, x0, q(2), 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', [0.15 0.15 0.15], 'MarkerSize', 7, ...
            'HandleVisibility', 'off');
    elseif numel(y) >= opts.min_n
        plot(ax, x0 + [-1 1] * opts.width, [y(1) y(1)], '-', 'Color', opts.color * 0.7, 'LineWidth', 2, 'HandleVisibility', 'off');
    end
    rs = RandStream('mt19937ar', 'Seed', 7);
    xj = x0 + (rand(rs, numel(y), 1) - 0.5) * 2 * opts.jitter;
    pc = opts.color * 0.6;
    h.points = gobjects(0);
    if any(~opts.hollow)
        h.points(end+1) = scatter(ax, xj(~opts.hollow), y(~opts.hollow), opts.marker_size, pc, 'filled', 'HandleVisibility', 'off');
    end
    if any(opts.hollow)
        h.points(end+1) = scatter(ax, xj(opts.hollow), y(opts.hollow), opts.marker_size, pc, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    for i = find(opts.show_label(:))'
        text(ax, xj(i) + 0.03, y(i), opts.labels{i}, 'FontSize', 7, 'VerticalAlignment', 'middle');
    end
end
