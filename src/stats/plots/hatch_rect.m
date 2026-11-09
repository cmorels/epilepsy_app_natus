function hatch_rect(ax, x0, x1, y0, y1, color, step, slant)
% HATCH_RECT  Diagonal hatching inside the rectangle [x0 x1] x [y0 y1]
% (data units): lines every `step` along x, each rising `slant` in x over
% the rectangle's height, clipped to the rectangle. MATLAB patches have no
% native hatch.

    if x1 <= x0, return; end
    xs = (x0 - slant):step:x1;
    X = NaN(3, numel(xs));
    Y = NaN(3, numel(xs));
    for i = 1:numel(xs)
        t = [max(0, (x0 - xs(i)) / slant), min(1, (x1 - xs(i)) / slant)];
        if t(2) <= t(1), continue; end
        X(1:2, i) = xs(i) + t * slant;
        Y(1:2, i) = y0 + t * (y1 - y0);
    end
    plot(ax, X(:), Y(:), '-', 'Color', color, 'LineWidth', 0.6, 'HandleVisibility', 'off');
end
