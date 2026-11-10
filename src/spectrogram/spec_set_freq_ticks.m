function spec_set_freq_ticks(ax, cfg, frange, scale, ticks)
% SPEC_SET_FREQ_TICKS  Frequency (y) axis of a time-frequency panel with
% plain Hz labels: cfg.spec.freq_ticks (default 10, 30, 100) inside the
% panel's range, plus the panel's maximum frequency as the last label
% (cfg.spec.freq_ticks_max), e.g. 10 30 100 150 for a 2-150 Hz map.
%
%   spec_set_freq_ticks(ax, cfg, frange)            % log axis
%   spec_set_freq_ticks(ax, cfg, frange, 'linear')  % linear axis
%   spec_set_freq_ticks(ax, cfg, frange, scale, ticks)   % these ticks instead of freq_ticks
%
%   cfg    : cfg.spec, or a cfg with a .spec field (fields missing in a .mat
%            saved by an older version take the defaults above)
%   frange : [fmin fmax] of that panel
% On a linear axis a tick below fmin (e.g. 0) extends the axis down to it
% (no data there: grey background).

    if nargin < 4 || isempty(scale)
        scale = 'log';
    end
    if isfield(cfg, 'spec')
        cfg = cfg.spec;
    end
    if nargin < 5 || isempty(ticks)
        ticks = [10 30 100];
        if isfield(cfg, 'freq_ticks') && ~isempty(cfg.freq_ticks)
            ticks = cfg.freq_ticks;
        end
    end
    add_max = ~isfield(cfg, 'freq_ticks_max') || cfg.freq_ticks_max;

    ticks = sort(ticks(:))';
    if strcmp(scale, 'linear')
        t = ticks(ticks <= frange(2));               % ticks below fmin extend the axis
    else
        t = ticks(ticks >= frange(1) & ticks <= frange(2));
    end
    fmax = frange(2);
    if add_max && (isempty(t) || abs(t(end) - fmax) > 1e-6 * fmax)
        t = [t, fmax];
    end
    labels = arrayfun(@(v) sprintf('%g', round(v, 1)), t, 'UniformOutput', false);
    if strcmp(scale, 'linear')
        lo = frange(1);
        if ~isempty(t)
            lo = min(t(1), frange(1));
        end
        ylim(ax, [lo fmax]);
        set(ax, 'YScale', 'linear', 'YMinorTick', 'off');
    else
        set(ax, 'YScale', 'log', 'YMinorTick', 'off');
    end
    if ~isempty(t)
        set(ax, 'YTick', t, 'YTickLabel', labels);
    end
end
