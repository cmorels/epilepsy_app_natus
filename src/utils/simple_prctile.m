function p = simple_prctile(x, pct)
% SIMPLE_PRCTILE  Linear-interpolation percentile (matches MATLAB's own
% prctile default method), without requiring the Statistics and Machine
% Learning Toolbox -- this project only depends on MATLAB base + Signal
% Processing Toolbox.
%
%   p = simple_prctile(x, pct)   % pct in [0, 100], NaNs in x ignored

    x = x(:);
    x = x(~isnan(x));
    x = sort(x);
    n = numel(x);

    if n == 0
        p = NaN;
        return;
    end
    if n == 1
        p = x(1);
        return;
    end

    rank = (pct / 100) * n + 0.5;
    rank = min(max(rank, 1), n);
    lo = floor(rank);
    hi = ceil(rank);
    if lo == hi
        p = x(lo);
    else
        p = x(lo) + (rank - lo) * (x(hi) - x(lo));
    end
end
