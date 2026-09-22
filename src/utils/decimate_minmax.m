function [t_dec, y_dec] = decimate_minmax(t, y, target_points)
% DECIMATE_MINMAX  Downsample (t, y) for plotting a full-recording
% overview without losing visible spikes/transients: each output bin
% keeps both the min AND max sample in that bin (in their real time
% order), instead of a single strided sample that could skip right over
% a narrow event. A no-op if y already has <= 2*target_points samples.
%
% Long recordings (e.g. a 20+ hour EDF at 2 kHz is 100+ million samples)
% make a single-line plot of every sample slow, memory-heavy, and prone
% to failing to save as .fig; this keeps panorama figures tractable while
% zoom figures (a few seconds of data) are unaffected and still plot
% every sample.
%
%   [t_dec, y_dec] = decimate_minmax(t, y, target_points)

    t = t(:);
    y = y(:);
    n = numel(y);

    if n <= 2 * target_points
        t_dec = t;
        y_dec = y;
        return;
    end

    block = ceil(n / target_points);
    n_blocks = ceil(n / block);
    t_dec = NaN(2 * n_blocks, 1);
    y_dec = NaN(2 * n_blocks, 1);

    for b = 1:n_blocks
        s = (b - 1) * block + 1;
        e = min(b * block, n);
        seg = y(s:e);
        [mn, mn_i] = min(seg);
        [mx, mx_i] = max(seg);
        if mn_i <= mx_i
            t_dec(2*b-1) = t(s + mn_i - 1); y_dec(2*b-1) = mn;
            t_dec(2*b)   = t(s + mx_i - 1); y_dec(2*b)   = mx;
        else
            t_dec(2*b-1) = t(s + mx_i - 1); y_dec(2*b-1) = mx;
            t_dec(2*b)   = t(s + mn_i - 1); y_dec(2*b)   = mn;
        end
    end
end
