function [x_tf, x_disp] = spec_prepare_signal(x, fs, cfg)
% SPEC_PREPARE_SIGNAL  The two versions of a loaded window.
%
%   [x_tf, x_disp] = spec_prepare_signal(x, fs, cfg)
%
%   x_tf   : signal for time-frequency: mean and linear trend removed per
%            gap-free block, NO band-pass (filter_by_blocks + detrend).
%   x_disp : display trace: zero-phase Butterworth band-pass
%            cfg.spec.trace_band (filtfilt) per gap-free block; a block too
%            short for filtfilt is NaN.
% Neither crosses a NaN gap, neither interpolates.

    x = x(:);
    blocks = mask_to_segments(~isnan(x));
    x_tf = filter_by_blocks(x, blocks, @(v) detrend(v));

    band = cfg.spec.trace_band;
    band(2) = min(band(2), 0.95 * fs / 2);
    n_half = max(1, round(cfg.spec.trace_filter_order / 2));   % butter band-pass of n -> order 2n
    [z, p, k] = butter(n_half, band / (fs / 2), 'bandpass');
    [sos, g] = zp2sos(z, p, k);
    min_len = 3 * (2 * size(sos, 1)) + 1;
    x_disp = filter_by_blocks(x, blocks, @(v) filt_or_nan(sos, g, v, min_len));
end

function y = filt_or_nan(sos, g, v, min_len)
    if numel(v) <= min_len
        y = NaN(size(v));
    else
        y = filtfilt(sos, g, v);
    end
end
