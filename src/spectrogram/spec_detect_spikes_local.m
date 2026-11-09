function [t_peaks, info] = spec_detect_spikes_local(x, fs, t_rel, include, p)
% SPEC_DETECT_SPIKES_LOCAL  Spike peaks of one interval with the same steps
% as detect_iid.m (which is not called nor modified), so threshold
% variants can be tried without re-running module 2.
%
%   [t_peaks, info] = spec_detect_spikes_local(x, fs, t_rel, include, p)
%
%   x       : clean signal of the loaded window, uV (NaN at gaps)
%   t_rel   : t_rel of each sample of x
%   include : logical, samples where a peak may be kept and that enter the
%             threshold (interval, not NaN, outside exclusion zones)
%   p       : cfg.spec.iid.detect
%
% Steps (detect_iid.m): threshold from |x| (NOT band-passed) over the
% included samples; x_mV = (x - median) / 1000; band-pass p.bandpass_band
% per gap-free block (bandpass, Steepness); |.|; findpeaks per block
% (MinPeakDistance, MinPeakProminence, MinPeakHeight = threshold,
% MaxPeakWidth); drop peaks above p.upper_uV and outside include.
%   threshold_mode 'pipeline_grid': baseline = detect_iid's search on the
%     grid_start:grid_step:grid_max uV grid for the grid_fraction (0.97)
%     quantile; threshold = threshold_factor * baseline
%   'p97'   : threshold = threshold_factor * 97th percentile of |x|
%   'fixed' : threshold = threshold_uV
% The band-pass sees only the loaded window: within its first/last ~0.5 s
% it can differ slightly from detect_iid on the whole recording.

    x = x(:);
    valid = ~isnan(x);
    absx_mV = abs(x(include)) / 1000;
    if isempty(absx_mV)
        error('spec_detect_spikes_local:NoSamples', 'No included samples to set the threshold.');
    end
    switch p.threshold_mode
        case 'pipeline_grid'
            sv = p.grid_start_uV;
            frac = 0;
            while frac < p.grid_fraction && sv <= p.grid_max_uV
                frac = mean(absx_mV < sv / 1000);
                if frac < p.grid_fraction
                    sv = sv + p.grid_step_uV;
                end
            end
            if sv > p.grid_max_uV
                base_uV = p.grid_max_uV - p.grid_step_uV;
            else
                base_uV = sv - p.grid_step_uV;
            end
            thr_uV = p.threshold_factor * base_uV;
        case 'p97'
            base_uV = 1000 * prctile(absx_mV, 97);
            thr_uV = p.threshold_factor * base_uV;
        case 'fixed'
            base_uV = NaN;
            thr_uV = p.threshold_uV;
            if isnan(thr_uV)
                error('spec_detect_spikes_local:NoThreshold', 'threshold_mode ''fixed'' needs threshold_uV.');
            end
        otherwise
            error('spec_detect_spikes_local:BadMode', 'threshold_mode must be pipeline_grid, p97 or fixed.');
    end

    blocks = mask_to_segments(valid);
    x_mV = (x - median(x(valid))) / 1000;
    bp = filter_by_blocks(x_mV, blocks, @(v) bp_or_nan(v, p, fs));
    a = abs(bp);
    min_dist = max(1, round(p.min_distance_s * fs));
    max_w = max(1, round(p.max_width_s * fs));
    pos = zeros(0, 1);
    amp = zeros(0, 1);
    for i = 1:size(blocks, 1)
        s = blocks(i, 1);
        e = blocks(i, 2);
        blk = a(s:e);
        if numel(blk) < 3 || all(isnan(blk))
            continue;
        end
        [pk, po] = findpeaks(blk, 'MinPeakDistance', min_dist, 'MinPeakProminence', p.min_prominence_mV, ...
            'MinPeakHeight', thr_uV / 1000, 'MaxPeakWidth', max_w);
        pos = [pos; s - 1 + po(:)]; %#ok<AGROW>
        amp = [amp; pk(:)]; %#ok<AGROW>
    end
    keep = amp <= p.upper_uV / 1000;
    pos = pos(keep);
    amp = amp(keep);
    keep = include(pos);
    pos = pos(keep);
    amp = amp(keep);
    t_peaks = t_rel(pos);
    info = struct('threshold_mode', p.threshold_mode, 'baseline_uV', base_uV, 'threshold_uV', thr_uV, ...
        'upper_uV', p.upper_uV, 'min_distance_s', p.min_distance_s, 'amplitude_mV', amp);
end

function y = bp_or_nan(v, p, fs)
    if numel(v) < 64
        y = NaN(size(v));
        return;
    end
    y = bandpass(v, p.bandpass_band, fs, 'Steepness', p.bandpass_steepness);
end
