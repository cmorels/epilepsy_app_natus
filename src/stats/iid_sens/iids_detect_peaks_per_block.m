function [spikes, pos, width, prom, block_id] = iids_detect_peaks_per_block(abs_mV_full, raw_blocks, fs, lower_mV, cfg)
% IIDS_DETECT_PEAKS_PER_BLOCK  Verbatim copy of the local function detect_peaks_per_block of src/detect_iid.m
% (git HEAD 3672166), renamed only so it does not shadow anything on the path.
% Used by cs_iid_sens_channel.m (threshold sensitivity, v2 section 5.2);
% the lower threshold is already a parameter where it applies.

    MinPeakDistance = max(1, round(cfg.iid.min_peak_distance_s * fs));
    MaxPeakWidth = max(1, round(cfg.iid.max_peak_width_s * fs));

    spikes = zeros(0, 1); pos = zeros(0, 1); width = zeros(0, 1); prom = zeros(0, 1); block_id = zeros(0, 1);
    for i = 1:size(raw_blocks, 1)
        s = raw_blocks(i, 1);
        e = raw_blocks(i, 2);
        blk = abs_mV_full(s:e);
        [sp, po, wi, pr] = findpeaks(blk, ...
            'MinPeakDistance', MinPeakDistance, ...
            'MinPeakProminence', cfg.iid.min_peak_prominence_mV, ...
            'MinPeakHeight', lower_mV, ...
            'MaxPeakWidth', MaxPeakWidth);
        if isempty(po)
            continue;
        end
        spikes = [spikes; sp(:)]; pos = [pos; s - 1 + po(:)]; width = [width; wi(:)]; prom = [prom; pr(:)]; %#ok<AGROW>
        block_id = [block_id; repmat(i, numel(po), 1)]; %#ok<AGROW>
    end
end
