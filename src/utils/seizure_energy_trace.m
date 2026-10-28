function trace = seizure_energy_trace(data, cfg)
% SEIZURE_ENERGY_TRACE  Stage-1 computation of detect_seizures.m (global
% min-max normalization -> per-block bandpass -> per-block Hilbert
% envelope^power -> per-block movmean -> global median*factor threshold),
% factored out into a shared utility WITHOUT modifying detect_seizures.m
% (which is required to stay byte-for-byte unchanged -- see
% src/detect_seizures.m and README.md "dos ramas").
%
% Used by:
%   - evaluate_detections.m, to reconstruct the energy trace for a missed
%     ground-truth window (peak energy, longest run above threshold),
%     since detect_seizures.m's own returned struct only carries scalar
%     .metrics, not the full trace.
%   - detect_seizures_robust.m, whose Stage 1 is required to be numerically
%     identical to detect_seizures.m's.
%
% Every parameter comes from cfg.seizure.* (the LEGACY config), regardless
% of which branch calls this -- Stage 1 is shared and unchanged; only what
% happens to the resulting trace afterward differs between branches.
%
%   trace = seizure_energy_trace(data, cfg)
%
% OUTPUT (struct trace): energy_full (NaN outside surviving-block extent),
% bp_full (band-passed signal, same padding), t_rel, fs, threshold,
% median_energy, pct_above, raw_blocks (all valid blocks), blocks (struct
% array of SURVIVING blocks: block_id, trimmed_start, trimmed_end, energy),
% rejected_blocks (table, same shape as detect_seizures.m's qc field).

    signal = data.signal(:);
    fs = data.fs;
    t_rel = data.t_rel(:);
    valid_mask = data.valid_mask(:);

    valid_signal = signal(valid_mask);
    if numel(valid_signal) < 2
        error('seizure_energy_trace:NotEnoughValidSamples', '%s has fewer than 2 valid (non-gap) samples.', data.file);
    end

    min_val = min(valid_signal);
    max_val = max(valid_signal);
    lfp_norm = (signal - min_val) / (max_val - min_val);

    raw_blocks = mask_to_segments(valid_mask);
    lfp_bp_full = filter_by_blocks(lfp_norm, raw_blocks, ...
        @(v) bandpass(v, cfg.seizure.bandpass_band, fs, 'Steepness', cfg.seizure.bandpass_steepness));

    trim_samples = round(cfg.seizure.edge_trim_s * fs);
    window_samples = round(cfg.seizure.window_sec * fs);
    min_len_samples = 2 * trim_samples + ceil(cfg.seizure.min_seizure_duration * fs);

    blocks = struct('block_id', {}, 'trimmed_start', {}, 'trimmed_end', {}, 'energy', {});
    rejected_rows = cell(0, 5);
    energy_full = NaN(size(signal));

    for i = 1:size(raw_blocks, 1)
        s = raw_blocks(i, 1);
        e = raw_blocks(i, 2);
        if (e - s + 1) < min_len_samples
            rejected_rows(end+1, :) = {i, t_rel(s), t_rel(e), (e - s + 1) / fs, 'too_short_after_trim'}; %#ok<AGROW>
            continue;
        end

        ts = s + trim_samples;
        te = e - trim_samples;

        bp_trimmed = lfp_bp_full(ts:te);
        envelope = abs(hilbert(bp_trimmed));
        power = envelope .^ cfg.seizure.power_exponent;
        energy = movmean(power, window_samples);

        idx = numel(blocks) + 1;
        blocks(idx).block_id = i;
        blocks(idx).trimmed_start = ts;
        blocks(idx).trimmed_end = te;
        blocks(idx).energy = energy(:);
        energy_full(ts:te) = energy(:);
    end

    rejected_blocks = cell2table(rejected_rows, ...
        'VariableNames', {'block_id', 'start_s', 'end_s', 'duration_s', 'reason'});

    if isempty(blocks)
        energy_all = zeros(0, 1);
    else
        energy_all = vertcat(blocks.energy);
    end

    if isempty(energy_all)
        median_energy = NaN;
        threshold = NaN;
        pct_above = 0;
    else
        median_energy = median(energy_all);
        threshold = median_energy * cfg.seizure.median_factor;
        pct_above = 100 * sum(energy_all > threshold) / numel(energy_all);
    end

    trace = struct();
    trace.energy_full = energy_full;
    trace.bp_full = lfp_bp_full;
    trace.t_rel = t_rel;
    trace.fs = fs;
    trace.threshold = threshold;
    trace.median_energy = median_energy;
    trace.pct_above = pct_above;
    trace.raw_blocks = raw_blocks;
    trace.blocks = blocks;
    trace.rejected_blocks = rejected_blocks;
end
