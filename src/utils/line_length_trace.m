function [ll_full, ll_median_global] = line_length_trace(x_bp, valid_blocks, fs, cfg)
% LINE_LENGTH_TRACE  Per-sample line length (moving sum of |diff(x)|,
% centered window) computed on the SAME band-passed, min-max normalized
% signal the energy trace uses (see seizure_energy_trace.m), over the same
% trimmed block extents, so both traces are sample-for-sample coherent.
%
% The min-max normalization cancels out in any ratio ll/ll_median_global
% (both scale by the same constant), so ll_ratio is scale-invariant --
% this is what lets detect_seizures_robust.m compare it against a fixed
% threshold regardless of a channel's absolute gain.
%
%   [ll_full, ll_median_global] = line_length_trace(x_bp, valid_blocks, fs, cfg)
%
%   x_bp         : band-passed signal (e.g. seizure_energy_trace's .bp_full)
%   valid_blocks : Nx2 [start end] index pairs -- the TRIMMED, surviving
%                  block extents (e.g. from seizure_energy_trace's .blocks,
%                  [.trimmed_start .trimmed_end]), NOT the raw untrimmed
%                  blocks -- never computed across a NaN or a block edge.
%   cfg          : struct from pipeline_config.m (uses cfg.seizure_robust.ll_window_s)
%
% OUTPUT: ll_full (NaN outside valid_blocks), ll_median_global (median of
% the trace over the concatenation of all blocks -- the normalization
% reference for ll_ratio).

    window_samples = round(cfg.seizure_robust.ll_window_s * fs);
    ll_full = NaN(size(x_bp));
    block_vals = cell(size(valid_blocks, 1), 1);

    for i = 1:size(valid_blocks, 1)
        s = valid_blocks(i, 1);
        e = valid_blocks(i, 2);
        seg = x_bp(s:e);
        seg = seg(:);

        if numel(seg) < 2
            ll_full(s:e) = 0;
            block_vals{i} = zeros(numel(seg), 1);
            continue;
        end

        d = abs(diff(seg));
        d_aligned = [d(1); d];  % length == numel(seg), one |diff| value per sample
        ll_seg = movsum(d_aligned, window_samples);

        ll_full(s:e) = ll_seg;
        block_vals{i} = ll_seg;
    end

    if isempty(block_vals)
        ll_median_global = NaN;
    else
        ll_median_global = median(vertcat(block_vals{:}));
    end
end
