function [y, n_blocks_skipped, freqs_used] = apply_notch_blocks(x, valid_mask, fs, cfg)
% APPLY_NOTCH_BLOCKS  Zero-phase mains notch (designfilt bandstopiir +
% filtfilt), applied per contiguous valid block, never across a NaN gap.
%
%   [y, n_blocks_skipped, freqs_used] = apply_notch_blocks(x, valid_mask, fs, cfg)
%
% filtfilt (zero-phase) is mandatory here, not a stylistic choice: shifting
% event times is exactly what this whole project measures against, so any
% phase distortion from the notch would corrupt the thing being validated.
%
% Notches cfg.precondition.notch_freqs plus any cfg.precondition.notch_harmonics
% multiples of them that land under Nyquist and inside the widest configured
% analysis band (the union of cfg.seizure.bandpass_band and
% cfg.iid.bandpass_band) -- notching a frequency neither detector looks at
% would be a no-op at best.
%
% A block shorter than 3 * cfg.precondition.notch_order samples is left
% UNFILTERED (original values kept, not zeroed/NaN-ed) and counted in
% n_blocks_skipped, since filtfilt needs enough samples to pad each stage
% internally. Filters are designed once per (fs, frequency) and cached
% (persistent) across calls, so repeated calls across channels/blocks in
% the same session don't re-run designfilt every time.

    x = x(:);
    blocks = mask_to_segments(valid_mask);

    freqs = resolve_notch_freqs(fs, cfg);
    freqs_used = freqs;
    min_len = 3 * cfg.precondition.notch_order;

    if isempty(blocks) || isempty(freqs)
        y = x;
        n_blocks_skipped = 0;
        return;
    end

    filters = cell(1, numel(freqs));
    for i = 1:numel(freqs)
        filters{i} = get_cached_filter(fs, freqs(i), cfg.precondition.notch_halfwidth_hz, cfg.precondition.notch_order);
    end

    block_lengths = blocks(:, 2) - blocks(:, 1) + 1;
    n_blocks_skipped = sum(block_lengths < min_len);

    y = filter_by_blocks(x, blocks, @(v) notch_one_block(v, min_len, filters));
end

%% ======================================================================
function v_out = notch_one_block(v, min_len, filters)
    if numel(v) < min_len
        v_out = v;
        return;
    end
    v_out = v;
    for i = 1:numel(filters)
        v_out = filtfilt(filters{i}, v_out);
    end
end

function freqs = resolve_notch_freqs(fs, cfg)
    base = cfg.precondition.notch_freqs(:)';
    all_freqs = base;
    for h = cfg.precondition.notch_harmonics(:)'
        all_freqs = [all_freqs, base * h]; %#ok<AGROW>
    end
    all_freqs = unique(all_freqs);

    bands = [cfg.seizure.bandpass_band; cfg.iid.bandpass_band];
    widest = [min(bands(:, 1)), max(bands(:, 2))];

    nyquist = fs / 2;
    freqs = all_freqs(all_freqs < nyquist & all_freqs >= widest(1) & all_freqs <= widest(2));
end

function d = get_cached_filter(fs, f0, halfwidth_hz, order)
    persistent cache
    if isempty(cache)
        cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end
    key = sprintf('fs=%.6f_f0=%.6f_hw=%.6f_ord=%d', fs, f0, halfwidth_hz, order);
    if isKey(cache, key)
        d = cache(key);
        return;
    end
    d = designfilt('bandstopiir', 'FilterOrder', order, ...
        'HalfPowerFrequency1', f0 - halfwidth_hz, 'HalfPowerFrequency2', f0 + halfwidth_hz, ...
        'SampleRate', fs);
    cache(key) = d;
end
