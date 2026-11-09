function [complex_mat, complex_block_id, spike_complex_ids] = iids_group_into_complexes(pos, pos_s, spikes, spike_block_id, fs, cfg)
% IIDS_GROUP_INTO_COMPLEXES  Verbatim copy of the local function group_into_complexes of src/detect_iid.m
% (git HEAD 3672166), renamed only so it does not shadow anything on the path.
% Used by cs_iid_sens_channel.m (threshold sensitivity, v2 section 5.2);
% the lower threshold is already a parameter where it applies.

    if isempty(pos)
        complex_mat = zeros(0, 8);
        complex_block_id = zeros(0, 1);
        spike_complex_ids = zeros(0, 1);
        return;
    end

    max_interspike_samples = round((cfg.iid.max_interspike_ms / 1000) * fs);
    ISI_samples = diff(pos);
    same_block = diff(spike_block_id) == 0;

    complex_id = 1;
    spike_complex_ids = ones(size(pos));
    for i = 1:numel(ISI_samples)
        if ISI_samples(i) > max_interspike_samples || ~same_block(i)
            complex_id = complex_id + 1;
        end
        spike_complex_ids(i+1) = complex_id;
    end

    num_complexes = max(spike_complex_ids);
    complex_mat = zeros(num_complexes, 8);
    complex_block_id = zeros(num_complexes, 1);

    for c = 1:num_complexes
        idx = find(spike_complex_ids == c);
        n_spikes = numel(idx);
        first_i = idx(1);
        last_i = idx(end);
        start_s = pos_s(first_i);
        end_s = pos_s(last_i);
        duration_ms = (end_s - start_s) * 1000;
        max_amp = max(spikes(idx));
        mean_amp = mean(spikes(idx));
        is_poly = n_spikes >= cfg.iid.polyspike_min_n;

        complex_mat(c, :) = [c, start_s, end_s, duration_ms, n_spikes, max_amp, mean_amp, is_poly];
        complex_block_id(c) = spike_block_id(first_i);
    end
end
