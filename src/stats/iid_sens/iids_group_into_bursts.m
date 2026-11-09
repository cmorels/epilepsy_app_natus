function burst_mat = iids_group_into_bursts(complex_mat, complex_block_id, cfg)
% IIDS_GROUP_INTO_BURSTS  Verbatim copy of the local function group_into_bursts of src/detect_iid.m
% (git HEAD 3672166), renamed only so it does not shadow anything on the path.
% Used by cs_iid_sens_channel.m (threshold sensitivity, v2 section 5.2);
% the lower threshold is already a parameter where it applies.

    if isempty(complex_mat)
        burst_mat = zeros(0, 5);
        return;
    end

    complex_times_s = complex_mat(:, 2);  % Start_s
    Interval_dt = diff(complex_times_s);
    same_block = diff(complex_block_id) == 0;
    Status = (Interval_dt <= cfg.iid.burst_max_gap_s) & same_block;

    D = zeros(0, 5);
    complex_tally = 1;
    dur = 0;
    t_ini = complex_times_s(1);
    blk = complex_block_id(1);

    for k = 1:numel(Interval_dt)
        if Status(k)
            complex_tally = complex_tally + 1;
            dur = dur + Interval_dt(k);
        else
            D = [D; complex_tally, dur, t_ini, t_ini + dur, blk]; %#ok<AGROW>
            complex_tally = 1;
            dur = 0;
            t_ini = complex_times_s(k+1);
            blk = complex_block_id(k+1);
        end
    end
    D = [D; complex_tally, dur, t_ini, t_ini + dur, blk];
    burst_mat = D;

    if ~isempty(burst_mat)
        mask = burst_mat(:, 1) > cfg.iid.burst_min_complexes & ...
               burst_mat(:, 2) > cfg.iid.burst_min_duration_s & ...
               burst_mat(:, 2) < cfg.iid.burst_max_duration_s;
        burst_mat = burst_mat(mask, :);
    end
end
