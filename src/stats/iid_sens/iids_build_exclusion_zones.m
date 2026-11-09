function [exclusion_table, inclusion_mask] = iids_build_exclusion_zones(seizures, valid_mask, t_rel, cfg)
% IIDS_BUILD_EXCLUSION_ZONES  Verbatim copy of the local function build_exclusion_zones of src/detect_iid.m
% (git HEAD 3672166), renamed only so it does not shadow anything on the path.
% Used by cs_iid_sens_channel.m (threshold sensitivity, v2 section 5.2);
% the lower threshold is already a parameter where it applies.

    rows = cell(0, 4);  % {reason, source_id, start_s, end_s}

    for i = 1:height(seizures)
        rows(end+1, :) = {'seizure', seizures.id(i), ...
            seizures.start_s(i) - cfg.iid.exclusion_buffer_s, seizures.end_s(i) + cfg.iid.exclusion_buffer_s}; %#ok<AGROW>
    end

    gap_idx = mask_to_segments(~valid_mask);
    for i = 1:size(gap_idx, 1)
        rows(end+1, :) = {'gap', i, t_rel(gap_idx(i, 1)), t_rel(gap_idx(i, 2))}; %#ok<AGROW>
    end

    manual = cfg.iid.exclusion_zones_manual;
    for i = 1:size(manual, 1)
        rows(end+1, :) = {'manual', i, manual(i, 1), manual(i, 2)}; %#ok<AGROW>
    end

    inclusion_mask = valid_mask;

    if isempty(rows)
        exclusion_table = table('Size', [0 8], ...
            'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double', 'datetime', 'datetime'}, ...
            'VariableNames', {'zone_id', 'reason', 'source_id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs'});
        return;
    end

    t_min = t_rel(1);
    t_max = t_rel(end);
    reason = string(rows(:, 1));
    source_id = cell2mat(rows(:, 2));
    start_s = max(cell2mat(rows(:, 3)), t_min);
    end_s = min(cell2mat(rows(:, 4)), t_max);

    keep = end_s > start_s;
    reason = reason(keep); source_id = source_id(keep); start_s = start_s(keep); end_s = end_s(keep);

    for i = 1:numel(start_s)
        inclusion_mask(t_rel >= start_s(i) & t_rel <= end_s(i)) = false;
    end

    zone_id = (1:numel(start_s))';
    exclusion_table = table(zone_id, reason, source_id, start_s, end_s, end_s - start_s, ...
        'VariableNames', {'zone_id', 'reason', 'source_id', 'start_s', 'end_s', 'duration_s'});
end
