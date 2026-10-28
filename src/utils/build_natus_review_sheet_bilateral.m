function T = build_natus_review_sheet_bilateral(seizures_events, iid_bursts, gaps_summary, file_subject_map, tz)
% BUILD_NATUS_REVIEW_SHEET_BILATERAL  natus_review_sheet.csv for a
% reconciled run (cfg.bilateral.rescue_mode ~= 'off'): same 7 columns as
% build_natus_review_sheet.m plus seizure_id, detection_status,
% accepted_in_regions and source_file.
%
% One row per (seizure_id, channel). Seizure rows are placed at their
% event's reference-window start, so every channel's row of one event sits
% together (ordered by region) and the review is done event by event;
% abs_time itself stays each row's OWN start. IID bursts and gaps are
% interleaved by their own start time. Events are keyed by
% (subject_id, source_file, seizure_id): seizure_id is only unique within
% one animal + session.

    n_s = height(seizures_events);
    n_b = height(iid_bursts);
    n_g = height(gaps_summary);

    if n_s > 0
        ref_abs = seizures_events.start_abs - seconds(seizures_events.start_s - seizures_events.ref_start_s);
    else
        ref_abs = seizures_events.start_abs;
    end
    sort_abs = [ref_abs; iid_bursts.start_abs; gaps_summary.start_abs];
    rows_abs = [seizures_events.start_abs; iid_bursts.start_abs; gaps_summary.start_abs];
    rows_dur = [seizures_events.duration_s; iid_bursts.duration_s; gaps_summary.duration_s];
    rows_type = [repmat({'seizure'}, n_s, 1); repmat({'iid_burst'}, n_b, 1); repmat({'gap'}, n_g, 1)];
    rows_region = [seizures_events.region; iid_bursts.region; repmat({''}, n_g, 1)];
    gap_subjects = cellfun(@(s) lookup_subject(s, file_subject_map), gaps_summary.source_file, 'UniformOutput', false);
    rows_subject = [seizures_events.subject_id; iid_bursts.subject_id; gap_subjects];
    rows_sid = [seizures_events.seizure_id; nan(n_b + n_g, 1)];
    rows_status = [seizures_events.detection_status; repmat({''}, n_b + n_g, 1)];
    rows_acc = [seizures_events.accepted_in_regions; repmat({''}, n_b + n_g, 1)];
    rows_src = [seizures_events.source_file; repmat({''}, n_b, 1); gaps_summary.source_file];

    n = numel(rows_abs);
    names = {'abs_time', 'clock_time', 'event_type', 'duration_s', 'region', 'subject_id', 'natus_confirmed', ...
        'seizure_id', 'detection_status', 'accepted_in_regions', 'source_file'};
    if n == 0
        T = table('Size', [0 11], 'VariableTypes', {'datetime', 'cell', 'cell', 'double', 'cell', 'cell', 'cell', ...
            'double', 'cell', 'cell', 'cell'}, 'VariableNames', names);
        T.abs_time.TimeZone = tz;
        return;
    end

    key = table(sort_abs, rows_subject, rows_src, rows_sid, rows_region, rows_abs);
    [~, order] = sortrows(key);
    clock_time = cellstr(string(rows_abs(order), 'HH:mm:ss'));
    T = table(rows_abs(order), clock_time, rows_type(order), rows_dur(order), rows_region(order), ...
        rows_subject(order), repmat({''}, n, 1), rows_sid(order), rows_status(order), rows_acc(order), rows_src(order), ...
        'VariableNames', names);
end

function subj = lookup_subject(source_file, file_subject_map)
    if isKey(file_subject_map, source_file)
        subj = file_subject_map(source_file);
    else
        subj = '';
    end
end
