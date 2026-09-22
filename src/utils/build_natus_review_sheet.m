function T = build_natus_review_sheet(seizures_events, iid_bursts, gaps_summary, file_subject_map, tz)
% BUILD_NATUS_REVIEW_SHEET  Merge seizures + IID bursts + gaps into one
% chronologically-sorted sheet for manual comparison against a Natus
% review. Used by both run_pipeline_edf.m (single run) and
% merge_pipeline_runs.m (several runs combined).

    rows_abs = [seizures_events.start_abs; iid_bursts.start_abs; gaps_summary.start_abs];
    rows_dur = [seizures_events.duration_s; iid_bursts.duration_s; gaps_summary.duration_s];
    rows_type = [repmat({'seizure'}, height(seizures_events), 1); ...
                 repmat({'iid_burst'}, height(iid_bursts), 1); ...
                 repmat({'gap'}, height(gaps_summary), 1)];
    rows_region = [seizures_events.region; iid_bursts.region; repmat({''}, height(gaps_summary), 1)];
    gap_subjects = cellfun(@(s) lookup_subject(s, file_subject_map), gaps_summary.source_file, 'UniformOutput', false);
    rows_subject = [seizures_events.subject_id; iid_bursts.subject_id; gap_subjects];

    n = numel(rows_abs);
    if n == 0
        T = table('Size', [0 7], ...
            'VariableTypes', {'datetime', 'cell', 'cell', 'double', 'cell', 'cell', 'cell'}, ...
            'VariableNames', {'abs_time', 'clock_time', 'event_type', 'duration_s', 'region', 'subject_id', 'natus_confirmed'});
        T.abs_time.TimeZone = tz;
        return;
    end

    [abs_sorted, order] = sort(rows_abs);
    clock_time = cellstr(string(abs_sorted, 'HH:mm:ss'));

    T = table(abs_sorted, clock_time, rows_type(order), rows_dur(order), rows_region(order), rows_subject(order), ...
        repmat({''}, n, 1), ...
        'VariableNames', {'abs_time', 'clock_time', 'event_type', 'duration_s', 'region', 'subject_id', 'natus_confirmed'});
end

function subj = lookup_subject(source_file, file_subject_map)
    if isKey(file_subject_map, source_file)
        subj = file_subject_map(source_file);
    else
        subj = '';
    end
end
