function T = build_natus_review_sheet_events(seizures_events, iid_bursts, gaps_summary, file_subject_map, tz)
% BUILD_NATUS_REVIEW_SHEET_EVENTS  natus_review_sheet.csv when the review
% band / event categories are on (bilateral_events.m): the 7 columns of
% build_natus_review_sheet.m plus event_columns('natus_review_sheet'):
% event_id, category, ll_status, ll_ratio, review_priority,
% figure_joint_path, source_file. natus_confirmed stays the (empty) column
% for the reviewer's verdict.
%
% One row per (event_id, channel) that detected it, placed at the event's
% reference-window start so all rows of one event sit together (keyed by
% subject_id + source_file + event_id: event_id is only unique per animal +
% session); IID bursts and gaps interleaved by their own start.
% review_priority: 1 = Candidates, 2 = Candidates_in_band, 3 = Crisis
% (the contested ones first); NaN for legacy / non-seizure rows.

    E = seizures_events;
    n_s = height(E); n_b = height(iid_bursts); n_g = height(gaps_summary);
    prio = nan(n_s, 1);
    prio(strcmp(E.category, 'Candidates')) = 1;
    prio(strcmp(E.category, 'Candidates_in_band')) = 2;
    prio(strcmp(E.category, 'Crisis')) = 3;
    ref_abs = E.start_abs;
    has_ev = ~isnan(E.event_id);
    ref_abs(has_ev) = E.start_abs(has_ev) - seconds(E.start_s(has_ev) - E.ref_start_s(has_ev));

    sort_abs = [ref_abs; iid_bursts.start_abs; gaps_summary.start_abs];
    rows_abs = [E.start_abs; iid_bursts.start_abs; gaps_summary.start_abs];
    rows_dur = [E.duration_s; iid_bursts.duration_s; gaps_summary.duration_s];
    rows_type = [repmat({'seizure'}, n_s, 1); repmat({'iid_burst'}, n_b, 1); repmat({'gap'}, n_g, 1)];
    rows_region = [E.region; iid_bursts.region; repmat({''}, n_g, 1)];
    gap_subjects = cellfun(@(s) lookup_subject(s, file_subject_map), gaps_summary.source_file, 'UniformOutput', false);
    rows_subject = [E.subject_id; iid_bursts.subject_id; gap_subjects];
    pad_c = repmat({''}, n_b + n_g, 1); pad_n = nan(n_b + n_g, 1);
    rows_eid = [E.event_id; pad_n];
    rows_cat = [E.category; pad_c];
    rows_st = [E.ll_status; pad_c];
    rows_ll = [E.ll_ratio; pad_n];
    rows_prio = [prio; pad_n];
    rows_fig = [E.figure_joint_path; pad_c];
    rows_src = [E.source_file; repmat({''}, n_b, 1); gaps_summary.source_file];

    names = [{'abs_time', 'clock_time', 'event_type', 'duration_s', 'region', 'subject_id', 'natus_confirmed'}, ...
        event_columns('natus_review_sheet')];
    n = numel(rows_abs);
    if n == 0
        T = table('Size', [0 numel(names)], 'VariableTypes', {'datetime', 'cell', 'cell', 'double', 'cell', 'cell', 'cell', ...
            'double', 'cell', 'cell', 'double', 'double', 'cell', 'cell'}, 'VariableNames', names);
        T.abs_time.TimeZone = tz;
        return;
    end
    key = table(sort_abs, rows_subject, rows_src, rows_eid, rows_region, rows_abs);
    [~, order] = sortrows(key);
    T = table(rows_abs(order), cellstr(string(rows_abs(order), 'HH:mm:ss')), rows_type(order), rows_dur(order), ...
        rows_region(order), rows_subject(order), repmat({''}, n, 1), rows_eid(order), rows_cat(order), rows_st(order), ...
        rows_ll(order), rows_prio(order), rows_fig(order), rows_src(order), 'VariableNames', names);
end

function subj = lookup_subject(source_file, file_subject_map)
    if isKey(file_subject_map, source_file)
        subj = file_subject_map(source_file);
    else
        subj = '';
    end
end
