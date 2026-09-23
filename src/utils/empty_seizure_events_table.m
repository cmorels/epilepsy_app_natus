function T = empty_seizure_events_table(tz)
% Correctly-typed 0-row template for seizures_events.csv. The schema here
% is the single source of truth shared by run_pipeline_edf.m (writer) and
% read_pipeline_csv.m (reader) -- keep both in sync with this file.
%
% seizure_mode + the 5 columns after it (over_max_duration, ll_ratio,
% peak_energy_ratio, hf_ratio_db, envelope_cv) exist for every row
% regardless of which branch produced it (see detect_seizures_robust.m
% and run_pipeline_edf.m's build_seizure_event_rows): legacy rows get
% seizure_mode='legacy' and false/NaN for the rest, since
% detect_seizures.m has no such concepts. read_pipeline_csv.m backfills
% the same defaults when reading a seizures_events.csv written before
% these columns existed.
    T = table('Size', [0 18], ...
        'VariableTypes', {'cell', 'cell', 'datetime', 'cell', 'double', 'double', 'double', 'double', 'datetime', 'datetime', 'double', 'logical', ...
        'cell', 'logical', 'double', 'double', 'double', 'double'}, ...
        'VariableNames', {'subject_id', 'region', 'session_start', 'source_file', 'seizure_id', 'start_s', 'end_s', ...
        'duration_s', 'start_abs', 'end_abs', 'block_id', 'adjacent_to_gap', ...
        'seizure_mode', 'over_max_duration', 'll_ratio', 'peak_energy_ratio', 'hf_ratio_db', 'envelope_cv'});
    T = apply_tz(T, tz);
end
