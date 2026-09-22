function T = empty_seizure_events_table(tz)
% Correctly-typed 0-row template for seizures_events.csv. The schema here
% is the single source of truth shared by run_pipeline_edf.m (writer) and
% read_pipeline_csv.m (reader) -- keep both in sync with this file.
    T = table('Size', [0 12], ...
        'VariableTypes', {'cell', 'cell', 'datetime', 'cell', 'double', 'double', 'double', 'double', 'datetime', 'datetime', 'double', 'logical'}, ...
        'VariableNames', {'subject_id', 'region', 'session_start', 'source_file', 'seizure_id', 'start_s', 'end_s', ...
        'duration_s', 'start_abs', 'end_abs', 'block_id', 'adjacent_to_gap'});
    T = apply_tz(T, tz);
end
