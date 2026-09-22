function T = empty_iid_bursts_table(tz)
% Correctly-typed 0-row template for iid_bursts.csv. Single source of
% truth shared by run_pipeline_edf.m (writer) and read_pipeline_csv.m
% (reader) -- keep both in sync with this file.
    T = table('Size', [0 8], ...
        'VariableTypes', {'cell', 'cell', 'double', 'double', 'datetime', 'datetime', 'double', 'double'}, ...
        'VariableNames', {'subject_id', 'region', 'start_s', 'end_s', 'start_abs', 'end_abs', 'duration_s', 'n_complexes'});
    T = apply_tz(T, tz);
end
