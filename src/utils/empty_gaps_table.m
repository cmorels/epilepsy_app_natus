function T = empty_gaps_table(tz)
% Correctly-typed 0-row template for gaps_summary.csv. Single source of
% truth shared by edf_import.m / run_pipeline_edf.m (writers) and
% read_pipeline_csv.m (reader) -- keep all in sync with this file.
    T = table('Size', [0 9], ...
        'VariableTypes', {'double', 'double', 'double', 'double', 'datetime', 'datetime', 'double', 'double', 'cell'}, ...
        'VariableNames', {'gap_id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', 'prev_record_idx', 'next_record_idx', 'source_file'});
    T = apply_tz(T, tz);
end
