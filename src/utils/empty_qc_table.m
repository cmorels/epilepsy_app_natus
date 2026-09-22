function T = empty_qc_table(tz)
% Correctly-typed 0-row template for qc_report.csv. Single source of
% truth shared by run_pipeline_edf.m (writer) and read_pipeline_csv.m
% (reader) -- keep both in sync with this file.
    T = table('Size', [0 12], ...
        'VariableTypes', {'cell', 'cell', 'cell', 'datetime', 'cell', 'double', 'cell', 'double', 'double', 'double', 'double', 'cell'}, ...
        'VariableNames', {'subject_id', 'region', 'source_file', 'session_start', 'stages_completed', ...
        'n_errors', 'error_messages', 'n_blocks_rejected_short', 'outlier_pct', 'nan_pct', 'n_warnings', 'warning_messages'});
    T = apply_tz(T, tz);
end
