function T = empty_seizure_summary_table(tz)
% Correctly-typed 0-row template for seizures_summary.csv. Single source
% of truth shared by run_pipeline_edf.m (writer) and read_pipeline_csv.m
% (reader) -- keep both in sync with this file.
    names = {'subject_id', 'region', 'session_start', 'source_file', 'total_duration_min', 'valid_duration_min', ...
        'n_gaps', 'gap_duration_min', 'n_seizures', 'total_seizure_time_s', 'pct_time_in_seizure', ...
        'mean_duration_s', 'min_duration_s', 'max_duration_s', 'median_energy', 'threshold_value', 'pct_above_thr', ...
        'n_segments', 'n_rejected', 'bandpass_low', 'bandpass_high', 'power_exponent', 'window_s', 'median_factor', 'min_seizure_duration'};
    types = [{'cell', 'cell', 'datetime', 'cell'}, repmat({'double'}, 1, numel(names) - 4)];
    T = table('Size', [0 numel(names)], 'VariableTypes', types, 'VariableNames', names);
    T = apply_tz(T, tz);
end
