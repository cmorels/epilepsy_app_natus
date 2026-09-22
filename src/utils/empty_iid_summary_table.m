function T = empty_iid_summary_table(tz)
% Correctly-typed 0-row template for iid_summary.csv. Single source of
% truth shared by run_pipeline_edf.m (writer) and read_pipeline_csv.m
% (reader) -- keep both in sync with this file.
    names = {'subject_id', 'region', 'session_start', 'source_file', 'total_duration_min', 'analyzed_duration_min', ...
        'excluded_duration_min', 'n_exclusion_zones', 'baseline_uV', 'lower_threshold_uV', 'upper_threshold_uV', ...
        'total_peaks', 'total_complexes', 'n_single', 'n_polyspike', 'pct_polyspike', 'complexes_per_min', 'single_per_min', ...
        'polyspikes_per_min', 'n_bursts', 'bursts_per_hour', 'mean_spikes_per_polyspike'};
    types = [{'cell', 'cell', 'datetime', 'cell'}, repmat({'double'}, 1, numel(names) - 4)];
    T = table('Size', [0 numel(names)], 'VariableTypes', types, 'VariableNames', names);
    T = apply_tz(T, tz);
end
