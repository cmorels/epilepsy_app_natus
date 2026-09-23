function T = empty_qc_table(tz)
% Correctly-typed 0-row template for qc_report.csv. Single source of
% truth shared by run_pipeline_edf.m (writer) and read_pipeline_csv.m
% (reader) -- keep both in sync with this file.
    names = qc_column_names();
    types = cell(1, numel(names));
    text_cols = {'subject_id', 'region', 'source_file', 'stages_completed', 'error_messages', 'warning_messages', ...
        'case_applied', 'case_source', 'suggested_case', 'quality_class', 'reference_source', 'gain_source', ...
        'seizure_threshold_mode', 'iid_threshold_mode'};
    for i = 1:numel(names)
        if strcmp(names{i}, 'session_start')
            types{i} = 'datetime';
        elseif strcmp(names{i}, 'notch_applied')
            types{i} = 'logical';
        elseif ismember(names{i}, text_cols)
            types{i} = 'cell';
        else
            types{i} = 'double';
        end
    end
    T = table('Size', [0 numel(names)], 'VariableTypes', types, 'VariableNames', names);
    T = apply_tz(T, tz);
end
