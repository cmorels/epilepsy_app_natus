function [names, types, defaults] = bilateral_columns(kind)
% BILATERAL_COLUMNS  Extra columns that bilateral_reconcile.m adds to the
% seizure CSVs (cfg.bilateral.rescue_mode ~= 'off'). Single source of
% truth for the writer (run_pipeline_edf.m), the reader
% (read_pipeline_csv.m) and merge_pipeline_runs.m.
%
%   [names, types, defaults] = bilateral_columns(kind)
%
%   kind : 'seizures_events' | 'seizures_summary' | 'natus_review_sheet'
%
% These columns are APPENDED to the existing schema and only when
% reconciliation is on, so a run with rescue_mode='off' writes exactly the
% same files as before this stage existed. defaults are the values used
% to backfill a table that was never reconciled (e.g. merging an 'off'
% run with a reconciled one): detection_status = 'not_reconciled', so an
% unreconciled row can never be mistaken for an accepted/rescued/imputed one.

    switch kind
        case 'seizures_events'
            names = {'channel_seizure_id', 'detection_status', 'accepted_in_n_channels', 'accepted_in_regions', ...
                'n_channels_total', 'is_bilateral_accepted', 'fragmented', 'ref_start_s', 'ref_end_s', 'rejected_by'};
            types = {'double', 'cell', 'double', 'cell', 'double', 'logical', 'logical', 'double', 'double', 'cell'};
            defaults = {NaN, {'not_reconciled'}, NaN, {''}, NaN, false, false, NaN, NaN, {''}};
        case 'seizures_summary'
            names = {'n_seizures_accepted', 'n_seizures_reported', 'n_events_reported', 'n_rescued', 'n_imputed', ...
                'total_seizure_time_s_accepted', 'pct_time_in_seizure_accepted', ...
                'total_seizure_time_s_reported', 'pct_time_in_seizure_reported'};
            types = repmat({'double'}, 1, numel(names));
            defaults = repmat({NaN}, 1, numel(names));
        case 'natus_review_sheet'
            names = {'seizure_id', 'detection_status', 'accepted_in_regions', 'source_file'};
            types = {'double', 'cell', 'cell', 'cell'};
            defaults = {NaN, {''}, {''}, {''}};
        otherwise
            error('bilateral_columns:BadKind', 'Unknown kind ''%s''.', kind);
    end
end
