function [names, types, defaults] = event_columns(kind)
% EVENT_COLUMNS  Columns added by the review band / event categories
% (resolve_ll_band.m, bilateral_events.m). Single source of truth for the
% writer (run_pipeline_edf.m), the reader (read_pipeline_csv.m) and
% merge_pipeline_runs.m. All are APPENDED to the existing schemas, and only
% when the band is on, so a band-off run writes exactly what it always did.
%
%   kind : 'll_status'            -> {'ll_status'} (seizures_events; present
%                                    whenever the band is on, legacy rows too)
%          'seizures_events'      -> event columns (only if the run has robust
%                                    channels)
%          'seizures_summary'     -> per-channel event counts / rates
%          'natus_review_sheet'   -> extra review-sheet columns
    switch kind
        case 'll_status'
            names = {'ll_status'};
            types = {'cell'};
            defaults = {{'accepted'}};
        case 'seizures_events'
            names = {'event_id', 'category', 'n_accepted', 'n_in_band', 'n_channels_in_group', 'fragmented', ...
                'hemisphere', 'accepted_in_regions', 'ref_start_s', 'ref_end_s', 'figure_individual_path', 'figure_joint_path'};
            types = {'double', 'cell', 'double', 'double', 'double', 'logical', 'cell', 'cell', 'double', 'double', 'cell', 'cell'};
            defaults = {NaN, {''}, NaN, NaN, NaN, false, {''}, {''}, NaN, NaN, {''}, {''}};
        case 'seizures_summary'
            names = {'n_events_crisis', 'n_events_candidates', 'n_events_candidates_in_band', ...
                'n_accepted_this_channel', 'n_in_band_this_channel', ...
                'total_seizure_time_s_crisis', 'pct_time_in_seizure_crisis', ...
                'total_seizure_time_s_crisis_candidates', 'pct_time_in_seizure_crisis_candidates'};
            types = repmat({'double'}, 1, numel(names));
            defaults = repmat({NaN}, 1, numel(names));
        case 'natus_review_sheet'
            names = {'event_id', 'category', 'll_status', 'll_ratio', 'review_priority', 'figure_joint_path', 'source_file'};
            types = {'double', 'cell', 'cell', 'double', 'double', 'cell', 'cell'};
            defaults = {NaN, {''}, {''}, NaN, NaN, {''}, {''}};
        otherwise
            error('event_columns:BadKind', 'Unknown kind ''%s''.', kind);
    end
end
