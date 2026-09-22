function result = merge_pipeline_runs(output_roots, merged_output_dir, tz)
% MERGE_PIPELINE_RUNS  Combine several separate run_pipeline_edf.m runs
% (e.g. one per subject, or one per recording day) into one consolidated
% set of summaries.
%
% run_pipeline_edf.m only consolidates the files it processes in a single
% call; a second call with a different input folder does not append to an
% existing 05_summaries/, it overwrites it. This function is the
% supported way to combine multiple such runs afterward, instead of
% re-running everything through one shared cfg.edf.subject_id (which
% would mislabel every file with the same subject).
%
%   result = merge_pipeline_runs(output_roots, merged_output_dir)
%   result = merge_pipeline_runs(output_roots, merged_output_dir, tz)
%
%   output_roots     : cell array of cfg.paths.output_root paths, one per
%                       previous run_pipeline_edf call (each must contain
%                       a 05_summaries/ folder)
%   merged_output_dir : where to write the combined 05_summaries/
%   tz                : TimeZone for datetime columns (default
%                       'Europe/Paris' -- must match what every run
%                       actually used; plain CSV text carries no
%                       timezone, see README.md)
%
% OUTPUT (struct result): the 9 merged tables (seizures_events,
% seizures_summary, iid_events, iid_summary, iid_bursts, gaps_summary,
% qc_report, natus_review_sheet) and .paths to the written files.

    if nargin < 3 || isempty(tz)
        tz = 'Europe/Paris';
    end
    if ~iscell(output_roots) || isempty(output_roots)
        error('merge_pipeline_runs:BadInput', 'output_roots must be a non-empty cell array of paths.');
    end

    seizure_event_parts = {}; seizure_summary_parts = {};
    iid_event_parts = {}; iid_summary_parts = {}; iid_burst_parts = {};
    gap_parts = {}; qc_parts = {};
    file_subject_map = containers.Map('KeyType', 'char', 'ValueType', 'char');

    for r = 1:numel(output_roots)
        summaries_dir = fullfile(output_roots{r}, '05_summaries');
        if ~isfolder(summaries_dir)
            error('merge_pipeline_runs:MissingRun', ...
                'No 05_summaries/ found under %s -- is this a run_pipeline_edf.m output_root?', output_roots{r});
        end
        fprintf('merge_pipeline_runs: reading %s\n', summaries_dir);

        seizure_event_parts{end+1} = read_pipeline_csv(fullfile(summaries_dir, 'seizures_events.csv'), 'seizures_events', tz); %#ok<AGROW>
        seizure_summary_parts{end+1} = read_pipeline_csv(fullfile(summaries_dir, 'seizures_summary.csv'), 'seizures_summary', tz); %#ok<AGROW>
        iid_event_parts{end+1} = read_pipeline_csv(fullfile(summaries_dir, 'iid_events.csv'), 'iid_events', tz); %#ok<AGROW>
        iid_summary_parts{end+1} = read_pipeline_csv(fullfile(summaries_dir, 'iid_summary.csv'), 'iid_summary', tz); %#ok<AGROW>
        iid_burst_parts{end+1} = read_pipeline_csv(fullfile(summaries_dir, 'iid_bursts.csv'), 'iid_bursts', tz); %#ok<AGROW>
        gaps_r = read_pipeline_csv(fullfile(summaries_dir, 'gaps_summary.csv'), 'gaps', tz);
        gap_parts{end+1} = gaps_r; %#ok<AGROW>
        qc_r = read_pipeline_csv(fullfile(summaries_dir, 'qc_report.csv'), 'qc', tz);
        qc_parts{end+1} = qc_r; %#ok<AGROW>

        for i = 1:height(qc_r)
            if ~isKey(file_subject_map, qc_r.source_file{i})
                file_subject_map(qc_r.source_file{i}) = qc_r.subject_id{i};
            end
        end
    end

    seizures_events = vertcat_or_empty(seizure_event_parts, @() empty_seizure_events_table(tz));
    seizures_summary = vertcat_or_empty(seizure_summary_parts, @() empty_seizure_summary_table(tz));
    iid_events = vertcat_or_empty(iid_event_parts, @() empty_iid_events_table(tz));
    iid_summary = vertcat_or_empty(iid_summary_parts, @() empty_iid_summary_table(tz));
    iid_bursts = vertcat_or_empty(iid_burst_parts, @() empty_iid_bursts_table(tz));
    gaps_summary = vertcat_or_empty(gap_parts, @() empty_gaps_table(tz));
    qc_report = vertcat_or_empty(qc_parts, @() empty_qc_table(tz));

    natus_review_sheet = build_natus_review_sheet(seizures_events, iid_bursts, gaps_summary, file_subject_map, tz);

    if ~isfolder(merged_output_dir)
        mkdir(merged_output_dir);
    end
    summaries_out = fullfile(merged_output_dir, '05_summaries');
    if ~isfolder(summaries_out)
        mkdir(summaries_out);
    end

    paths = write_all_summaries(summaries_out, seizures_events, seizures_summary, ...
        iid_events, iid_summary, iid_bursts, gaps_summary, qc_report, natus_review_sheet);

    fprintf('merge_pipeline_runs: %d run(s) -> %d seizure(s), %d IID complex(es), %d burst(s)\n', ...
        numel(output_roots), height(seizures_events), height(iid_events), height(iid_bursts));

    result = struct();
    result.seizures_events = seizures_events;
    result.seizures_summary = seizures_summary;
    result.iid_events = iid_events;
    result.iid_summary = iid_summary;
    result.iid_bursts = iid_bursts;
    result.gaps_summary = gaps_summary;
    result.qc_report = qc_report;
    result.natus_review_sheet = natus_review_sheet;
    result.paths = paths;
end
