function paths = write_all_summaries(summaries_dir, seizures_events, seizures_summary, ...
        iid_events, iid_summary, iid_bursts, gaps_summary, qc_report, natus_review_sheet)
% WRITE_ALL_SUMMARIES  Writes the 8 CSVs + pipeline_summary.xlsx that make
% up 05_summaries/. Used by both run_pipeline_edf.m (single run) and
% merge_pipeline_runs.m (several runs combined) so the two never drift
% into different output formats.
%
% Every datetime column is given an explicit millisecond-precision Format
% before writing -- MATLAB's default text format for datetime
% (dd-MMM-yyyy HH:mm:ss) silently truncates to whole seconds, which would
% otherwise throw away real sub-second precision that start_s/end_s (the
% authoritative relative-time columns) still carry.

    seizures_events = set_dt_format(seizures_events);
    seizures_summary = set_dt_format(seizures_summary);
    iid_events = set_dt_format(iid_events);
    iid_summary = set_dt_format(iid_summary);
    iid_bursts = set_dt_format(iid_bursts);
    gaps_summary = set_dt_format(gaps_summary);
    qc_report = set_dt_format(qc_report);
    natus_review_sheet = set_dt_format(natus_review_sheet);

    paths = struct();
    paths.seizures_events = fullfile(summaries_dir, 'seizures_events.csv');
    paths.seizures_summary = fullfile(summaries_dir, 'seizures_summary.csv');
    paths.iid_events = fullfile(summaries_dir, 'iid_events.csv');
    paths.iid_summary = fullfile(summaries_dir, 'iid_summary.csv');
    paths.iid_bursts = fullfile(summaries_dir, 'iid_bursts.csv');
    paths.gaps_summary = fullfile(summaries_dir, 'gaps_summary.csv');
    paths.qc_report = fullfile(summaries_dir, 'qc_report.csv');
    paths.natus_review_sheet = fullfile(summaries_dir, 'natus_review_sheet.csv');
    paths.pipeline_summary_xlsx = fullfile(summaries_dir, 'pipeline_summary.xlsx');

    writetable(seizures_events, paths.seizures_events);
    writetable(seizures_summary, paths.seizures_summary);
    writetable(iid_events, paths.iid_events);
    writetable(iid_summary, paths.iid_summary);
    writetable(iid_bursts, paths.iid_bursts);
    writetable(gaps_summary, paths.gaps_summary);
    writetable(qc_report, paths.qc_report);
    writetable(natus_review_sheet, paths.natus_review_sheet);

    if exist(paths.pipeline_summary_xlsx, 'file') == 2
        delete(paths.pipeline_summary_xlsx);
    end
    writetable(seizures_events, paths.pipeline_summary_xlsx, 'Sheet', 'seizures_events');
    writetable(seizures_summary, paths.pipeline_summary_xlsx, 'Sheet', 'seizures_summary');
    writetable(iid_events, paths.pipeline_summary_xlsx, 'Sheet', 'iid_events');
    writetable(iid_summary, paths.pipeline_summary_xlsx, 'Sheet', 'iid_summary');
    writetable(iid_bursts, paths.pipeline_summary_xlsx, 'Sheet', 'iid_bursts');
    writetable(gaps_summary, paths.pipeline_summary_xlsx, 'Sheet', 'gaps');
    writetable(qc_report, paths.pipeline_summary_xlsx, 'Sheet', 'qc');
end

function T = set_dt_format(T)
    for name = T.Properties.VariableNames
        if isdatetime(T.(name{1}))
            T.(name{1}).Format = 'dd-MMM-yyyy HH:mm:ss.SSS';
        end
    end
end
