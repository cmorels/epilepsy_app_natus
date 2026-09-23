function run_all_subjects(data_root, output_root, subjects, cases_file, recording_log)
% RUN_ALL_SUBJECTS  Run run_pipeline_edf.m once per subject folder and
% merge every run into one set of summaries.
%
%   run_all_subjects()                                  % defaults below
%   run_all_subjects(data_root, output_root)
%   run_all_subjects(data_root, output_root, {'001-s','005-s'})
%   run_all_subjects(data_root, output_root, {}, 'cases.csv')
%   run_all_subjects(data_root, output_root, {}, '', 'EEG_recording_log.xlsx')
%
% data_root   folder holding one subfolder per subject ('001-s', '005-s', ...)
%             containing .EDF files. Default: parent folder of this repo.
% output_root where results go: <output_root>/<subject>/... per subject and
%             <output_root>/merged/05_summaries for the combined CSVs.
%             Default: <data_root>/pipeline_output.
% subjects    cell array of subfolder names to process; {} = all subfolders
%             matching '*-s' that contain at least one .edf/.EDF.
% cases_file  optional cases CSV (see src/load_cases.m); '' = none.
% recording_log  xlsx giving, per (animal, EDF), the Port and the
%             HPCr/HPCl channels (cfg.edf.channels.mode = 'log').
%             Default: <data_root>/EEG_recording_log.xlsx.
%
% subject_id is taken from the folder name ('001-s' -> '001'), because the
% Natus EDF filenames do not carry it. A subject that fails is logged and
% the loop continues with the next one; so is a single EDF missing from the
% recording log or whose logged channels are not in the file (see each
% subject's 05_summaries/qc_report.csv).
%
% From a terminal (see run_all_subjects.ps1):
%   matlab -batch "run_all_subjects"

    repo_root = fileparts(mfilename('fullpath'));
    addpath(fullfile(repo_root, 'src'));
    addpath(fullfile(repo_root, 'src', 'utils'));

    if nargin < 1 || isempty(data_root),   data_root = fileparts(repo_root); end
    if nargin < 2 || isempty(output_root), output_root = fullfile(data_root, 'pipeline_output'); end
    if nargin < 3 || isempty(subjects),    subjects = find_subject_folders(data_root); end
    if nargin < 4, cases_file = ''; end
    if nargin < 5 || isempty(recording_log), recording_log = fullfile(data_root, 'EEG_recording_log.xlsx'); end
    if ischar(subjects) || isstring(subjects), subjects = cellstr(subjects); end
    if exist(recording_log, 'file') ~= 2
        error('run_all_subjects:NoRecordingLog', 'Recording log not found: %s', recording_log);
    end

    fprintf('Data root:   %s\nOutput root: %s\nChannel log: %s\nSubjects:    %d\n\n', ...
        data_root, output_root, recording_log, numel(subjects));

    done_dirs = {};
    failed = {};
    t_all = tic;
    for k = 1:numel(subjects)
        name = subjects{k};
        in_dir = fullfile(data_root, name);
        out_dir = fullfile(output_root, name);
        fprintf('[%d/%d] %s ...\n', k, numel(subjects), name);
        t = tic;
        try
            cfg = pipeline_config();
            cfg.edf.subject_id = regexprep(name, '-s$', '');
            cfg.paths.output_root = out_dir;
            cfg.edf.channels.mode = 'log';
            cfg.edf.channels.log_file = recording_log;
            if ~isempty(cases_file)
                cfg.cases.file = cases_file;
            end
            run_pipeline_edf(in_dir, cfg);
            done_dirs{end+1} = out_dir; %#ok<AGROW>
            fprintf('[%d/%d] %s OK (%.1f min)\n\n', k, numel(subjects), name, toc(t)/60);
        catch err
            failed{end+1} = sprintf('%s: %s', name, err.message); %#ok<AGROW>
            fprintf(2, '[%d/%d] %s FAILED: %s\n\n', k, numel(subjects), name, err.message);
        end
    end

    if ~isempty(done_dirs)
        fprintf('Merging %d run(s) into %s ...\n', numel(done_dirs), fullfile(output_root, 'merged'));
        try
            merge_pipeline_runs(done_dirs, fullfile(output_root, 'merged'));
        catch err
            failed{end+1} = sprintf('merge: %s', err.message);
            fprintf(2, 'Merge FAILED: %s\n', err.message);
        end
    end

    fprintf('\nFinished in %.1f min: %d OK, %d failed.\n', toc(t_all)/60, numel(done_dirs), numel(failed));
    for k = 1:numel(failed)
        fprintf(2, '  %s\n', failed{k});
    end
end

function names = find_subject_folders(data_root)
    d = dir(fullfile(data_root, '*-s'));
    d = d([d.isdir]);
    names = {};
    for k = 1:numel(d)
        f = fullfile(data_root, d(k).name);
        if ~isempty(dir(fullfile(f, '*.edf'))) || ~isempty(dir(fullfile(f, '*.EDF')))
            names{end+1} = d(k).name; %#ok<AGROW>
        else
            fprintf('Skipping %s (no .EDF files)\n', d(k).name);
        end
    end
end
