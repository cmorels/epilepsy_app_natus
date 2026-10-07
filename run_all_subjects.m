function run_all_subjects(data_root, output_root, subjects, cases_file, recording_log, files_csv)
% RUN_ALL_SUBJECTS  Run the EDF campaign (src/campaign_config.m): one
% run_pipeline_edf.m call per animal, then merge every run into one set of
% summaries.
%
%   run_all_subjects()                                   % every animal folder
%   run_all_subjects(data_root, output_root)
%   run_all_subjects(data_root, output_root, {'004-s','005-s'})
%   run_all_subjects(data_root, output_root, {}, '', '', 'edf_list.csv')   % only the listed EDFs
%
% data_root     folder holding one subfolder per animal containing .EDF
%               files. Default: parent folder of this repo.
% output_root   <output_root>/<animal>/... per animal,
%               <output_root>/merged/05_summaries for the combined CSVs,
%               <output_root>/reference_*.csv for the gain reference.
%               Default: <data_root>/pipeline_output_campaign_<yyyyMMdd>
%               (pipeline_output_campaign_preliminar_<yyyyMMdd> with a
%               files_csv).
% subjects      cell array of animal folder names to process; {} = all.
% cases_file    optional cases CSV (src/load_cases.m); '' = none. When
%               given it REPLACES the Excel-driven case (cfg.cases.from_excel
%               = false).
% recording_log the Excel (EEG_recording_log.xlsx) giving, per (animal,
%               EDF), Port, HPCr/HPCl channels and Attenuation.
%               Default: <data_root>/EEG_recording_log.xlsx.
% files_csv     optional CSV listing exactly which EDFs to process, with
%               columns animal_id, edf_file (the EDF file name, with or
%               without .EDF). '' = every EDF of every animal folder.
%
% Animal ID = folder name, used exactly as written ('005-s' and '005' are
% different animals). With no files_csv, a subfolder is an animal folder
% when its name appears exactly in the Excel's Animal ID column and it
% contains .EDF files; anything else (pipeline_output*, the repo, ...) is
% ignored.
%
% Steps: (1) list the EDFs; (2) compute_run_reference.m imports them all
% once and measures the gain reference from the NON-attenuated channels of
% this run; (3) run_pipeline_edf.m per animal (reusing that import);
% (4) merge_pipeline_runs.m. A listed EDF that does not exist, an animal
% or a file that fails, is logged and the run continues.
%
% From a terminal (see run_all_subjects.ps1):
%   matlab -batch "run_all_subjects"

    repo_root = fileparts(mfilename('fullpath'));
    addpath(fullfile(repo_root, 'src'));
    addpath(fullfile(repo_root, 'src', 'utils'));

    if nargin < 1 || isempty(data_root),   data_root = fileparts(repo_root); end
    if nargin < 6, files_csv = ''; end
    if nargin < 2 || isempty(output_root)
        if isempty(files_csv)
            output_root = fullfile(data_root, ['pipeline_output_campaign_' datestr(now, 'yyyymmdd')]); %#ok<TNOW1,DATST>
        else
            output_root = fullfile(data_root, ['pipeline_output_campaign_preliminar_' datestr(now, 'yyyymmdd')]); %#ok<TNOW1,DATST>
        end
    end
    if nargin < 3, subjects = {}; end
    if nargin < 4, cases_file = ''; end
    if nargin < 5 || isempty(recording_log), recording_log = fullfile(data_root, 'EEG_recording_log.xlsx'); end
    if ischar(subjects) || isstring(subjects), subjects = cellstr(subjects); end
    if exist(recording_log, 'file') ~= 2
        error('run_all_subjects:NoRecordingLog', 'Excel not found: %s', recording_log);
    end
    if ~isfolder(output_root)
        mkdir(output_root);
    end

    cfg_base = campaign_config();
    cfg_base.edf.channels.log_file = recording_log;
    if ~isempty(cases_file)
        cfg_base.cases.from_excel = false;
        cfg_base.cases.file = cases_file;
    end

    excel = load_recording_log(recording_log);
    excel_ids = unique(excel.animal_id);

    %% (1) EDFs of this run
    if isempty(files_csv)
        if isempty(subjects)
            subjects = find_subject_folders(data_root, excel_ids);
        end
        [jobs, skipped] = jobs_from_folders(data_root, subjects, output_root);
    else
        [jobs, skipped] = jobs_from_csv(files_csv, data_root, output_root);
        if ~isempty(subjects)
            jobs = jobs(ismember({jobs.subject_id}, subjects));
        end
    end
    animals = unique({jobs.subject_id}, 'stable');

    fprintf('Data root:   %s\nOutput root: %s\nExcel:       %s\nEDF list:    %s\nAnimals:     %d\nEDF files:   %d\n\n', ...
        data_root, output_root, recording_log, pick(files_csv, '(every EDF of every animal folder)'), ...
        numel(animals), numel(jobs));
    write_skipped(skipped, output_root);

    t_all = tic;

    %% (2) gain reference from the non-attenuated channels of this run
    if ~isempty(jobs) && cfg_base.cases.from_excel
        fprintf('Step 1/3: importing every EDF and measuring the gain reference ...\n');
        ref = compute_run_reference(jobs, cfg_base, output_root);
        if ref.map.Count > 0
            cfg_base.quality.reference_sigma_uV = ref.map;
        end
        cfg_base.quality.reference_sigma_fallback_uV = ref.fallback;
        fprintf('\n');
    end

    %% (3) one run per animal
    fprintf('Step 2/3: processing %d animal(s) ...\n\n', numel(animals));
    done_dirs = {};
    failed = {};
    for k = 1:numel(animals)
        name = animals{k};
        mine = jobs(strcmp({jobs.subject_id}, name));
        out_dir = mine(1).output_root;
        fprintf('[%d/%d] %s (%d EDF) ...\n', k, numel(animals), name, numel(mine));
        t = tic;
        try
            cfg = cfg_base;
            cfg.edf.subject_id = name;
            cfg.paths.output_root = out_dir;
            run_pipeline_edf({mine.edf_path}, cfg);
            done_dirs{end+1} = out_dir; %#ok<AGROW>
            fprintf('[%d/%d] %s OK (%.1f min)\n\n', k, numel(animals), name, toc(t)/60);
        catch err
            failed{end+1} = sprintf('%s: %s', name, err.message); %#ok<AGROW>
            fprintf(2, '[%d/%d] %s FAILED: %s\n\n', k, numel(animals), name, err.message);
        end
    end

    %% (4) merge
    if ~isempty(done_dirs)
        fprintf('Step 3/3: merging %d run(s) into %s ...\n', numel(done_dirs), fullfile(output_root, 'merged'));
        try
            merge_pipeline_runs(done_dirs, fullfile(output_root, 'merged'), cfg_base.general.timezone);
        catch err
            failed{end+1} = sprintf('merge: %s', err.message);
            fprintf(2, 'Merge FAILED: %s\n', err.message);
        end
    end

    fprintf('\nFinished in %.1f min: %d animal(s) OK, %d failed, %d listed EDF(s) skipped.\n', ...
        toc(t_all)/60, numel(done_dirs), numel(failed), numel(skipped));
    for k = 1:numel(failed)
        fprintf(2, '  %s\n', failed{k});
    end
end

%% ======================================================================
function names = find_subject_folders(data_root, excel_ids)
% Animal folder = subfolder whose name is EXACTLY an Animal ID of the Excel
% and that contains .EDF files.
    d = dir(data_root);
    d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
    names = {};
    for k = 1:numel(d)
        if ~ismember(d(k).name, excel_ids)
            continue;
        end
        f = fullfile(data_root, d(k).name);
        if ~isempty(dir(fullfile(f, '*.edf'))) || ~isempty(dir(fullfile(f, '*.EDF')))
            names{end+1} = d(k).name; %#ok<AGROW>
        else
            fprintf('Skipping %s (no .EDF files)\n', d(k).name);
        end
    end
end

function [jobs, skipped] = jobs_from_folders(data_root, subjects, output_root)
    jobs = empty_jobs();
    skipped = empty_skipped();
    for k = 1:numel(subjects)
        folder = fullfile(data_root, subjects{k});
        files = [dir(fullfile(folder, '*.edf')); dir(fullfile(folder, '*.EDF'))];
        if isempty(files)
            skipped(end+1) = struct('animal_id', subjects{k}, 'edf_file', '', 'reason', 'animal folder not found or without .EDF files'); %#ok<AGROW>
            continue;
        end
        paths = unique(fullfile({files.folder}, {files.name}), 'stable');
        for i = 1:numel(paths)
            jobs(end+1) = struct('subject_id', subjects{k}, 'edf_path', paths{i}, ...
                'output_root', fullfile(output_root, subjects{k})); %#ok<AGROW>
        end
    end
end

function [jobs, skipped] = jobs_from_csv(files_csv, data_root, output_root)
% Reads files_csv (columns animal_id, edf_file; everything as text so an
% ID like '005' is never turned into 5). A row whose EDF does not exist is
% skipped and reported.
    if exist(files_csv, 'file') ~= 2
        error('run_all_subjects:NoFilesCsv', 'EDF list not found: %s', files_csv);
    end
    opts = detectImportOptions(files_csv, 'VariableNamingRule', 'preserve');  % delimiter (, ; tab) detected automatically
    need = {'animal_id', 'edf_file'};
    missing = setdiff(need, opts.VariableNames);
    if ~isempty(missing)
        error('run_all_subjects:BadFilesCsv', 'EDF list %s is missing column(s): %s (expected: animal_id, edf_file)', ...
            files_csv, strjoin(missing, ', '));
    end
    opts.SelectedVariableNames = need;
    opts = setvartype(opts, need, 'char');
    T = readtable(files_csv, opts);

    jobs = empty_jobs();
    skipped = empty_skipped();
    for i = 1:height(T)
        animal = strtrim(T.animal_id{i});
        edf = strtrim(T.edf_file{i});
        if isempty(animal) && isempty(edf)
            continue;
        end
        path = resolve_edf_path(fullfile(data_root, animal), edf);
        if isempty(path)
            skipped(end+1) = struct('animal_id', animal, 'edf_file', edf, 'reason', 'EDF not found in the animal folder'); %#ok<AGROW>
            fprintf(2, 'Skipping %s | %s: EDF not found in %s\n', animal, edf, fullfile(data_root, animal));
            continue;
        end
        if any(strcmp({jobs.subject_id}, animal) & strcmp({jobs.edf_path}, path))
            continue;  % listed twice
        end
        jobs(end+1) = struct('subject_id', animal, 'edf_path', path, 'output_root', fullfile(output_root, animal)); %#ok<AGROW>
    end
end

function path = resolve_edf_path(folder, edf)
    path = '';
    [~, ~, ext] = fileparts(edf);
    if strcmpi(ext, '.edf')
        candidates = {edf};
    else
        candidates = {[edf '.EDF'], [edf '.edf']};
    end
    for c = 1:numel(candidates)
        p = fullfile(folder, candidates{c});
        if exist(p, 'file') == 2
            path = p;
            return;
        end
    end
end

function jobs = empty_jobs()
    jobs = struct('subject_id', {}, 'edf_path', {}, 'output_root', {});
end

function s = empty_skipped()
    s = struct('animal_id', {}, 'edf_file', {}, 'reason', {});
end

function write_skipped(skipped, output_root)
    if isempty(skipped)
        return;
    end
    try
        T = table({skipped.animal_id}', {skipped.edf_file}', {skipped.reason}', ...
            'VariableNames', {'animal_id', 'edf_file', 'reason'});
        writetable(T, fullfile(output_root, 'skipped_inputs.csv'));
        fprintf('%d listed EDF(s) skipped -> %s\n\n', height(T), fullfile(output_root, 'skipped_inputs.csv'));
    catch ME
        warning('run_all_subjects:SkippedWrite', 'Could not write skipped_inputs.csv: %s', ME.message);
    end
end

function v = pick(a, b)
    if isempty(a)
        v = b;
    else
        v = a;
    end
end
