function cases = load_cases(csv_path, cfg)
% LOAD_CASES  Read and validate a case-assignment CSV.
%
%   cases = load_cases(csv_path, cfg)
%
% Minimal required format (see README.md "sistema de casos" for the full
% write-up and resolve_case.m for how a row is matched and applied):
%
%   subject_id,source_file,region,case,gain,notch,notes
%   39920,39920_20250505_HPCright.txt,,attenuated,,,bajo a 10 uV/mm en Natus
%   39921,39921_20250612_HPCleft.txt,HPCleft,both,20,,5 uV/mm + filtro de red
%   39922,39922_20250613_HPCright.txt,,line,,on,zumbido desde las 02:00
%
% Only 'source_file' is required; every other expected column
% (subject_id, region, case, gain, notch, notes, suggested_case) is added
% as empty/NaN if the file doesn't have it -- a hand-written minimal CSV
% and the full CSV make_cases_template.m produces (which also carries
% suggested_case and diagnostic columns) both load into the same shape.
%
% Validates every non-empty 'case' value against cfg.cases.profiles right
% away (fail fast, one clear error naming the bad row's source_file and
% listing the valid case names) rather than only failing when that
% specific row happens to get resolved during a run.
%
% Does NOT check that every row matches a file that actually gets
% processed -- run_pipeline_edf.m does that after the run, since it's the
% only place that knows what was actually processed.

    if exist(csv_path, 'file') ~= 2
        error('load_cases:FileNotFound', 'Cases CSV not found: %s', csv_path);
    end

    % VariableNamingRule must be passed to detectImportOptions itself, not
    % set on the returned object afterward (verified empirically -- the
    % latter does not retroactively un-rename anything): 'case' is a
    % MATLAB reserved word, and without this, detectImportOptions renames
    % that column to 'xCase' at detection time, before there's any object
    % left to fix up.
    opts = detectImportOptions(csv_path, 'VariableNamingRule', 'preserve');
    if ~ismember('source_file', opts.VariableNames)
        error('load_cases:MissingColumn', 'Cases CSV %s has no "source_file" column.', csv_path);
    end

    text_cols = {'subject_id', 'source_file', 'region', 'case', 'notch', 'notes', 'suggested_case'};
    for i = 1:numel(text_cols)
        if ismember(text_cols{i}, opts.VariableNames)
            opts = setvartype(opts, text_cols{i}, 'char');
        end
    end
    if ismember('gain', opts.VariableNames)
        opts = setvartype(opts, 'gain', 'double');
    end

    T = readtable(csv_path, opts);

    all_cols = [text_cols, {'gain'}];
    n = height(T);
    for i = 1:numel(all_cols)
        col = all_cols{i};
        if ~ismember(col, T.Properties.VariableNames)
            if strcmp(col, 'gain')
                T.(col) = nan(n, 1);
            else
                T.(col) = repmat({''}, n, 1);
            end
        end
    end

    % readtable can hand back an empty text cell as {0x0 char} instead of
    % {''}; normalize so downstream isempty()/strcmpi() checks are uniform.
    for i = 1:numel(text_cols)
        col = text_cols{i};
        T.(col) = cellfun(@(v) char(v), T.(col), 'UniformOutput', false);
    end

    valid_cases = fieldnames(cfg.cases.profiles);
    for i = 1:n
        c = strtrim(T.('case'){i});
        if ~isempty(c) && ~ismember(c, valid_cases)
            error('load_cases:BadCaseName', ...
                'Row %d (source_file="%s"): case "%s" is not defined in cfg.cases.profiles. Valid cases: %s', ...
                i, T.source_file{i}, c, strjoin(valid_cases, ', '));
        end
    end

    cases = T;
end
