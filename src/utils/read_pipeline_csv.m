function T = read_pipeline_csv(csv_path, kind, tz)
% READ_PIPELINE_CSV  Safe reader for this pipeline's own summary CSVs.
%
%   T = read_pipeline_csv(csv_path, kind, tz)
%
%   kind : one of 'seizures_events', 'seizures_summary', 'iid_events',
%          'iid_summary', 'iid_bursts', 'gaps', 'qc' -- selects the
%          expected schema (same empty_*_table.m functions run_pipeline_edf.m
%          uses to write these files, so reader and writer never drift).
%   tz   : TimeZone to apply to datetime columns (default 'Europe/Paris';
%          plain CSV text carries no timezone info, see README.md).
%
% Plain readtable() is NOT safe on these files: MATLAB auto-detects a
% column like subject_id="097" as the number 97, and a 0-row CSV (a file
% with e.g. zero seizures) comes back with every column typed 'double',
% including datetime ones. This function forces every column to the type
% the pipeline actually uses, regardless of row count.

    if nargin < 3 || isempty(tz)
        tz = 'Europe/Paris';
    end

    template = empty_table_for(kind, tz);
    names = template.Properties.VariableNames;

    opts = detectImportOptions(csv_path);
    missing = setdiff(names, opts.VariableNames);
    optional_defaults = optional_columns_for(kind);
    required_missing = setdiff(missing, keys(optional_defaults));
    if ~isempty(required_missing)
        error('read_pipeline_csv:MissingColumn', ...
            'Expected column(s) not found in %s: %s', csv_path, strjoin(required_missing, ', '));
    end
    % `missing` (a strict subset with a registered default, if non-empty here) covers
    % files written before a schema extension (e.g. seizures_events.csv before the
    % seizure_mode/confidence columns) -- see optional_columns_for below.
    present_names = setdiff(names, missing, 'stable');

    text_cols = {};
    for i = 1:numel(present_names)
        name = present_names{i};
        col = template.(name);
        if iscell(col) || isdatetime(col)
            opts = setvartype(opts, name, 'char');  % read as text; convert explicitly below (sidesteps readtable's 0-row datetime bug)
            text_cols{end+1} = name; %#ok<AGROW>
        elseif islogical(col)
            opts = setvartype(opts, name, 'logical');
        else
            opts = setvartype(opts, name, 'double');
        end
    end

    Traw = readtable(csv_path, opts);
    Traw = Traw(:, present_names);
    n = height(Traw);

    cols = cell(1, numel(names));
    for i = 1:numel(names)
        name = names{i};
        if ismember(name, missing)
            cols{i} = repmat(optional_defaults(name), n, 1);
            continue;
        end
        if isdatetime(template.(name))
            if n == 0
                cols{i} = datetime.empty(0, 1);
            else
                cols{i} = datetime(Traw.(name), 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS');
            end
            cols{i}.TimeZone = tz;
        else
            cols{i} = Traw.(name);
        end
    end

    T = table(cols{:}, 'VariableNames', names);
end

%% ======================================================================
function m = optional_columns_for(kind)
% Columns a schema has grown since some existing summary CSVs were
% written, with the default value to backfill when reading an
% older-shaped file (never an error -- see the missing-column check
% above). Every file predating a given column was necessarily written
% by whatever this default represents (e.g. every seizures_events.csv
% written before seizure_mode existed came only from the legacy branch).
    m = containers.Map('KeyType', 'char', 'ValueType', 'any');
    if strcmp(kind, 'seizures_events')
        m('seizure_mode') = {'legacy'};
        m('over_max_duration') = false;
        m('ll_ratio') = NaN;
        m('peak_energy_ratio') = NaN;
        m('hf_ratio_db') = NaN;
        m('envelope_cv') = NaN;
    end
end

function T = empty_table_for(kind, tz)
    switch kind
        case 'seizures_events'
            T = empty_seizure_events_table(tz);
        case 'seizures_summary'
            T = empty_seizure_summary_table(tz);
        case 'iid_events'
            T = empty_iid_events_table(tz);
        case 'iid_summary'
            T = empty_iid_summary_table(tz);
        case 'iid_bursts'
            T = empty_iid_bursts_table(tz);
        case 'gaps'
            T = empty_gaps_table(tz);
        case 'qc'
            T = empty_qc_table(tz);
        otherwise
            error('read_pipeline_csv:BadKind', ...
                'kind must be one of seizures_events, seizures_summary, iid_events, iid_summary, iid_bursts, gaps, qc (got ''%s'').', kind);
    end
end
