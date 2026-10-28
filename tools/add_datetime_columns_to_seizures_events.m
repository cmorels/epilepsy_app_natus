function T = add_datetime_columns_to_seizures_events(output_root, tz)
% ADD_DATETIME_COLUMNS_TO_SEIZURES_EVENTS  Adds start_datetime and
% end_datetime (right after seizure_id, 'dd/MM/yyyy HH:mm:ss' local time,
% e.g. 29/03/2026 17:35:45) to seizures_events of outputs written before
% write_all_summaries.m did it itself (2026-10-04). Nothing is re-processed.
%
%   T = add_datetime_columns_to_seizures_events(output_root)
%   T = add_datetime_columns_to_seizures_events(output_root, tz)    % default 'Europe/Paris'
%
% output_root : a campaign output (one subfolder per animal + merged/) or a
%               single run_pipeline_edf.m output_root. Every 05_summaries/
%               found is updated: seizures_events.csv and the
%               seizures_events sheet of pipeline_summary.xlsx. No other
%               file, sheet or column is touched; existing columns are
%               written back exactly as the pipeline writes them.
% tz          : TimeZone the run used (cfg.general.timezone).
%
% Safe to run twice (the two columns are recomputed, never duplicated).
% A folder that fails (e.g. the xlsx is open in Excel) is reported and the
% others are still updated.
%
% OUTPUT T: table folder, n_rows, csv, xlsx (status of each write).

    if nargin < 2 || isempty(tz)
        tz = 'Europe/Paris';
    end
    repo = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(repo, 'src'), fullfile(repo, 'src', 'utils'));

    dirs = summaries_dirs(output_root);
    rows = cell(0, 4);
    for i = 1:numel(dirs)
        d = dirs{i};
        csv = fullfile(d, 'seizures_events.csv');
        xlsx = fullfile(d, 'pipeline_summary.xlsx');
        if exist(csv, 'file') ~= 2
            continue;
        end
        try
            S = read_pipeline_csv(csv, 'seizures_events', tz);
            S = add_seizure_datetime_columns(S);
            S = set_dt_format(S);
        catch ME
            rows(end+1, :) = {d, NaN, ['read failed: ' ME.message], 'not touched'}; %#ok<AGROW>
            continue;
        end
        try
            writetable(S, csv);
            csv_status = 'updated';
        catch ME
            csv_status = ['write failed: ' ME.message];
        end
        if exist(xlsx, 'file') == 2
            try
                writetable(S, xlsx, 'Sheet', 'seizures_events', 'WriteMode', 'overwritesheet');
                xlsx_status = 'updated';
            catch ME
                xlsx_status = ['write failed: ' ME.message];
            end
        else
            xlsx_status = 'no xlsx';
        end
        rows(end+1, :) = {d, height(S), csv_status, xlsx_status}; %#ok<AGROW>
        fprintf('%s: %d row(s) | csv %s | xlsx %s\n', d, height(S), csv_status, xlsx_status);
    end
    T = cell2table(rows, 'VariableNames', {'folder', 'n_rows', 'csv', 'xlsx'});
    fprintf('add_datetime_columns_to_seizures_events: %d folder(s) with seizures_events.csv\n', height(T));
end

function dirs = summaries_dirs(output_root)
    dirs = {};
    if isfolder(fullfile(output_root, '05_summaries'))
        dirs{end+1} = fullfile(output_root, '05_summaries');
    end
    d = dir(output_root);
    d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
    for k = 1:numel(d)
        p = fullfile(output_root, d(k).name, '05_summaries');
        if isfolder(p)
            dirs{end+1} = p; %#ok<AGROW>
        end
    end
end

function T = set_dt_format(T)
% Same as write_all_summaries.m, so start_abs etc. are written unchanged.
    for name = T.Properties.VariableNames
        if isdatetime(T.(name{1}))
            T.(name{1}).Format = 'dd-MMM-yyyy HH:mm:ss.SSS';
        end
    end
end
