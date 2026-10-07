function T = add_seizure_datetime_columns(T)
% ADD_SEIZURE_DATETIME_COLUMNS  Adds start_datetime and end_datetime to a
% seizures_events table, right after seizure_id: the seizure's start / end
% as local clock time, 'dd/MM/yyyy HH:mm:ss' (e.g. 29/03/2026 17:35:45),
% the form used to look an event up in Natus.
%
%   T = add_seizure_datetime_columns(T)
%
% Computed from start_abs / end_abs (already in the recording's local time,
% DST included; seconds truncated, no milliseconds). Stored as TEXT so
% Excel shows it exactly like that instead of reformatting it. Existing
% columns are never changed; if the two columns are already there they are
% recomputed in place. Used by write_all_summaries.m (every run and every
% merge) and tools/add_datetime_columns_to_seizures_events.m (outputs
% written before 2026-10-04). read_pipeline_csv.m ignores both columns.

    fmt = 'dd/MM/yyyy HH:mm:ss';
    names = {'start_datetime', 'end_datetime'};
    sources = {'start_abs', 'end_abs'};
    if width(T) == 0
        return;
    end
    for i = 1:numel(names)
        if ismember(names{i}, T.Properties.VariableNames)
            T.(names{i}) = [];
        end
    end

    n = height(T);
    vals = cell(1, numel(names));
    for i = 1:numel(names)
        if ismember(sources{i}, T.Properties.VariableNames)
            vals{i} = to_text(T.(sources{i}), fmt, n);
        else
            vals{i} = repmat({''}, n, 1);
        end
    end

    if ismember('seizure_id', T.Properties.VariableNames)
        anchor = 'seizure_id';
    else
        anchor = T.Properties.VariableNames{end};
    end
    T = addvars(T, vals{1}, vals{2}, 'NewVariableNames', names, 'After', anchor);
end

function c = to_text(v, fmt, n)
    try
        if ~isdatetime(v)
            % text fallback of robust_vertcat.m: the writer's own format
            v = datetime(v, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS', 'Locale', 'en_US');
        end
        s = string(v, fmt);
        s(ismissing(s)) = "";
        c = cellstr(s(:));
    catch
        c = repmat({''}, n, 1);
    end
    if numel(c) ~= n
        c = repmat({''}, n, 1);
    end
end
