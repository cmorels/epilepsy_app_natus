function T = robust_vertcat(parts, empty_fn, tz, label)
% ROBUST_VERTCAT  vertcat(parts{:}) that never throws, for the final
% consolidation of a run (run_pipeline_edf.m, merge_pipeline_runs.m).
%
%   T = robust_vertcat(parts, empty_fn, tz, label)
%
% 1. Every datetime column of every part is stamped with TimeZone tz
%    (apply_tz.m), so a part carrying an unzoned NaT (e.g. the qc row of an
%    EDF that failed to import) can no longer break the concatenation --
%    the failure that used to throw away a whole subject's results.
% 2. If vertcat still fails (any other column type/schema mismatch), the
%    parts are concatenated AS TEXT instead: union of all columns, every
%    value converted to text, missing columns left empty. The CSV written
%    from it holds the same information; a warning names the table.
% empty_fn() is returned for no parts (as vertcat_or_empty.m). label only
% names the table in the warning.

    if nargin < 4
        label = 'table';
    end
    if isempty(parts)
        T = empty_fn();
        return;
    end
    for i = 1:numel(parts)
        try
            parts{i} = apply_tz(parts{i}, tz);
        catch
            % leave this part as is; the text fallback below still handles it
        end
    end
    try
        T = vertcat(parts{:});
        return;
    catch ME
        warning('robust_vertcat:TextFallback', '%s: %s -- concatenated as text instead.', label, ME.message);
    end
    T = vertcat_as_text(parts);
end

function T = vertcat_as_text(parts)
    names = {};
    for i = 1:numel(parts)
        if istable(parts{i})
            names = [names, setdiff(parts{i}.Properties.VariableNames, names, 'stable')]; %#ok<AGROW>
        end
    end
    cols = cell(1, numel(names));
    for j = 1:numel(names)
        col = cell(0, 1);
        for i = 1:numel(parts)
            p = parts{i};
            if ~istable(p)
                continue;
            end
            if ismember(names{j}, p.Properties.VariableNames)
                col = [col; to_text(p.(names{j}), height(p))]; %#ok<AGROW>
            else
                col = [col; repmat({''}, height(p), 1)]; %#ok<AGROW>
            end
        end
        cols{j} = col;
    end
    T = table(cols{:}, 'VariableNames', names);
end

function c = to_text(v, n)
    try
        if isdatetime(v)
            v.Format = 'dd-MMM-yyyy HH:mm:ss.SSS';
        end
        s = string(v);
        if size(s, 2) > 1
            s = join(s, ' ', 2);
        end
        s(ismissing(s)) = "";
        c = cellstr(s(:));
        if numel(c) ~= n
            c = repmat({''}, n, 1);
        end
    catch
        c = repmat({'?'}, n, 1);
    end
end
