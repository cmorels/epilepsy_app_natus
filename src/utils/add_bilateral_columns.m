function T = add_bilateral_columns(T, kind)
% ADD_BILATERAL_COLUMNS  Append any bilateral_columns(kind) column that T
% lacks, filled with that column's default (see bilateral_columns.m).
% Used for 0-row templates and to backfill never-reconciled tables so they
% can be vertcat-ed with reconciled ones.
    [names, types, defaults] = bilateral_columns(kind);
    n = height(T);
    for i = 1:numel(names)
        if ismember(names{i}, T.Properties.VariableNames)
            continue;
        end
        if n == 0
            switch types{i}
                case 'cell',    T.(names{i}) = cell(0, 1);
                case 'logical', T.(names{i}) = false(0, 1);
                otherwise,      T.(names{i}) = zeros(0, 1);
            end
        else
            T.(names{i}) = repmat(defaults{i}, n, 1);
        end
    end
end
