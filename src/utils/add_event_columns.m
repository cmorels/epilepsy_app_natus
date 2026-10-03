function T = add_event_columns(T, kind)
% ADD_EVENT_COLUMNS  Append any event_columns(kind) column T lacks, filled
% with its default (see event_columns.m). For 0-row templates and to make
% tables from different files/runs stackable.
    [names, types, defaults] = event_columns(kind);
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
