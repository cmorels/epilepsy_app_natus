function T = vertcat_or_empty(parts, empty_fn)
% VERTCAT_OR_EMPTY  vertcat(parts{:}), or empty_fn() if parts is empty.
% empty_fn takes no arguments (wrap with an anonymous function to pass a
% timezone or other context, e.g. @() empty_seizure_events_table(tz)).
    if isempty(parts)
        T = empty_fn();
    else
        T = vertcat(parts{:});
    end
end
