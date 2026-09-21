function abs_dt = rel_to_abs_time(t_rel_s, session_start)
% REL_TO_ABS_TIME  Convert time relative to session start (seconds) into
% absolute datetime, honoring session_start's explicit TimeZone.
%
%   abs_dt = rel_to_abs_time(t_rel_s, session_start)
%
%   t_rel_s      : numeric, any shape (seconds since session_start)
%   session_start: scalar datetime with a non-empty TimeZone

    if ~isa(session_start, 'datetime') || ~isscalar(session_start)
        error('rel_to_abs_time:BadSessionStart', 'session_start must be a scalar datetime.');
    end
    if isempty(session_start.TimeZone)
        error('rel_to_abs_time:NoTimeZone', 'session_start must have an explicit TimeZone.');
    end
    if ~isnumeric(t_rel_s)
        error('rel_to_abs_time:BadInput', 't_rel_s must be numeric.');
    end

    abs_dt = session_start + seconds(t_rel_s);
end
