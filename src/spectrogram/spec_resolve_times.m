function varargout = spec_resolve_times(action, varargin)
% SPEC_RESOLVE_TIMES  The three clocks of an event, for one session.
%
%   t_rel = spec_resolve_times('to_t_rel', value, time_ref, r)
%   [natus, pipeline] = spec_resolve_times('clocks', t_rel, r)
%   s = spec_resolve_times('text', dt)              % 'dd/MM/yyyy HH:mm:ss.SSS'
%
%   r : struct from spec_resolve_files.m (session_start, wall_start)
%
% t_rel    : seconds since the first sample of the txt, t_rel = (i-1)/fs.
% pipeline : session_start (with its TimeZone) + seconds(t_rel), i.e.
%            rel_to_abs_time.m -- the clock of start_abs/end_abs in
%            seizures_events. It follows daylight-saving changes.
% natus    : wall clock at session start + seconds(t_rel), an UNZONED
%            datetime. Natus shows the start clock time plus elapsed time
%            and does not jump at a daylight-saving change, so after
%            29/03/2026 (spring forward) pipeline = natus + 1 h for a
%            session that started before it.
%
% value (to_t_rel): number (already t_rel, time_ref ignored), datetime,
% or text 'dd/MM/yyyy HH:mm:ss.SSS' (also without .SSS, or yyyy-MM-dd).
% time_ref: 'natus' | 'pipeline' | 't_rel'. A pipeline-clock text that
% falls in the repeated hour of an autumn change is ambiguous and errors.

    switch action
        case 'to_t_rel'
            varargout{1} = to_t_rel(varargin{:});
        case 'clocks'
            t_rel = varargin{1};
            r = varargin{2};
            natus = r.wall_start + seconds(t_rel);
            natus.Format = 'dd/MM/yyyy HH:mm:ss.SSS';
            pipeline = rel_to_abs_time(t_rel, r.session_start);
            pipeline.Format = 'dd/MM/yyyy HH:mm:ss.SSS';
            varargout = {natus, pipeline};
        case 'text'
            dt = varargin{1};
            s = string(dt, 'dd/MM/yyyy HH:mm:ss.SSS');
            s(ismissing(s)) = "";
            if isscalar(s)
                varargout{1} = char(s);
            else
                varargout{1} = cellstr(s);
            end
        case 'parse'
            varargout{1} = parse_text(varargin{1});
        otherwise
            error('spec_resolve_times:BadAction', 'Unknown action ''%s''.', action);
    end
end

%% ======================================================================
function t = to_t_rel(value, time_ref, r)
    if isnumeric(value)
        t = double(value);
        return;
    end
    if iscell(value) && isscalar(value)
        value = value{1};
    end
    if isstring(value)
        value = char(value);
    end
    if ischar(value)
        num = str2double(value);
        if ~isnan(num) && ~isempty(regexp(strtrim(value), '^[-+]?[\d.]+([eE][-+]?\d+)?$', 'once'))
            t = num;
            return;
        end
        dt = parse_text(value);
    elseif isdatetime(value)
        dt = value;
        if ~isempty(dt.TimeZone)
            % a zoned datetime is an absolute instant whatever time_ref says
            t = seconds(dt - r.session_start);
            return;
        end
    else
        error('spec_resolve_times:BadValue', 'Time must be a number (t_rel), a datetime or text.');
    end

    switch lower(time_ref)
        case 'natus'
            t = seconds(dt - r.wall_start);
        case 'pipeline'
            tz = r.session_start.TimeZone;
            z = datetime(dt.Year, dt.Month, dt.Day, dt.Hour, dt.Minute, dt.Second, 'TimeZone', tz);
            if ~same_wall(z, dt)
                error('spec_resolve_times:NonexistentTime', ...
                    '%s does not exist in %s (skipped by the spring daylight-saving change).', char(dt), tz);
            end
            if same_wall(z - hours(1), dt) || same_wall(z + hours(1), dt)
                error('spec_resolve_times:AmbiguousTime', ...
                    '%s occurs twice in %s (autumn daylight-saving change); give it as t_rel or Natus time.', ...
                    char(dt), tz);
            end
            t = seconds(z - r.session_start);
        case 't_rel'
            error('spec_resolve_times:BadValue', 'time_ref ''t_rel'' needs a number, got "%s".', char(dt));
        otherwise
            error('spec_resolve_times:BadTimeRef', 'time_ref must be ''natus'', ''pipeline'' or ''t_rel'' (got ''%s'').', time_ref);
    end
end

function tf = same_wall(z, dt)
    tf = z.Year == dt.Year && z.Month == dt.Month && z.Day == dt.Day && z.Hour == dt.Hour && ...
        z.Minute == dt.Minute && abs(z.Second - dt.Second) < 1e-6;
end

function dt = parse_text(s)
    s = strtrim(char(s));
    fmts = {'dd/MM/yyyy HH:mm:ss.SSS', 'dd/MM/yyyy HH:mm:ss', 'dd/MM/yyyy HH:mm', ...
        'yyyy-MM-dd HH:mm:ss.SSS', 'yyyy-MM-dd HH:mm:ss', 'dd-MMM-yyyy HH:mm:ss.SSS', 'dd-MMM-yyyy HH:mm:ss'};
    for k = 1:numel(fmts)
        try
            dt = datetime(s, 'InputFormat', fmts{k}, 'Locale', 'en_US');
            if ~isnat(dt)
                return;
            end
        catch
        end
    end
    error('spec_resolve_times:BadText', ...
        'Could not read the time "%s" (expected dd/MM/yyyy HH:mm:ss.SSS).', s);
end
