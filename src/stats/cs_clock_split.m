function [bin_s, day_tbl] = cs_clock_split(iv, tz, edges_min)
% CS_CLOCK_SPLIT  Seconds of the intervals iv (Nx2 POSIX s) falling in each
% local wall-clock bin of the day. Bin edges are minutes after local
% midnight (e.g. 0:60:1440 for hours, [0 480 1170 1440] for dark/light/dark),
% built per calendar day with datetime(y,m,d,0,min,0,'TimeZone',tz) so they
% stay on the wall clock across DST changes.
%
%   bin_s   : 1 x nbins seconds
%   day_tbl : table(date, seconds_per_bin [1 x nbins], total_s), one row per
%             local calendar day touched by iv

    iv = cs_iv('union', iv);
    nb = numel(edges_min) - 1;
    bin_s = zeros(1, nb);
    dates = NaT(0, 1);
    vals = zeros(0, nb);
    for i = 1:size(iv, 1)
        t0 = datetime(iv(i, 1), 'ConvertFrom', 'posixtime', 'TimeZone', tz);
        d = dateshift(t0, 'start', 'day');
        while posixtime(d) < iv(i, 2)
            row = zeros(1, nb);
            for k = 1:nb
                bs = posixtime(datetime(year(d), month(d), day(d), 0, edges_min(k), 0, 'TimeZone', tz));
                be = posixtime(datetime(year(d), month(d), day(d), 0, edges_min(k+1), 0, 'TimeZone', tz));
                lo = max(bs, iv(i, 1));
                hi = min(be, iv(i, 2));
                if hi > lo
                    row(k) = hi - lo;
                end
            end
            dd = datetime(year(d), month(d), day(d));
            j = find(dates == dd, 1);
            if isempty(j)
                dates(end+1, 1) = dd; %#ok<AGROW>
                vals(end+1, :) = row; %#ok<AGROW>
            else
                vals(j, :) = vals(j, :) + row;
            end
            bin_s = bin_s + row;
            d = datetime(year(d), month(d), day(d) + 1, 'TimeZone', tz);
        end
    end
    [dates, o] = sort(dates);
    vals = vals(o, :);
    day_tbl = table(dates, vals, sum(vals, 2), 'VariableNames', {'date', 'seconds_per_bin', 'total_s'});
end
