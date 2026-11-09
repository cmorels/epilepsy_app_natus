function t = rp_start(E, tz)
% RP_START  Time of the first sample, as edf_import.m sets session_start:
% header start date/time + onset of the first record (EDF+ time keeping).
    t = datetime([E.hdr.startdate ' ' E.hdr.starttime], 'InputFormat', 'dd.MM.yy HH.mm.ss', 'TimeZone', tz) + ...
        seconds(E.tal0_s);
    t.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
end
