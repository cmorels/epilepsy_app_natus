function [t0, t1] = cs_edf_header_span(file, tz)
% CS_EDF_HEADER_SPAN  Start and end (POSIX s) of an EDF from its fixed
% 256-byte header only: startdate (bytes 169-176, dd.mm.yy), starttime
% (177-184, hh.mm.ss, local wall clock -> TimeZone tz), number of data
% records (237-244) x duration of a record in s (245-252). Years 85-99 ->
% 19xx, otherwise 20xx (EDF spec). NaN when the file cannot be read.

    t0 = NaN; t1 = NaN;
    fid = fopen(file, 'r');
    if fid < 0, return; end
    h = fread(fid, 256, '*char')';
    fclose(fid);
    if numel(h) < 256, return; end
    d = sscanf(h(169:176), '%d.%d.%d');
    t = sscanf(h(177:184), '%d.%d.%d');
    nrec = str2double(strtrim(h(237:244)));
    drec = str2double(strtrim(h(245:252)));
    if numel(d) ~= 3 || numel(t) ~= 3, return; end
    y = d(3) + 2000 - 100 * (d(3) >= 85);
    t0 = posixtime(datetime(y, d(2), d(1), t(1), t(2), t(3), 'TimeZone', tz));
    t1 = t0 + nrec * drec;
end
