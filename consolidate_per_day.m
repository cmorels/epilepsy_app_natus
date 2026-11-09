%% CONSOLIDATE_PER_DAY  Seizures per day of each mouse, from the consolidation
% (consolidate_results.m), in a new workbook seizures_per_day.xlsx with two
% definitions of a day (the user, 2026-11-05), one sheet each:
%   D_24h_desde_inicio  D1 = the first 24 h from the start of the mouse's
%                       recording (its first analyzed, not excluded EDF),
%                       D2 = the next 24 h, ...
%   D_calendario        calendar days (00:00-24:00 local time): D1 = the date
%                       the mouse's recording starts, D2 = the next date, ...
% Columns: mouse, day (D<k>), seizures (onset in the day), hours (observed
% seizure time in the day: per_mouse's observed_h, T_seizure_obs),
% seizures_per_24h = seizures / hours * 24, under_12h = 1 if the day has fewer
% than low_h observed hours (kept, only marked; the user, 2026-11-05), else 0.
% Days without observed time are left out; D<k> follows elapsed time, so the
% numbering skips them.
% Seizures = the rows of seizures_all.xlsx (sheet seizures), with their exact
% onset from _work/state.mat. Checks: every seizure falls in a day with
% observed time; per mouse, the days add up to per_mouse's n_seizures and
% observed_h. Days below low_h observed hours are also printed.
% Reads only consolidated results; writes only seizures_per_day.xlsx.
%
% Run:  matlab -batch "run('<repo>/consolidate_per_day.m')"
% Optional before run(): consolidate_out (default <batch 2026>/consolidated_20261102).

repo = fileparts(mfilename('fullpath'));
batch = fileparts(repo);
if ~exist('consolidate_out', 'var'), consolidate_out = fullfile(batch, 'consolidated_20261102'); end
restoredefaultpath;
addpath(genpath(fullfile(repo, 'src')));
low_h = 12;   % marked below this many observed hours (run_campaign_stats' cfg.min_valid_h_day)

%% ============================ INPUTS =====================================
q = load(fullfile(consolidate_out, '_work', 'state.mat'), 'cfg', 'O', 'SZ');
cfg = q.cfg; O = q.O; SZ = q.SZ; tz = cfg.tz;
clear q
xa = fullfile(consolidate_out, 'seizures_all.xlsx');
opts = detectImportOptions(xa, 'Sheet', 'seizures', 'VariableNamingRule', 'preserve');
opts = setvartype(opts, {'mouse', 'edf_file', 'original_id', 'onset'}, 'char');
SA = readtable(xa, opts);
opts = detectImportOptions(fullfile(consolidate_out, 'seizures_summary.xlsx'), 'Sheet', 'per_mouse', 'VariableNamingRule', 'preserve');
opts = setvartype(opts, {'mouse', 'group'}, 'char');
SM = readtable(fullfile(consolidate_out, 'seizures_summary.xlsx'), opts);
% the exact onset of each seizure of seizures_all, by mouse | EDF | original_id
key = @(m, e, o) strcat(m, '|', erase(e, '.EDF'), '|', o);
[tf, j] = ismember(key(SA.mouse, SA.edf_file, SA.original_id), key(SZ.mouse, SZ.edf, SZ.original_id));
assert(all(tf) && numel(unique(j)) == height(SA) && height(SA) == height(SZ), ...
    'seizures_all.xlsx and _work/state.mat do not hold the same seizures');
assert(isequal(SA.onset, SZ.onset(j)), 'seizures_all.xlsx: onsets differ from _work/state.mat');
SA.start_posix = SZ.start_posix(j);

%% ============================ DAYS =======================================
om = O.mouse(strcmp({O.mouse.scope}, 'merged'));
R = {cell(0, 4), cell(0, 4)};   % {mouse, day, seizures, hours} per sheet
for i = 1:height(SM)
    m = SM.mouse{i};
    o = om(strcmp({om.mouse}, m));
    T = o.T_sz;
    t0 = o.rec(1, 1);   % start of the mouse's recording
    s = SA.start_posix(strcmp(SA.mouse, m));
    assert(all(arrayfun(@(t) any(t >= T(:, 1) & t < T(:, 2)), s)), '%s: a seizure outside the observed time', m);
    % 24 h blocks from t0
    n_in = 0;
    for k = 1:ceil((T(end, 2) - t0) / 86400)
        b = t0 + [k - 1, k] * 86400;
        h = cs_iv('length', cs_iv('intersect', T, b)) / 3600;
        n = sum(s >= b(1) & s < b(2));
        if h > 0
            R{1}(end+1, :) = {m, sprintf('D%d', k), n, h};
            n_in = n_in + n;
        end
    end
    assert(n_in == numel(s), '%s: seizures outside the 24 h blocks with observed time', m);
    % calendar days (local), numbered from the date of t0
    [~, dd] = cs_clock_split(T, tz, [0 1440]);
    loc = @(t) datetime(t, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
    day0 = loc(t0); day0 = datetime(year(day0), month(day0), day(day0));
    ls = loc(s); ls = datetime(year(ls), month(ls), day(ls));
    n_in = 0;
    for k = 1:height(dd)
        if dd.total_s(k) <= 0, continue; end
        n = sum(ls == dd.date(k));
        R{2}(end+1, :) = {m, sprintf('D%d', days(dd.date(k) - day0) + 1), n, dd.total_s(k) / 3600};
        n_in = n_in + n;
    end
    assert(n_in == numel(s), '%s: seizures outside the calendar days with observed time', m);
end
names = {'D_24h_desde_inicio', 'D_calendario'};
T2 = cell(1, 2);
for p = 1:2
    t = cell2table(R{p}, 'VariableNames', {'mouse', 'day', 'seizures', 'hours'});
    t.seizures_per_24h = t.seizures ./ t.hours * 24;
    t.(sprintf('under_%gh', low_h)) = double(t.hours < low_h);
    % per mouse the days give per_mouse's n_seizures and observed_h
    for i = 1:height(SM)
        k = strcmp(t.mouse, SM.mouse{i});
        assert(sum(t.seizures(k)) == SM.n_seizures(i) && abs(sum(t.hours(k)) - SM.observed_h(i)) < 1e-6, ...
            '%s, %s: the days do not add up to per_mouse', names{p}, SM.mouse{i});
    end
    T2{p} = t;
end

%% ============================ WRITE ======================================
xo = fullfile(consolidate_out, 'seizures_per_day.xlsx');
if isfile(xo), delete(xo); end   % this script's own file
for p = 1:2, writetable(T2{p}, xo, 'Sheet', names{p}); end

%% ============================ REPORT =====================================
fprintf('\n==== %s\n', xo);
for p = 1:2
    t = T2{p};
    lo = t(t.hours < low_h, :);
    fprintf('\n%s: %d days (%d mice); %d with < %g h observed (%d of them with seizures)\n', names{p}, height(t), ...
        numel(unique(t.mouse)), height(lo), low_h, sum(lo.seizures > 0));
    for i = 1:height(lo)
        fprintf('  LOW|%s|%s|%s|%.2f|%d|%.2f\n', names{p}, lo.mouse{i}, lo.day{i}, lo.hours(i), lo.seizures(i), lo.seizures_per_24h(i));
    end
end
