%% IID_PERIICTAL_EXCLUSION  Basal IID rates per mouse and hemisphere without
% peri-ictal windows: an additive analysis of the consolidation
% (consolidated_20261102). Nothing is re-detected and no existing file is
% modified (the user's spec, 2026-11-09).
%
% For each configuration {side, window_min} (CONFIG below), every seizure of
% the mouse (seizures_all.xlsx, sheet seizures, any source) gets an exclusion
% zone: both = [onset - W, offset + W], pre = [onset - W, onset], post =
% [offset, offset + W], none = no zone. A mouse's zones are unioned in absolute
% time (they may cross EDFs) and apply to both hemispheres; each channel (mouse
% x EDF x hemisphere) loses the part of its valid IID mask they cover and the
% complexes whose onset falls in them (zone bounds included). The 4-OHT sheet
% uses the window of consolidate_4OHT.m, with the zones of all the mouse's
% seizures (one just outside the window can reach into it).
%
% Inputs (read only): the stored re-detection of every channel of
% iid_summary.xlsx (_work/state.mat iid_C, the 2026-11-01 run's
% phase3_results.mat C or _work/stats/iid_rerun_final, as consolidate_4OHT.m;
% each must give its per_mouse_EDF_hemisphere row), its valid mask and complex
% onsets (cs_iid_channel_intervals: channel span minus detect_iid's exclusion
% intervals), the seizures' exact onset / offset (_work/state.mat). observed_min:
% per_mouse = detect_iid's analyzed time minus the excluded time (side none
% gives iid_summary exactly; the interval mask differs from detect_iid's sample
% count by < 0.5 s per channel); per_mouse_4-OHT and the light / dark split =
% the interval mask, as iid_summary.xlsx and consolidate_4OHT.m.
% pct_excluded = excluded_min / (observed_min + excluded_min) * 100.
%
% Output: <batch 2026>/_iid_periictal_exclusion_<yyyyMMdd>/iid_periictal_exclusion.xlsx
% (per_mouse, per_mouse_4-OHT, README) and per_mouse.csv, per_mouse_4-OHT.csv.
% Checks (nothing is written if one fails): side none reproduces iid_summary
% (n_IID, observed_min within tol_repro_min), excluded_min does not decrease
% with W and both lies between max(pre, post) and pre + post, the controls are
% the same in every configuration, none = kept + excluded (IIDs and time), no
% negative time.
%
% Run:  matlab -batch "run('<repo>/iid_periictal_exclusion.m')"
% Optional before run(): periictal_out (output folder; must not exist yet).

%% ============================ CONFIG =====================================
% {side, window_min}; side none = the reference (no zone). More windows (60, 120 ...) can be added.
CONFIG = {'none', 0; ...
          'both', 10; 'both', 20; 'both', 30; ...
          'pre',  10; 'pre',  20; 'pre',  30; ...
          'post', 10; 'post', 20; 'post', 30};
tol_repro_min = 0.01;   % side none vs iid_summary.xlsx, observed_min (min)
pre_h = 2; post_h = 48; % 4-OHT window (consolidate_4OHT.m)

repo = fileparts(mfilename('fullpath'));
batch = fileparts(repo);
consolidate_out = fullfile(batch, 'consolidated_20261102');
rec_log = '\\iws\epic\users\camila.morel\1-DATA\Video_ECoG\2026\natus\EDF export\decimated\EEG_recording_log.xlsx';
run_0111 = fullfile(batch, '_stats_revision_20261101_final', 'phase3_results.mat');
if ~exist('periictal_out', 'var'), periictal_out = fullfile(batch, ['_iid_periictal_exclusion_' char(datetime('now', 'Format', 'yyyyMMdd'))]); end
if isfolder(periictal_out), error('iid_periictal_exclusion:Out', '%s already exists -- nothing written', periictal_out); end
restoredefaultpath;
addpath(genpath(fullfile(repo, 'src')));
git = cs_git_info(repo, {[mfilename('fullpath') '.m']});
[~, br] = system(sprintf('git -C "%s" rev-parse --abbrev-ref HEAD', repo));
[~, sst] = system(sprintf('git -C "%s" status --porcelain -- "%s.m"', repo, mfilename('fullpath')));
started = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));

%% ============================ INPUTS =====================================
sf = fullfile(consolidate_out, '_work', 'state.mat');
q = load(sf, 'cfg', 'S', 'Z', 'rv', 'edf_tbl', 'SZ', 'iid_C');
cfg = q.cfg; S = q.S; Z = q.Z; rv = q.rv; edf_tbl = q.edf_tbl; SZ = q.SZ; tz = cfg.tz;
U = q.iid_C(strcmp(q.iid_C.variant, 'principal'), :);
clear q
q = load(run_0111, 'C');
C0 = q.C(strcmp(q.C.variant, 'principal'), :);
clear q
xi = fullfile(consolidate_out, 'iid_summary.xlsx');
rd = @(f, s) readtable(f, 'Sheet', s, 'TextType', 'char', 'VariableNamingRule', 'preserve');
IE = rd(xi, 'per_mouse_EDF_hemisphere'); IM = rd(xi, 'per_mouse'); I4 = rd(xi, 'per_mouse_4-OHT');
xa = fullfile(consolidate_out, 'seizures_all.xlsx');
opts = detectImportOptions(xa, 'Sheet', 'seizures', 'VariableNamingRule', 'preserve');
opts = setvartype(opts, {'mouse', 'edf_file', 'original_id', 'onset', 'offset'}, 'char');
SA = readtable(xa, opts);
% the seizures of seizures_all.xlsx with their exact onset / offset (state), by mouse | EDF | original_id
key = @(m, e, o) strcat(m, '|', erase(e, '.EDF'), '|', o);
[tf, j] = ismember(key(SA.mouse, SA.edf_file, SA.original_id), key(SZ.mouse, SZ.edf, SZ.original_id));
assert(all(tf) && numel(unique(j)) == height(SA) && height(SA) == height(SZ), 'seizures_all.xlsx and _work/state.mat do not hold the same seizures');
assert(isequal(SA.onset, SZ.onset(j)) && isequal(SA.offset, SZ.offset(j)), 'seizures_all.xlsx: onset / offset differ from _work/state.mat');
SA.on = SZ.start_posix(j); SA.off = SZ.end_posix(j);

% the stored re-detection of each channel (row of per_mouse_EDF_hemisphere)
[C, src] = stored_channels(IE, Z, U, C0, consolidate_out);
K = cs_iid_channel_intervals(C, S);
on = hm(cfg.lights_on); off = hm(cfg.lights_off);
dark = cell(height(C), 1);   % per complex: onset in the dark phase
for i = 1:height(C)
    dt = datetime(K.on_abs{i}, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
    mm = hour(dt) * 60 + minute(dt) + second(dt) / 60;
    dark{i} = mm < on | mm >= off;
end
W4 = windows_4oht(rec_log, unique(C.mouse), edf_tbl, rv, tz, pre_h, post_h);
mice = unique(C.mouse);
grp = cellfun(@(m) S.units.group{find(strcmp(S.units.mouse, m), 1)}, mice, 'UniformOutput', false);
n_sz = cellfun(@(m) sum(strcmp(SA.mouse, m)), mice);

%% ============================ CONFIGURATIONS =============================
nv = {'observed_min', 'excluded_min', 'n_IID', 'n_single', 'n_polyspike', 'observed_min_light', 'n_IID_light', ...
    'observed_min_dark', 'n_IID_dark', 'n_dropped'};
R = {[], []};   % per sheet: one row per mouse x hemisphere x configuration
meta = cell(0, 5);   % mouse, group, hemisphere, side, window_min
for c = 1:size(CONFIG, 1)
    side = CONFIG{c, 1}; W = CONFIG{c, 2} * 60;
    for k = 1:numel(mice)
        m = mice{k};
        s = strcmp(SA.mouse, m);
        Zm = zones(SA.on(s), SA.off(s), side, W);
        for h = {'HPCl', 'HPCr'}
            ch = find(strcmp(C.mouse, m) & strcmp(C.region, h{1}));
            if isempty(ch), continue; end
            a = zeros(2, numel(nv));
            for i = ch'
                for p = 1:2   % 1 = whole recording, 2 = 4-OHT window
                    inc = K.inc{i};
                    t = K.on_abs{i};
                    if p == 2
                        if ~isfield(W4, matlab.lang.makeValidName(m)), continue; end
                        w = W4.(matlab.lang.makeValidName(m));
                        inc = cs_iv('intersect', inc, w);
                        sel = t >= w(1) & t < w(2);
                    else
                        sel = true(size(t));
                    end
                    ex = cs_iv('intersect', inc, Zm);
                    kept = cs_iv('subtract', inc, Zm);
                    drop = sel & in_iv(t, Zm);
                    keep = sel & ~drop;
                    poly = C.complexes{i}(:, 4) ~= 0;
                    lp = cs_clock_split(kept, tz, [0 on off 1440]);
                    obs = cs_iv('length', kept) / 60;
                    if p == 1, obs = C.analyzed_min(i) - cs_iv('length', ex) / 60; end
                    a(p, :) = a(p, :) + [obs, cs_iv('length', ex) / 60, sum(keep), sum(keep & ~poly), sum(keep & poly), lp(2) / 60, ...
                        sum(keep & ~dark{i}), (lp(1) + lp(3)) / 60, sum(keep & dark{i}), sum(sel & in_iv(t, cs_iv('intersect', Zm, K.inc{i})))];
                end
            end
            meta(end+1, :) = {m, grp{k}, h{1}, side, CONFIG{c, 2}}; %#ok<SAGROW>
            R{1}(end+1, :) = a(1, :);
            R{2}(end+1, :) = a(2, :);
        end
    end
end
T = cell(1, 2);
for p = 1:2
    t = cell2table(meta, 'VariableNames', {'mouse', 'group', 'hemisphere', 'side', 'window_min'});
    for v = 1:numel(nv), t.(nv{v}) = R{p}(:, v); end
    if p == 2, t = t(ismember(matlab.lang.makeValidName(t.mouse), fieldnames(W4)), :); end   % mice with a 4-OHT window
    T{p} = t;
end

%% ============================ CHECKS =====================================
chk = {};
sheets = {'per_mouse', 'per_mouse_4-OHT'};
ref = {IM, I4};
for p = 1:2
    t = T{p};
    n0 = t(strcmp(t.side, 'none'), :);
    % 1) side none reproduces iid_summary.xlsx
    [tf, j] = ismember(strcat(n0.mouse, '|', n0.hemisphere), strcat(ref{p}.mouse, '|', ref{p}.hemisphere));
    assert(all(tf) && height(n0) == height(ref{p}), '%s: side none has not the mouse x hemisphere rows of iid_summary.xlsx', sheets{p});
    dn = n0.n_IID - ref{p}.n_IID(j); dobs = n0.observed_min - ref{p}.observed_min(j);
    bad = dn ~= 0 | abs(dobs) > tol_repro_min;
    if any(bad)
        disp(table(n0.mouse(bad), n0.hemisphere(bad), dn(bad), dobs(bad), 'VariableNames', {'mouse', 'hemisphere', 'd_n_IID', 'd_observed_min'}));
        error('iid_periictal_exclusion:Repro', '%s: side none does not reproduce iid_summary.xlsx (%d rows) -- nothing written', sheets{p}, sum(bad));
    end
    oth = {'n_single', 'n_polyspike', 'observed_min_light', 'n_IID_light', 'observed_min_dark', 'n_IID_dark'};
    dmax = cellfun(@(v) max(abs(n0.(v) - ref{p}.(v)(j))), oth);
    chk(end+1, :) = {sprintf('1. reproducción (%s)', sheets{p}), sprintf(['PASA: %d ratón x hemisferio; n_IID idéntico, observed_min máx. |dif| %.2g min ' ...
        '(tolerancia %g); además máx. |dif| %s'], height(n0), max(abs(dobs)), tol_repro_min, ...
        strjoin(strcat(oth, {' '}, arrayfun(@(x) sprintf('%.2g', x), dmax, 'UniformOutput', false)), ', '))}; %#ok<SAGROW>
    % 2) monotony in W; both between max(pre, post) and pre + post
    g = strcat(t.mouse, '|', t.hemisphere);
    nm = 0; nb = 0;
    for u = unique(g)'
        r = t(strcmp(g, u{1}), :);
        for sd = {'both', 'pre', 'post'}
            x = sortrows(r(strcmp(r.side, sd{1}), :), 'window_min');
            nm = nm + sum(diff(x.excluded_min) < -1e-9);
        end
        for Wm = unique(r.window_min(~strcmp(r.side, 'none')))'
            e = @(sd) r.excluded_min(strcmp(r.side, sd) & r.window_min == Wm);
            if isscalar(e('both')) && isscalar(e('pre')) && isscalar(e('post'))
                nb = nb + (e('both') < max(e('pre'), e('post')) - 1e-9 || e('both') > e('pre') + e('post') + 1e-9);
            end
        end
    end
    assert(nm == 0 && nb == 0, '%s: monotony failed (%d decreasing steps, %d both outside [max(pre, post), pre + post]) -- nothing written', sheets{p}, nm, nb);
    chk(end+1, :) = {sprintf('2. monotonía (%s)', sheets{p}), 'PASA: excluded_min no decrece de W en W en ningún side; both entre max(pre, post) y pre + post en todas las filas'}; %#ok<SAGROW>
    % 3) controls unchanged across configurations
    ctl = unique(t.mouse(strcmp(t.group, 'control')));
    nc = 0;
    for u = ctl'
        for h = {'HPCl', 'HPCr'}
            r = t(strcmp(t.mouse, u{1}) & strcmp(t.hemisphere, h{1}), nv);
            nc = nc + ~isequaln(r{:, :}, repmat(r{1, :}, height(r), 1));
        end
    end
    assert(nc == 0, '%s: %d control mouse x hemisphere change across configurations -- nothing written', sheets{p}, nc);
    chk(end+1, :) = {sprintf('3. controles sin cambios (%s)', sheets{p}), sprintf('PASA: %d controles (%s), valores idénticos en las %d configuraciones', ...
        numel(ctl), strjoin(ctl', ', '), size(CONFIG, 1))}; %#ok<SAGROW>
    % 4) balance with side none (IIDs: kept + those with onset in the zones cut to the valid mask)
    [~, j0] = ismember(g, strcat(n0.mouse, '|', n0.hemisphere));
    bi = t.n_IID + t.n_dropped - n0.n_IID(j0); bo = t.observed_min + t.excluded_min - n0.observed_min(j0);
    assert(all(bi == 0) && all(abs(bo) < 1e-6), '%s: balance with side none fails (%d IID rows, %d time rows) -- nothing written', sheets{p}, ...
        sum(bi ~= 0), sum(abs(bo) >= 1e-6));
    chk(end+1, :) = {sprintf('4. cuadratura (%s)', sheets{p}), sprintf(['PASA: en las %d filas n_IID(none) = n_IID + IIDs excluidos y observed_min(none) = ' ...
        'observed_min + excluded_min (máx. |dif| %.1g min)'], height(t), max(abs(bo)))}; %#ok<SAGROW>
    % 5) no negative time
    assert(all(t.observed_min >= 0) && all(t.excluded_min >= 0), '%s: negative time -- nothing written', sheets{p});
    chk(end+1, :) = {sprintf('5. rangos (%s)', sheets{p}), sprintf('PASA: observed_min mínimo %.1f min, ninguna fila < 0', min(t.observed_min))}; %#ok<SAGROW>
end

%% ============================ WRITE ======================================
cols = {'mouse', 'group', 'hemisphere', 'side', 'window_min', 'observed_min', 'excluded_min', 'pct_excluded', 'n_IID', 'IID_per_min', ...
    'n_single', 'n_polyspike', 'observed_min_light', 'n_IID_light', 'IID_per_min_light', 'observed_min_dark', 'n_IID_dark', 'IID_per_min_dark'};
for p = 1:2
    t = T{p};
    t.pct_excluded = 100 * t.excluded_min ./ (t.observed_min + t.excluded_min);
    t.IID_per_min = t.n_IID ./ t.observed_min;
    t.IID_per_min_light = t.n_IID_light ./ t.observed_min_light;
    t.IID_per_min_dark = t.n_IID_dark ./ t.observed_min_dark;
    T{p} = sortrows(t(:, cols), {'mouse', 'hemisphere'});   % configurations in CONFIG order within each mouse x hemisphere (stable)
end
L = {'fecha', started; 'script', [mfilename('fullpath') '.m']; 'rama git', strtrim(br); 'git HEAD', git.head; ...
    'script en HEAD sin cambios', mat2str(isempty(strtrim(sst))); 'script SHA-1', git.files{1, 2}; ...
    'otros cambios sin commit en el repositorio', strrep(git.status, newline, ' | '); ...
    'código de la consolidación (state.mat)', sprintf('git HEAD %s, cambios sin commit: %s', cfg.git.head, mat2str(cfg.git.dirty))};
for c = 1:size(CONFIG, 1), L(end+1, :) = {sprintf('CONFIG %d', c), sprintf('side = %s, window_min = %g', CONFIG{c, 1}, CONFIG{c, 2})}; end %#ok<SAGROW>
L = [L; {'tol_repro_min', num2str(tol_repro_min); 'ventana 4-OHT', sprintf('de %d h antes de la inyección a %d h después o al fin del registro (consolidate_4OHT.m)', pre_h, post_h); ...
    'zonas', ['both = [onset - W, offset + W], pre = [onset - W, onset], post = [offset, offset + W]; unidas por ratón en tiempo absoluto, ' ...
    'iguales para HPCl y HPCr, recortadas a la máscara válida de IIDs; un IID se excluye si su onset cae en una zona (bordes incluidos)']; ...
    'observed_min', ['per_mouse: tiempo analizado de detect_iid - tiempo excluido (side none = iid_summary); per_mouse_4-OHT y luz/oscuridad: ' ...
    'máscara de intervalos (como iid_summary.xlsx y consolidate_4OHT.m)']; ...
    'pct_excluded', 'excluded_min / (observed_min + excluded_min) * 100'; ...
    'entrada: estado', sf; 'entrada: iid_summary', xi; 'entrada: crisis', [xa ' (hoja seizures)']; 'entrada: corrida 01/11', run_0111; ...
    'entrada: re-detección guardada', fullfile(consolidate_out, '_work', 'stats', 'iid_rerun_final'); 'entrada: registro 4-OHT', rec_log}];
[u, ~, gs] = unique(src);
for k = 1:numel(u), L(end+1, :) = {'canales IID', sprintf('%d de %s', sum(gs == k), u{k})}; end %#ok<SAGROW>
L(end+1, :) = {'crisis usadas (total)', sprintf('%d (%s)', height(SA), strjoin(strcat(unique(SA.source), '=', ...
    arrayfun(@(s) num2str(sum(strcmp(SA.source, s))), unique(SA.source), 'UniformOutput', false))', ', '))};
for k = 1:numel(mice), L(end+1, :) = {sprintf('crisis usadas: %s (%s)', mice{k}, grp{k}), num2str(n_sz(k))}; end %#ok<SAGROW>
L = [L; chk];
mkdir(periictal_out);
xo = fullfile(periictal_out, 'iid_periictal_exclusion.xlsx');
for p = 1:2
    writetable(T{p}, xo, 'Sheet', sheets{p});
    writetable(T{p}, fullfile(periictal_out, [sheets{p} '.csv']), 'Encoding', 'UTF-8');
end
writetable(cell2table(L, 'VariableNames', {'item', 'value'}), xo, 'Sheet', 'README');

%% ============================ REPORT =====================================
fprintf('\n==== %s\n', xo);
for i = 1:size(chk, 1), fprintf('  %s: %s\n', chk{i, 1}, chk{i, 2}); end
fprintf('\nboth 30 min, per mouse (HPCl = HPCr for pct_excluded unless shown): pct_excluded | IID/min mean(HPCl, HPCr) none -> both30 (relative change)\n');
for p = 1:2
    t = T{p};
    for k = 1:numel(mice)
        a = t(strcmp(t.mouse, mice{k}) & strcmp(t.side, 'both') & t.window_min == 30, :);
        b = t(strcmp(t.mouse, mice{k}) & strcmp(t.side, 'none'), :);
        if isempty(a), continue; end
        r0 = mean(b.IID_per_min); r1 = mean(a.IID_per_min);
        fprintf('REP|%s|%s|%s|%s|%.2f|%.4f|%.4f|%.2f|%d\n', sheets{p}, mice{k}, grp{k}, strjoin(arrayfun(@(x) sprintf('%.2f', x), unique(round(a.pct_excluded, 4))', ...
            'UniformOutput', false), '/'), mean(a.excluded_min), r0, r1, 100 * (r1 / r0 - 1), n_sz(k));
    end
end

%% ======================================================================
function Zn = zones(on, off, side, W)
% the mouse's exclusion zones (POSIX s), unioned
    switch side
        case 'none', Zn = zeros(0, 2);
        case 'both', Zn = [on - W, off + W];
        case 'pre',  Zn = [on - W, on];
        case 'post', Zn = [off, off + W];
        otherwise, error('iid_periictal_exclusion:Side', 'unknown side %s', side);
    end
    Zn = cs_iv('union', Zn);
end

function b = in_iv(t, Zn)
% onsets t inside any zone, bounds included
    b = false(size(t));
    for k = 1:size(Zn, 1), b = b | (t >= Zn(k, 1) & t <= Zn(k, 2)); end
end

function m = hm(s)
    v = sscanf(s, '%d:%d'); m = v(1) * 60 + v(2);
end

function [C, src] = stored_channels(IE, Z, U, C0, out)
% one row per row of per_mouse_EDF_hemisphere: the state's iid_C, the 2026-11-01 run (phase3_results.mat)
% or _work/stats/iid_rerun_final, the first that holds the channel (as consolidate_results.m /
% consolidate_4OHT.m); it must give the row's n_IID and observed_min
    P = Z.channels(Z.channels.run, :);
    key = @(m, e, h) strcat(m, '|', erase(e, '.EDF'), '|', h);
    kU = key(U.mouse, U.edf, U.region); k0 = key(C0.mouse, C0.edf, C0.region); kP = key(P.mouse, P.edf, P.region);
    kE = key(IE.mouse, IE.edf_file, IE.hemisphere);
    cols = {'campaign', 'mouse', 'edf', 'region', 'variant', 'start_posix', 'duration_s', 'analyzed_min', 'n_complexes', 'complexes', 'excl'};
    rows = cell(height(IE), 1); src = cell(height(IE), 1);
    for i = 1:height(IE)
        jp = find(strcmp(kP, kE{i})); ju = find(strcmp(kU, kE{i})); j0 = find(strcmp(k0, kE{i}));
        if ~isscalar(jp), error('iid_periictal_exclusion:Plan', '%s: not once in the IID plan', kE{i}); end
        if isscalar(ju)
            r = U(ju, cols); src{i} = '_work/state.mat (iid_C)';
        elseif isscalar(j0)
            r = C0(j0, cols); src{i} = '_stats_revision_20261101_final/phase3_results.mat (C)';
        else
            f = dir(fullfile(out, '_work', 'stats', 'iid_rerun_final', IE.mouse{i}, erase(IE.edf_file{i}, '.EDF'), ...
                sprintf('*_%s_IID_rerun.mat', IE.hemisphere{i})));
            if ~isscalar(f), error('iid_periictal_exclusion:Stored', '%s: no stored re-detection', kE{i}); end
            w = load(fullfile(f.folder, f.name), 'rr');
            s = w.rr.summary; cx = w.rr.spike_complex_table; e = w.rr.exclusion_table;
            r = table(P.campaign(jp), P.mouse(jp), P.edf(jp), P.region(jp), {'principal'}, P.start_posix(jp), P.duration_s(jp), ...
                s.analyzed_duration_min, s.total_complexes, {[cx.Start_s cx.End_s cx.N_spikes double(cx.Is_polyspike) cx.Max_amplitude_mV]}, ...
                {[e.start_s e.end_s]}, 'VariableNames', cols);
            src{i} = '_work/stats/iid_rerun_final';
        end
        if r.n_complexes ~= IE.n_IID(i) || abs(r.analyzed_min - IE.observed_min(i)) > 1e-6
            error('iid_periictal_exclusion:Stored', '%s: the stored re-detection differs from per_mouse_EDF_hemisphere', kE{i});
        end
        rows{i} = r;
    end
    C = vertcat(rows{:});
end

function W = windows_4oht(rec_log, mice, edf_tbl, rv, tz, pre_h, post_h)
% the 4-OHT window [w0, w1) of each mouse, as consolidate_4OHT.m: from pre_h before the injection (sheet
% "4-OHT" of the recording log) to post_h after, or to the end of the recording (last analyzed, not excluded
% EDF overlapping it), extended through a disconnection reviewed on video that continues the recording.
% Field = matlab.lang.makeValidName(mouse); mice without a window have no field
    opts = detectImportOptions(rec_log, 'Sheet', '4-OHT', 'VariableNamingRule', 'preserve');
    opts = setvartype(opts, opts.VariableNames{1}, 'char');
    X = readtable(rec_log, opts);
    X = X(:, 1:3);
    X.Properties.VariableNames = {'mouse', 'rec_start', 'inj'};
    X.mouse = strtrim(X.mouse);
    X.inj.TimeZone = tz;
    W = struct();
    for i = 1:numel(mice)
        m = mice{i};
        k = find(strcmp(X.mouse, m));
        if numel(k) ~= 1 || isnat(X.inj(k)), continue; end
        inj = posixtime(X.inj(k));
        e = edf_tbl(strcmp(edf_tbl.mouse, m), :);
        ex = ismember(strcat(e.mouse, '|', e.edf), strcat(rv.excluded.mouse, '|', rv.excluded.edf));
        ov = ~ex & e.end_posix > inj - pre_h * 3600 & e.start_posix < inj + post_h * 3600;
        if ~any(ov), continue; end
        rec_end = max(e.end_posix(ov));
        sz_end = rec_end;
        dx = rv.dx(strcmp(rv.dx.mouse, m) & rv.dx.video_reviewed & rv.dx.start_posix <= rec_end + 1 & rv.dx.end_posix > rec_end, :);
        if height(dx), sz_end = max(sz_end, max(dx.end_posix)); end
        W.(matlab.lang.makeValidName(m)) = [inj - pre_h * 3600, min(inj + post_h * 3600, sz_end)];
    end
end
