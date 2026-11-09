%% CONSOLIDATE_4OHT  Statistics of the 4-OHT recombination window, added to
% the consolidation (consolidate_results.m) as new sheets:
%   seizures_all.xlsx      seizures_4-OHT     the seizures with onset in the window, + 4-OHT_delay
%   seizures_summary.xlsx  per_mouse_4-OHT    per_mouse's columns and criteria, inside the window
%                          epileptic_groups   SAN / IAN (a seizure with delay in [-2, +12] h)
%   iid_summary.xlsx       per_mouse_4-OHT    per_mouse's columns and criteria, inside the window
% plus a short section at the end of README.md. Existing sheets and
% figures/ are not touched (checked cell by cell after writing).
%
% Window: from 2 h before the injection (sheet "4-OHT" of the recording log)
% to 48 h after, or to the end of the recording if it ends earlier (last
% analyzed, not excluded EDF of the mouse that overlaps it). For the seizures
% only, a disconnection reviewed on video that continues the recording extends
% that end to its own (055-s, DX-23: the EEG was cut there; 2026-11-04; since
% 2026-11-05 it ends where the video review ends, and that time is observed).
% A seizure is in the window if its onset is; delay = (onset - injection) in h;
% seizures that straddle a window edge are listed in the README. Observed time
% for seizures = per_mouse's (T_seizure_obs: recording minus disconnections not
% reviewed on video, plus 055-s's DX-23) cut to the window; for IIDs = each
% channel's analyzed time (span minus the exclusions the re-detection applied)
% cut to the window. Columns the user adds to epileptic_groups (after its own)
% are kept, by mouse (2026-11-05).
%
% Reads only consolidated results: _work/state.mat, seizures_all.xlsx, and the
% re-detected IID complexes of every consolidated channel -- those of the
% state (iid_C, the 2026-11-04 update's mice), the 2026-11-01 run's
% phase3_results.mat (C) or _work/stats/iid_rerun_final -- each checked
% against iid_summary.xlsx (n_IID, observed_min). Nothing is re-detected, no
% individual IID is saved. Before writing, the same computation over an
% unbounded window must reproduce the existing per_mouse sheets.
%
% Run:  matlab -batch "run('<repo>/consolidate_4OHT.m')"
% Optional before run(): consolidate_out (default <batch 2026>/consolidated_20261102).

repo = fileparts(mfilename('fullpath'));
batch = fileparts(repo);
if ~exist('consolidate_out', 'var'), consolidate_out = fullfile(batch, 'consolidated_20261102'); end
restoredefaultpath;
addpath(genpath(fullfile(repo, 'src')));
rec_log = '\\iws\epic\users\camila.morel\1-DATA\Video_ECoG\2026\natus\EDF export\decimated\EEG_recording_log.xlsx';
run_0111 = fullfile(batch, '_stats_revision_20261101_final', 'phase3_results.mat');
% mice left out of the 4-OHT sheets (reported): {mouse, reason}. 113-s (out on 2026-11-03: the EDF of
% its window was missing) is in since its EDFs were added (2026-11-04)
left_out = cell(0, 2);
pre_h = 2; post_h = 48; san_h = [-2 12];

%% ============================ INPUTS =====================================
sf = fullfile(consolidate_out, '_work', 'state.mat');
q = load(sf, 'cfg', 'S', 'O', 'Z', 'rv', 'edf_tbl', 'SZ');
cfg = q.cfg; S = q.S; O = q.O; Z = q.Z; rv = q.rv; edf_tbl = q.edf_tbl; SZ = q.SZ;
iid_C = table();   % complexes of the channels consolidated by the last update (consolidate_results.m)
if ~isempty(whos('-file', sf, 'iid_C'))
    q = load(sf, 'iid_C'); iid_C = q.iid_C(strcmp(q.iid_C.variant, 'principal'), :);
end
clear q
tz = cfg.tz;
xs = fullfile(consolidate_out, {'seizures_all.xlsx', 'seizures_summary.xlsx', 'iid_summary.xlsx'});
txt = {'mouse', 'group', 'edf_file', 'original_id', 'event_id', 'hemispheres', 'onset', 'offset', 'source', 'light_phase', ...
    'comments', 'hemisphere'};
SA = read_sheet(xs{1}, 'seizures', txt);
SM = read_sheet(xs{2}, 'per_mouse', txt);
IE = read_sheet(xs{3}, 'per_mouse_EDF_hemisphere', txt);
IM = read_sheet(xs{3}, 'per_mouse', txt);
[tf, j] = ismember(SA.event_id, SZ.event_id);
assert(all(tf) && height(SA) == height(SZ), 'seizures_all.xlsx and _work/state.mat do not hold the same seizures');
SA.start_posix = SZ.start_posix(j); SA.in_dx = SZ.in_dx(j);

opts = detectImportOptions(rec_log, 'Sheet', '4-OHT', 'VariableNamingRule', 'preserve');
opts = setvartype(opts, opts.VariableNames{1}, 'char');
X = readtable(rec_log, opts);
X = X(:, 1:3);
X.Properties.VariableNames = {'mouse', 'rec_start', 'inj'};
X.mouse = strtrim(X.mouse);
X.inj.TimeZone = tz; X.rec_start.TimeZone = tz;   % Natus local wall clock

mice = unique(S.units.mouse(S.units.included));
notes = {};
W = table();   % one row per mouse with a window
for i = 1:numel(mice)
    m = mice{i};
    if any(strcmp(left_out(:, 1), m))
        notes{end+1} = sprintf('%s: fuera de las hojas 4-OHT: %s', m, left_out{strcmp(left_out(:, 1), m), 2}); %#ok<SAGROW>
        continue;
    end
    k = find(strcmp(X.mouse, m));
    if numel(k) ~= 1 || isnat(X.inj(k))
        notes{end+1} = sprintf('%s: %d fila(s) con inyección en la hoja 4-OHT; fuera de las hojas 4-OHT', m, numel(k)); %#ok<SAGROW>
        continue;
    end
    inj = posixtime(X.inj(k));
    if X.inj(k) < X.rec_start(k)
        notes{end+1} = sprintf('%s: la inyección (%s) es anterior al inicio de registro de la hoja (%s)', m, ...
            char(X.inj(k), 'dd/MM/yyyy HH:mm'), char(X.rec_start(k), 'dd/MM/yyyy HH:mm')); %#ok<SAGROW>
    end
    e = edf_tbl(strcmp(edf_tbl.mouse, m), :);
    excl = ismember(strcat(e.mouse, '|', e.edf), strcat(rv.excluded.mouse, '|', rv.excluded.edf));
    ov = ~excl & e.end_posix > inj - pre_h * 3600 & e.start_posix < inj + post_h * 3600;
    if ~any(ov)
        why = 'ningún EDF analizado cubre la ventana';
        exov = excl & e.end_posix > inj - pre_h * 3600 & e.start_posix < inj + post_h * 3600;
        if any(exov), why = sprintf('el EDF de la ventana está excluido (%s)', strjoin(e.edf(exov)', ', ')); end
        notes{end+1} = sprintf('%s: sin datos en la ventana (inyección %s): %s; fuera de las hojas 4-OHT', m, ...
            char(X.inj(k), 'dd/MM/yyyy HH:mm'), why); %#ok<SAGROW>
        continue;
    end
    rec_end = max(e.end_posix(ov));
    % a disconnection reviewed on video that continues the recording: the window goes on through it for the
    % seizures (055-s, DX-23); its time is observed for seizures (consolidate_results.m), not for IIDs
    sz_end = rec_end;
    dx = rv.dx(strcmp(rv.dx.mouse, m) & rv.dx.video_reviewed & rv.dx.start_posix <= rec_end + 1 & rv.dx.end_posix > rec_end, :);
    for k2 = 1:height(dx)
        sz_end = max(sz_end, dx.end_posix(k2));
        notes{end+1} = sprintf(['%s: el EEG termina a %+.1f h de la inyección (desconexión %s); la ventana sigue por esa desconexión, ' ...
            'revisada en video, hasta %+.1f h: sus crisis y su tiempo cuentan para crisis; los IIDs son los del EEG'], m, ...
            (rec_end - inj) / 3600, dx.id(k2), (dx.end_posix(k2) - inj) / 3600); %#ok<SAGROW>
    end
    W = [W; table({m}, S.units.group(find(strcmp(S.units.mouse, m), 1)), inj, inj - pre_h * 3600, ...
        min(inj + post_h * 3600, sz_end), sz_end >= inj + post_h * 3600, sum(ov), rec_end, sz_end > rec_end, ...
        'VariableNames', {'mouse', 'group', 'inj', 'w0', 'w1', 'reached_48h', 'n_edfs', 'eeg_end', 'video_end'})]; %#ok<AGROW>
end
% seizures across a window edge (only the onset decides): listed
for i = 1:height(W)
    for e2 = [W.w0(i) W.w1(i)]
        for j = find(strcmp(SA.mouse, W.mouse{i}) & SA.start_posix < e2 & SA.start_posix + SA.duration_s > e2)'
            inside = 'no entra'; if SA.start_posix(j) >= W.w0(i), inside = 'entra'; end
            notes{end+1} = sprintf(['%s: la crisis %d de %s (%s; onset %s, offset %s) cruza el borde de la ventana (%s): %s en ' ...
                '`seizures_4-OHT` (cuenta el onset)'], W.mouse{i}, SA.seizure_id(j), SA.edf_file{j}, SA.original_id{j}, SA.onset{j}, ...
                SA.offset{j}, char(datetime(e2, 'ConvertFrom', 'posixtime', 'TimeZone', tz, 'Format', 'dd/MM/yyyy HH:mm:ss')), inside); %#ok<SAGROW>
        end
    end
end
om = O.mouse(strcmp({O.mouse.scope}, 'merged'));

%% ============================ SEIZURES ===================================
% self-test: an unbounded window must give per_mouse as it is
full = W; full.w0(:) = -Inf; full.w1(:) = Inf;
full.n_edfs = SM.n_edfs(cellfun(@(m) find(strcmp(SM.mouse, m)), full.mouse));
chk = seizure_rows(full, SA, S, om);
ref = SM(cellfun(@(m) find(strcmp(SM.mouse, m)), full.mouse), :);
assert(isequal(chk.Properties.VariableNames, SM.Properties.VariableNames), 'per_mouse_4-OHT columns differ from per_mouse');
assert(same_table(chk, ref), 'the window computation does not reproduce per_mouse over the whole recording');

SW = seizure_rows(W, SA, S, om);
in_w = false(height(SA), 1); delay = NaN(height(SA), 1);
for i = 1:height(W)
    s = strcmp(SA.mouse, W.mouse{i}) & SA.start_posix >= W.w0(i) & SA.start_posix < W.w1(i);
    in_w(s) = true;
    delay(s) = (SA.start_posix(s) - W.inj(i)) / 3600;
end
T4 = SA(in_w, :);
T4.delay = delay(in_w);
T4 = table(T4.mouse, T4.edf_file, T4.seizure_id, T4.original_id, T4.hemispheres, T4.onset, T4.offset, T4.duration_s, T4.delay, ...
    T4.peak_energy_ratio_HPCl, T4.peak_energy_ratio_HPCr, T4.mean_to_thr_HPCl, T4.mean_to_thr_HPCr, T4.auc_to_thr_s_HPCl, ...
    T4.auc_to_thr_s_HPCr, T4.source, T4.light_phase, 'VariableNames', {'mouse', 'edf', 'seizure_id', 'original_id', 'hemispheres', ...
    'onset', 'offset', 'duration', '4-OHT_delay', 'peak_energy_ratio_HPCl', 'peak_energy_ratio_HPCr', 'mean_to_thr_HPCl', ...
    'mean_to_thr_HPCr', 'auc_to_thr_s_HPCl', 'auc_to_thr_s_HPCr', 'source', 'light_phase'});
assert(height(T4) == sum(SW.n_seizures), 'seizures_4-OHT rows ~= sum of n_seizures');

ep = W(strcmp(W.group, 'epileptic'), :);
G = table(ep.mouse, repmat({''}, height(ep), 1), zeros(height(ep), 1), NaN(height(ep), 1), ...
    'VariableNames', {'mouse', 'group', 'n_seizures', 'closest_seizure_delay_h'});
for i = 1:height(ep)
    d = T4.('4-OHT_delay')(strcmp(T4.mouse, ep.mouse{i}));
    G.n_seizures(i) = numel(d);
    G.group{i} = 'IAN';
    if any(d >= san_h(1) & d <= san_h(2)), G.group{i} = 'SAN'; end
    if ~isempty(d)
        [~, k] = min(abs(d));
        G.closest_seizure_delay_h(i) = d(k);
    end
end

%% ============================ IIDs =======================================
% re-detected complexes of every consolidated channel, checked against iid_summary
q = load(run_0111, 'C');
C = q.C(strcmp(q.C.variant, 'principal'), :);
clear q
P = Z.channels(Z.channels.run, :);
kC = strcat(C.mouse, '|', C.edf, '|', C.region);
kP = strcat(P.mouse, '|', P.edf, '|', P.region);
kU = {};
if ~isempty(iid_C), kU = strcat(iid_C.mouse, '|', iid_C.edf, '|', iid_C.region); end
IE.edf = erase(IE.edf_file, '.EDF');
kE = strcat(IE.mouse, '|', IE.edf, '|', IE.hemisphere);
CH = struct('mouse', IE.mouse, 'hemisphere', IE.hemisphere, 'start_posix', NaN, 'duration_s', NaN, 'complexes', [], 'bursts', [], ...
    'excl', [], 'src', '');
for i = 1:height(IE)
    j = find(strcmp(kC, kE{i}));
    ju = find(strcmp(kU, kE{i}));
    if isscalar(ju)   % a mouse of the last update: its own results
        r = struct('start', iid_C.start_posix(ju), 'dur', iid_C.duration_s(ju), 'cx', iid_C.complexes{ju}, 'b', iid_C.bursts{ju}, ...
            'ex', iid_C.excl{ju}, 'n', iid_C.n_complexes(ju), 'amin', iid_C.analyzed_min(ju), 'src', 'update');
    elseif isscalar(j)
        r = struct('start', C.start_posix(j), 'dur', C.duration_s(j), 'cx', C.complexes{j}, 'b', C.bursts{j}, 'ex', C.excl{j}, ...
            'n', C.n_complexes(j), 'amin', C.analyzed_min(j), 'src', '20261101');
    else
        f = dir(fullfile(consolidate_out, '_work', 'stats', 'iid_rerun_final', IE.mouse{i}, IE.edf{i}, ...
            sprintf('*_%s_IID_rerun.mat', IE.hemisphere{i})));
        jp = find(strcmp(kP, kE{i}));
        if ~isscalar(f) || ~isscalar(jp)
            error('consolidate_4OHT:IidMissing', '%s: no re-detected IIDs in the 2026-11-01 run nor in _work -- ask the user', kE{i});
        end
        w = load(fullfile(f.folder, f.name), 'rr');
        cx = w.rr.spike_complex_table; b = w.rr.burst_table; ex = w.rr.exclusion_table;
        r = struct('start', P.start_posix(jp), 'dur', P.duration_s(jp), ...
            'cx', [cx.Start_s cx.End_s cx.N_spikes double(cx.Is_polyspike) cx.Max_amplitude_mV], ...
            'b', [b.Start_s b.End_s b.N_complexes], 'ex', [ex.start_s ex.end_s], ...
            'n', w.rr.summary.total_complexes, 'amin', w.rr.summary.analyzed_duration_min, 'src', '_work');
    end
    if r.n ~= IE.n_IID(i) || abs(r.amin - IE.observed_min(i)) > 1e-6
        error('consolidate_4OHT:IidDiffers', '%s: re-detected IIDs (%s) differ from iid_summary.xlsx', kE{i}, r.src);
    end
    CH(i).start_posix = r.start; CH(i).duration_s = r.dur; CH(i).complexes = r.cx;
    CH(i).bursts = r.b; CH(i).excl = r.ex; CH(i).src = r.src;
end
% self-test: an unbounded window must give per_mouse as it is (all mice of iid_summary; one row per
% mouse x hemisphere since 2026-11-04)
fullI = table(IM.mouse, IM.group, IM.hemisphere, -Inf(height(IM), 1), Inf(height(IM), 1), ...
    'VariableNames', {'mouse', 'group', 'hemisphere', 'w0', 'w1'});
chk = iid_rows(fullI, CH, cfg);
assert(isequal(chk.Properties.VariableNames, IM.Properties.VariableNames), 'IID per_mouse_4-OHT columns differ from per_mouse');
% per_mouse's observed_min is the re-detection's analyzed_duration, which differs
% from span - exclusions (used here, as for its light / dark split) by < 1 s per
% channel (run_campaign_stats check); hence a 1e-5 relative tolerance on it and its rates
assert(same_table(chk, IM, {'observed_min', 'IID_per_min', 'bursts_per_h'}, 1e-5), ...
    'the window computation does not reproduce the IID per_mouse over the whole recording');
hs = unique(IE(:, {'mouse', 'hemisphere'}));   % mouse, then HPCl before HPCr
[tf, j] = ismember(hs.mouse, W.mouse);
WH = [W(j(tf), {'mouse', 'group'}), hs(tf, 'hemisphere'), W(j(tf), {'w0', 'w1'})];
IW = iid_rows(WH, CH, cfg);

%% ============================ WRITE ======================================
% into copies first; the originals are replaced only if every existing sheet is unchanged
tmp = fullfile(tempdir, 'consolidate_4OHT');
if ~isfolder(tmp), mkdir(tmp); end
xt = cellfun(@(f) fullfile(tmp, [char(java.util.UUID.randomUUID) '.xlsx']), xs, 'UniformOutput', false);
for f = 1:numel(xs), copyfile(xs{f}, xt{f}); end
writetable(T4, xt{1}, 'Sheet', 'seizures_4-OHT', 'WriteMode', 'overwritesheet');
writetable(SW, xt{2}, 'Sheet', 'per_mouse_4-OHT', 'WriteMode', 'overwritesheet');
writetable(G, xt{2}, 'Sheet', 'epileptic_groups', 'WriteMode', 'overwritesheet');
% the user's columns of epileptic_groups (after the script's; 2026-11-05: '**' on 058-s and 085-s), back in their mouse's row
UG = user_columns(xs{2}, 'epileptic_groups', width(G));
put_user_columns(xt{2}, 'epileptic_groups', UG, G.mouse);
assert(isequal(user_columns(xt{2}, 'epileptic_groups', width(G)), UG), 'epileptic_groups: the user''s columns were not kept');
for i = 1:size(UG, 1)
    who = 'encabezado'; if ~isempty(UG{i, 1}), who = UG{i, 1}; end
    notes{end+1} = sprintf('`epileptic_groups`: celda del usuario conservada (columna %c, %s): "%s"', 'A' + UG{i, 2} - 1, who, val_text(UG{i, 3})); %#ok<SAGROW>
end
writetable(IW, xt{3}, 'Sheet', 'per_mouse_4-OHT', 'WriteMode', 'overwritesheet');
new4 = {'seizures_4-OHT', 'per_mouse_4-OHT', 'epileptic_groups'};
for f = 1:numel(xs)
    before = sheets_content(xs{f});
    after = sheets_content(xt{f});
    for s = fieldnames(before)'
        nm = before.(s{1}).name;
        if any(strcmp(nm, new4)), continue; end
        assert(isfield(after, s{1}) && strcmp(after.(s{1}).name, nm) && isequaln(before.(s{1}).cells, after.(s{1}).cells), ...
            '%s: sheet %s would change -- nothing written', xs{f}, nm);
    end
end
for f = 1:numel(xs), copyfile(xt{f}, xs{f}); delete(xt{f}); end
write_readme_4oht(fullfile(consolidate_out, 'README.md'), W, SW, G, IW, notes, pre_h, post_h, san_h);

%% ============================ REPORT =====================================
fprintf('\n==== 4-OHT: %s\n', consolidate_out);
for i = 1:height(W)
    k = strcmp(SW.mouse, W.mouse{i});
    fin = 'llega a +48 h';
    if ~W.reached_48h(i), fin = sprintf('recortada en %+.1f h', (W.w1(i) - W.inj(i)) / 3600); end
    if W.eeg_end(i) < W.w1(i), fin = sprintf('%s (EEG hasta %+.1f h)', fin, (W.eeg_end(i) - W.inj(i)) / 3600); end
    grp = 'control';
    if any(strcmp(G.mouse, W.mouse{i})), grp = G.group{strcmp(G.mouse, W.mouse{i})}; end
    fprintf('  %s (%s): observado %.1f h en la ventana; %s; %d crisis; %s\n', W.mouse{i}, W.group{i}, SW.observed_h(k), fin, ...
        SW.n_seizures(k), grp);
end
fprintf('  SAN %d | IAN %d (epilépticos con ventana: %d)\n', sum(strcmp(G.group, 'SAN')), sum(strcmp(G.group, 'IAN')), height(G));
for i = 1:numel(notes), fprintf('  nota: %s\n', notes{i}); end

%% ======================================================================
function T = seizure_rows(W, SA, S, om)
% per_mouse of seizures_summary.xlsx (consolidate_results.m: seizure_summary /
% summary_rows) restricted to each mouse's window [w0, w1)
    n = height(W);
    z = NaN(n, 1);
    T = table(W.mouse, W.group, W.n_edfs, zeros(n, 1), zeros(n, 1), z, z, z, z, z, zeros(n, 1), zeros(n, 1), zeros(n, 1), zeros(n, 1), ...
        'VariableNames', {'mouse', 'group', 'n_edfs', 'observed_h', 'n_seizures', 'seizures_per_24h', 'dur_mean_s', 'dur_median_s', ...
        'dur_min_s', 'dur_max_s', 'n_TP', 'n_FP', 'n_missed', 'n_missed_DX'});
    fp = S.fp_events;
    for i = 1:n
        win = [W.w0(i) W.w1(i)];
        Tsz = om(strcmp({om.mouse}, W.mouse{i})).T_sz;
        if all(isfinite(win)), Tsz = cs_iv('intersect', Tsz, win); end
        T.observed_h(i) = cs_iv('length', Tsz) / 3600;
        s = strcmp(SA.mouse, W.mouse{i}) & SA.start_posix >= win(1) & SA.start_posix < win(2);
        in_obs = arrayfun(@(t) any(t >= Tsz(:, 1) & t < Tsz(:, 2)), SA.start_posix(s));
        d = SA.duration_s(s);
        T.n_seizures(i) = sum(s);
        if T.observed_h(i) > 0, T.seizures_per_24h(i) = sum(in_obs) / T.observed_h(i) * 24; end
        if ~isempty(d)
            T.dur_mean_s(i) = mean(d); T.dur_median_s(i) = median(d); T.dur_min_s(i) = min(d); T.dur_max_s(i) = max(d);
        end
        missed = strcmp(SA.source(s), 'missed_by_pipeline');
        dx = ~cellfun(@isempty, SA.in_dx(s));
        T.n_TP(i) = sum(strcmp(SA.source(s), 'pipeline_TP'));
        T.n_FP(i) = sum(strcmp(fp.mouse, W.mouse{i}) & fp.start_posix >= win(1) & fp.start_posix < win(2));
        T.n_missed(i) = sum(missed & ~dx);
        T.n_missed_DX(i) = sum(missed & dx);
    end
    assert(all(T.n_seizures == T.n_TP + T.n_missed + T.n_missed_DX), 'n_seizures ~= n_TP + n_missed + n_missed_DX');
end

function T = iid_rows(W, CH, cfg)
% per_mouse of iid_summary.xlsx (consolidate_results.m: iid_summary /
% add_rates, one row per mouse x hemisphere, pooled over that hemisphere's
% EDFs) restricted to each row's window [w0, w1)
    on = hm(cfg.lights_on); off = hm(cfg.lights_off);
    n = height(W);
    A = zeros(n, 9);   % observed_s, n, light_s, n_light, dark_s, n_dark, n_poly, n_bursts, n_edfs
    amp = cell(n, 1);
    for i = 1:n
        win = [W.w0(i) W.w1(i)];
        for c = find(strcmp({CH.mouse}, W.mouse{i}) & strcmp({CH.hemisphere}, W.hemisphere{i}))
            h = CH(c);
            inc = cs_iv('subtract', [0 h.duration_s], h.excl) + h.start_posix;
            if all(isfinite(win)), inc = cs_iv('intersect', inc, win); end
            if cs_iv('length', inc) <= 0, continue; end
            t = h.complexes(:, 1) + h.start_posix;
            k = t >= win(1) & t < win(2);
            p = cs_clock_split(inc, cfg.tz, [0 on off 1440]);
            nd = 0;
            if any(k)
                dt = datetime(t(k), 'ConvertFrom', 'posixtime', 'TimeZone', cfg.tz);
                mm = hour(dt) * 60 + minute(dt) + second(dt) / 60;
                nd = sum(mm < on | mm >= off);
            end
            tb = h.bursts(:, 1) + h.start_posix;
            A(i, :) = A(i, :) + [cs_iv('length', inc), sum(k), p(2), sum(k) - nd, p(1) + p(3), nd, ...
                sum(h.complexes(k, 4) ~= 0), sum(tb >= win(1) & tb < win(2)), 1];
            amp{i} = [amp{i}; h.complexes(k, 5)];
        end
    end
    T = table(W.mouse, W.group, W.hemisphere, A(:, 9), A(:, 1) / 60, A(:, 2), ...
        'VariableNames', {'mouse', 'group', 'hemisphere', 'n_edfs', 'observed_min', 'n_IID'});
    T.IID_per_min = T.n_IID ./ T.observed_min;
    T.observed_min_light = A(:, 3) / 60; T.n_IID_light = A(:, 4); T.IID_per_min_light = A(:, 4) ./ (A(:, 3) / 60);
    T.observed_min_dark = A(:, 5) / 60; T.n_IID_dark = A(:, 6); T.IID_per_min_dark = A(:, 6) ./ (A(:, 5) / 60);
    T.n_single = A(:, 2) - A(:, 7); T.n_polyspike = A(:, 7); T.pct_polyspike = 100 * A(:, 7) ./ A(:, 2);
    T.n_bursts = A(:, 8); T.bursts_per_h = A(:, 8) ./ (A(:, 1) / 3600);
    T.median_amplitude_mV = cellfun(@(a) median(a, 'omitnan'), amp);
end

function m = hm(s)
    v = sscanf(s, '%d:%d'); m = v(1) * 60 + v(2);
end

function ok = same_table(A, B, loose, loose_tol)
% same rows and values (numbers to 1e-9 relative, the columns in loose to
% loose_tol; NaN == NaN); prints the first column that differs
    if nargin < 3, loose = {}; loose_tol = 0; end
    ok = height(A) == height(B);
    for v = A.Properties.VariableNames
        if ~ok, return; end
        a = A.(v{1}); b = B.(v{1});
        if isnumeric(a)
            tol = 1e-9;
            if any(strcmp(v{1}, loose)), tol = loose_tol; end
            ok = isequal(isnan(a), isnan(b)) && all(abs(a(~isnan(a)) - b(~isnan(b))) <= tol * max(1, abs(b(~isnan(b)))));
        else
            ok = isequal(a, b);
        end
        if ~ok, fprintf(2, 'same_table: column %s differs\n', v{1}); end
    end
end

function T = read_sheet(f, s, txt)
% a sheet of the consolidated workbooks, its text columns kept as text
    opts = detectImportOptions(f, 'Sheet', s, 'VariableNamingRule', 'preserve');
    t = intersect(opts.VariableNames, txt);
    opts = setvartype(opts, t, 'char');
    T = readtable(f, opts);
end

function U = user_columns(f, s, nc)
% the non-empty cells after the first nc columns of sheet s (the user's): {mouse ('' = header row), column, value}
    U = cell(0, 3);
    if ~any(strcmp(sheetnames(f), s)), return; end
    c = readcell(f, 'Sheet', s);
    for r = 1:size(c, 1)
        for k = nc + 1:size(c, 2)
            if ~isa(c{r, k}, 'missing')
                key = '';
                if r > 1, key = c{r, 1}; end
                U(end+1, :) = {key, k, c{r, k}}; %#ok<AGROW>
            end
        end
    end
end

function put_user_columns(f, s, U, mice)
% U (user_columns) into sheet s of f: header cells in row 1, the others in their mouse's row (mice = rows 2..)
    for i = 1:size(U, 1)
        r = 1;
        if ~isempty(U{i, 1})
            r = find(strcmp(mice, U{i, 1})) + 1;
            assert(isscalar(r), 'consolidate_4OHT: %s has a user cell in %s but no row any more -- ask the user', U{i, 1}, s);
        end
        assert(U{i, 2} <= 26, 'consolidate_4OHT: user cell of %s beyond column Z', s);
        writecell(U(i, 3), f, 'Sheet', s, 'Range', sprintf('%c%d', 'A' + U{i, 2} - 1, r));
    end
end

function s = val_text(v)
    if ischar(v) || isstring(v), s = char(v); else, s = num2str(v); end
end

function C = sheets_content(f)
% every sheet's cells (to check that writing the new sheets left the others as they were)
    C = struct();
    sh = sheetnames(f);
    for i = 1:numel(sh)
        C.(sprintf('s%d', i)) = struct('name', sh{i}, 'cells', {readcell(f, 'Sheet', sh{i})});
    end
end

function write_readme_4oht(file, W, SW, G, IW, notes, pre_h, post_h, san_h)
% replaces (or appends) the 4-OHT section at the end of README.md
    head = '## Ventana de recombinación 4-OHT (hojas `*_4-OHT` y `epileptic_groups`)';
    txt = fileread(file);
    k = strfind(txt, head);
    if ~isempty(k), txt = txt(1:k(1) - 1); end
    txt = strtrim(txt);
    L = {'', '', head, ''};
    L{end+1} = sprintf(['Generado por `epilepsy_app_natus/consolidate_4OHT.m` el %s a partir de los resultados consolidados (no se re-detectó nada). ' ...
        'Ventana: de %d h antes de la inyección de 4-OHT (hoja "4-OHT" de `EEG_recording_log.xlsx`) a %d h después, o hasta el fin del registro si ' ...
        'termina antes. Una crisis entra si su onset está en la ventana; `4-OHT_delay` = (onset − inyección) en horas.'], ...
        char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), pre_h, post_h);
    L{end+1} = '';
    L{end+1} = ['- `seizures_all.xlsx` → `seizures_4-OHT`: las crisis de la ventana con sus datos de `seizures` (`edf` = `edf_file`, ' ...
        '`duration` = `duration_s`, en s) y `4-OHT_delay`.'];
    L{end+1} = ['- `seizures_summary.xlsx` → `per_mouse_4-OHT`: mismas columnas y criterios que `per_mouse`, solo dentro de la ventana ' ...
        '(tiempo observado recortado a la ventana, sin las desconexiones). `epileptic_groups`: `SAN` si hubo al menos una crisis con delay ' ...
        sprintf('entre %+d y %+d h (incluidos), `IAN` si no; `n_seizures` en toda la ventana; `closest_seizure_delay_h` = delay de la crisis ', san_h) ...
        'de la ventana más cercana a la inyección (con signo; vacío sin crisis).'];
    L{end+1} = ['- `iid_summary.xlsx` → `per_mouse_4-OHT`: mismas columnas y criterios que `per_mouse` (una fila por ratón y hemisferio, HPCl y ' ...
        'HPCr, sin sumar hemisferios), solo dentro de la ventana. Tiempos de ' ...
        'los complejos de la re-detección: `_stats_revision_20261101_final/phase3_results.mat`, `_work/stats/iid_rerun_final` y, para los ratones ' ...
        'de las actualizaciones desde el 04/11, los de cada actualización (`iid_C` en `_work/state.mat`), verificados canal por canal contra ' ...
        '`per_mouse_EDF_hemisphere`; con una ventana sin límites el cálculo reproduce exactamente las hojas `per_mouse`.'];
    L{end+1} = ['- Si el EEG termina en una desconexión revisada en video (055-s, DX-23: video revisado hasta las 20:31:39), la ventana sigue ' ...
        'por esa desconexión para las crisis: sus crisis (en `n_seizures`, `n_missed_DX` y SAN/IAN) y su tiempo (`observed_h`) cuentan; los ' ...
        'IIDs son los del EEG.'];
    L{end+1} = ['- `epileptic_groups`: las columnas que el usuario agrega después de las del script (p. ej. la marca `**`) se conservan en la ' ...
        'fila de su ratón al reescribir la hoja.'];
    L{end+1} = sprintf('- SAN: %d, IAN: %d (de %d ratones epilépticos con ventana).', sum(strcmp(G.group, 'SAN')), sum(strcmp(G.group, 'IAN')), height(G));
    L{end+1} = '';
    L{end+1} = '| ratón | grupo | inyección | fin de la ventana | observado crisis (h) | observado IIDs HPCl (min) | observado IIDs HPCr (min) | crisis |';
    L{end+1} = '|---|---|---|---|---|---|---|---|';
    for i = 1:height(W)
        k = strcmp(SW.mouse, W.mouse{i});
        io = NaN(1, 2);
        hn = {'HPCl', 'HPCr'};
        for h = 1:2
            ki = strcmp(IW.mouse, W.mouse{i}) & strcmp(IW.hemisphere, hn{h});
            if any(ki), io(h) = IW.observed_min(ki); end
        end
        fin = '+48 h';
        lab = 'fin del registro';
        if W.video_end(i), lab = 'fin del video revisado'; end
        if ~W.reached_48h(i), fin = sprintf('%+.1f h (%s)', (W.w1(i) - W.inj(i)) / 3600, lab); end
        if W.eeg_end(i) < W.w1(i), fin = sprintf('%s; EEG hasta %+.1f h', fin, (W.eeg_end(i) - W.inj(i)) / 3600); end
        L{end+1} = sprintf('| %s | %s | %s | %s | %.1f | %.0f | %.0f | %d |', W.mouse{i}, W.group{i}, ...
            char(datetime(W.inj(i), 'ConvertFrom', 'posixtime', 'TimeZone', 'Europe/Paris', 'Format', 'dd/MM/yyyy HH:mm')), fin, ...
            SW.observed_h(k), io(1), io(2), SW.n_seizures(k)); %#ok<AGROW>
    end
    if ~isempty(notes)
        L{end+1} = '';
        L{end+1} = 'Ratones fuera de estas hojas o con observaciones:';
        for i = 1:numel(notes), L{end+1} = ['- ' notes{i}]; end %#ok<AGROW>
    end
    fid = fopen(file, 'w', 'n', 'UTF-8');
    fprintf(fid, '%s\n', txt, L{:});
    fclose(fid);
end
