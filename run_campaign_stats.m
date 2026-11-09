%% RUN_CAMPAIGN_STATS  Seizure and IID statistics of campaigns C1/C2 from
% EXISTING outputs plus the user's review (seizures_events.csv AE/AF columns
% and revision_campaigns_crisis.xlsx). Nothing upstream is re-run:
% extraction, decimation, precondition/clean_lfp, reference calibration and
% seizure detection are never called; only detect_iid is re-run (phase 2/3,
% through src/stats). Existing outputs are never modified; everything new
% goes to cfg.out_dir.
%
% v1 phases (cfg.phase):
%   1  final seizure set, observation time, seizure stats, detector
%      performance, seizure-threshold diagnostic (sections 5.1, 5.2, 6.1-6.3)
%   2  IID exclusion zones, detect_iid wrapper, reproducibility test (5.3)
%   3  IID re-detection on every channel, sections 6.4-6.5, checks (7-8)
% v2 (cfg.run_tag = 'v2', prompt_claude_code_v2_iid_graficos.md):
%   A  label snapshot, excluded EDFs, checks vs v1, channels to re-detect;
%      stops when one of its STOP conditions holds (until cfg.phaseA_stop_ack)
%   B  the v1 phases 1-3 from the snapshot, IID re-detection incremental on v1
%   C-E  flat segments + IID threshold sensitivity, 24-h days, figures
% final (cfg.run_tag = 'final'): the v2 flow with committed code; caches and
% IID runs of v2/v1 reused only where code and zones are identical
%
% Run:  matlab -batch "run('<repo>/run_campaign_stats.m')"

%% ============================== CONFIG ==================================
cfg = struct();
cfg.repo_root        = fileparts(mfilename('fullpath'));
cfg.review_xlsx      = fullfile(cfg.repo_root, 'revision_campaigns_crisis.xlsx');
cfg.controls         = {'051-s','056-s','070-s','083-s','109-s'};   % confirmed by the user (phase 0)
cfg.epileptic        = {};          % empty = every in-scope non-control mouse (listed in README)
cfg.lights_on        = '08:00';     % light phase 08:00-19:30 (11.5 h)
cfg.lights_off       = '19:30';     % dark phase 19:30-08:00 (12.5 h) -> always normalized by exposure
cfg.min_valid_h_day  = 12;          % a day enters daily analyses with >= 12 h of valid observation
cfg.periictal_win_min= [-60 -30; -30 0; 0 30; 30 60];   % rel. to onset (pre) and to end (post)
cfg.baseline_min_dist_h = 2;        % interictal baseline: >= 2 h from any seizure
cfg.cluster_gap_h    = 24;          % seizures <= 24 h apart = same cluster (6 h and 12 h also reported)
cfg.merge_review_gap_s = 30;        % TPs of a mouse closer than 30 s -> listed for review (not merged)
cfg.iid_excl_buffer_s  = 5;         % +/- s for ALL IID exclusion zones (= exclusion_buffer_s)
cfg.iid_campaign_for_duplicates = 'C2';   % EDF analyzed in both campaigns -> IIDs of C2 (none exist, phase 0)
cfg.sz_thr_diag_min_frac = 0.02;    % seizure-threshold diagnostic when disconnected fraction > 2 %
cfg.sz_thr_change_flag   = 0.05;    % threshold change >= 5 % -> flag
cfg.n_boot           = 10000;
cfg.rng_seed         = 1;

%% ---- paths and phase-0 decisions (2026-10-29) --------------------------
cfg.phase            = 3;           % run phases 1..cfg.phase
cfg.run_threshold_diag = true;      % section 6.3 (reads the large *_clean.txt; cached per channel)
cfg.batch_root       = fileparts(cfg.repo_root);
cfg.campaigns        = struct('name', {'C1', 'C2'}, ...
    'root', {fullfile(cfg.batch_root, 'pipeline_output_campaign_20261006'), ...
             fullfile(cfg.batch_root, 'pipeline_output_campaign_20261014')});
cfg.tz               = 'Europe/Paris';      % EDF header times are Natus local wall clock
% 077-s (C2): its seizures_events.csv was rewritten by the 2026-10-19 rerun;
% the labels live in this copy (same rows as merged/05_summaries).
cfg.label_file_override = {'C2', '077-s', 'seizures_events_excel_20261019.csv'};
cfg.comment_map      = {'digging', 'other'};   % unlisted AF values mapped by decision (else R6)
% crisette_focal: {obs id, mouse, campaign, EDF, seizure_id} -- outside the seizure total
cfg.crisette = { ...
    'O-07', '088-s', 'C1', '20260709A~ TRA_9298ff0f-2a8a-4829-a1dc-94016f3d4f25', 1; ...
    'O-08', '088-s', 'C2', '20260708~ TRAP_9762f514-4176-4494-bbf9-af131b8824b4', 2; ...
    'O-13', '084-s', 'C2', '20260713~ TRAP_5c69fdc2-d56a-4917-832b-f840a83ad5fa', 1};
% manual merges: {obs id, mouse, campaign, EDF, seizure_ids} -> one seizure from the first onset to the last end
% O-02: 110-s crisis 15 = one seizure (its two fragments share seizure_id 15)
% O-33..O-36 (2026-10-30): consecutive crisis that are one seizure (124-s, 062-s)
cfg.merge_events = { ...
    'O-02', '110-s', 'C1', '20260803~ TRAP_2bbce7f4-2ccd-4764-b26b-42267b0d33f6', 15; ...
    'O-33', '124-s', 'C2', '20260729~ TRAP_d5154696-610f-41eb-91af-d50e5658bcd3', [8 9]; ...
    'O-34', '062-s', 'C1', '20260528~ TRAP_8362c089-4ccf-4eee-a218-d3152509e3cf', [1 2]; ...
    'O-35', '062-s', 'C1', '20260528~ TRAP_8362c089-4ccf-4eee-a218-d3152509e3cf', [3 4]; ...
    'O-36', '062-s', 'C1', '20260603~ TRAP_c74056d4-3d8d-4a43-a124-e89056f37721', [8 9]};
cfg.nd_source = {'ND-01', 'missed_eeg'; 'ND-02', 'missed_eeg'; 'ND-03', 'missed_eeg'; 'ND-04', 'missed_eeg'; ...
                 'ND-05', 'video'; 'ND-06', 'video'};   % ND-07 has no times -> exclusion_log
cfg.disconnect_run_max_gap_s = 60;  % R5: consecutive 'disconnect' FPs join into one zone only if <= 60 s apart
cfg.cluster_gaps_h   = [6 12 24];
cfg.transition_window_h = 2;        % rate in the 2 h after each light transition vs rest of its phase
% IID re-detection (phases 2-3)
cfg.iid_rerun_scope  = 'in_scope';  % 'in_scope' = channels of reviewed units | 'all' = every channel of C1+C2
cfg.iid_rerun_figures = false;      % detect_iid's panorama/polyspike figures (~4 min/channel); results identical
% reproducibility test (5.3): {campaign, mouse, EDF prefix, region} -- C2, not 077-s (stale .mat)
cfg.iid_repro_channels = { ...
    'C2', '020-s', '20260406~',  'HPCr'; ...   % accepted + rescued + imputed rows
    'C2', '083-s', '20260722A~', 'HPCl'; ...   % accepted + imputed rows, 55 h
    'C2', '076-s', '20260713~',  'HPCr'};      % no detections (empty list)
% the user's answers to the flags of the phase reports (2026-10-29): {topic prefix, decision}
cfg.flag_decisions = { ...
    'fuera de alcance (sin revision)',      'decidido (fase 0): unidades sin revision fuera de las estadisticas'; ...
    'R5: que es "consecutivo"',             'decidido (fase 1): limite de 60 s aceptado'; ...
    'ND pendiente en el Excel',             'decidido (fase 1): ND-02..04 se suman'; ...
    'TP separados < 30 s',                  'decidido (fase 1): se cuentan como 2 crisis'; ...
    'R7 TP con texto en AF',                'decidido (fase 1): se mantienen como crisis'; ...
    'AF no reconocido (R6/R7)',             'decidido (fase 1): TP con texto en AF se mantienen como crisis (R7)'; ...
    'umbral de crisis cambia >= 5 %',       'informativo (6.3, sin re-detectar): revisado en fase 1'; ...
    'ND solapa una deteccion',              'revisado (fase 1): ND-05 dentro de DX-06, esperable'; ...
    'FP con AF mixto (IIDs + artefacto)',   'decidido (fase 2): excluir como R4 en ambos hemisferios'; ...
    'canal IID con < 1 h analizable',       'decidido (fase 2): se mantiene'; ...
    'alcance de la re-deteccion IID',       'decidido (fase 2): las unidades sin user_check/comments no se corren'; ...
    'figuras IID por canal',                'decidido (fase 2): sin figuras por ahora'; ...
    'A.3 048-s',                            'decidido (fase A v2): se excluye solo 20260622~ TRAP_fd13b728 (DX-20/ED-16)'; ...
    'EDF excluido por AE incompleto',       'decidido (fase A v2): regla A.2, solo esos EDF quedan fuera'; ...
    'tramos planos',                        'decidido (v2, 2026-10-30): por ahora no se excluyen; solo se listan'};

%% ---- v2 (2026-10-29): corrections, incremental IID re-detection --------
cfg.run_tag              = 'final';   % 2026-10-30: committed code, after v2 (all flags decided)
cfg.v1_dir               = fullfile(cfg.batch_root, '_stats_revision_20261029');   % reference of v1_vs_v2
cfg.prev_dir             = fullfile(cfg.batch_root, '_stats_revision_20261029_v2'); % caches (6.3, repro, sens) copied from here
cfg.out_dir              = fullfile(cfg.batch_root, ['_stats_revision_20261030_' cfg.run_tag]);
if exist('out_dir_override', 'var')   % set before run() to write a check run elsewhere (IID runs included)
    cfg.out_dir = out_dir_override;
end
cfg.stop_after           = 'E';     % last phase this run executes: 'A' .. 'E'
cfg.edf_root             = '\\iws\epic\users\camila.morel\1-DATA\Video_ECoG\2026\natus\EDF export\decimated';   % EDF headers only (day t0)
cfg.phaseA_stop_ack      = true;    % the user answered phase A's STOP (2026-10-29): 048-s only fd13b728; 084-s C2 by rule A.2
cfg.labels_snapshot_refresh = false;   % false: reuse labels_snapshot/ if it exists (one snapshot per run)
cfg.recording_log        = '\\iws\epic\users\camila.morel\1-DATA\Video_ECoG\2026\natus\EDF export\decimated\EEG_recording_log.xlsx';
% EDFs out of everything: {ED id (3_EDF_estado), mouse, EDF, kind}; DX rows with 'EDF completo' are added automatically
cfg.edf_excluded = { ...
    'ED-16', '048-s', '20260622~ TRAP_fd13b728-0df8-45ad-8f8c-2b79a894d3ef', 'descartado'; ...
    'ED-17', '128-s', '20260729~ TRAP_d5154696-610f-41eb-91af-d50e5658bcd3', 'pendiente de reprocesamiento'};
% 128-s 20260729 reprocessed as two EDFs, A1 only and A2 only (batch 2026\128-s\reprocess_20260729_A1A2,
% src/stats/reprocess/rp128_*.m). true: they replace the exclusion ED-17 (DX-21 stays), their IIDs are
% re-detected with the wrapper and everything is recomputed. Keep false until their seizures_events.csv is labeled.
cfg.include_128s_0729_reprocessed = false;
if exist('include_128s_0729_reprocessed_override', 'var')
    cfg.include_128s_0729_reprocessed = include_128s_0729_reprocessed_override;
end
% units kept in scope EDF by EDF: an EDF with any event lacking AE is left out (flag), not the whole unit
cfg.edf_level_units      = {'C2', '084-s'};
% IID results: {this run's root, earlier runs' roots}; a channel with the same zones and code is reused
cfg.iid_rerun_dirs = struct( ...
    'principal',    {{fullfile(cfg.out_dir, ['iid_rerun_' cfg.run_tag]), fullfile(cfg.prev_dir, 'iid_rerun_v2'), ...
                      fullfile(cfg.v1_dir, 'iid_rerun')}}, ...
    'sin_crisette', {{fullfile(cfg.out_dir, ['iid_rerun_' cfg.run_tag '_sin_crisette']), fullfile(cfg.prev_dir, 'iid_rerun_v2_sin_crisette'), ...
                      fullfile(cfg.v1_dir, 'iid_rerun_sin_crisette')}});
% Days
cfg.day_t0_rule          = 'first_edf_start';   % day 1 starts at the first sample of the mouse's first EDF (any campaign)
cfg.violin_min_h         = 20;      % mouse-days entering the per-day violin (counts)
cfg.violin_min_n         = 3;       % violin only if n >= 3, else points
cfg.raster_show_excluded_days = true;   % days with recording but excluded are shown shaded, not removed
% Post-ictal and peri-ictal
cfg.postictal_excl_min   = [60 120];    % 60 = main; 120 = sensitivity
cfg.lead_min_gap_h       = 4;       % lead seizure: no seizure in the previous 4 h
cfg.lead_min_coverage    = 0.8;     % ... and those 4 h observed >= 80 %
cfg.periictal_min_coverage = 0.8;   % peri-ictal window valid if >= 80 % observed
% Flat segments
cfg.flat_thr_uV          = 5;
cfg.flat_min_s           = 60;
% IID threshold sensitivity (descriptive)
cfg.sens_fixed_uV        = [312.5 400 500 625];
cfg.sens_pct_bp          = [99 99.5 99.9];   % percentiles of |x_bp| (15-70 Hz, rectified) per channel

if exist('run_threshold_diag_override', 'var')   % set before run() for a quick pass without 6.3
    cfg.run_threshold_diag = run_threshold_diag_override;
end

%% ---- derived ------------------------------------------------------------
on = sscanf(cfg.lights_on, '%d:%d'); on = on(1) * 60 + on(2);
off = sscanf(cfg.lights_off, '%d:%d'); off = off(1) * 60 + off(2);
w = cfg.transition_window_h * 60;
cfg.phase_edges_min = [0 on on + w off off + w 1440];   % dark | light-post | light | dark-post | dark
cfg.phase_dark_bins = [1 4 5];
cfg.phase_bin_after_on = 2;   cfg.phase_bins_light_rest = 3;
cfg.phase_bin_after_off = 4;  cfg.phase_bins_dark_rest = [1 5];
clear on off w

restoredefaultpath;                         % the user path shadows toolbox functions (e.g. chronux findpeaks)
addpath(genpath(fullfile(cfg.repo_root, 'src')));
if cfg.include_128s_0729_reprocessed
    cfg = rp128_integrate('config', cfg);   % drops ED-17 from cfg.edf_excluded
end
if ~isfolder(cfg.out_dir), mkdir(cfg.out_dir); end
cfg.run_started = char(datetime('now', 'TimeZone', cfg.tz, 'Format', 'yyyy-MM-dd HH:mm:ss Z'));
[cfg.git, cfg.provenance] = cs_provenance(cfg);   % IID code must be HEAD's, unmodified
fprintf('run_campaign_stats %s: HEAD %s (dirty=%d, IID code unmodified=%d), started %s\n', cfg.run_tag, cfg.git.head, ...
    cfg.git.dirty, cfg.git.iid_code_unmodified, cfg.run_started);
t_phase = struct();

%% ============================ PHASE A (v2) ================================
% snapshot, excluded EDFs, checks vs v1, channels to re-detect (section 3)
tA = tic;
code_state = cs_code_state(cfg);                       % git_head.txt + code_state.patch
for d = {cfg.prev_dir, cfg.v1_dir}   % 6.3 / repro / sens caches (keyed by code + inputs) of the latest earlier run
    if ~isfolder(fullfile(cfg.out_dir, 'cache')) && isfolder(fullfile(d{1}, 'cache'))
        copyfile(fullfile(d{1}, 'cache'), fullfile(cfg.out_dir, 'cache'));
    end
end
[snap, cfg.labels_snapshot_dir] = cs_label_snapshot(cfg);   % every label is read from here
if cfg.include_128s_0729_reprocessed
    snap = rp128_integrate('snapshot', cfg, snap);
end
fprintf('labels: %d files in %s (reused=%d)\n', height(snap), cfg.labels_snapshot_dir, snap.reused(1));

[edf_tbl, ch_tbl] = cs_load_edf_table(cfg);
[det_rows, review_rows, det_log] = cs_load_detections(cfg, edf_tbl);
if cfg.include_128s_0729_reprocessed   % original 20260729 out, A1 / A2 in (02_clean, detections + labels)
    [edf_tbl, ch_tbl] = rp128_integrate('tables', cfg, edf_tbl, ch_tbl);
    [det_rows, review_rows, rp_log] = rp128_integrate('detections', cfg, edf_tbl, det_rows, review_rows);
    det_log = [det_log, rp_log];
end
fprintf('%s\n', det_log{:});
events = cs_build_events(det_rows);
rv = cs_load_review(cfg, edf_tbl);
[rv.excluded, f_excl] = cs_excluded_edfs(cfg, rv, edf_tbl, events);
disp(rv.excluded(:, {'id', 'mouse', 'campaign', 'edf', 'kind'}));

edf_lists = table();   % listed for analysis but absent from 02_clean
for c = 1:numel(cfg.campaigns)
    L = readtable(fullfile(cfg.campaigns(c).root, 'edf_list.csv'), 'TextType', 'char', 'Delimiter', ',');
    L.campaign = repmat({cfg.campaigns(c).name}, height(L), 1);
    L.reason = repmat({'error en el analisis (ver log de la campaign)'}, height(L), 1);
    sk = fullfile(cfg.campaigns(c).root, 'skipped_inputs.csv');
    if isfile(sk)
        K = readtable(sk, 'TextType', 'char', 'Delimiter', ',');
        [tf, j] = ismember(strcat(L.animal_id, '|', L.edf_file), strcat(K.animal_id, '|', K.edf_file));
        L.reason(tf) = K.reason(j(tf));
    end
    edf_lists = [edf_lists; table(L.campaign, L.animal_id, L.edf_file, L.reason, ...
        'VariableNames', {'campaign', 'mouse', 'edf', 'reason'})]; %#ok<AGROW>
end
if cfg.include_128s_0729_reprocessed
    edf_lists = rp128_integrate('edf_lists', cfg, edf_lists);
end

S = cs_build_seizures(cfg, events, rv, edf_tbl, edf_lists);
O = cs_observation_time(cfg, S, rv, edf_tbl);

% A.1 labels vs v1, A.5 IID zones vs v1, checks + STOP conditions
LC = cs_label_changes(cfg, snap, det_rows);
q = load(fullfile(cfg.v1_dir, 'phase1_results.mat'), 'S');
v1_units = strcat(q.S.units.campaign(q.S.units.included), '|', q.S.units.mouse(q.S.units.included));
clear q
if ~cfg.git.iid_code_unmodified
    error('run_campaign_stats:IidCode', 'detect_iid / load_lfp_txt / pipeline_config (or their helpers) differ from HEAD');
end
Z = cs_iid_zones(cfg, S, rv, ch_tbl, det_rows);
chg = LC(LC.read_by_analysis & LC.changed_vs_v1, :);
RD = cs_iid_compare_v1(cfg, Z, strcat(chg.campaign, '|', chg.mouse), v1_units);
if cfg.include_128s_0729_reprocessed
    RD = rp128_integrate('iid_plan', cfg, RD);
end
[A, stopA, f_A] = cs_phaseA_checks(cfg, edf_tbl, S.events, rv, S, LC, RD);
disp(A(:, {'punto', 'resultado', 'stop'}));
t_phase.A = toc(tA);
if stopA && ~cfg.phaseA_stop_ack
    xf = fullfile(cfg.out_dir, ['stats_campaigns_' cfg.run_tag '.xlsx']);
    if isfile(xf), delete(xf); end
    cs_write_sheet(xf, 'README', cell2table(readme_v2(cfg, code_state, snap, t_phase), 'VariableNames', {'clave', 'valor'}));
    cs_write_sheet(xf, 'checks_phaseA', A);
    cs_write_sheet(xf, 'labels_snapshot', LC);
    cs_write_sheet(xf, 'excluded_edfs', removevars(rv.excluded, {'start_posix', 'end_posix'}));
    cs_write_sheet(xf, 'iid_redetect_v2', RD);
    cs_write_sheet(xf, 'review_scope', S.units);
    fl = cell2table([f_A; f_excl; S.flags], 'VariableNames', {'tipo', 'tema', 'detalle'});
    cs_write_sheet(xf, 'flags', fl);
    save(fullfile(cfg.out_dir, 'phaseA_results.mat'), 'cfg', 'snap', 'LC', 'A', 'RD', 'rv', 'S', 'O', 'Z', 'edf_tbl', 'ch_tbl', '-v7.3');
    fprintf('PHASE A STOP (%s): %s\n', xf, strjoin(A.punto(A.stop)', ' | '));
    return;
end

%% ============================ PHASE B (v2) = v1 phases 1-3 ================
tB = tic;
R = cs_seizure_stats(cfg, S, O);
P = cs_detector_performance(cfg, S, O, rv);
if cfg.run_threshold_diag
    edf_diag = edf_tbl;
    if cfg.include_128s_0729_reprocessed   % their *_seizures.mat are not in the C2 folder
        [edf_diag, f_rp] = rp128_integrate('threshold_diag_edfs', cfg, edf_tbl);
        S.flags = [S.flags; f_rp];
    end
    [TD, TDnd] = cs_threshold_diag(cfg, S, rv, edf_diag);
else
    TD = table(); TDnd = table();
end

%% ---- flags -------------------------------------------------------------
flags = [f_A; f_excl; S.flags; O.flags];
u = S.units(~S.units.included, :);
for i = 1:height(u)
    flags(end+1, :) = {'decidir', 'fuera de alcance (sin revision)', sprintf('%s %s: %s, %d/%d eventos sin etiqueta', ...
        u.campaign{i}, u.mouse{i}, u.status{i}, u.n_events_unlabeled(i), u.n_events(i))}; %#ok<SAGROW>
end
ins = ismember(strcat(det_rows.campaign, '|', det_rows.mouse), strcat(S.units.campaign(S.units.included), '|', S.units.mouse(S.units.included)));
x = det_rows(ins & startsWith(det_rows.comment, 'unrecognized:'), :);
for i = 1:height(x)
    flags(end+1, :) = {'revisar', 'AF no reconocido (R6/R7)', sprintf('%s %s %s #%d %s: "%s" (AE=%s)', x.campaign{i}, ...
        x.mouse{i}, x.edf{i}, x.seizure_id(i), x.region{i}, x.af_raw{i}, x.label{i})}; %#ok<SAGROW>
end
x = det_rows(ins & ~cellfun(@isempty, det_rows.comment_mapped), :);
if ~isempty(x)
    flags(end+1, :) = {'info', 'AF mapeado por decision', sprintf('%d filas "%s" -> other (R4)', height(x), strjoin(unique(x.comment_mapped)', ','))};
end
gaps_r5 = [];   % gaps between consecutive 'disconnect' FPs of an EDF (rule R5 bound)
ein = S.events(S.events.in_scope, :);
[gk, ~] = findgroups(ein(:, {'campaign', 'mouse', 'edf'}));
for j = 1:max([gk; 0])
    e = sortrows(ein(gk == j, :), 'start_posix');
    isd = strcmp(e.label, 'FP') & contains(e.comments, 'disconnect');
    k = find(isd(1:end-1) & isd(2:end));
    gaps_r5 = [gaps_r5; e.start_posix(k + 1) - e.end_posix(k)]; %#ok<AGROW>
end
flags(end+1, :) = {'decidir', 'R5: que es "consecutivo"', sprintf(['%d pares de FP disconnect consecutivos: %d a <= %g s (unidos en un tramo), ' ...
    '%d entre %g s y 10 min, %d > 10 min (max %.1f h). Sin limite, 070-s 20260807 quedaba 84 %% desconectado.'], numel(gaps_r5), ...
    sum(gaps_r5 <= cfg.disconnect_run_max_gap_s), cfg.disconnect_run_max_gap_s, sum(gaps_r5 > cfg.disconnect_run_max_gap_s & gaps_r5 <= 600), ...
    cfg.disconnect_run_max_gap_s, sum(gaps_r5 > 600), max([gaps_r5; 0]) / 3600)};
if ~isempty(TD)
    for i = find(TD.flag_change_ge_5pct)'
        flags(end+1, :) = {'revisar', 'umbral de crisis cambia >= 5 %', sprintf('%s %s %s %s: desconectado %.1f %%, umbral energia %+.1f %%, mediana LL %+.1f %%', ...
            TD.campaign{i}, TD.mouse{i}, TD.edf{i}, TD.region{i}, 100 * TD.disconnected_frac(i), ...
            TD.energy_thr_change_pct(i), TD.ll_median_change_pct(i))}; %#ok<SAGROW>
    end
end
%% ============================ PHASE 2 ===================================
% IID exclusion zones (5.3) + reproducibility test of the detect_iid wrapper
if cfg.phase >= 2
    disp(Z.rule_summary);   % Z: phase A
    [RT, RTcfg] = cs_iid_repro_test(cfg, ch_tbl, cfg.provenance);
    disp(RT(:, {'mouse', 'region', 'variante', 'baseline_uV_orig', 'baseline_uV_wrapper', 'n_complejos_orig', ...
        'n_complejos_wrapper', 'analyzed_min_orig', 'analyzed_min_wrapper', 'reproduce', 'carga_s', 'deteccion_s'}));
    disp(RTcfg);
    flags = [flags; Z.flags];
    p = Z.channels;
    nr = ~p.run & p.unit_in_scope;
    for i = find(nr)'
        flags(end+1, :) = {'info', 'canal en alcance no re-detectado', sprintf('%s %s %s %s: %s', p.campaign{i}, p.mouse{i}, ...
            p.edf{i}, p.region{i}, p.not_run_reason{i})}; %#ok<SAGROW>
    end
    if strcmp(cfg.iid_rerun_scope, 'in_scope')
        flags(end+1, :) = {'decidir', 'alcance de la re-deteccion IID', sprintf(['%d de %d canales (unidades revisadas, sin EDF descartado). ' ...
            'Las %d unidades sin revision (%d canales) no tienen etiquetas: todas sus detecciones serian R6. ' ...
            'Con cfg.iid_rerun_scope = ''all'' se re-detectan tambien.'], sum(p.run), height(p), sum(~S.units.included), ...
            sum(~p.unit_in_scope))};
    end
    if ~cfg.iid_rerun_figures
        flags(end+1, :) = {'info', 'figuras IID por canal', ['detect_iid no genera su panorama/polyspike por canal ' ...
            '(cfg.iid_rerun_figures = false; ~4 min/canal). Resultados identicos; se pueden generar para canales puntuales.']};
    end
    for i = find(~strcmp(RTcfg.diferencias_cfg_iid, 'identico'))'
        flags(end+1, :) = {'revisar', 'cfg.iid de la campaign distinto del actual', sprintf('%s %s: %s', RTcfg.campaign{i}, ...
            RTcfg.mouse{i}, RTcfg.diferencias_cfg_iid{i})}; %#ok<SAGROW>
    end
end

%% ============================ PHASE 3 ===================================
% IID re-detection on every planned channel (5.3), sections 6.4-6.5, checks (8)
if cfg.phase >= 3
    if ~all(RT.reproduce)
        error('run_campaign_stats:Repro', 'the reproducibility test (5.3) failed: no re-detection');
    end
    C = cs_iid_rerun_all(cfg, Z, 0);   % reuses v1 / finished channels (same zones + code), runs the rest
    fprintf('IID re-detection: %d runs (%d reused from an earlier run, %d reused, %d run now)\n', height(C), ...
        sum(strcmp(C.status, 'reutilizado de una corrida anterior')), sum(strcmp(C.status, 'reutilizado')), sum(strcmp(C.status, 'corrido')));
    if cfg.include_128s_0729_reprocessed   % the reprocessed channels' own 04_iid results
        [Corig, f_orig] = rp128_integrate('iid_original', cfg, Z.channels(Z.channels.run, :));
    else
        [Corig, f_orig] = cs_iid_original(cfg, Z.channels(Z.channels.run, :));
    end
    I = cs_iid_stats(cfg, S, R, Z, C, Corig);
    I.mouse.seizure_rate_per_24h = NaN(height(I.mouse), 1);
    for i = 1:height(I.mouse)
        k = strcmp(R.mouse.mouse, I.mouse.mouse{i}) & strcmp(R.mouse.scope, I.mouse.scope{i});
        if any(k), I.mouse.seizure_rate_per_24h(i) = R.mouse.rate_per_24h(k); end
    end
    [CR, f_corr] = cs_iid_correlations(cfg, I, R);
    [chk3, f_chk, CRI, fpw] = cs_iid_checks(cfg, Z, C, I, O);
    flags = [flags; f_orig; f_corr; f_chk];
    fl = I.channels(I.channels.p97_real_uV < 20, :);
    for i = 1:height(fl)
        flags(end+1, :) = {'decidir', 'senal plana en canal re-detectado', sprintf(['%s %s %s %s: P97 real %.1f uV en %.1f h ' ...
            'analizadas (%d complejos). Sin senal: cuenta como tiempo observado de crisis (T_seizure_obs) y de IIDs (T_iid). ' ...
            'Agregarlo como desconexion?'], fl.campaign{i}, fl.mouse{i}, fl.edf{i}, fl.region{i}, fl.p97_real_uV(i), ...
            fl.T_iid_h(i), fl.n_complexes(i))}; %#ok<SAGROW>
    end
    g = I.group(strcmp(I.group.metrica, 'complexes_per_min'), :);
    for i = 1:height(g)
        flags(end+1, :) = {'info', 'criterio >= 10x (control vs epileptico)', sprintf(['%s, %s: fold change de medianas %.1f ' ...
            '(IC95 %.1f-%.1f); criterio %s'], g.datos{i}, g.conjunto{i}, g.fold_change_medianas(i), g.fc_ci95_lo(i), ...
            g.fc_ci95_hi(i), string(g.criterio_10x_fold(i)))}; %#ok<SAGROW>
    end
end
t_phase.B = toc(tB);
run_phase = @(p) find('ABCDE' == p) <= find('ABCDE' == cfg.stop_after);

%% ============================ PHASE C (v2) ================================
% one pass per channel: flat segments + IID threshold sensitivity (section 5)
if run_phase('C')
    tC = tic;
    [SensT, SensV, SensF] = cs_iid_sens_all(cfg, Z, C, 0);   % background workers' cache; runs what is missing
    [BC, SM, SG, FL] = cs_iid_sens_stats(cfg, SensT, SensV, SensF, S, R, C, Z);
    for i = find(~SensT.actual_reproduce)'
        flags(end+1, :) = {'revisar', 'sensibilidad: variante actual no reproduce', sprintf('%s %s %s %s', SensT.campaign{i}, ...
            SensT.mouse{i}, SensT.edf{i}, SensT.region{i})}; %#ok<SAGROW>
    end
    nfl = sum(~strcmp(FL.cubierto, 'si'));
    flags(end+1, :) = {'decidir', 'tramos planos (|x| < 5 uV > 60 s)', sprintf(['%d tramos en %d canales (%.1f h); %d no cubiertos por ' ...
        'una exclusion (%.1f h). Nada se excluyo: hoja flatline_scan'], height(FL), height(unique(FL(:, {'mouse', 'edf', 'region'}))), ...
        sum(FL.duracion_s) / 3600, nfl, sum(FL.duracion_s(~strcmp(FL.cubierto, 'si'))) / 3600)};
    t_phase.C = toc(tC);
end

%% ============================ PHASE D (v2) ================================
% days of 24 h, daily association with post-ictal exclusion, peri-ictal with lead seizures (section 6)
if run_phase('D')
    tD = tic;
    K = cs_iid_channel_intervals(C, S);
    [DW, Mday] = cs_day_windows(cfg, S, O, rv, K, edf_tbl, edf_lists);
    [DA, DAd] = cs_iid_daily_assoc(cfg, S, K, DW, Mday);
    [PMa, PSa] = cs_periictal(cfg, S, O, K, false);
    [PMl, PSl, no_lead] = cs_periictal(cfg, S, O, K, true);
    I.periictal = [PMa; PMl];
    PL = outerjoin(PSa(:, 2:end), PSl(:, 2:end), 'Keys', 'ventana', 'MergeKeys', true);
    PL.Properties.VariableNames = regexprep(PL.Properties.VariableNames, '_left$', '_todas');
    PL.Properties.VariableNames = regexprep(PL.Properties.VariableNames, '_right$', '_lideres');
    PL.ratones_sin_lideres = repmat({strjoin(no_lead, ', ')}, height(PL), 1);
    CR = CR(strcmp(CR.analisis, 'entre animales'), :);   % within-animal: iid_daily_assoc (24-h days)
    flags(end+1, :) = {'info', 'peri-ictal: ratones sin crisis lideres', sprintf('%d: %s (fuera del analisis de lideres)', ...
        numel(no_lead), strjoin(no_lead, ', '))};
    t_phase.D = toc(tD);
end

%% ============================ PHASE E (v2) ================================
% section 7 figures (+ the v1 ones, regenerated at the end)
if run_phase('E')
    tE = tic;
    cfg.raster_v2 = true;
    [figs_v2, figchk] = cs_figures_v2(cfg, S, R, rv, DW, Mday, fullfile(cfg.out_dir, 'figures'));
    figs_v2{end+1} = cs_figure_sensitivity(cfg, SG, SM, fullfile(cfg.out_dir, 'figures'));
    t_phase.E = toc(tE);
end

flag_tbl = cell2table(flags, 'VariableNames', {'tipo', 'tema', 'detalle'});
% decisions the user already took on flags (phase reports, 2026-10-29)
flag_tbl.estado = repmat({'abierto'}, height(flag_tbl), 1);
for i = 1:size(cfg.flag_decisions, 1)
    k = startsWith(flag_tbl.tema, cfg.flag_decisions{i, 1});
    flag_tbl.estado(k) = cfg.flag_decisions(i, 2);
end
flag_tbl.estado(strcmp(flag_tbl.tipo, 'info')) = {'informativo'};

%% ---- consistency checks (section 8, phase-1 part) ------------------------
chk = {};
k110 = 'C1|110-s|20260803~ TRAP_2bbce7f4-2ccd-4764-b26b-42267b0d33f6|15';
n110 = sum(contains(S.seizures.origin_ids, k110));
chk(end+1, :) = {'110-s crisis 15 aparece una sola vez', n110 == 1, sprintf('%d fila(s)', n110)};
e114 = S.events(strcmp(S.events.mouse, '114-s') & startsWith(S.events.edf, '20260804~'), :);
d11 = rv.dx(rv.dx.id == "DX-11", :);
in11 = e114.start_posix >= d11.start_posix & e114.end_posix <= d11.end_posix;
chk(end+1, :) = {'114-s EDF 20260804: 0 TP y todas sus detecciones en DX-11', ~any(strcmp(e114.label, 'TP')) && all(in11), ...
    sprintf('%d eventos, %d TP, %d dentro de DX-11', height(e114), sum(strcmp(e114.label, 'TP')), sum(in11))};
s72 = S.seizures(strcmp(S.seizures.mouse, '072-s'), :);
o72 = O.mouse(strcmp({O.mouse.mouse}, '072-s') & strcmp({O.mouse.scope}, 'merged'));
okv = all(ismember({'ND-05', 'ND-06'}, s72.origin_ids(strcmp(s72.source, 'video'))));
chk(end+1, :) = {'072-s: ND-05/ND-06 video y T_seizure_obs > T_eeg', okv && cs_iv('length', o72.T_sz) > cs_iv('length', o72.T_eeg), ...
    sprintf('T_eeg %.2f h, T_seizure_obs %.2f h', cs_iv('length', o72.T_eeg) / 3600, cs_iv('length', o72.T_sz) / 3600)};
chk(end+1, :) = {'001-s incluye ND-01; 073-s incluye ND-02..04', any(strcmp(S.seizures.origin_ids, 'ND-01') & strcmp(S.seizures.mouse, '001-s')) && ...
    all(ismember({'ND-02', 'ND-03', 'ND-04'}, S.seizures.origin_ids(strcmp(S.seizures.mouse, '073-s')))), ''};
bad = O.edf_day.T_eeg_h > O.edf_day.rec_h + 1e-9 | O.edf_day.T_seizure_obs_h > O.edf_day.rec_h + 1e-9;
chk(end+1, :) = {'T_eeg <= T_seizure_obs <= registro, por raton x EDF x dia', ~any(bad) && all(O.edf_day.T_eeg_h <= O.edf_day.T_seizure_obs_h + 1e-9), ...
    sprintf('%d filas fuera', sum(bad))};
chk(end+1, :) = {'crisis por video = 2 (ND-05, ND-06)', sum(strcmp(S.seizures.source, 'video')) == 2, ''};
c13 = cfg.crisette(strcmp(cfg.crisette(:, 1), 'O-13'), :);
n13 = sum(strcmp(S.events.campaign, c13{3}) & strcmp(S.events.mouse, c13{2}) & strcmp(S.events.edf, c13{4}));
chk(end+1, :) = {'O-13 (084-s C2 20260713): una sola deteccion candidata', n13 == 1, sprintf('%d deteccion(es) en el EDF; clasificada #%d', n13, c13{5})};
for k = 1:height(rv.excluded)
    x = rv.excluded(k, :);
    n_t = sum(strcmp(O.edf_day.mouse, x.mouse{1}) & strcmp(O.edf_day.edf, x.edf{1}));
    n_s = sum(strcmp(S.seizures.mouse, x.mouse{1}) & contains(S.seizures.edf, x.edf{1}));
    n_i = sum(Z.channels.run & strcmp(Z.channels.mouse, x.mouse{1}) & strcmp(Z.channels.edf, x.edf{1}));
    n_p = sum(S.events.usable & strcmp(S.events.mouse, x.mouse{1}) & strcmp(S.events.edf, x.edf{1}));
    chk(end+1, :) = {sprintf('EDF excluido %s (%s %s) no aporta tiempo, crisis, IIDs ni desempeno', x.id{1}, x.mouse{1}, x.edf{1}(1:10)), ...
        n_t + n_s + n_i + n_p == 0, sprintf('%d filas de tiempo, %d crisis, %d canales IID, %d eventos en desempeno', n_t, n_s, n_i, n_p)}; %#ok<SAGROW>
end
if cfg.phase >= 2
    chk(end+1, :) = {'prueba de reproducibilidad 5.3 (wrapper vs 04_iid C2)', all(RT.reproduce), sprintf('%d/%d corridas (%d canales x variantes A/B) identicas', ...
        sum(RT.reproduce), height(RT), height(RT) / 2)};
end
if cfg.phase >= 3
    chk = [chk; chk3];
end
if run_phase('C')
    chk(end+1, :) = {'sensibilidad: la variante actual reproduce exactamente la re-deteccion en todos los canales', ...
        all(SensT.actual_reproduce), sprintf('%d/%d canales', sum(SensT.actual_reproduce), height(SensT))};
end
if run_phase('D')
    om = O.mouse(strcmp({O.mouse.scope}, 'merged'));
    dmax = 0;
    for i = 1:numel(om)
        dmax = max(dmax, abs(sum(DW.T_seizure_obs_h(strcmp(DW.mouse, om(i).mouse))) * 3600 - cs_iv('length', om(i).T_sz)));
    end
    chk(end+1, :) = {'day_windows: suma de T_seizure_obs por raton = total (tolerancia 1 s)', dmax <= 1, sprintf('max |diferencia| = %.3f s', dmax)};
    for k = find(~strcmp(rv.excluded.kind, 'sin AE completo'))'
        x = rv.excluded(k, :);
        sp = [x.start_posix x.end_posix];
        m = Mday(strcmp({Mday.mouse}, x.mouse{1}));
        n_sz = 0; h_sz = 0; h_ex = 0;
        if ~isempty(m)
            h_sz = cs_iv('length', cs_iv('intersect', m.T_sz, sp));
            h_ex = cs_iv('length', cs_iv('intersect', [m.excluded; m.pending], sp));
            n_sz = sum(m.seizures.start_posix >= sp(1) & m.seizures.start_posix < sp(2));
        end
        n_iid = sum(strcmp(K.mouse, x.mouse{1}) & strcmp(K.edf, x.edf{1}));
        chk(end+1, :) = {sprintf('%s %s: sin tiempo, crisis ni IIDs en day_windows; el raster lo dibuja excluido/pendiente', x.mouse{1}, ...
            x.edf{1}(1:10)), h_sz == 0 && n_sz == 0 && n_iid == 0 && abs(h_ex - diff(sp)) < 1, sprintf(['%.0f s observados, %d crisis, ' ...
            '%d canales IID; %.1f de %.1f h dibujadas como excluidas/pendientes'], h_sz, n_sz, n_iid, h_ex / 3600, diff(sp) / 3600)}; %#ok<SAGROW>
    end
end
if run_phase('E')
    nod = cellfun(@isempty, figchk(:, 2));
    chk(end+1, :) = {'ninguna figura de la seccion 7 contiene fechas calendario', all(nod), sprintf('%d/%d figuras sin fechas; %s', ...
        sum(nod), numel(nod), strjoin(strcat(figchk(~nod, 1), ': ', figchk(~nod, 2))', ' || '))};
    nv = sum(DW.enters_violin & strcmp(DW.group, 'epileptic'));
    chk(end+1, :) = {'raton-dias del violin por dia <= filas de day_windows con datos', nv <= sum(DW.rec_h > 0), ...
        sprintf('%d raton-dias en el violin, %d filas con registro', nv, sum(DW.rec_h > 0))};
end
chk_tbl = cell2table(chk, 'VariableNames', {'control', 'pasa', 'detalle'});
disp(chk_tbl);

%% ---- outputs -------------------------------------------------------------
xf = fullfile(cfg.out_dir, ['stats_campaigns_' cfg.run_tag '.xlsx']);
if isfile(xf), delete(xf); end
grp = S.units(S.units.included, :);
readme = [readme_v2(cfg, code_state, snap, t_phase); { ...
    'C1', cfg.campaigns(1).root; 'C2', cfg.campaigns(2).root; ...
    'unidades raton x campaign en alcance', strjoin(strcat(grp.campaign, ':', grp.mouse)', ', '); ...
    'n ratones epilepticos', sprintf('%d: %s', numel(unique(grp.mouse(strcmp(grp.group, 'epileptic')))), strjoin(unique(grp.mouse(strcmp(grp.group, 'epileptic')))', ', ')); ...
    'n ratones control', sprintf('%d: %s', numel(unique(grp.mouse(strcmp(grp.group, 'control')))), strjoin(unique(grp.mouse(strcmp(grp.group, 'control')))', ', ')); ...
    'chequeo tiempo absoluto', det_log{1}}];
if run_phase('C')
    readme(end+1, :) = {'nota: sensibilidad del umbral IID', ['La prominencia minima de findpeaks (200 uV) es la misma en todas ' ...
        'las variantes. Con umbrales bajos (por debajo de ~200 uV) ese criterio es el que limita: en varios canales las variantes ' ...
        'P99 / P99.5 de |x_bp| dan exactamente el mismo resultado.']};
    readme(end+1, :) = {'tiempo de ejecucion fase C (pasada por canal)', sprintf(['%.1f h sumando los %d canales (corridos en 4 ' ...
        'procesos MATLAB en paralelo, cs_iid_sens_worker.m); "fase C" arriba = solo la lectura del cache (%s, reutilizado solo con ' ...
        'el mismo codigo y las mismas zonas)'], sum(SensT.tiempo_s) / 3600, height(SensT), fullfile(cfg.out_dir, 'cache', 'sens'))};
end
if run_phase('D')
    readme(end+1, :) = {'nota: dias (v2)', sprintf(['Dia k = [t0 + 24(k-1) h, t0 + 24k h), t0 = primera muestra del primer EDF del ' ...
        'raton (cualquier campaign y estado, tambien no analizados: header del EDF). %d raton-dias de epilepticos entran en la ' ...
        'asociacion (>= %g h de T_iid y de T_seizure_obs; v1: 47 dias calendario)'], sum(DW.enters_assoc & strcmp(DW.group, 'epileptic')), ...
        cfg.min_valid_h_day)};
end
if cfg.phase >= 2
    for i = 1:height(RTcfg)
        readme(end+1, :) = {sprintf('cfg.iid de %s %s (config_used.mat) vs pipeline_config actual', RTcfg.campaign{i}, RTcfg.mouse{i}), ...
            RTcfg.diferencias_cfg_iid{i}}; %#ok<SAGROW>
    end
end
cs_write_sheet(xf, 'README', cell2table(readme, 'VariableNames', {'clave', 'valor'}));
cs_write_sheet(xf, 'checks_phaseA', A);
cs_write_sheet(xf, 'labels_snapshot', LC);
cs_write_sheet(xf, 'excluded_edfs', removevars(rv.excluded, {'start_posix', 'end_posix'}));
cs_write_sheet(xf, 'iid_redetect_v2', RD);
cs_write_sheet(xf, 'seizures_final', removevars(S.seizures, {'start_posix', 'end_posix'}));
cs_write_sheet(xf, 'crisette_focal', removevars(S.crisette, 'start_posix'));
cs_write_sheet(xf, 'observation_time', O.edf_day);
cs_write_sheet(xf, 'seizure_stats_mouse', [R.mouse; outerjoin_fill(R.group, R.mouse)]);
cs_write_sheet(xf, 'seizure_daily', R.daily);
cs_write_sheet(xf, 'detector_performance', P);
cs_write_sheet(xf, 'seizure_threshold_diag', TD);
cs_write_sheet(xf, 'seizure_threshold_ND', TDnd);
cs_write_sheet(xf, 'circadian', R.circ_tests);
cs_write_sheet(xf, 'circadian_hourly', R.circ_hourly);
cs_write_sheet(xf, 'review_scope', S.units);
cs_write_sheet(xf, 'close_pairs', S.close_pairs);
cs_write_sheet(xf, 'exclusion_log', S.exclusion_log);
if cfg.phase >= 2
    cs_write_sheet(xf, 'iid_exclusion_zones', removevars(Z.zones, {'start_posix', 'end_posix'}));
    cs_write_sheet(xf, 'iid_zone_summary', Z.rule_summary);
    p = Z.channels;
    p.duration_h = p.duration_s / 3600;
    cs_write_sheet(xf, 'iid_channel_plan', p(:, {'campaign', 'mouse', 'edf', 'region', 'unit_in_scope', 'run', 'not_run_reason', ...
        'duration_h', 'n_zones_seizures', 'n_zones_manual', 'excluded_h', 'expected_T_iid_h'}));
    if cfg.phase < 3
        cs_write_sheet(xf, 'iid_fp_iids_windows', removevars(Z.fp_iids, {'start_posix', 'end_posix'}));
    end
    cs_write_sheet(xf, 'iid_repro_test', RT);
end
if cfg.phase >= 3
    cs_write_sheet(xf, 'iid_by_mouse_edf', I.channels);
    cs_write_sheet(xf, 'iid_before_after', I.before_after);
    cs_write_sheet(xf, 'iid_stats', I.mouse);
    cs_write_sheet(xf, 'group_comparison', I.group);
    cs_write_sheet(xf, 'correlations', CR);
    cs_write_sheet(xf, 'periictal', I.periictal);
    cs_write_sheet(xf, 'iid_circadian', I.circ_tests);
    cs_write_sheet(xf, 'iid_circadian_hourly', I.circ_hourly);
    cs_write_sheet(xf, 'iid_daily', removevars(I.day, 'daykey'));
    cs_write_sheet(xf, 'iid_crisette', CRI);
    cs_write_sheet(xf, 'iid_fp_iids_windows', fpw);
end
if run_phase('D')
    cs_write_sheet(xf, 'v1_vs_v2', cs_v1_vs_v2(cfg, S, R, I, CR, DA, PSa, PSl, DW, rv));
    cs_write_sheet(xf, 'day_windows', DW);
    cs_write_sheet(xf, 'iid_daily_assoc', DA);
    cs_write_sheet(xf, 'periictal_lead', PL);
end
if run_phase('C')
    cs_write_sheet(xf, 'flatline_scan', FL);
    cs_write_sheet(xf, 'iid_threshold_by_channel', BC);
    cs_write_sheet(xf, 'iid_threshold_sensitivity', SM);
    cs_write_sheet(xf, 'iid_threshold_summary', SG);
end
cs_write_sheet(xf, 'flags', flag_tbl);
cs_write_sheet(xf, 'checks', chk_tbl);
figs = cs_figures_seizures(cfg, S, O, R, fullfile(cfg.out_dir, 'figures'));
save(fullfile(cfg.out_dir, 'phase1_results.mat'), 'cfg', 'edf_tbl', 'ch_tbl', 'det_rows', 'review_rows', 'events', ...
    'rv', 'S', 'O', 'R', 'P', 'TD', 'TDnd', 'flag_tbl', 'chk_tbl', '-v7.3');
if cfg.phase >= 2
    save(fullfile(cfg.out_dir, 'phase2_results.mat'), 'cfg', 'Z', 'RT', 'RTcfg', '-v7.3');
end
if cfg.phase >= 3
    figs = [figs, cs_figures_iid(cfg, I, CR, fullfile(cfg.out_dir, 'figures'))];
    save(fullfile(cfg.out_dir, 'phase3_results.mat'), 'cfg', 'C', 'Corig', 'I', 'CR', 'CRI', 'fpw', 'chk3', '-v7.3');
end
if run_phase('E')
    figs = [figs, figs_v2];
end
if run_phase('D')
    save(fullfile(cfg.out_dir, 'phaseCDE_results.mat'), 'cfg', 'SensT', 'SensV', 'SensF', 'BC', 'SM', 'SG', 'FL', 'K', 'DW', 'Mday', ...
        'DA', 'DAd', 'PMa', 'PSa', 'PMl', 'PSl', 'no_lead', 'PL', '-v7.3');
end
fprintf('%s written (phases through %s): %s (+ %d figures)\n', cfg.run_tag, cfg.stop_after, xf, numel(figs));

function r = readme_v2(cfg, code_state, snap, t_phase)
% README rows common to every v2 output: run, code, label snapshot, CONFIG, timings.
    r = {'fecha de corrida', cfg.run_started; 'corrida', cfg.run_tag; 'resultados v1 (referencia)', cfg.v1_dir; ...
        'git HEAD', code_state.head; 'arbol con cambios sin commit', mat2str(cfg.git.dirty); ...
        'codigo usado', sprintf('%s (git diff HEAD + %d archivos sin rastrear de src/stats y run_campaign_stats.m)', ...
            code_state.patch_file, numel(code_state.untracked)); ...
        'snapshot de etiquetas', sprintf('%s: %d archivos, tomado %s', cfg.labels_snapshot_dir, height(snap), snap.snapshot_taken{1}); ...
        'codigo IID (detect_iid y dependencias) sin cambios respecto de HEAD', mat2str(cfg.git.iid_code_unmodified)};
    for i = 1:size(cfg.git.iid_code, 1)
        r(end+1, :) = {['sha1 (IID) ' cfg.git.iid_code{i, 1}], cfg.git.iid_code{i, 2}}; %#ok<AGROW>
    end
    for i = 1:size(cfg.git.files, 1)
        r(end+1, :) = {['sha1 ' cfg.git.files{i, 1}], cfg.git.files{i, 2}}; %#ok<AGROW>
    end
    for f = fieldnames(t_phase)'
        r(end+1, :) = {['tiempo de ejecucion fase ' f{1}], sprintf('%.1f min', t_phase.(f{1}) / 60)}; %#ok<AGROW>
    end
    fn = fieldnames(cfg);
    for i = 1:numel(fn)
        if any(strcmp(fn{i}, {'git', 'campaigns', 'provenance', 'iid_rerun_dirs'})), continue; end
        v = cfg.(fn{i});
        if iscell(v), v = strjoin(cellfun(@(z) strtrim(evalc('disp(z)')), v(:)', 'UniformOutput', false), ' | ');
        elseif islogical(v) || isnumeric(v), v = mat2str(v);
        elseif ~ischar(v), v = strtrim(evalc('disp(v)')); end
        r(end+1, :) = {['cfg.' fn{i}], v}; %#ok<AGROW>
    end
end

function T = outerjoin_fill(G, M)
% Group rows under the per-mouse column layout (missing columns -> NaN).
    T = M([], :);
    for i = 1:height(G)
        r = M(1, :);
        for v = M.Properties.VariableNames
            if ismember(v{1}, G.Properties.VariableNames)
                r.(v{1}) = G.(v{1})(i);
            elseif isnumeric(r.(v{1}))
                r.(v{1}) = NaN;
            else
                r.(v{1}) = {''};
            end
        end
        r.rate_per_24h = G.pooled_rate_per_24h(i);
        T = [T; r]; %#ok<AGROW>
    end
end
