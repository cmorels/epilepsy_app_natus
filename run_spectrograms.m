%% RUN_SPECTROGRAMS  Spectrogram module (module 3 -> 06_spectrograms).
% The ONLY file to edit to choose what is processed: paths, events,
% parameters. Runs on outputs already on disk (02_clean, optionally
% 05_summaries / 04_iid); it never runs detection nor reads an EDF.
%
% 1. Fill A / B / C below (any combination).
% 2. Run with cfg.spec.dry_run = true: prints and writes
%    <output_root>\06_spectrograms_summary\_selection_preview.csv
%    (files, t_rel, Natus / pipeline clock, NaN %, peak source, status).
% 3. Check it, set dry_run = false and run again.
% See README.md, "Spectrogram module".

repo = fileparts(mfilename('fullpath'));
addpath(fullfile(repo, 'src'), fullfile(repo, 'src', 'utils'), fullfile(repo, 'src', 'spectrogram'));

cfg = spectrogram_config();

% Folder holding one subfolder per animal (<root>\<subject>\01_txt, 02_clean, ... or one folder per EDF).
% cfg.spec.pipeline_root = '\\iws\epic\users\camila.morel\2-Analysis\ECoG\batch 2026';
cfg.spec.pipeline_root = repo;              % local test animal: <repo>\073-s
cfg.spec.output_root   = '';                % '' = pipeline_root
cfg.spec.dry_run       = true;              % first run: preview only
cfg.spec.signal_source = 'clean';           % 'clean' | 'raw' (01_txt x gain, no outlier interpolation)
% Optional: seizures_events used for exclusion zones (IID) and baseline checks of manual events
cfg.spec.seizures_file = fullfile(repo, '073-s', 'seizures_events.csv');

events = struct([]);

%% A. Manual table
% tipo:       'seizure' | 'iid'
% subject:    exact ('005-s' ~= '005')
% sesion:     'yyyyMMdd_HHmmss' or the EDF file name
% region:     'HPCr' | 'HPCl' | 'both'
% inicio/fin: 'dd/MM/yyyy HH:mm:ss.SSS' in time_ref, or a number = t_rel (s)
% time_ref:   'natus' (as read in Natus) | 'pipeline' (as in seizures_events) | 't_rel'
% optional 9th column: per-event overrides, struct('win_s', 4) or 'win_s=4; fmax=60'
manual = {
% id        tipo       subject  sesion              region  inicio                     fin                        time_ref
 'zona1'   'iid'      '073-s'  '20260709B~ TRA_323846a4-5a4c-42a4-95ee-512ff8c1c068'  'both'  '09/07/2026 17:49:03.000'  '09/07/2026 22:30:00.000'  'natus'  struct('spike_source', 'detect_local')
 'zona2'   'iid'      '073-s'  '20260709B~ TRA_323846a4-5a4c-42a4-95ee-512ff8c1c068'  'both'  '10/07/2026 09:38:32.000'  '10/07/2026 15:43:55.000'  'natus'  struct('spike_source', 'detect_local')
};
% IID intervals longer than cfg.spec.iid.max_interval_s (600 s) run in 'zone' mode: spike-aligned
% maps over every spike of the zone + spike-rate overview (no continuous scalogram).
events = [events, spec_read_event_list(manual)];

% Same columns in an Excel / CSV (extra columns by name = overrides):
% events = [events, spec_read_event_list(fullfile(repo, 'spectrogram_events.xlsx'))];

%% B. Detected seizures by ID (onset/offset = start_s/end_s of seizures_events)
selB = struct();
selB.seizures_file = fullfile(repo, '073-s', 'seizures_events.csv');   % your reviewed copy (user_check / comments)
selB.region = 'both';                                                     % one figure with HPCr + HPCl
selB.ids = {'073-s', '20260706_130048', [1 2]                             % {subject, sesion, [seizure_id]; ...}
            '073-s', '20260709B~ TRA_323846a4-5a4c-42a4-95ee-512ff8c1c068.EDF', 1};   % EDF name; its txt is not in 073-s
events = [events, spec_select_seizures(selB, cfg)];

%% C. Detected seizures by filter
selC = struct();
selC.seizures_file = fullfile(repo, '073-s', 'seizures_events.csv');
selC.region = 'both';
selC.filter.subject    = {'073-s'};
selC.filter.sesion     = {'20260708_134229'};
selC.filter.user_check = 'TP';
selC.filter.max_n      = 1;
% selC.filter.date_range = {'08/07/2026 00:00:00', '09/07/2026 00:00:00'};   % pipeline clock
% selC.filter.duration_s = [10 60];
events = [events, spec_select_seizures(selC, cfg)];

%% Parameters (defaults: src/spectrogram/spectrogram_config.m)
% cfg.spec.norm = 'db_baseline';           % 'db' | 'db_baseline' | 'zscore' | 'percent'
% cfg.spec.clim = [-10 10];                % fixed colour limits to compare animals
% cfg.spec.freq_scale = 'log';
% cfg.spec.seizure.baseline = [-30 -5];    % or {'06/07/2026 13:20:00.000', '06/07/2026 13:20:25.000'}
% cfg.spec.iid.spike_source = 'pipeline';  % 'pipeline' | 'detect_local' | 'manual'

out = spec_run(cfg, events);
