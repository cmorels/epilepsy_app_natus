function cfg = pipeline_config()
% PIPELINE_CONFIG  Single source of truth for every parameter used by the
% EDF-based pipeline (src/edf_import.m, load_lfp_txt.m, clean_lfp.m,
% detect_seizures.m, detect_iid.m, run_pipeline_edf.m).
%
% Values marked "ported from X" are copied verbatim from the original
% Intan-era scripts and must NOT be changed without a separate, explicit
% methodology decision (see KNOWN_ISSUES.md).

cfg = struct();

%% ---- General / time --------------------------------------------------
cfg.general.timezone = 'Europe/Paris';                    % applied to every absolute timestamp
cfg.general.session_start_input_format = 'dd.MM.yy HH.mm.ss'; % Natus EDF StartDate+StartTime format
cfg.general.iso8601_format = 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX'; % header datetime format (ISO 8601 + offset)
cfg.general.overwrite = false;                             % shared "skip stage if output exists" switch

%% ---- Stage 1: EDF import (src/edf_import.m) ---------------------------
% Channel selection. mode = 'list' (region name = EDF label, no renaming),
% 'map' (containers.Map from EDF label -> region name), or 'all' (every
% signal in the file, region name = EDF label).
cfg.edf.channels.mode = 'list';
cfg.edf.channels.labels = {'A7C1', 'A7C3'};   % labels found in 097-s test EDF; edit per study
cfg.edf.channels.map = containers.Map({'EEG1', 'EEG2'}, {'HPCleft', 'HPCright'}); % TEMPLATE for mode='map' -- edit to real labels/regions before using

cfg.edf.gap_tol_s = 1e-3;   % |actual_dt - DataRecordDuration| beyond this is a gap candidate (matches edf_txt_conversion_Samara.m tol)
cfg.edf.gap_min_s = 0.5;    % candidates shorter than this are logged as jitter, not a real gap
cfg.edf.subject_id = '';    % '' -> derive from EDF filename (a warning is logged; set explicitly for real runs)
cfg.edf.output_dir = '';    % '' -> caller decides (defaults to pwd when edf_import is run standalone)

% Source-unit (EDF PhysicalDimensions, lower-cased) -> multiplier to microvolts.
cfg.edf.unit_aliases = containers.Map({'uv', 'microv', 'mv', 'v'}, {1, 1, 1000, 1e6});

%% ---- Stage 2: generic txt reader (src/load_lfp_txt.m) -----------------
% No tunable parameters: the reader is fully generic and driven by
% whatever "# key = value" header the txt file carries.

%% ---- Stage 3: outlier removal (src/clean_lfp.m) ------------------------
% Ported from complete_pipeline_seizures.m (STEP 2).
cfg.clean.use_fixed_thresholds = true;
cfg.clean.fixed_lower_uV = -2500;
cfg.clean.fixed_upper_uV = 2500;
cfg.clean.k_factor = 10;
cfg.clean.auto_switch_pct = 0.1;   % if use_fixed_thresholds=false and prelim. outliers > this %, switch to fixed
cfg.clean.output_dir = '';         % '' -> caller decides (defaults to pwd when clean_lfp is run standalone)

%% ---- Stage 4: seizure detection (src/detect_seizures.m) ---------------
% Ported from complete_pipeline_seizures.m (STEP 3).
cfg.seizure.bandpass_band = [5 75];
cfg.seizure.bandpass_steepness = 0.5;
cfg.seizure.edge_trim_s = 1;        % trim_samples = round(fs) in the original script
cfg.seizure.power_exponent = 4;
cfg.seizure.window_sec = 2;
cfg.seizure.median_factor = 10;
cfg.seizure.min_seizure_duration = 15;
cfg.seizure.zoom_margin_s = 5;      % per-seizure zoom figure margin
cfg.seizure.output_dir = '';        % '' -> caller decides (defaults to pwd when detect_seizures is run standalone)

%% ---- Stage 5: IID detection (src/detect_iid.m) -------------------------
% Ported from IID_detection_FINAL.m.
cfg.iid.baseline.start_uV = 100;
cfg.iid.baseline.step_uV = 5;
cfg.iid.baseline.target_fraction = 0.97;
cfg.iid.baseline.max_uV = 130;
cfg.iid.lower_threshold_factor = 2.5;   % lower_threshold_uV = factor * baseline_uV
cfg.iid.upper_threshold_uV = 2000;
cfg.iid.bandpass_band = [15 70];
cfg.iid.bandpass_steepness = 0.5;
cfg.iid.min_peak_distance_s = 0.030;
cfg.iid.max_peak_width_s = 0.060;
cfg.iid.min_peak_prominence_mV = 0.2;
cfg.iid.max_interspike_ms = 150;        % ISI cutoff to group spikes into one complex
cfg.iid.polyspike_min_n = 2;            % complex with >= this many spikes is a polyspike
cfg.iid.burst_max_gap_s = 5;            % max inter-complex gap to stay in the same burst
cfg.iid.burst_min_complexes = 3;        % burst kept if n_complexes > this
cfg.iid.burst_min_duration_s = 4;       % burst kept if duration_s > this
cfg.iid.burst_max_duration_s = 40;      % burst kept if duration_s < this
cfg.iid.exclusion_buffer_s = 5;         % buffer added around each exclusion zone
cfg.iid.exclusion_zones_manual = zeros(0, 2); % extra [start end] zones, unioned with seizures+gaps
cfg.iid.output_dir = '';                % '' -> caller decides (defaults to pwd when detect_iid is run standalone)

%% ---- Stage 6: batch orchestration (src/run_pipeline_edf.m) ------------
cfg.paths.output_root = fullfile(pwd, 'pipeline_output');

end
