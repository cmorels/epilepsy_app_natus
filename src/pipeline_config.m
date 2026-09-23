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
cfg.seizure.threshold_mode = 'median_factor'; % DEFAULT (unchanged behavior). Alternatives 'log_mad'/'moving_baseline' land with detect_seizures.m's threshold-mode support; see KNOWN_ISSUES.md
cfg.seizure.log_mad_k = NaN;        % k for threshold_mode='log_mad'; NaN until calibrated via calibrate_log_mad_k.m -- not used while threshold_mode='median_factor'
cfg.seizure.baseline_window_s = 300; % moving baseline window for threshold_mode='moving_baseline' -- not used while threshold_mode='median_factor'

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
cfg.iid.threshold_mode = 'absolute';    % DEFAULT (unchanged behavior: fixed 100:5:130 baseline grid). Alternative 'relative_mad' lands with detect_iid.m's threshold-mode support; see KNOWN_ISSUES.md
cfg.iid.prominence_k = NaN;             % MinPeakProminence = this * mad(x_bp,1) for threshold_mode='relative_mad' -- not used while threshold_mode='absolute'

%% ---- Signal quality assessment (src/signal_quality.m) ------------------
% Runs on every channel, always, on the RAW (pre-gain) signal -- measuring
% is free and never changes any result on its own (see KNOWN_ISSUES.md /
% README.md "sistema de casos"). Only precondition_lfp.m (gated by the
% case system) ever acts on what this measures.
cfg.quality.welch_window_s = 4;         % pwelch analysis window, with 50% overlap between windows
cfg.quality.scale_band = [15 70];       % band sigma_band_uV is measured in -- matches what detect_iid actually analyzes
cfg.quality.line_exclusion_band = [48 52]; % excluded from sigma_band_uV so mains contamination can't masquerade as signal
cfg.quality.line_freq_hz = 50;          % mains frequency; line_ratio_db bands are derived from this (49-51 / 45-48 / 52-55 at 50 Hz)
cfg.quality.line_ratio_thr_db = 10;     % line_ratio_db above this -> line-contaminated
cfg.quality.line_intermittent_pct = 10; % pct_time_line_high above this -> line-contaminated even if the median looks clean
cfg.quality.low_amp_snr_db = 26;        % snr_quantization_db below this -> low_amplitude (absolute criterion, no reference needed)
cfg.quality.min_snr_quant_db = 12;      % snr_quantization_db below this -> unusable regardless of everything else
cfg.quality.attenuation_ratio_thr = 3;  % reference/sigma_band_uV above this -> suggested_case='attenuated' in make_cases_template.m
cfg.quality.reference_sigma_uV = NaN;   % NaN (uncalibrated) by default -- scalar or containers.Map keyed by region; see estimate_reference_sigma.m
cfg.quality.recursive_scan = false;     % make_cases_template.m / estimate_reference_sigma.m: search input folder(s) recursively for *.txt
cfg.quality.min_reference_channels = 5; % warn when a region has fewer channels than this to calibrate a reference from

%% ---- Case system (src/load_cases.m, src/make_cases_template.m) --------
% Declares, per file/channel, whether precondition_lfp.m/clean_lfp.m may
% ACT on what signal_quality.m measured. Default run: every channel is
% 'normal' (gain off, notch off) -- byte-identical to the pipeline before
% this system existed. See README.md "sistema de casos".
cfg.cases.default = 'normal';    % applied when nothing else says otherwise
cfg.cases.force = '';            % '' or a case name: forces EVERY channel in the run to this case
cfg.cases.file = '';             % '' or a path to a cases CSV (see load_cases.m for the format)
% seizure_mode: 'legacy' -> src/detect_seizures.m (untouched, unmodified);
% 'robust' -> src/detect_seizures_robust.m (see README.md "dos ramas" --
% built for the case where the legacy energy-threshold-only method fails
% on low-quality signal: real seizures near threshold get cut by a
% duration filter applied before merging, while IID trains produce
% sustained plateaus that read as seizures; the robust branch merges
% first, then filters by duration, then filters by line-length ratio).
% Editable: setting seizure_mode='robust' on the 'normal' profile (or
% 'legacy' on any other) works without touching any code.
cfg.cases.profiles = struct( ...
    'normal',     struct('gain_mode', 'off',  'notch_mode', 'off', 'seizure_mode', 'legacy'), ...
    'attenuated', struct('gain_mode', 'auto', 'notch_mode', 'off', 'seizure_mode', 'robust'), ...
    'line',       struct('gain_mode', 'off',  'notch_mode', 'on',  'seizure_mode', 'robust'), ...
    'both',       struct('gain_mode', 'auto', 'notch_mode', 'on',  'seizure_mode', 'robust'));

%% ---- Preconditioning: gain + notch (src/precondition_lfp.m, --------
%%      utils/apply_notch_blocks.m), gated by the case system above -----
cfg.precondition.gain_deadband = 1.5;    % |gain_estimate| within [1/this, this] -> forced to exactly 1 (no-op)
cfg.precondition.gain_max = 50;          % clamp on any resolved gain (explicit or auto), with a warning
cfg.precondition.notch_freqs = 50;       % Hz; notch center frequency/ies
cfg.precondition.notch_harmonics = [];   % extra multiples of notch_freqs to also notch, if below Nyquist and the widest analysis band
cfg.precondition.notch_halfwidth_hz = 2; % each notch stopband is [f0-this, f0+this]
cfg.precondition.notch_order = 4;        % designfilt bandstopiir FilterOrder

%% ---- Robust seizure branch (src/detect_seizures_robust.m), used only ---
%%      when the resolved case's seizure_mode='robust' -- see README.md
%%      "dos ramas". Every value below comes from the measured analysis of
%%      animal 005 (20.25 h, 2 channels, 4 video-confirmed seizures) --
%%      see KNOWN_ISSUES.md for that validation's size and limits.
cfg.seizure_robust.ll_window_s = 2;       % line-length window; same time scale as the energy trace's movmean
cfg.seizure_robust.ll_threshold = 1.75;   % x ll_median_global; stable operating point measured between 1.7 and 1.8 on animal 005
cfg.seizure_robust.merge_gap_s = 2;       % merge energy-threshold crossings separated by this much or less, BEFORE the duration filter
cfg.seizure_robust.min_duration_s = 10;   % applied AFTER merging (applying it before, like the legacy branch does, is what lost the animal-005 GT1/GT4 seizures)
cfg.seizure_robust.max_duration_s = 60;   % bounds the merge so an IID train can't chain indefinitely; also reflects clinical experience that a seizure rarely runs longer. Events past this are KEPT and flagged over_max_duration, never silently dropped
cfg.seizure_robust.bilateral_tol_s = 5;   % overlap tolerance for the cross-channel bilateral-concordance label (confidence tag only, never a filter)
cfg.seizure_robust.hf_band = [80 250];    % Hz; informative hf_ratio_db band (movement-related high-frequency power) -- never a filter, see KNOWN_ISSUES.md non-convulsive-seizure caveat
cfg.seizure_robust.hf_reference_band = [5 40]; % Hz; hf_ratio_db's reference/denominator band
cfg.seizure_robust.dual_report = false;   % when true, run_pipeline_edf.m runs BOTH branches on every file and writes the non-primary one to seizures_events_alt.csv

%% ---- Detector evaluation against ground truth (src/evaluate_detections.m)
cfg.eval.match_tol_s = 10;                % seconds of allowed edge slack when matching a detection to a ground-truth event

%% ---- Stage 6: batch orchestration (src/run_pipeline_edf.m) ------------
cfg.paths.output_root = fullfile(pwd, 'pipeline_output');

end
