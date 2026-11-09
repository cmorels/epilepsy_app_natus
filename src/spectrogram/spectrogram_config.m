function cfg = spectrogram_config()
% SPECTROGRAM_CONFIG  Defaults of the spectrogram module (module 3,
% output stage 06_spectrograms). Independent of pipeline_config.m and
% campaign_config.m: nothing here is read by the detectors and nothing
% there is read by this module. Every field is under cfg.spec.
%
% Edit the copy returned in run_spectrograms.m, never this file, unless a
% default itself should change. See README.md, "Spectrogram module".

    s = struct();

    %% ---- paths ---------------------------------------------------------
    s.pipeline_root = '';      % folder holding one subfolder per animal (<root>\<subject>\...)
    s.output_root = '';        % '' -> pipeline_root. Writes <output_root>\<subject>\06_spectrograms\
                               % and <output_root>\06_spectrograms_summary\ (nothing else)
    s.seizures_file = '';      % optional seizures_events (.csv/.xlsx), used for exclusion zones and
                               % baseline checks of events that did not come from it

    %% ---- run control ---------------------------------------------------
    s.dry_run = true;          % true: resolve + write _selection_preview.csv only, compute nothing
    s.tz = 'Europe/Paris';     % only used when a txt header carries no timezone
    s.signal_source = 'clean'; % 'clean' (02_clean *_clean.txt) | 'raw' (01_txt *.txt x gain_applied,
                               % + the same 50 Hz notch if the clean header says notch_applied)
    s.max_nan_frac = 0.05;     % skip an event when more than this fraction of its window is NaN
    s.count_outliers = true;   % count samples clean_lfp interpolated in the window (needs the 01_txt)
    s.cache_disk = true;       % keep the binary signal cache in <out>\<subject>\06_spectrograms\_cache
                               % (false: cache in tempdir, deleted at the end of the run)
    s.cache_chunk = 5e6;       % samples per textscan chunk when building the cache
    s.trace_band = [1 100];    % display trace: zero-phase band-pass per gap-free block
    s.trace_filter_order = 4;

    %% ---- normalization / figure ---------------------------------------
    s.norm = 'db_baseline';    % 'db' | 'db_baseline' | 'zscore' | 'percent'
    s.clim = 'auto';           % 'auto' or [min max] (same units as norm), fixed across animals
    s.clim_prctile = [2 98];   % 'auto' = these percentiles of the displayed values
    s.cmap = 'parula';
    s.freq_scale = 'linear';   % 'linear' | 'log' (y axis of every TF panel)
    s.visible = false;         % figures hidden while running (saved .fig opens visible)
    s.save_fig = true;
    s.save_png = true;
    s.png_dpi = 300;

    %% ---- seizures ------------------------------------------------------
    z = struct();
    z.pad_pre_s = 30;          % context before onset (t = 0 at onset)
    z.pad_post_s = 30;         % context after offset
    z.baseline = [-30 -5];     % [a b] s relative to onset; or an absolute interval:
                               % {'dd/MM/yyyy HH:mm:ss.SSS', '...'} in baseline_time_ref, or
                               % [t1 t2] t_rel with baseline_is_t_rel = true
    z.baseline_is_t_rel = false;
    z.baseline_time_ref = 'natus';   % clock of a text baseline: 'natus' | 'pipeline'
    z.postictal_s = 30;        % PSD of [offset, offset + this]; 0 = no post-ictal PSD
    z.mt = struct('win_s', 2, 'step_s', 0.1, 'TW', 3, 'K', 5, 'fmin', 0.5, 'fmax', 100);
    z.morlet = struct('fmin', 1, 'fmax', 100, 'nf', 80, 'cycles', [3 10], ...
        'cycles_scale', 'linear', 't_step_s', 0.01);
    z.bands = struct('name', {'delta', 'theta', 'beta', 'gamma'}, ...
        'range', {[1 4], [4 12], [12 30], [30 100]});
    z.feature_fmin = 1;        % peak / dominant frequency searched in [feature_fmin, mt.fmax] (notch band excluded)
    z.domfreq_method = 'aperiodic';   % dominant frequency of each multitaper window:
                                      % 'aperiodic'      peak above the 1/f background (log-log line fit)
                                      % 'baseline_ratio' max of the power in dB re baseline (broadband
                                      %                  seizures make it jump across the band)
                                      % 'absolute_band'  max of the absolute power within domfreq_band
    z.domfreq_band = [2 30];
    s.seizure = z;

    %% ---- IIDs / spikes -------------------------------------------------
    d = struct();
    d.spike_source = 'pipeline';   % 'pipeline' (04_iid *_IID_results.mat) | 'detect_local' | 'manual'
    d.polarity = 'auto';           % 'auto' | 'neg' | 'pos'
    d.align_ms = 20;               % re-align each peak to the extremum within +-align_ms
    d.align_on = 'analysis';       % 'analysis' (detrended signal used for TF) | 'display' (1-100 Hz trace)
    d.epoch_s = [-0.5 0.5];
    d.baseline = [-0.5 -0.3];      % s relative to each peak, or an absolute clean interval (as seizures)
    d.baseline_is_t_rel = false;
    d.baseline_time_ref = 'natus';
    d.exclusion_buffer_s = 5;      % around seizures, when zones are built here (not 'pipeline')
    d.p2p_window_ms = 100;         % peak-to-peak and FWHM searched within +-this
    d.maxfreq_window_ms = 50;      % frequency of max induced power searched within +-this
    d.itc = true;
    d.max_interval_s = 600;        % up to this: 'interval' mode (continuous scalogram of the interval);
                                   % longer (e.g. hours of a spiking period): 'zone' mode -- spike-aligned
                                   % maps over every spike + spike-rate overview, no continuous scalogram
    d.zone_chunk_s = 60;           % zone mode: CWT computed in blocks of this length (RAM)
    d.rate_bin_s = 60;             % zone mode: spike-rate bins
    d.overview_points = 20000;     % zone mode: min/max points of the trace overview
    d.max_traces_saved = 300;      % individual spike traces kept in the .mat / figure (mean / SEM use all)
    d.morlet = struct('fmin', 2, 'fmax', 150, 'nf', 80, 'cycles', [3 7], ...
        'cycles_scale', 'linear', 't_step_s', 0.005);
    % 'detect_local': same steps as detect_iid.m (|band-pass 15-70 Hz| + findpeaks),
    % restricted to the interval. Values copied from pipeline_config.m cfg.iid.* (2026-11).
    d.detect = struct( ...
        'threshold_mode', 'pipeline_grid', ... % 'pipeline_grid' (detect_iid: P97 of |signal| on the
        ...                                    %  100:5:130 uV grid) | 'p97' (P97 of the interval,
        ...                                    %  no grid) | 'fixed' (threshold_uV)
        'threshold_factor', 2.5, ...
        'threshold_uV', NaN, ...
        'upper_uV', 2000, ...
        'min_distance_s', 0.030, ...
        'max_width_s', 0.060, ...
        'min_prominence_mV', 0.2, ...
        'bandpass_band', [15 70], ...
        'bandpass_steepness', 0.5, ...
        'grid_start_uV', 100, 'grid_step_uV', 5, 'grid_max_uV', 130, 'grid_fraction', 0.97);
    s.iid = d;

    cfg = struct('spec', s);
end
