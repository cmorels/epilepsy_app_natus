function K = rp128_checks()
% RP128_CHECKS  Phase 4: controls after the pipeline run of the A1 / A2
% EDFs, plus the review figures (<out>/figures_review):
%   1. reference / configuration actually used (config_used.mat, qc_report)
%   2. seizure thresholds per hemisphere vs 128-s 20260730 (C2)
%   3. seizure detections: number, duration, hemisphere, none > 60 s, none
%      in 17:18:39-17:23:55, rate per hour vs the other EDFs of 128-s
%   4. IIDs (provisional): outputs of both hemispheres, baseline, threshold, complexes/min
%   5. flat segments (|x| < 5 uV for > 60 s) in the two EDFs
%   6. figures: panorama of both hemispheres with the A1 -> A2 junction,
%      raster of this recording (v2 format: rows Day k, no dates)

    o = rp128_paths();
    addpath(fullfile(o.repo, 'src'), fullfile(o.repo, 'src', 'utils'), fullfile(o.repo, 'src', 'stats'));
    fig_dir = fullfile(o.out, 'figures_review');
    if ~isfolder(fig_dir), mkdir(fig_dir); end
    K = struct();
    fs = 512;

    % ---- 1. reference and configuration -------------------------------------
    used = load(fullfile(o.out, 'config_used.mat'), 'cfg');
    c2 = load(o.c2_cfg_mat, 'cfg');
    R = readtable(o.c2_ref_csv, 'TextType', 'char');
    ref = used.cfg.quality.reference_sigma_uV;
    K.ref = table(R.region, R.reference_uV, [cellfun(@(r) ref(r), R.region(1:end-1)); used.cfg.quality.reference_sigma_fallback_uV], ...
        'VariableNames', {'region', 'C2_reference_uV', 'used_uV'});
    K.ref.identical = K.ref.C2_reference_uV == K.ref.used_uV;
    K.cfg_diff = rp_cfg_diff(used.cfg, c2.cfg, {'paths.output_root', 'edf.channels.log_file'});
    Q = readtable(fullfile(o.out, '05_summaries', 'qc_report.csv'), 'TextType', 'char', 'VariableNamingRule', 'preserve');
    K.qc = Q(:, {'source_file', 'region', 'case_applied', 'gain_applied', 'notch_applied', 'reference_used', 'sigma_band_uV', ...
        'line_ratio_db', 'quality_class', 'n_errors'});
    fprintf('1. Reference used = C2''s: %d (%s); config differs from the C2 in: %s\n', all(K.ref.identical), ...
        strjoin(compose('%s %.10g', string(K.ref.region), K.ref.used_uV), ', '), strjoin([K.cfg_diff, {'nothing'}], ', '));
    disp(K.ref); disp(K.qc);

    % ---- 2. seizure thresholds vs 20260730 -----------------------------------
    % (seizures_summary.csv: *_seizures.mat is saved only for channels with >= 1 detection)
    rd = @(f) readtable(f, 'TextType', 'char', 'VariableNamingRule', 'preserve');
    S30 = rd(fullfile(o.c2_subject, '05_summaries', 'seizures_summary.csv'));
    S30 = S30(contains(S30.source_file, o.next_stem), :);
    SN = rd(fullfile(o.out, '05_summaries', 'seizures_summary.csv'));
    rows = {};
    for r = o.regions
        b = S30(strcmp(S30.region, r{1}), :);
        for i = 1:2
            a = SN(strcmp(SN.region, r{1}) & contains(SN.source_file, o.seg(i).stem), :);
            rows(end+1, :) = {r{1}, o.seg(i).name, a.threshold_value, b.threshold_value, a.threshold_value / b.threshold_value, ...
                a.n_segments, a.n_segments / (a.valid_duration_min / 60), b.n_segments / (b.valid_duration_min / 60), ...
                a.n_seizures}; %#ok<AGROW>
        end
    end
    K.thr = cell2table(rows, 'VariableNames', {'region', 'edf', 'energy_threshold', 'threshold_0730', 'ratio_vs_0730', ...
        'n_candidates', 'candidates_per_h', 'candidates_per_h_0730', 'n_seizures'});
    fprintf(['2. Seizure thresholds vs 20260730 (C2). Energy of the min-max normalized signal (each file scaled by its own ' ...
        'range), threshold = 10 x its median: the absolute value is per file, not comparable across files of different length/artifacts\n']);
    disp(K.thr);

    % ---- 3. detections -------------------------------------------------------------
    D = read_events(fullfile(o.out, '05_summaries', 'seizures_events_pipeline.csv'), fullfile(o.out, '05_summaries', 'seizures_events.csv'));
    if height(D) > 0
        D.start = datetime(D.start_abs, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS', 'TimeZone', o.tz, 'Locale', 'en_US');
    else
        D.start = NaT(0, 1, 'TimeZone', o.tz);
    end
    D.duration_s = D.end_s - D.start_s;
    K.events = D(:, {'source_file', 'seizure_id', 'region', 'detection_status', 'start', 'duration_s', 'll_ratio'});
    K.events.start.Format = 'yyyy-MM-dd HH:mm:ss';
    K.n_rows = height(D);
    K.n_events = height(unique(D(:, {'source_file', 'seizure_id'})));
    K.max_duration_s = max([D.duration_s; 0]);
    K.in_dx21 = sum(D.start >= o.cut1 & D.start < o.cut2);
    fprintf('3. Detections: %d event(s), %d row(s); max duration %.1f s (> 60 s: %d); onsets in DX-21: %d\n', K.n_events, ...
        K.n_rows, K.max_duration_s, sum(D.duration_s > 60), K.in_dx21);
    disp(K.events);
    K.rates = rates(o, D);
    disp(K.rates);

    % ---- 4. IIDs -------------------------------------------------------------------
    I = readtable(fullfile(o.out, '05_summaries', 'iid_summary.csv'), 'TextType', 'char', 'VariableNamingRule', 'preserve');
    I.n_files_04_iid = zeros(height(I), 1);
    for i = 1:height(I)
        stem = regexprep(I.source_file{i}, '\.EDF$', '', 'ignorecase');
        I.n_files_04_iid(i) = numel(dir(fullfile(o.out, '04_iid', stem, sprintf('*_%s_*', I.region{i}))));
    end
    K.iid = I(:, {'source_file', 'region', 'analyzed_duration_min', 'baseline_uV', 'lower_threshold_uV', 'total_complexes', ...
        'complexes_per_min', 'n_files_04_iid'});
    I30 = readtable(fullfile(o.c2_subject, '05_summaries', 'iid_summary.csv'), 'TextType', 'char', 'VariableNamingRule', 'preserve');
    K.iid_0730 = I30(contains(I30.source_file, o.next_stem), {'region', 'baseline_uV', 'lower_threshold_uV', 'complexes_per_min'});
    fprintf('4. IIDs (provisional, pipeline zones):\n'); disp(K.iid);
    fprintf('   20260730 (C2), for reference:\n'); disp(K.iid_0730);
    % continuity: complexes per hour, last 6 h of A2 vs first 6 h of 20260730 (same port, consecutive)
    J = rd(fullfile(o.out, '05_summaries', 'iid_events.csv'));
    J30 = rd(fullfile(o.c2_subject, '05_summaries', 'iid_events.csv'));
    J30 = J30(contains(J30.source_file, o.next_stem), :);
    t = @(T) datetime(T.start_abs, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS', 'TimeZone', o.tz, 'Locale', 'en_US');
    a2_end = datetime(2026, 7, 30, 13, 2, 11, 'TimeZone', o.tz);
    n30 = datetime(2026, 7, 30, 13, 3, 36, 'TimeZone', o.tz);
    rows = {};
    for r = o.regions
        ta = t(J(strcmp(J.region, r{1}) & contains(J.source_file, o.seg(2).stem), :));
        tb = t(J30(strcmp(J30.region, r{1}), :));
        for hh = -6:5
            if hh < 0, w = a2_end + hours([hh hh + 1]); x = ta; src = 'A2 (29-30/07)';
            else, w = n30 + hours([hh hh + 1]); x = tb; src = '20260730 (C2)'; end
            rows(end+1, :) = {r{1}, src, hh, sum(x >= w(1) & x < w(2)) / 60}; %#ok<AGROW>
        end
    end
    K.iid_continuity = cell2table(rows, 'VariableNames', {'region', 'edf', 'hour_from_junction', 'complexes_per_min'});
    fprintf('   continuity A2 -> 20260730 (complexes/min per hour, hour 0 = 13:03 on 30/07):\n');
    disp(unstack(K.iid_continuity(:, {'region', 'hour_from_junction', 'complexes_per_min'}), 'complexes_per_min', 'region'));

    % ---- 5. flat segments + 6. figures ----------------------------------------------
    t0 = []; seg = struct();
    flat = cell(0, 4);
    for i = 1:2
        E = rp_edf_read(o.seg(i).edf, o.seg(i).labels);
        st = rp_start(E, o.tz);
        if i == 1, t0 = st; end
        seg(i).t_h = hours(st - t0);
        for c = 1:2
            x = rp_edf_phys(E.x(:, c), E.sig(c));
            runs = rp_runs(abs(x) < 5);
            runs = runs((runs(:, 2) - runs(:, 1) + 1) / fs > 60, :);
            for k = 1:size(runs, 1)
                flat(end+1, :) = {o.seg(i).name, o.regions{c}, st + seconds((runs(k, 1) - 1) / fs), (runs(k, 2) - runs(k, 1) + 1) / fs}; %#ok<AGROW>
            end
            [seg(i).lo(:, c), seg(i).hi(:, c)] = envelope_(x, 2 * fs);
            k0 = rp_first_sample(st, o.cut1 - minutes(20), fs);   % zoom window around the junction
            k1 = rp_first_sample(st, o.cut2 + minutes(20), fs);
            k0 = max(k0, 1); k1 = min(k1, numel(x));
            if k1 > k0
                [seg(i).zlo{c}, seg(i).zhi{c}] = envelope_(x(k0:k1), fs / 4);
                seg(i).z_t0 = st + seconds((k0 - 1) / fs);
            end
        end
        seg(i).n = size(E.x, 1);
    end
    K.flat = cell2table(flat, 'VariableNames', {'edf', 'region', 'start', 'duration_s'});
    fprintf('5. Flat segments (|x| < 5 uV for > 60 s): %d\n', height(K.flat));
    disp(K.flat);
    K.fig_panorama = panorama(o, seg, D, t0, fs, fig_dir);
    K.fig_raster = raster(o, D, t0, fig_dir);
    K.fig_detections = {};
    for i = 1:2
        d = dir(fullfile(o.out, '03_seizures', o.seg(i).stem, '*seizure*.png'));
        K.fig_detections = [K.fig_detections; fullfile({d.folder}, {d.name})'];
    end
    fprintf('6. Figures: %s | %s | %d per-detection PNG(s) in 03_seizures\n', K.fig_panorama, K.fig_raster, numel(K.fig_detections));
    save(fullfile(o.out, 'phase4_checks.mat'), 'K');
end

function D = read_events(pristine, labeled)
    f = pristine;
    if ~isfile(f), f = labeled; end
    opts = detectImportOptions(f, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
    opts = setvartype(opts, opts.VariableNames, 'char');
    D = readtable(f, opts);
    for v = {'seizure_id', 'start_s', 'end_s', 'll_ratio'}
        D.(v{1}) = str2double(D.(v{1}));
    end
end

function T = rates(o, D)
% detection events per hour of recording, every EDF of 128-s (C1, C2, this run)
    rows = {};
    roots = {fullfile(o.batch, 'pipeline_output_campaign_20261006', o.subject), o.c2_subject};
    camp = {'C1', 'C2'};
    for c = 1:2
        E = read_events('', fullfile(roots{c}, '05_summaries', 'seizures_events.csv'));
        hasl = ismember('user_check', E.Properties.VariableNames);
        ff = dir(fullfile(roots{c}, '02_clean', '*_HPCr_clean.txt'));
        for k = 1:numel(ff)
            h = rp_txt_header(fullfile(ff(k).folder, ff(k).name));
            x = E(strcmp(E.source_file, h.source_file), :);
            ntp = NaN;
            if hasl, ntp = height(unique(x(strcmp(strtrim(x.user_check), 'TP'), {'seizure_id'}))); end
            rows(end+1, :) = {camp{c}, regexprep(h.source_file, '\.EDF$', ''), str2double(h.total_duration_s) / 3600, ...
                height(unique(x(:, {'seizure_id'}))), ntp}; %#ok<AGROW>
        end
    end
    for i = 1:2
        h = rp_txt_header(fullfile(o.out, '02_clean', sprintf('%s_HPCr_clean.txt', base_of(o, i))));
        x = D(strcmp(D.source_file, [o.seg(i).stem '.EDF']), :);
        rows(end+1, :) = {'reprocess', o.seg(i).stem, str2double(h.total_duration_s) / 3600, height(unique(x(:, {'seizure_id'}))), NaN}; %#ok<AGROW>
    end
    T = cell2table(rows, 'VariableNames', {'source', 'edf', 'hours', 'n_detection_events', 'n_TP_events'});
    T.detections_per_h = T.n_detection_events ./ T.hours;
    T.TP_per_h = T.n_TP_events ./ T.hours;
end

function b = base_of(o, i)
    d = dir(fullfile(o.out, '02_clean', '*_HPCr_clean.txt'));
    for k = 1:numel(d)
        h = rp_txt_header(fullfile(d(k).folder, d(k).name));
        if strcmp(h.source_file, [o.seg(i).stem '.EDF']), b = erase(d(k).name, '_HPCr_clean.txt'); return; end
    end
    error('rp128_checks:Clean', 'no clean file for %s', o.seg(i).stem);
end

function r = rp_runs(m)
    d = diff([0; m(:); 0]);
    r = [find(d == 1), find(d == -1) - 1];
end

function [lo, hi] = envelope_(x, w)
    n = floor(numel(x) / w) * w;
    X = reshape(x(1:n), w, []);
    lo = min(X, [], 1)'; hi = max(X, [], 1)';
end

function f = panorama(o, seg, D, t0, fs, fig_dir)
    fig = figure('Visible', 'off', 'Position', [50 50 1500 900], 'Color', 'w');
    col = {[0.20 0.45 0.70], [0.85 0.40 0.10]};
    dx = hours([o.cut1 o.cut2] - t0);
    allv = abs([seg(1).lo(:); seg(1).hi(:); seg(2).lo(:); seg(2).hi(:)]);
    yl = prctile(allv, 99.9) * [-1.1 1.1];
    for c = 1:2
        ax = subplot(4, 1, c, 'Parent', fig); hold(ax, 'on');
        patch(ax, dx([1 2 2 1]), yl([1 1 2 2]), [0.88 0.25 0.25], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        for i = 1:2
            t = seg(i).t_h + ((1:numel(seg(i).lo(:, c)))' - 0.5) * 2 / 3600;
            plot(ax, [t; flipud(t)], [seg(i).lo(:, c); flipud(seg(i).hi(:, c))], '-', 'Color', col{i});
        end
        xline(ax, dx(2), 'k--', 'A1 \rightarrow A2', 'LabelVerticalAlignment', 'bottom');
        x = D(strcmp(D.region, o.regions{c}), :);
        plot(ax, hours(x.start - t0), repmat(yl(2) * 0.9, height(x), 1), 'v', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        ylim(ax, yl); xlim(ax, [0, seg(2).t_h + seg(2).n / fs / 3600]);
        ylabel(ax, sprintf('%s (\\muV)', o.regions{c}));
        title(ax, sprintf('%s: A1 (%s) until 17:18:39, A2 (%s) from 17:23:55; red = DX-21; \\nabla = detection', ...
            o.regions{c}, o.seg(1).labels{c}, o.seg(2).labels{c}));
        xlabel(ax, 'hours since the first sample (min/max per 2 s)');
    end
    for c = 1:2
        ax = subplot(4, 1, 2 + c, 'Parent', fig); hold(ax, 'on');
        z = hours([o.cut1 o.cut2] - t0) * 60;
        patch(ax, z([1 2 2 1]), yl([1 1 2 2]), [0.88 0.25 0.25], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
        for i = 1:2
            if ~isfield(seg(i), 'zlo') || numel(seg(i).zlo) < c || isempty(seg(i).zlo{c}), continue; end
            t = minutes(seg(i).z_t0 - t0) + ((1:numel(seg(i).zlo{c}))' - 0.5) * 0.25 / 60;
            plot(ax, [t; flipud(t)], [seg(i).zlo{c}; flipud(seg(i).zhi{c})], '-', 'Color', col{i});
        end
        xlim(ax, z + [-20 20]); ylim(ax, yl);
        ylabel(ax, sprintf('%s (\\muV)', o.regions{c}));
        xlabel(ax, 'minutes since the first sample (zoom +/- 20 min around DX-21, min/max per 0.25 s)');
        title(ax, sprintf('%s: A1 -> A2 junction', o.regions{c}));
    end
    f = fullfile(fig_dir, 'panorama_A1_A2.png');
    exportgraphics(fig, f, 'Resolution', 200);
    close(fig);
end

function f = raster(o, D, t0, fig_dir)
% same drawing as the English v2 rasters (raster_en.m): rows Day k from the
% first sample, hours since start of day, dark phase 19:30-08:00
    t0p = posixtime(t0);
    obs = zeros(0, 2);
    for i = 1:2
        E = rp_edf_read(o.seg(i).edf);
        st = posixtime(rp_start(E, o.tz));
        obs(end+1, :) = [st, st + E.hdr.nrec * E.hdr.rec_dur_s]; %#ok<AGROW>
    end
    dx = posixtime([o.cut1 o.cut2]);
    nr = ceil((obs(end, 2) - t0p) / 86400);
    fig = figure('Visible', 'off', 'Position', [100 100 1150 max(320, 90 + 26 * nr)], 'Color', 'w');
    ax = axes(fig, 'Position', [0.08 0.10 0.68 0.80]); hold(ax, 'on');
    c = struct('eeg', [0.45 0.75 0.45], 'dx', [0.88 0.25 0.25], 'dark', [0.86 0.86 0.90]);
    for d = 1:nr
        w = t0p + [(d - 1) d] * 86400;
        band(ax, dark_iv(w, o.tz), w, d, 0.5, c.dark);
        band(ax, obs, w, d, 0.32, c.eeg);
        band(ax, dx, w, d, 0.32, c.dx);
    end
    ts = posixtime(D.start);
    k = floor((ts - t0p) / 86400) + 1;
    plot(ax, (ts - t0p - (k - 1) * 86400) / 3600, k, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'MarkerSize', 6, ...
        'HandleVisibility', 'off');
    L = [patch(ax, NaN, NaN, c.eeg, 'EdgeColor', 'none', 'DisplayName', 'EEG observed'), ...
        patch(ax, NaN, NaN, c.dx, 'EdgeColor', 'none', 'DisplayName', 'disconnection (DX-21)'), ...
        patch(ax, NaN, NaN, c.dark, 'EdgeColor', 'none', 'DisplayName', 'dark phase (19:30-08:00)'), ...
        plot(ax, NaN, NaN, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'DisplayName', 'detection (not labeled yet)')];
    legend(ax, L, 'Location', 'eastoutside', 'FontSize', 8);
    set(ax, 'YDir', 'reverse', 'YTick', 1:nr, 'YTickLabel', arrayfun(@(d) sprintf('Day %d', d), 1:nr, 'UniformOutput', false), ...
        'XLim', [0 24], 'XTick', 0:2:24, 'YLim', [0.4 nr + 0.6], 'FontSize', 9, 'Layer', 'top');
    xlabel(ax, 'Hours since start of day');
    title(ax, '128-s (epileptic): 20260729 reprocessed, A1 + A2');
    f = fullfile(fig_dir, 'raster_128s_20260729_A1A2.png');
    exportgraphics(fig, f, 'Resolution', 300);
    close(fig);
end

function band(ax, X, w, row, hh, face)
    x = cs_iv('intersect', X, w);
    for k = 1:size(x, 1)
        a = (x(k, 1) - w(1)) / 3600; b = (x(k, 2) - w(1)) / 3600;
        patch(ax, [a b b a], row + [-hh -hh hh hh], face, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

function iv = dark_iv(w, tz)
    d0 = dateshift(datetime(w(1), 'ConvertFrom', 'posixtime', 'TimeZone', tz), 'start', 'day');
    iv = zeros(0, 2);
    for d = [d0 - caldays(1), d0, d0 + caldays(1)]
        [y, mo, dd] = ymd(d);
        iv(end+1, :) = posixtime(datetime(y, mo, dd, [0 8], 0, 0, 'TimeZone', tz)); %#ok<AGROW>
        iv(end+1, :) = [posixtime(datetime(y, mo, dd, 19, 30, 0, 'TimeZone', tz)), posixtime(datetime(y, mo, dd + 1, 0, 0, 0, 'TimeZone', tz))]; %#ok<AGROW>
    end
    iv = cs_iv('intersect', iv, w);
end
