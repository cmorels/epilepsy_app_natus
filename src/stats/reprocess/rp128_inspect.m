function P = rp128_inspect()
% RP128_INSPECT  Phase 1 (read only): header of the new EDF, case A/B,
% times (start = original's, no overlap with 20260730), fs, hemisphere
% mapping, A1 vs A2 compatibility per hemisphere (60 min before 17:18:39
% vs 60 min after 17:23:55), the C2 configuration of 128-s and whether the
% current campaign_config reproduces it. P.stop lists every STOP reason.

    o = rp128_paths();
    addpath(fullfile(o.repo, 'src'), fullfile(o.repo, 'src', 'utils'));
    P = struct('stop', {{}});

    % ---- 1. header --------------------------------------------------------
    lab = {'EEG A1C2', 'EEG A1C4', 'EEG A2C2', 'EEG A2C4'};
    E = rp_edf_read(o.src_edf, lab);
    fs = E.sig(1).fs;
    start = rp_start(E, o.tz);
    n = size(E.x, 1);
    P.header = E.hdr;
    P.labels = {E.sig.label};
    P.fs = unique([E.sig(1:4).fs]);
    P.start = start;
    P.end = start + seconds(n / fs);
    P.n_samples = n;
    P.units = unique({E.sig(1:4).dim});
    fprintf('EDF %s\n  %s, ns=%d, nrec=%d x %g s, %d samples/ch, fs=%g Hz, units=%s\n  start %s, end %s (%.4f h)\n', ...
        o.src_edf, E.hdr.reserved, E.hdr.ns, E.hdr.nrec, E.hdr.rec_dur_s, n, fs, strjoin(P.units, ','), ...
        char(P.start, 'yyyy-MM-dd HH:mm:ss.SSS'), char(P.end, 'yyyy-MM-dd HH:mm:ss.SSS'), n / fs / 3600);
    fprintf('  signals: %s\n', strjoin(P.labels, ' | '));
    if ~isequal(P.fs, 512), P.stop{end+1} = sprintf('fs = %s, not 512 Hz', mat2str(P.fs)); end

    % ---- 2. case A / B ----------------------------------------------------
    k1 = rp_first_sample(start, o.cut1, fs);   % first sample at/after 17:18:39 (1-based)
    k2 = rp_first_sample(start, o.cut2, fs);   % first sample at/after 17:23:55
    seg = {1:k1 - 1, k2:n};
    P.activity = table(lab', zeros(4, 1), zeros(4, 1), 'VariableNames', {'label', 'mad_before_cut_uV', 'mad_after_cut_uV'});
    for c = 1:4
        y = rp_edf_phys(E.x(:, c), E.sig(c));
        P.activity.mad_before_cut_uV(c) = mad(y(seg{1}), 1);
        P.activity.mad_after_cut_uV(c) = mad(y(seg{2}), 1);
    end
    P.case = 'B';
    P.case_why = sprintf(['4 EEG signals (A1C2, A1C4, A2C2, A2C4) over the whole file (%d samples each), not one ' ...
        'combined pair; MAD (uV) before 17:18:39 / after 17:23:55: %s'], n, strjoin(compose('%s %.1f / %.1f', ...
        string(lab'), P.activity.mad_before_cut_uV, P.activity.mad_after_cut_uV), '; '));
    fprintf('Case %s: %s\n', P.case, P.case_why);
    disp(P.activity);

    % ---- 3. times -----------------------------------------------------------
    c2h = rp_txt_header(o.c2_txt.HPCr);
    P.c2_start = datetime(c2h.session_start_datetime, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSSXXX', 'TimeZone', o.tz);
    P.c2_n_samples = str2double(c2h.n_samples);
    P.c2_end = P.c2_start + seconds(P.c2_n_samples / fs);
    En = rp_edf_read(o.next_edf);
    P.next_start = rp_start(En, o.tz);
    P.start_matches = abs(milliseconds(P.start - P.c2_start)) < 0.5;
    P.overlap_s = max(0, seconds(P.end - P.next_start));
    fprintf('C2 original: start %s, %d samples (%.0f s), end %s\n', char(P.c2_start, 'yyyy-MM-dd HH:mm:ss.SSS'), ...
        P.c2_n_samples, P.c2_n_samples / fs, char(P.c2_end, 'yyyy-MM-dd HH:mm:ss.SSS'));
    fprintf('Next EDF (%s): start %s; overlap with the new EDF: %.3f s\n', o.next_stem, ...
        char(P.next_start, 'yyyy-MM-dd HH:mm:ss.SSS'), P.overlap_s);
    if ~P.start_matches, P.stop{end+1} = 'start of the new EDF differs from the original'; end
    if P.overlap_s > 0, P.stop{end+1} = sprintf('overlaps 20260730 by %.3f s', P.overlap_s); end
    % tail beyond the original's last sample (record padding of the 12 s export)
    tail = P.c2_n_samples + 1:n;
    P.tail = struct('n', numel(tail), 'distinct_digital', arrayfun(@(c) numel(unique(E.x(tail, c))), 1:4));
    fprintf('Samples beyond the original end: %d (%.3f s); distinct digital values per channel: %s\n', numel(tail), ...
        numel(tail) / fs, mat2str(P.tail.distinct_digital));

    % ---- 3b. A1 of the new EDF == the C2 import of the original EDF -------
    for r = 1:2
        y0 = rp_read_txt(o.c2_txt.(o.regions{r}), k1 - 1);
        y1 = rp_edf_phys(E.x(1:k1 - 1, r), E.sig(r));
        P.a1_vs_c2.(o.regions{r}) = max(abs(y1 - y0));
        fprintf('A1 %s (new EDF) vs C2 01_txt, samples 1..%d: max |diff| = %.2e uV\n', o.regions{r}, k1 - 1, P.a1_vs_c2.(o.regions{r}));
    end

    % ---- 4. mapping and C2 case from the master log ------------------------
    T = load_recording_log(o.master_log);
    row = T(strcmp(T.filename, o.stem) & strcmp(T.animal_id, o.subject), :);
    P.log_row = row;
    c2 = load(o.c2_cfg_mat, 'cfg');
    cfg_c2 = c2.cfg;
    [case_spec, ~] = resolve_case_excel(row.attenuation{1}, cfg_c2);
    P.c2_case = case_spec.case_applied;
    fprintf('Master log: %s port %s, HPCr=%s, HPCl=%s, Attenuation="%s" -> C2 case %s\n', o.subject, row.port{1}, ...
        row.hpcr{1}, row.hpcl{1}, row.attenuation{1}, P.c2_case);
    ch = struct('HPCr', row.hpcr{1}, 'HPCl', row.hpcl{1});

    % ---- 5. reference and configuration of the C2 --------------------------
    R = readtable(o.c2_ref_csv, 'TextType', 'char');
    P.ref = R;
    disp(R);
    cfg = rp128_cfg(o, o.master_log);
    P.cfg_diff = rp_cfg_diff(cfg, cfg_c2);
    fprintf('campaign_config (HEAD) + C2 reference vs C2 config_used.mat of 128-s: %s\n', ...
        ifelse(isempty(P.cfg_diff), 'identical', ['differs in ' strjoin(P.cfg_diff, ', ')]));
    if ~isempty(P.cfg_diff), P.stop{end+1} = ['configuration differs from the C2: ' strjoin(P.cfg_diff, ', ')]; end

    % ---- 6. A1 vs A2 compatibility ------------------------------------------
    w = o.compat_window_s * fs;
    rows = {};
    for r = 1:2
        reg = o.regions{r};
        ia = find(strcmp(lab, ['EEG A1' ch.(reg)]));
        ib = find(strcmp(lab, ['EEG A2' ch.(reg)]));
        a = rp_edf_phys(E.x(k1 - w:k1 - 1, ia), E.sig(ia));
        b = rp_edf_phys(E.x(k2:k2 + w - 1, ib), E.sig(ib));
        ma = rp_metrics(a, fs, E.sig(ia), reg, cfg);
        mb = rp_metrics(b, fs, E.sig(ib), reg, cfg);
        rows(end+1, :) = {reg, lab{ia}, ma.med_abs, ma.mad, ma.p97, ma.rel50, ma.line_ratio_db, ma.sigma_band, ma.suggested}; %#ok<AGROW>
        rows(end+1, :) = {reg, lab{ib}, mb.med_abs, mb.mad, mb.p97, mb.rel50, mb.line_ratio_db, mb.sigma_band, mb.suggested}; %#ok<AGROW>
        rows(end+1, :) = {reg, 'A2/A1', mb.med_abs / ma.med_abs, mb.mad / ma.mad, mb.p97 / ma.p97, mb.rel50 / ma.rel50, ...
            mb.line_ratio_db - ma.line_ratio_db, mb.sigma_band / ma.sigma_band, ''}; %#ok<AGROW>
        for m = {'mad', 'p97'}
            rr = mb.(m{1}) / ma.(m{1});
            if rr < o.compat_ratio(1) || rr > o.compat_ratio(2)
                P.stop{end+1} = sprintf('%s: A2/A1 %s = %.3f outside [%.2f, %.2f]', reg, upper(m{1}), rr, o.compat_ratio); %#ok<AGROW>
            end
        end
        if ~strcmp(ma.suggested, mb.suggested)
            P.stop{end+1} = sprintf('%s: suggested case A1 = %s, A2 = %s', reg, ma.suggested, mb.suggested); %#ok<AGROW>
        end
    end
    P.compat = cell2table(rows, 'VariableNames', {'region', 'channel', 'median_abs_uV', 'MAD_uV', 'P97_abs_uV', ...
        'rel_power_49_51Hz', 'line_ratio_dB', 'sigma_band_uV', 'suggested_case'});
    fprintf('A1: %s -> %s (60 min); A2: %s -> %s (60 min)\n', char(start + seconds((k1 - w - 1) / fs), 'HH:mm:ss.SSS'), ...
        char(start + seconds((k1 - 1) / fs), 'HH:mm:ss.SSS'), char(start + seconds((k2 - 1) / fs), 'HH:mm:ss.SSS'), ...
        char(start + seconds((k2 + w - 1) / fs), 'HH:mm:ss.SSS'));
    disp(P.compat);

    % ---- 7. annotations other than time keeping ----------------------------
    P.annotations = rp_annotations(E);
    fprintf('Annotations (besides time keeping): %d\n', height(P.annotations));
    disp(P.annotations);

    P.k1 = k1; P.k2 = k2;
    if isempty(P.stop)
        fprintf('\nPHASE 1: all checks passed\n');
    else
        fprintf(2, '\nPHASE 1 STOP:\n  %s\n', strjoin(P.stop, '\n  '));
    end
end

function m = rp_metrics(x, fs, s, region, cfg)
    m.med_abs = median(abs(x));
    m.mad = mad(x, 1);
    m.p97 = prctile(abs(x), 97);
    [pxx, f] = pwelch(x, hamming(4 * fs), 2 * fs, [], fs);
    m.rel50 = sum(pxx(f >= 49 & f <= 51)) / sum(pxx(f >= 1 & f <= 250));
    step = abs(s.pmax - s.pmin) / (s.dmax - s.dmin);
    data = struct('signal', x, 'fs', fs, 't_rel', (0:numel(x) - 1)' / fs, 'valid_mask', true(numel(x), 1), ...
        'meta', struct('quantization_step_uV', num2str(step, 17), 'physical_min_uV', num2str(s.pmin), ...
        'physical_max_uV', num2str(s.pmax)));
    q = signal_quality(data, cfg);
    m.line_ratio_db = q.line_ratio_db;
    m.sigma_band = q.sigma_band_uV;
    % make_cases_template.m's rule, with the C2 reference of this region
    ref = cfg.quality.reference_sigma_uV(region);
    att = q.sigma_band_uV < ref / cfg.quality.attenuation_ratio_thr || q.snr_quantization_db < cfg.quality.low_amp_snr_db;
    line = q.line_ratio_db > cfg.quality.line_ratio_thr_db || q.pct_time_line_high > cfg.quality.line_intermittent_pct;
    if strcmp(q.quality_class, 'unusable'), m.suggested = 'unusable';
    elseif att && line, m.suggested = 'both';
    elseif att, m.suggested = 'attenuated';
    elseif line, m.suggested = 'line';
    else, m.suggested = 'normal';
    end
end

function A = rp_annotations(E)
    rows = {};
    for r = 1:size(E.ann, 2)
        tals = strsplit(char(E.ann(:, r)'), char(0));
        for t = 1:numel(tals)
            parts = strsplit(tals{t}, char(20));
            if numel(parts) < 2, continue; end
            txt = strjoin(parts(2:end), ' | ');
            if isempty(strtrim(strrep(txt, '|', ''))), continue; end
            rows(end+1, :) = {r, parts{1}, txt}; %#ok<AGROW>
        end
    end
    A = cell2table(rows, 'VariableNames', {'record', 'onset', 'text'});
end

function y = rp_read_txt(p, n)
    fid = fopen(p, 'r'); cl = onCleanup(@() fclose(fid));
    c = textscan(fid, '%f', n, 'CommentStyle', '#');
    y = c{1};
end

function s = ifelse(c, a, b)
    if c, s = a; else, s = b; end
end
