function [D, ND] = cs_threshold_diag(cfg, S, rv, edf_tbl)
% CS_THRESHOLD_DIAG  Section 6.3 -- INFORMATIVE ONLY, seizure detection is
% never run. The robust detector's thresholds are per channel x EDF file:
%   energy threshold = median(energy over the file's trimmed valid blocks)
%                      * cfg.seizure.median_factor   (seizure_energy_trace.m)
%   ll_median_global = median(line length over the same samples)
%                      (line_length_trace.m); kept iff ll_ratio >= 1.75
% Both are recomputed here from the existing *_clean.txt with the same
% functions and the run's own config_used.mat, once over all samples
% (= original, checked against the saved metrics.median_energy) and once
% without the disconnected samples (DX intervals + FP 'disconnect' runs).
%
% Channels: in-scope EDFs whose disconnected fraction > cfg.sz_thr_diag_min_frac,
% plus every EDF holding a missed seizure with EEG (ND windows).
% Crossings use the saved candidates (Stage-3 survivors) of *_seizures.mat:
% ll_ratio' = ll_ratio * ll_med_orig / ll_med_corr (event bounds held
% fixed, hence approximate); peak_energy_ratio' = peak * thr_orig / thr_corr.
% Results are cached per channel in <out_dir>/cache.

    cache_dir = fullfile(cfg.out_dir, 'cache');
    if ~isfolder(cache_dir), mkdir(cache_dir); end
    units = S.units(S.units.included, :);
    E = edf_tbl(edf_tbl.analyzed & ismember(strcat(edf_tbl.campaign, '|', edf_tbl.mouse), ...
        strcat(units.campaign, '|', units.mouse)), :);
    E = E(~ismember(strcat(E.mouse, '|', E.edf), strcat(rv.excluded.mouse, '|', rv.excluded.edf)), :);
    nd = S.seizures(strcmp(S.seizures.source, 'missed_eeg'), :);
    ev = S.events;
    D = {}; ND = {};
    for i = 1:height(E)
        span = [E.start_posix(i) E.end_posix(i)];
        dxm = rv.dx(strcmp(rv.dx.mouse, E.mouse{i}), :);
        dxz = cs_iv('intersect', [dxm.start_posix dxm.end_posix], span);
        dce = ev(strcmp(ev.campaign, E.campaign{i}) & strcmp(ev.mouse, E.mouse{i}) & strcmp(ev.edf, E.edf{i}) & ...
            strcmp(ev.label, 'FP') & contains(ev.comments, 'disconnect'), :);
        runs = cs_disconnect_runs(ev(strcmp(ev.campaign, E.campaign{i}) & strcmp(ev.mouse, E.mouse{i}) & strcmp(ev.edf, E.edf{i}), :), ...
            cfg.disconnect_run_max_gap_s);
        zone = cs_iv('union', dxz, runs);
        frac = cs_iv('length', zone) / diff(span);
        has_nd = any(nd.start_posix >= span(1) & nd.start_posix <= span(2) & strcmp(nd.mouse, E.mouse{i}));
        if frac <= cfg.sz_thr_diag_min_frac && ~has_nd
            continue;
        end
        dx_ids = strjoin(cellstr(dxm.id(cs_iv('overlap', [dxm.start_posix dxm.end_posix], span) > 0))', ',');
        for reg = {'HPCr', 'HPCl'}
            file = E.(['clean_' reg{1}]){i};
            key = regexprep(sprintf('%s_%s_%s_%s', E.campaign{i}, E.mouse{i}, E.edf{i}, reg{1}), '[^\w-]', '_');
            cf = fullfile(cache_dir, ['thrdiag_' key '.mat']);
            r = [];
            if isfile(cf)
                r = load(cf); r = r.r;
                if ~isfield(r, 'zone') || ~isequal(size(r.zone), size(zone)) || any(abs(r.zone(:) - zone(:)) > 1e-3)
                    r = [];  % cached with another exclusion zone -> recompute
                end
            end
            if isempty(r)
                r = diag_channel(cfg, E(i, :), reg{1}, file, zone, nd(strcmp(nd.mouse, E.mouse{i}), :));
                r.zone = zone;
                save(cf, 'r');
            end
            chg_thr = 100 * (r.thr_corr / r.thr_orig - 1);
            chg_ll = 100 * (r.ll_corr / r.ll_orig - 1);
            D(end+1, :) = {E.campaign{i}, E.mouse{i}, E.edf{i}, reg{1}, frac, dx_ids, height(dce), ...
                r.saved_median_energy, r.med_orig, r.med_corr, r.thr_orig, r.thr_corr, chg_thr, r.ll_orig, r.ll_corr, chg_ll, ...
                abs(chg_thr) >= 100 * cfg.sz_thr_change_flag || abs(chg_ll) >= 100 * cfg.sz_thr_change_flag, ...
                r.pct_above_orig, r.pct_above_corr, r.n_cand, r.n_kept, r.n_kept_to_reject, r.n_reject_to_keep, ...
                r.n_cand_energy_lost, r.note}; %#ok<AGROW>
            for k = 1:numel(r.nd)
                w = r.nd(k);
                ND(end+1, :) = {w.id, E.campaign{i}, E.mouse{i}, E.edf{i}, reg{1}, w.dur_s, w.maxE_over_thr_orig, ...
                    w.maxE_over_thr_corr, w.longest_run_orig_s, w.longest_run_corr_s, w.ll_ratio_orig, w.ll_ratio_corr, ...
                    w.criterion_orig, w.criterion_corr, w.overlapping_candidates}; %#ok<AGROW>
            end
        end
    end
    D = cell2table(D, 'VariableNames', {'campaign', 'mouse', 'edf', 'region', 'disconnected_frac', 'dx_ids', ...
        'n_FP_disconnect_events', 'saved_median_energy', 'median_energy_orig', 'median_energy_corr', 'energy_thr_orig', ...
        'energy_thr_corr', 'energy_thr_change_pct', 'll_median_orig', 'll_median_corr', 'll_median_change_pct', ...
        'flag_change_ge_5pct', 'pct_time_above_thr_orig', 'pct_time_above_thr_corr', 'n_candidates', 'n_kept', ...
        'n_kept_would_be_rejected_ll', 'n_rejected_would_be_kept_ll', 'n_candidates_peak_below_new_thr', 'note'});
    ND = cell2table(ND, 'VariableNames', {'nd_id', 'campaign', 'mouse', 'edf', 'region', 'window_s', ...
        'maxE_over_thr_orig', 'maxE_over_thr_corr', 'longest_supra_run_orig_s', 'longest_supra_run_corr_s', ...
        'll_ratio_orig', 'll_ratio_corr', 'failing_criterion_orig', 'failing_criterion_corr', 'overlapping_candidates'});
end

function r = diag_channel(cfg, Erow, region, file, zone, nd)
    t0 = tic;
    root = cfg.campaigns(strcmp({cfg.campaigns.name}, Erow.campaign{1})).root;
    C = load(fullfile(root, Erow.mouse{1}, 'config_used.mat'));
    C = C.cfg;
    data = load_lfp_txt(file);
    fs = data.fs;
    n = numel(data.signal);
    s0 = data.session_start;
    s0 = posixtime(s0);
    trace = seizure_energy_trace(data, C);
    tb = [vertcat(trace.blocks.trimmed_start), vertcat(trace.blocks.trimmed_end)];
    [ll_full, ll0] = line_length_trace(trace.bp_full, tb, fs, C);
    inb = false(n, 1);
    for b = 1:size(tb, 1)
        inb(tb(b, 1):tb(b, 2)) = true;
    end
    zm = false(n, 1);
    for z = 1:size(zone, 1)
        a = max(1, floor((zone(z, 1) - s0) * fs) + 1);
        b = min(n, ceil((zone(z, 2) - s0) * fs) + 1);
        if b >= a, zm(a:b) = true; end
    end
    En = trace.energy_full;
    r.med_orig = median(En(inb));
    r.med_corr = median(En(inb & ~zm));
    r.thr_orig = r.med_orig * C.seizure.median_factor;
    r.thr_corr = r.med_corr * C.seizure.median_factor;
    r.ll_orig = ll0;
    r.ll_corr = median(ll_full(inb & ~zm));
    r.pct_above_orig = 100 * mean(En(inb) > r.thr_orig);
    r.pct_above_corr = 100 * mean(En(inb & ~zm) > r.thr_corr);
    r.saved_median_energy = NaN;
    r.n_cand = NaN; r.n_kept = NaN; r.n_kept_to_reject = NaN; r.n_reject_to_keep = NaN; r.n_cand_energy_lost = NaN;
    r.note = '';
    [~, base] = fileparts(file);
    base = regexprep(base, '_clean$', '');
    mf = fullfile(root, Erow.mouse{1}, '03_seizures', Erow.edf{1}, [base '_seizures.mat']);
    cand = table();
    if isfile(mf)
        sr = load(mf, 'seizure_results');
        sr = sr.seizure_results;
        r.saved_median_energy = sr.metrics.median_energy;
        if abs(r.saved_median_energy / r.med_orig - 1) > 1e-6
            r.note = sprintf('median_energy recomputed differs from saved by %.3g%%', 100 * (r.med_orig / r.saved_median_energy - 1));
        end
        cand = sr.candidates;
        llc = cand.ll_ratio * r.ll_orig / r.ll_corr;
        pkc = cand.peak_energy_ratio * r.thr_orig / r.thr_corr;
        kept = logical(cand.kept);
        r.n_cand = height(cand); r.n_kept = sum(kept);
        r.n_kept_to_reject = sum(kept & llc < C.seizure_robust.ll_threshold);
        r.n_reject_to_keep = sum(~kept & llc >= C.seizure_robust.ll_threshold);
        r.n_cand_energy_lost = sum(pkc < 1);
    else
        r.note = 'sin *_seizures.mat (0 detecciones): sin candidatos guardados';
    end
    % missed-seizure windows in this file
    r.nd = struct('id', {}, 'dur_s', {}, 'maxE_over_thr_orig', {}, 'maxE_over_thr_corr', {}, 'longest_run_orig_s', {}, ...
        'longest_run_corr_s', {}, 'll_ratio_orig', {}, 'll_ratio_corr', {}, 'criterion_orig', {}, 'criterion_corr', {}, ...
        'overlapping_candidates', {});
    for k = 1:height(nd)
        a = floor((nd.start_posix(k) - s0) * fs) + 1;
        b = ceil((nd.end_posix(k) - s0) * fs) + 1;
        if a < 1 || b > n, continue; end
        w = struct();
        w.id = nd.origin_ids{k};
        w.dur_s = nd.duration_s(k);
        mE = max(En(a:b));
        w.maxE_over_thr_orig = mE / r.thr_orig;
        w.maxE_over_thr_corr = mE / r.thr_corr;
        w.longest_run_orig_s = longest_run(En, a, b, r.thr_orig, fs, C);
        w.longest_run_corr_s = longest_run(En, a, b, r.thr_corr, fs, C);
        w.ll_ratio_orig = median(ll_full(a:b), 'omitnan') / r.ll_orig;
        w.ll_ratio_corr = median(ll_full(a:b), 'omitnan') / r.ll_corr;
        w.criterion_orig = criterion(w.maxE_over_thr_orig, w.longest_run_orig_s, w.ll_ratio_orig, C);
        w.criterion_corr = criterion(w.maxE_over_thr_corr, w.longest_run_corr_s, w.ll_ratio_corr, C);
        w.overlapping_candidates = '';
        if ~isempty(cand)
            ws = (a - 1) / fs; we = (b - 1) / fs;
            o = find(cand.start_s < we & cand.end_s > ws);
            w.overlapping_candidates = strjoin(arrayfun(@(q) sprintf('cand %d kept=%d ll=%.2f', cand.candidate_id(q), ...
                cand.kept(q), cand.ll_ratio(q)), o, 'UniformOutput', false)', '; ');
        end
        r.nd(end+1) = w;
    end
    fprintf('thr-diag %s %s %s %s: %.0f s\n', Erow.campaign{1}, Erow.mouse{1}, Erow.edf{1}, region, toc(t0));
end

function d = longest_run(En, a, b, thr, fs, C)
% Longest merged supra-threshold segment touching [a b], built like the
% detector's Stage 3 (merge gaps <= merge_gap_s) over a +/-120 s context.
    pad = round(120 * fs);
    lo = max(1, a - pad); hi = min(numel(En), b + pad);
    seg = mask_to_segments(En(lo:hi) > thr);
    d = 0;
    if isempty(seg), return; end
    gap = round(C.seizure_robust.merge_gap_s * fs);
    m = seg(1, :);
    for k = 2:size(seg, 1)
        if seg(k, 1) - m(end, 2) - 1 <= gap
            m(end, 2) = seg(k, 2);
        else
            m(end+1, :) = seg(k, :); %#ok<AGROW>
        end
    end
    m = m + lo - 1;
    t = m(m(:, 1) <= b & m(:, 2) >= a, :);
    if ~isempty(t)
        d = max(t(:, 2) - t(:, 1) + 1) / fs;
    end
end

function c = criterion(maxE, run_s, llr, C)
    if maxE < 1
        c = 'energia bajo umbral';
    elseif run_s < C.seizure_robust.min_duration_s
        c = sprintf('duracion supra-umbral < %g s', C.seizure_robust.min_duration_s);
    elseif llr < C.seizure_robust.ll_threshold
        c = sprintf('ll_ratio < %g', C.seizure_robust.ll_threshold);
    else
        c = 'cumple los 3 criterios (revisar)';
    end
end
