function report = evaluate_detections(seizures_csv, ground_truth_csv, cfg, opts)
% EVALUATE_DETECTIONS  Score a seizures_events.csv against expert-reviewed
% ground truth. Without this, any future threshold change is made blind --
% this is what turns "I think it looks better" into a measured comparison.
%
%   report = evaluate_detections(seizures_csv, ground_truth_csv, cfg)
%   report = evaluate_detections(seizures_csv, ground_truth_csv, cfg, opts)
%
%   seizures_csv     : path to a seizures_events.csv (any pipeline run --
%                      legacy or robust, this function only reads the
%                      table, see seizure_mode column note below)
%   ground_truth_csv : subject_id,region,start_abs,end_abs,notes
%                      (region empty = must be found in ANY channel of
%                      that subject; a non-empty region means the expert
%                      localized that event to that channel specifically)
%   cfg              : struct from pipeline_config.m (uses cfg.eval.match_tol_s)
%   opts (optional) struct:
%     .clean_txt_lookup : containers.Map region -> *_clean.txt path.
%                          Enables the "why was this missed" / FP ll_ratio
%                          diagnostics (peak energy, longest run above
%                          threshold, ll_ratio), by reconstructing the
%                          energy/line-length traces for that region via
%                          seizure_energy_trace.m / line_length_trace.m.
%                          Without it, those columns are NaN.
%     .candidates_table  : table of ALL candidate events (matched AND
%                          line-length-rejected, with an ll_ratio column)
%                          from detect_seizures_robust.m, across whatever
%                          ll_threshold values you want swept. Enables the
%                          precision/recall-vs-ll_threshold curve. Without
%                          it, .pr_curve is empty (the legacy branch has no
%                          ll_ratio-tagged candidates to sweep).
%
% OUTPUT (struct report):
%   .per_region   : table (region, n_gt, n_detected, TP, FP, FN,
%                    fragmentation, precision, recall, f1)
%   .per_subject  : same, aggregated so a GT event found in ANY applicable
%                    channel counts once (answers "was this crisis caught
%                    at all", separate from the per-channel question)
%   .missed       : table of (gt_id, region) misses with peak_energy_ratio,
%                    longest_run_above_thr_s, ll_ratio (NaN if opts.clean_txt_lookup
%                    was not given)
%   .false_positives : table of unmatched detections, with ll_ratio if available,
%                    sorted by ll_ratio descending (closest to a line-length
%                    cutoff first -- the ones that most need a human look)
%   .fragments    : table of extra detections matched to an already-matched
%                    GT event (NOT counted as false positives)
%   .pr_curve     : table (ll_threshold, precision, recall, f1) if
%                    opts.candidates_table was given, else empty
%   .detections, .ground_truth : the tables actually used, annotated with
%                    match_status ('tp'|'fragment'|'fp') / found (bool)

    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'clean_txt_lookup')
        opts.clean_txt_lookup = containers.Map('KeyType', 'char', 'ValueType', 'char');
    end
    if ~isfield(opts, 'candidates_table')
        opts.candidates_table = table();
    end

    tz = cfg.general.timezone;
    det = read_pipeline_csv(seizures_csv, 'seizures_events', tz);
    gt = read_ground_truth(ground_truth_csv, tz);

    tol = seconds(cfg.eval.match_tol_s);

    det.match_status = repmat({'fp'}, height(det), 1);
    det.gt_id = nan(height(det), 1);
    gt.found_any = false(height(gt), 1);
    gt.found_regions = repmat({{}}, height(gt), 1);

    % Regions to evaluate = union of {regions with at least one detection}
    % and {regions the caller supplied a clean-txt lookup for}: a channel
    % that legitimately detected ZERO events still has to be evaluated
    % (every one of its applicable GT events is a miss), and it would
    % otherwise silently disappear since it contributes no rows to
    % seizures_events.csv at all.
    regions_present = union(unique(det.region, 'stable'), keys(opts.clean_txt_lookup));
    subjects_present = unique(gt.subject_id, 'stable');
    per_region_rows = cell(0, 10);
    missed_rows = cell(0, 6);

    for r = 1:numel(regions_present)
        region = regions_present{r};

        for s = 1:numel(subjects_present)
            subject_id = subjects_present{s};
            det_mask = strcmp(det.region, region) & strcmp(det.subject_id, subject_id);
            gt_here = strcmp(gt.subject_id, subject_id) & (strcmp(gt.region, region) | strcmp(gt.region, ''));
            if ~any(det_mask) && ~any(gt_here)
                continue;  % nothing to evaluate for this (region, subject) pair
            end
            gt_idx = find(gt_here);

            tp = 0; frag = 0;
            for g = gt_idx(:)'
                cand = find(det_mask & (det.start_abs <= gt.end_abs(g) + tol) & (det.end_abs >= gt.start_abs(g) - tol));
                cand = cand(isnan(det.gt_id(cand)));  % don't re-claim a detection already matched to a different GT
                if isempty(cand)
                    missed_rows(end+1, :) = {gt.gt_id(g), region, subject_id, gt.start_abs(g), gt.end_abs(g), NaN}; %#ok<AGROW>
                    continue;
                end
                [~, order] = sort(det.start_abs(cand));
                cand = cand(order);
                det.match_status{cand(1)} = 'tp';
                det.gt_id(cand(1)) = gt.gt_id(g);
                tp = tp + 1;
                gt.found_any(g) = true;
                gt.found_regions{g} = [gt.found_regions{g}, {region}];
                if numel(cand) > 1
                    det.match_status(cand(2:end)) = {'fragment'};
                    det.gt_id(cand(2:end)) = gt.gt_id(g);
                    frag = frag + numel(cand) - 1;
                end
            end

            fp = sum(det_mask & strcmp(det.match_status, 'fp'));
            fn = numel(gt_idx) - tp;
            precision = safe_div(tp, tp + fp);
            recall = safe_div(tp, tp + fn);
            f1 = safe_div(2 * precision * recall, precision + recall);
            per_region_rows(end+1, :) = {region, subject_id, numel(gt_idx), sum(det_mask), tp, fp, fn, frag, precision, recall}; %#ok<AGROW>
        end
    end

    per_region = cell2table(per_region_rows, 'VariableNames', ...
        {'region', 'subject_id', 'n_gt', 'n_detected', 'TP', 'FP', 'FN', 'fragmentation', 'precision', 'recall'});
    per_region.f1 = safe_div(2 .* per_region.precision .* per_region.recall, per_region.precision + per_region.recall);

    missed = cell2table(missed_rows, 'VariableNames', {'gt_id', 'region', 'subject_id', 'start_abs', 'end_abs', 'unused'});
    missed.unused = [];
    missed = add_miss_diagnostics(missed, opts.clean_txt_lookup, cfg);

    fp_mask = strcmp(det.match_status, 'fp');
    false_positives = det(fp_mask, {'subject_id', 'region', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs'});
    false_positives = add_fp_diagnostics(false_positives, opts.clean_txt_lookup, cfg);
    if height(false_positives) > 0
        [~, order] = sort(false_positives.ll_ratio, 'descend', 'MissingPlacement', 'last');
        false_positives = false_positives(order, :);
    end

    frag_mask = strcmp(det.match_status, 'fragment');
    fragments = det(frag_mask, {'subject_id', 'region', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', 'gt_id'});

    per_subject = build_per_subject(gt, det, regions_present);

    pr_curve = table();
    if height(opts.candidates_table) > 0
        pr_curve = sweep_ll_threshold(opts.candidates_table, gt, cfg);
    end

    report = struct();
    report.per_region = per_region;
    report.per_subject = per_subject;
    report.missed = missed;
    report.false_positives = false_positives;
    report.fragments = fragments;
    report.pr_curve = pr_curve;
    report.detections = det;
    report.ground_truth = gt;

    fprintf('evaluate_detections: %d GT event(s), %d detection(s) -> TP=%d FP=%d FN=%d fragmentation=%d\n', ...
        height(gt), height(det), sum(strcmp(det.match_status, 'tp')), sum(fp_mask), ...
        sum(per_region.FN), sum(frag_mask));
end

%% ======================================================================
function gt = read_ground_truth(path, tz)
    if exist(path, 'file') ~= 2
        error('evaluate_detections:GTNotFound', 'Ground truth CSV not found: %s', path);
    end
    opts = detectImportOptions(path, 'VariableNamingRule', 'preserve');
    text_cols = {'subject_id', 'region', 'start_abs', 'end_abs', 'notes'};
    for i = 1:numel(text_cols)
        if ismember(text_cols{i}, opts.VariableNames)
            opts = setvartype(opts, text_cols{i}, 'char');
        end
    end
    T = readtable(path, opts);
    n = height(T);
    for i = 1:numel(text_cols)
        col = text_cols{i};
        if ~ismember(col, T.Properties.VariableNames)
            T.(col) = repmat({''}, n, 1);
        else
            T.(col) = cellfun(@(v) strtrim(char(v)), T.(col), 'UniformOutput', false);
        end
    end
    start_dt = datetime(T.start_abs, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', tz);
    end_dt = datetime(T.end_abs, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', tz);
    gt = table((1:n)', T.subject_id, T.region, start_dt, end_dt, T.notes, ...
        'VariableNames', {'gt_id', 'subject_id', 'region', 'start_abs', 'end_abs', 'notes'});
end

function v = safe_div(a, b)
    v = a ./ b;
    v(b == 0) = 0;
end

function per_subject = build_per_subject(gt, det, regions_present)
    subjects = unique(gt.subject_id, 'stable');
    rows = cell(0, 7);
    for s = 1:numel(subjects)
        subject_id = subjects{s};
        mask = strcmp(gt.subject_id, subject_id);
        n_gt = sum(mask);
        tp = sum(gt.found_any(mask));
        fn = n_gt - tp;
        fp = sum(strcmp(det.subject_id, subject_id) & strcmp(det.match_status, 'fp'));
        precision = safe_div(tp, tp + fp);
        recall = safe_div(tp, tp + fn);
        f1 = safe_div(2 * precision * recall, precision + recall);
        rows(end+1, :) = {subject_id, n_gt, tp, fp, fn, precision, recall}; %#ok<AGROW>
    end
    per_subject = cell2table(rows, 'VariableNames', {'subject_id', 'n_gt', 'TP', 'FP', 'FN', 'precision', 'recall'});
    per_subject.f1 = safe_div(2 .* per_subject.precision .* per_subject.recall, per_subject.precision + per_subject.recall);
end

%% ======================================================================
function missed = add_miss_diagnostics(missed, lookup, cfg)
    n = height(missed);
    peak_energy_ratio = nan(n, 1);
    longest_run_s = nan(n, 1);
    ll_ratio = nan(n, 1);

    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:n
        region = missed.region{i};
        if ~isKey(lookup, region)
            continue;
        end
        tr = get_cached_traces(cache, lookup, region, cfg);
        if isempty(tr)
            continue;
        end
        [s_idx, e_idx] = window_to_indices(missed.start_abs(i), missed.end_abs(i), tr.session_start, tr.fs, tr.n, cfg.eval.match_tol_s);
        if isempty(s_idx)
            continue;
        end
        peak_energy_ratio(i) = max(tr.energy(s_idx:e_idx), [], 'omitnan') / tr.threshold;
        longest_run_s(i) = longest_run_above(tr.energy(s_idx:e_idx), tr.threshold) / tr.fs;

        [s_idx2, e_idx2] = window_to_indices(missed.start_abs(i), missed.end_abs(i), tr.session_start, tr.fs, tr.n, 0);
        if ~isempty(s_idx2)
            ll_ratio(i) = median(tr.ll(s_idx2:e_idx2), 'omitnan') / tr.ll_median_global;
        end
    end

    missed.peak_energy_ratio = peak_energy_ratio;
    missed.longest_run_above_thr_s = longest_run_s;
    missed.ll_ratio = ll_ratio;
end

function fp = add_fp_diagnostics(fp, lookup, cfg)
    n = height(fp);
    ll_ratio = nan(n, 1);
    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for i = 1:n
        region = fp.region{i};
        if ~isKey(lookup, region)
            continue;
        end
        tr = get_cached_traces(cache, lookup, region, cfg);
        if isempty(tr)
            continue;
        end
        [s_idx, e_idx] = window_to_indices(fp.start_abs(i), fp.end_abs(i), tr.session_start, tr.fs, tr.n, 0);
        if isempty(s_idx)
            continue;
        end
        ll_ratio(i) = median(tr.ll(s_idx:e_idx), 'omitnan') / tr.ll_median_global;
    end
    fp.ll_ratio = ll_ratio;
end

function tr = get_cached_traces(cache, lookup, region, cfg)
    if isKey(cache, region)
        tr = cache(region);
        return;
    end
    try
        data = load_lfp_txt(lookup(region));
        etrace = seizure_energy_trace(data, cfg);
        trimmed_blocks = zeros(numel(etrace.blocks), 2);
        for b = 1:numel(etrace.blocks)
            trimmed_blocks(b, :) = [etrace.blocks(b).trimmed_start, etrace.blocks(b).trimmed_end];
        end
        [ll_full, ll_med] = line_length_trace(etrace.bp_full, trimmed_blocks, etrace.fs, cfg);
        tr = struct('energy', etrace.energy_full, 'threshold', etrace.threshold, 'll', ll_full, ...
            'll_median_global', ll_med, 'fs', etrace.fs, 'session_start', data.session_start, 'n', numel(data.signal));
    catch ME
        warning('evaluate_detections:TraceReconstructionFailed', 'Could not reconstruct traces for region "%s": %s', region, ME.message);
        tr = [];
    end
    cache(region) = tr;
end

function [s_idx, e_idx] = window_to_indices(start_abs, end_abs, session_start, fs, n, pad_s)
    if isnat(session_start)
        s_idx = []; e_idx = [];
        return;
    end
    start_s = seconds(start_abs - session_start) - pad_s;
    end_s = seconds(end_abs - session_start) + pad_s;
    s_idx = max(1, round(start_s * fs) + 1);
    e_idx = min(n, round(end_s * fs) + 1);
    if s_idx > e_idx || s_idx > n || e_idx < 1
        s_idx = []; e_idx = [];
    end
end

function run_len = longest_run_above(energy_window, threshold)
    above = energy_window > threshold;
    above(isnan(energy_window)) = false;
    d = diff([0; above(:); 0]);
    starts = find(d == 1);
    ends = find(d == -1) - 1;
    if isempty(starts)
        run_len = 0;
    else
        run_len = max(ends - starts + 1);
    end
end

%% ======================================================================
function pr_curve = sweep_ll_threshold(candidates, gt, cfg)
% candidates: table with at least region, subject_id, start_abs, end_abs, ll_ratio.
% Mirrors the main matching semantics exactly: a kept candidate counts as
% matched (not a false positive) if it overlaps ANY applicable GT event,
% even if another kept candidate already matched that same GT event (that
% would be fragmentation in the main report, not a reason to inflate FP
% here) -- and a GT event is a true positive if ANY kept candidate
% overlaps it, regardless of how many do.
    thresholds = 1.0:0.05:3.0;
    rows = cell(numel(thresholds), 4);
    tol = seconds(cfg.eval.match_tol_s);
    for k = 1:numel(thresholds)
        thr = thresholds(k);
        kept = candidates(candidates.ll_ratio >= thr, :);
        kept_matched = false(height(kept), 1);
        gt_found = false(height(gt), 1);
        for g = 1:height(gt)
            cand_mask = strcmp(kept.subject_id, gt.subject_id{g}) & ...
                (strcmp(gt.region{g}, '') | strcmp(kept.region, gt.region{g})) & ...
                (kept.start_abs <= gt.end_abs(g) + tol) & (kept.end_abs >= gt.start_abs(g) - tol);
            if any(cand_mask)
                gt_found(g) = true;
                kept_matched(cand_mask) = true;
            end
        end
        tp = sum(gt_found);
        fn = height(gt) - tp;
        fp = sum(~kept_matched);
        precision = safe_div(tp, tp + fp);
        recall = safe_div(tp, tp + fn);
        rows(k, :) = {thr, precision, recall, safe_div(2 * precision * recall, precision + recall)};
    end
    pr_curve = cell2table(rows, 'VariableNames', {'ll_threshold', 'precision', 'recall', 'f1'});
end
