function results = bilateral_reconcile(per_channel_results, cfg)
% BILATERAL_RECONCILE  Cross-channel seizure reconciliation: a seizure is
% dropped only if EVERY channel of the animal rejected it; if at least one
% channel accepted it, it is reported in all of them, under one shared
% seizure_id, with each row saying whether that channel accepted it on its
% own merit ('accepted'), had it as a rejected candidate ('rescued'), or
% had nothing there at all ('imputed').
%
%   results = bilateral_reconcile(per_channel_results, cfg)
%
% Runs AFTER seizure detection on every channel of one recording and
% BEFORE any seizure CSV row, figure or IID exclusion zone is built. It
% never re-detects and never changes a threshold: it only regroups and
% measures what detect_seizures_robust.m returned.
%
% ROBUST BRANCH ONLY. Only channels whose case resolved to
% seizure_mode='robust' take part; a legacy channel (detect_seizures.m) is
% returned untouched (trace = []) and is neither reconciled nor used as a
% reference for other channels -- run_pipeline_edf.m tags its rows
% 'not_reconciled'.
%
% INPUT per_channel_results: struct array, one element per channel, with
%   .region, .subject_id, .source_file, .session_start
%   .seizure_mode      'legacy' | 'robust'
%   .seizure_results   the detector's output ([] if the clean or detection
%                      stage failed -- such a robust channel contributes no
%                      events but still gets one 'imputed' row per event of
%                      its session, with NO metrics and
%                      rejected_by='detection_failed', and counts in
%                      n_channels_total, so the failure is reported)
%   .data              the *_clean.txt struct the detector ran on
%
% OUTPUT results: the same array, where each reconciled channel's
%   .seizure_results.seizures           is REPLACED by the reconciled table
%   .seizure_results.seizures_accepted  keeps the detector's own table
%   .seizure_results.bilateral          counters (n_accepted, n_reported,
%                                       n_events, n_rescued, n_imputed, ...)
%   .trace                              energy / band-passed / line-length
%                                       traces of that channel (for figures)
%
% Reconciled table columns: id (the detector's own per-channel index, NaN
% on rescued/imputed rows), seizure_id (shared), start_s, end_s,
% duration_s, start_abs, end_abs, block_id, over_max_duration, ll_ratio,
% peak_energy_ratio, hf_ratio_db, envelope_cv, detection_status,
% accepted_in_n_channels, accepted_in_regions, n_channels_total,
% is_bilateral_accepted, fragmented, ref_start_s, ref_end_s, rejected_by.
%
% ALGORITHM
%  1. Group: accepted detections of DIFFERENT channels belong to the same
%     event if they overlap or are less than cfg.bilateral.match_tol_s
%     apart; grouping is transitive (union-find). Two detections of the
%     SAME channel are never linked directly; if they end up in one group
%     through another channel they stay separate rows with the same
%     seizure_id and fragmented=true. Groups are numbered 1..K by the
%     earliest start in the group. Reference window = [earliest start,
%     latest end] of the group's accepted rows.
%  2. For each event and each channel without an accepted row in it:
%     rescue = the channel's rejected candidate closest in spirit (the one
%     with the highest ll_ratio) among those within match_tol_s of the
%     reference window, WITH ITS OWN limits; else, in 'rescue_and_impute',
%     impute = the reference window, clipped to the single valid block of
%     that channel it overlaps most (a row never crosses a gap or a block
%     edge). In 'rescue' mode an event with no candidate gets no row.
%     A candidate is rescued into at most one event (the one whose
%     reference window it is closest to), never reused for two.
%  3. Metrics (ll_ratio, peak_energy_ratio, hf_ratio_db, envelope_cv,
%     duration) are computed on each row's own window in its own channel
%     -- never copied across channels. Accepted and rescued rows keep the
%     detector's own values (they were computed on exactly that window);
%     imputed rows get them computed here with the same formulas
%     (see detect_seizures_robust.m).
%
% "Rejected candidate" = the detector's .candidates rows with kept=false
% (rejected by ll_ratio); rejected_by = 'll_ratio'. Crossings the robust
% branch dropped for being shorter than cfg.seizure_robust.min_duration_s
% after merging are not exposed by the detector (detect_seizures_robust.m
% is not modified), so they are not rescue candidates: such an event is
% imputed instead.
%
% Channels are grouped only within the same (subject_id, source_file,
% session_start) -- or (source_file, session_start) if
% cfg.bilateral.require_same_subject is false. A single-channel recording
% passes through unchanged apart from receiving its seizure_id.

    results = per_channel_results;
    if isempty(results)
        return;
    end
    [results.trace] = deal([]);

    mode = cfg.bilateral.rescue_mode;
    if strcmp(mode, 'off')
        return;
    end
    if ~ismember(mode, {'rescue', 'rescue_and_impute'})
        error('bilateral_reconcile:BadMode', ...
            'cfg.bilateral.rescue_mode must be ''off'', ''rescue'' or ''rescue_and_impute'' (got ''%s'').', mode);
    end

    robust = arrayfun(@(p) strcmp(p.seizure_mode, 'robust'), results);
    usable = robust & arrayfun(@(p) ~isempty(p.seizure_results) && ~isempty(p.data), results);
    failed = robust & ~usable;  % robust channel whose clean/detection stage failed
    keys = arrayfun(@(p) session_key(p, cfg), results, 'UniformOutput', false);
    ukeys = unique(keys(usable), 'stable');
    for k = 1:numel(ukeys)
        in_session = strcmp(keys, ukeys{k});
        results = reconcile_session(results, find(usable & in_session), find(failed & in_session), cfg);
    end
end

%% ======================================================================
function key = session_key(p, cfg)
    if isnat(p.session_start)
        ss = 'NaT';
    else
        ss = char(string(p.session_start, 'yyyy-MM-dd''T''HH:mm:ss.SSS'));
    end
    if cfg.bilateral.require_same_subject
        key = sprintf('%s|%s|%s', p.subject_id, p.source_file, ss);
    else
        key = sprintf('%s|%s', p.source_file, ss);
    end
end

function results = reconcile_session(results, idx, idx_failed, cfg)
    mode = cfg.bilateral.rescue_mode;
    tol = cfg.bilateral.match_tol_s;
    nch = numel(idx);
    n_total = nch + numel(idx_failed);
    regions = arrayfun(@(p) p.region, results(idx), 'UniformOutput', false);

    %% per-channel inputs: accepted table, traces, rejected candidates
    ch = struct('acc', {}, 'tr', {}, 'rej', {});
    for j = 1:nch
        p = results(idx(j));
        ch(j).acc = p.seizure_results.seizures;
        ch(j).tr = channel_trace(p.data, cfg);
        ch(j).rej = rejected_candidates(p);
    end

    %% step 1: group accepted detections across channels
    A = zeros(0, 4);  % [channel j, row r, start_s, end_s]
    for j = 1:nch
        a = ch(j).acc;
        A = [A; repmat(j, height(a), 1), (1:height(a))', a.start_s, a.end_s]; %#ok<AGROW>
    end
    nA = size(A, 1);
    parent = 1:nA;
    for a1 = 1:nA
        for a2 = a1+1:nA
            if A(a1, 1) ~= A(a2, 1) && interval_gap(A(a1, 3:4), A(a2, 3:4)) < tol
                parent = uf_union(parent, a1, a2);
            end
        end
    end
    roots = arrayfun(@(i) uf_find(parent, i), 1:nA);
    [ur, ~, gidx] = unique(roots);
    nG = numel(ur);
    g_start = accumarray(gidx(:), A(:, 3), [nG 1], @min);
    g_end = accumarray(gidx(:), A(:, 4), [nG 1], @max);
    [~, order] = sortrows([g_start, g_end]);
    seizure_id_of_group = zeros(nG, 1);
    seizure_id_of_group(order) = 1:nG;
    ref_s = zeros(nG, 1); ref_e = zeros(nG, 1);
    ref_s(seizure_id_of_group) = g_start;
    ref_e(seizure_id_of_group) = g_end;
    A_sid = seizure_id_of_group(gidx(:));  % shared id of each accepted detection

    acc_chan = false(nG, nch);
    for i = 1:nA
        acc_chan(A_sid(i), A(i, 1)) = true;
    end

    %% step 2+3: rows per channel
    for j = 1:nch
        c = idx(j);
        p = results(c);
        tr = ch(j).tr;
        acc = ch(j).acc;
        rej = ch(j).rej;
        rows = {};
        n_acc = 0; n_res = 0; n_imp = 0; n_unfilled = 0;

        % A rejected candidate is rescued into at most ONE event: the one
        % whose reference window it is closest to (ties -> earliest), so one
        % candidate never shows up as two different seizures.
        cand_event = zeros(height(rej), 1);
        for q = 1:height(rej)
            g = max(rej.start_s(q), ref_s) - min(rej.end_s(q), ref_e);  % gap to every event's reference window
            g(g >= tol) = Inf;
            [gmin, best_sid] = min(g);
            if ~isempty(gmin) && isfinite(gmin)
                cand_event(q) = best_sid;
            end
        end

        for sid = 1:nG
            ev = struct('sid', sid, 'ref_s', ref_s(sid), 'ref_e', ref_e(sid), ...
                'n_acc', nnz(acc_chan(sid, :)), 'regions', strjoin(regions(acc_chan(sid, :)), ';'), ...
                'n_total', n_total);
            mine = find(A(:, 1) == j & A_sid == sid);
            if ~isempty(mine)
                for m = mine'
                    r = A(m, 2);
                    rows{end+1} = accepted_row(acc(r, :), p, ev, numel(mine) > 1); %#ok<AGROW>
                    n_acc = n_acc + 1;
                end
                continue;
            end

            near = find(cand_event == sid);
            if ~isempty(near)
                ll = rej.ll_ratio(near);
                ll(isnan(ll)) = -Inf;
                [~, best] = max(ll);
                q = near(best);
                rows{end+1} = rescued_row(rej(q, :), p, cfg, ev); %#ok<AGROW>
                n_res = n_res + 1;
            elseif strcmp(mode, 'rescue_and_impute')
                [row, ok] = imputed_row(p, tr, cfg, ev);
                rows{end+1} = row; %#ok<AGROW>
                n_imp = n_imp + 1;
                if ~ok
                    n_unfilled = n_unfilled + 1;
                end
            end
        end

        if isempty(rows)
            T = empty_reconciled_table(p.session_start);
        else
            T = vertcat(rows{:});
            T = sortrows(T, {'seizure_id', 'start_s'});
        end

        sr = p.seizure_results;
        sr.seizures_accepted = acc;
        sr.seizures = T;
        sr.bilateral = struct('rescue_mode', mode, 'match_tol_s', tol, 'n_channels_total', n_total, ...
            'n_accepted', n_acc, 'n_reported', height(T), 'n_events', numel(unique(T.seizure_id)), ...
            'n_rescued', n_res, 'n_imputed', n_imp, 'n_imputed_no_valid_signal', n_unfilled, ...
            'n_events_session', nG);
        if n_unfilled > 0
            warning('bilateral_reconcile:NoValidSignal', ...
                '%s: %d imputed event(s) fall outside every valid block of this channel; their metrics are NaN.', ...
                p.data.file, n_unfilled);
        end
        results(c).seizure_results = sr;
        results(c).trace = tr;
    end

    %% robust channels whose clean/detection stage failed: every event of
    %% the session is still reported there, as 'imputed' with the reference
    %% window and NO metrics (there is no usable trace to measure), tagged
    %% rejected_by='detection_failed', so the failure is visible per event.
    for c = idx_failed(:)'
        p = results(c);
        rows = cell(1, nG);
        for sid = 1:nG
            ev = struct('sid', sid, 'ref_s', ref_s(sid), 'ref_e', ref_e(sid), ...
                'n_acc', nnz(acc_chan(sid, :)), 'regions', strjoin(regions(acc_chan(sid, :)), ';'), ...
                'n_total', n_total);
            if ~isempty(p.data)
                dur = (idx_of(ev.ref_e, p.data.fs) - idx_of(ev.ref_s, p.data.fs) + 1) / p.data.fs;
            else
                dur = ev.ref_e - ev.ref_s;
            end
            m = struct('ll_ratio', NaN, 'peak_energy_ratio', NaN, 'hf_ratio_db', NaN, 'envelope_cv', NaN, ...
                'over_max_duration', false);
            rows{sid} = make_row(NaN, ev, ev.ref_s, ev.ref_e, dur, NaN, m, 'imputed', 'detection_failed', false, p.session_start);
        end
        empty_acc = empty_reconciled_table(p.session_start);
        if nG == 0
            T = empty_acc;
        else
            T = vertcat(rows{:});
        end
        sr = struct('seizures', T, 'seizures_accepted', empty_acc, 'candidates', [], ...
            'metrics', struct('median_energy', NaN, 'threshold', NaN, 'pct_above', NaN, 'n_segments', NaN, 'n_rejected', NaN), ...
            'qc', struct('n_blocks_total', NaN, 'n_blocks_rejected_short', NaN, 'rejected_blocks', table()), ...
            'detection_failed', true, 'mat_file', '', 'fig_file', '', 'png_file', '');
        sr.bilateral = struct('rescue_mode', mode, 'match_tol_s', tol, 'n_channels_total', n_total, ...
            'n_accepted', 0, 'n_reported', height(T), 'n_events', nG, 'n_rescued', 0, 'n_imputed', nG, ...
            'n_imputed_no_valid_signal', nG, 'n_events_session', nG);
        results(c).seizure_results = sr;
        results(c).trace = [];
    end
end

%% ======================================================================
function tr = channel_trace(data, cfg)
    t = seizure_energy_trace(data, cfg);
    if isempty(t.blocks)
        trimmed = zeros(0, 2);
        block_ids = zeros(0, 1);
    else
        trimmed = [vertcat(t.blocks.trimmed_start), vertcat(t.blocks.trimmed_end)];
        block_ids = vertcat(t.blocks.block_id);
    end
    [ll_full, ll_med] = line_length_trace(t.bp_full, trimmed, data.fs, cfg);
    tr = struct('energy_full', t.energy_full, 'bp_full', t.bp_full, 'll_full', ll_full, ...
        'll_median_global', ll_med, 'threshold', t.threshold, 'trimmed', trimmed, ...
        'block_ids', block_ids, 'blocks', t.blocks);
end

function rej = rejected_candidates(p)
% detect_seizures_robust.m's candidates rejected by the line-length filter.
    cand = p.seizure_results.candidates;
    cand = cand(~cand.kept, :);
    rej = table(cand.start_s, cand.end_s, cand.duration_s, cand.block_id, cand.ll_ratio, ...
        cand.peak_energy_ratio, cand.hf_ratio_db, cand.envelope_cv, repmat({'ll_ratio'}, height(cand), 1), ...
        'VariableNames', {'start_s', 'end_s', 'duration_s', 'block_id', 'll_ratio', ...
        'peak_energy_ratio', 'hf_ratio_db', 'envelope_cv', 'rejected_by'});
end

function m = window_metrics(tr, data, gs, ge, cfg)
% Same formulas as detect_seizures_robust.m's Stage-3/4 confidence metrics.
    fs = data.fs;
    energy_seg = tr.energy_full(gs:ge);
    ll_seg = tr.ll_full(gs:ge);
    m.ll_ratio = median(ll_seg, 'omitnan') / tr.ll_median_global;
    m.peak_energy_ratio = max(energy_seg, [], 'omitnan') / tr.threshold;
    m.envelope_cv = std(energy_seg, 0, 'omitnan') / mean(energy_seg, 'omitnan');
    if cfg.seizure_robust.hf_band(2) < fs / 2
        raw_seg = data.signal(gs:ge);
        raw_seg = raw_seg(~isnan(raw_seg));
        m.hf_ratio_db = 10 * log10(bandpower(raw_seg, fs, cfg.seizure_robust.hf_band) / ...
            bandpower(raw_seg, fs, cfg.seizure_robust.hf_reference_band));
    else
        m.hf_ratio_db = NaN;
    end
    m.over_max_duration = (ge - gs + 1) / fs > cfg.seizure_robust.max_duration_s;
end

function row = accepted_row(a, p, ev, fragmented)
    m = struct('ll_ratio', a.ll_ratio, 'peak_energy_ratio', a.peak_energy_ratio, ...
        'hf_ratio_db', a.hf_ratio_db, 'envelope_cv', a.envelope_cv, 'over_max_duration', a.over_max_duration);
    row = make_row(a.id, ev, a.start_s, a.end_s, a.duration_s, a.block_id, m, 'accepted', '', fragmented, p.session_start);
end

function row = rescued_row(r, p, cfg, ev)
    m = struct('ll_ratio', r.ll_ratio, 'peak_energy_ratio', r.peak_energy_ratio, ...
        'hf_ratio_db', r.hf_ratio_db, 'envelope_cv', r.envelope_cv, ...
        'over_max_duration', r.duration_s > cfg.seizure_robust.max_duration_s);
    row = make_row(NaN, ev, r.start_s, r.end_s, r.duration_s, r.block_id, m, 'rescued', r.rejected_by{1}, false, p.session_start);
end

function [row, ok] = imputed_row(p, tr, cfg, ev)
    fs = p.data.fs;
    gs = idx_of(ev.ref_s, fs);
    ge = idx_of(ev.ref_e, fs);
    ov = max(0, min(ge, tr.trimmed(:, 2)) - max(gs, tr.trimmed(:, 1)) + 1);
    [best_ov, b] = max(ov);
    ok = ~isempty(best_ov) && best_ov > 0;
    if ok
        gs = max(gs, tr.trimmed(b, 1));
        ge = min(ge, tr.trimmed(b, 2));
        m = window_metrics(tr, p.data, gs, ge, cfg);
        block_id = tr.block_ids(b);
        reason = 'no_candidate';
    else
        gs = min(max(gs, 1), numel(p.data.t_rel));
        ge = min(max(ge, gs), numel(p.data.t_rel));
        m = struct('ll_ratio', NaN, 'peak_energy_ratio', NaN, 'hf_ratio_db', NaN, 'envelope_cv', NaN, ...
            'over_max_duration', false);
        block_id = NaN;
        reason = 'no_valid_signal';
    end
    row = make_row(NaN, ev, p.data.t_rel(gs), p.data.t_rel(ge), (ge - gs + 1) / fs, block_id, m, ...
        'imputed', reason, false, p.session_start);
end

function row = make_row(local_id, ev, start_s, end_s, duration_s, block_id, m, status, rejected_by, fragmented, session_start)
    row = table(local_id, ev.sid, start_s, end_s, duration_s, ...
        abs_or_nat(start_s, session_start), abs_or_nat(end_s, session_start), block_id, ...
        logical(m.over_max_duration), m.ll_ratio, m.peak_energy_ratio, m.hf_ratio_db, m.envelope_cv, ...
        {status}, ev.n_acc, {ev.regions}, ev.n_total, ev.n_acc >= 2, fragmented, ev.ref_s, ev.ref_e, {rejected_by}, ...
        'VariableNames', reconciled_names());
end

function names = reconciled_names()
    names = {'id', 'seizure_id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', 'block_id', ...
        'over_max_duration', 'll_ratio', 'peak_energy_ratio', 'hf_ratio_db', 'envelope_cv', ...
        'detection_status', 'accepted_in_n_channels', 'accepted_in_regions', 'n_channels_total', ...
        'is_bilateral_accepted', 'fragmented', 'ref_start_s', 'ref_end_s', 'rejected_by'};
end

function T = empty_reconciled_table(session_start)
    T = table('Size', [0 22], 'VariableTypes', {'double', 'double', 'double', 'double', 'double', 'datetime', ...
        'datetime', 'double', 'logical', 'double', 'double', 'double', 'double', 'cell', 'double', 'cell', ...
        'double', 'logical', 'logical', 'double', 'double', 'cell'}, 'VariableNames', reconciled_names());
    T.start_abs.TimeZone = session_start.TimeZone;
    T.end_abs.TimeZone = session_start.TimeZone;
end

function t = abs_or_nat(t_rel_s, session_start)
    if isnat(session_start)
        t = NaT;
        t.TimeZone = session_start.TimeZone;
    else
        t = rel_to_abs_time(t_rel_s, session_start);
    end
end

function i = idx_of(t_s, fs)
    i = round(t_s * fs) + 1;  % inverse of t_rel = (i-1)/fs
end

function g = interval_gap(a, b)
% < 0 when the intervals overlap, otherwise the empty time between them.
    g = max(a(1), b(1)) - min(a(2), b(2));
end

function r = uf_find(parent, i)
    while parent(i) ~= i
        i = parent(i);
    end
    r = i;
end

function parent = uf_union(parent, a, b)
    ra = uf_find(parent, a);
    rb = uf_find(parent, b);
    if ra ~= rb
        parent(max(ra, rb)) = min(ra, rb);
    end
end
