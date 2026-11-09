function out = spec_run(cfg, events)
% SPEC_RUN  Entry point of the spectrogram module: validate, resolve every
% requested event, dry-run preview or compute.
%
%   out = spec_run(cfg, events)
%
%   cfg    : spectrogram_config() edited in run_spectrograms.m
%   events : spec_event structs from spec_read_event_list.m (manual table
%            or Excel/CSV) and/or spec_select_seizures.m (by ID / filter),
%            concatenated: [ev_manual, ev_ids, ev_filter]
%
% Steps
%   1. resolve each event x region: files (spec_resolve_files), t_rel of
%      inicio / fin (spec_resolve_times), NaN % (from the txt header /
%      gaps CSV, without reading the signal), source of the spike peaks
%   2. duplicates (same subject, session, region, tipo and interval) are
%      skipped with a warning
%   3. writes <output_root>\06_spectrograms_summary\_selection_preview.csv
%      (one row per event x region, status ok / error / duplicate)
%   4. cfg.spec.dry_run = true -> stops here (nothing computed)
%   5. otherwise each event runs in its own try/catch (an error is logged,
%      the batch goes on): spec_seizure_event / spec_iid_interval, figure,
%      spec_write_outputs -> <output_root>\<subject>\06_spectrograms\
%      <subject>_<session>\<id>_<region>.{mat,fig,png}
%   6. rebuilds spec_features_seizures.csv / spec_features_iid.csv and
%      appends to spec_run_log.txt (same summary folder)
%
% Per-event overrides (events(k).overrides, from the 9th column of the
% manual table or extra columns of the Excel), names case-insensitive:
%   signal_source norm clim cmap freq_scale max_nan_frac
%   seizure: win_s step_s tw k pad_pre_s pad_post_s postictal_s fmin fmax
%            cycles baseline baseline_is_t_rel baseline_time_ref
%   iid:     spike_source polarity align_ms epoch_s fmin fmax cycles baseline
%            baseline_is_t_rel baseline_time_ref threshold_mode
%            threshold_factor threshold_uV upper_uV min_distance_s
%
% out: preview (table), results (struct array: id, region, status, paths,
%      seconds), log_file, summary_dir

    s0 = cfg.spec;
    if isempty(s0.pipeline_root) || ~isfolder(s0.pipeline_root)
        error('spec_run:NoRoot', 'cfg.spec.pipeline_root not found: "%s"', s0.pipeline_root);
    end
    if isempty(s0.output_root)
        s0.output_root = s0.pipeline_root;
    end
    summary_dir = fullfile(s0.output_root, '06_spectrograms_summary');
    if ~isfolder(summary_dir)
        mkdir(summary_dir);
    end
    log_file = fullfile(summary_dir, 'spec_run_log.txt');
    spec_log('open', log_file);
    closer = onCleanup(@() spec_log('close'));
    spec_resolve_files('reset');
    spec_signal_cache('reset');
    git = spec_git_info();
    spec_log('==== spec_run %s   dry_run=%d   %d event(s)   git %s%s', char(datetime('now')), s0.dry_run, ...
        numel(events), git.head, tern(git.dirty, ' (uncommitted changes)', ''));
    spec_log('pipeline_root = %s   output_root = %s', s0.pipeline_root, s0.output_root);

    ctx = struct('seizures', [], 'cache_dir', '');
    if ~isempty(s0.seizures_file)
        [ctx.seizures, notes] = spec_read_seizures_table(s0.seizures_file, s0.tz);
        spec_log('seizures_file %s: %d row(s); %s', s0.seizures_file, height(ctx.seizures), strjoin(notes, '; '));
    end

    % ---- 1. resolve --------------------------------------------------
    rows = {};
    plan = struct('ev', {}, 'cfg', {}, 'R', {}, 'ok', {});
    seen = {};
    for k = 1:numel(events)
        ev = events(k);
        try
            s = apply_overrides(s0, ev.overrides, ev.tipo);
        catch ME
            rows(end+1, :) = preview_row(ev, ev.region, [], [], sprintf('error: %s', ME.message), s0); %#ok<AGROW>
            spec_log('ERROR %s: %s', ev.id, ME.message);
            continue;
        end
        regions = expand_regions(ev);
        R = struct([]);
        for j = 1:numel(regions)
            reg = regions{j};
            try
                r = spec_resolve_files(s.pipeline_root, ev.subject, ev.sesion, reg, s.signal_source, s.tz);
                [r.on, r.off] = event_times(ev, reg, r);
                dur = r.n_samples / r.fs;
                if ~(r.off > r.on)
                    error('spec_run:BadInterval', 'fin (%.3f s) must be after inicio (%.3f s) in t_rel.', r.off, r.on);
                end
                if r.on < 0 || r.off > dur
                    error('spec_run:OutsideRecording', 'interval [%.3f, %.3f] s t_rel is outside the recording (0-%.3f s).', ...
                        r.on, r.off, dur);
                end
                r.excel_attenuation = qc_attenuation(s.pipeline_root, r);
                key = sprintf('%s|%s|%s|%s|%d|%d', r.subject, r.session_key, reg, ev.tipo, round(r.on * 1000), round(r.off * 1000));
                if ismember(key, seen)
                    rows(end+1, :) = preview_row(ev, reg, r, ctx, 'duplicate (skipped)', s); %#ok<AGROW>
                    spec_log('WARNING %s %s: duplicate of an earlier event (same subject/session/region/interval), skipped', ev.id, reg);
                    continue;
                end
                seen{end+1} = key; %#ok<AGROW>
                rows(end+1, :) = preview_row(ev, reg, r, ctx, 'ok', s); %#ok<AGROW>
                R = [R, r]; %#ok<AGROW>
            catch ME
                rows(end+1, :) = preview_row(ev, reg, [], [], sprintf('error: %s', ME.message), s); %#ok<AGROW>
                spec_log('ERROR %s %s: %s', ev.id, reg, ME.message);
            end
        end
        if ~isempty(R)
            plan(end+1) = struct('ev', ev, 'cfg', s, 'R', R, 'ok', true); %#ok<AGROW>
        end
    end
    preview = cell2table(rows, 'VariableNames', preview_names());
    pfile = fullfile(summary_dir, '_selection_preview.csv');
    writetable(preview, pfile);
    spec_log('preview: %d row(s), %d ok -> %s', height(preview), sum(strcmp(preview.status, 'ok')), pfile);
    print_preview(preview);

    out = struct('preview', preview, 'results', struct('id', {}, 'region', {}, 'status', {}, 'paths', {}, 'seconds', {}), ...
        'log_file', log_file, 'summary_dir', summary_dir);
    if s0.dry_run
        spec_log('dry_run = true: nothing computed. Check %s, then set cfg.spec.dry_run = false.', pfile);
        return;
    end

    % ---- 5. compute ------------------------------------------------------
    results = out.results;
    for k = 1:numel(plan)
        ev = plan(k).ev;
        s = plan(k).cfg;
        R = plan(k).R;
        cfgE = struct('spec', s);
        sub = R(1).subject;
        if s.cache_disk
            ctx.cache_dir = fullfile(s.output_root, sub, '06_spectrograms', '_cache');
        else
            ctx.cache_dir = '';
        end
        out_dir = fullfile(s.output_root, sub, '06_spectrograms', sprintf('%s_%s', sub, R(1).session_key));
        if strcmp(ev.tipo, 'seizure')
            jobs = {R};
        else
            jobs = num2cell(R);
        end
        for j = 1:numel(jobs)
            Rj = jobs{j};
            tag = strjoin({Rj.region}, '+');
            t0 = tic;
            try
                if strcmp(ev.tipo, 'seizure')
                    res = spec_seizure_event(ev, Rj, ctx, cfgE);
                    fig = spec_plot_seizure(res);
                else
                    res = spec_iid_interval(ev, Rj, ctx, cfgE);
                    fig = spec_plot_iid(res);
                end
                paths = spec_write_outputs('event', res, fig, out_dir, git);
                close(fig);
                for c = cache_lines(res)
                    spec_log('  %s', c{1});
                end
                w = all_warnings(res);
                for i = 1:numel(w)
                    spec_log('  WARNING %s %s: %s', ev.id, tag, w{i});
                end
                results(end+1) = struct('id', ev.id, 'region', tag, 'status', 'ok', 'paths', paths, 'seconds', toc(t0)); %#ok<AGROW>
                spec_log('OK %s %s (%.1f s) -> %s', ev.id, tag, toc(t0), paths.mat);
            catch ME
                results(end+1) = struct('id', ev.id, 'region', tag, 'status', ['error: ' ME.message], ...
                    'paths', struct('mat', '', 'fig', '', 'png', ''), 'seconds', toc(t0)); %#ok<AGROW>
                spec_log('ERROR %s %s: %s (%s)', ev.id, tag, ME.message, where(ME));
                close(findobj('Type', 'figure', 'Visible', 'off'));
            end
        end
    end
    out.results = results;
    copyfile(pfile, fullfile(summary_dir, sprintf('_selection_preview_%s.csv', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')))));
    T = spec_write_outputs('summary', s0.output_root, summary_dir);
    spec_log('summary: %d seizure row(s), %d IID row(s) in %s', height(T.seizure), height(T.iid), summary_dir);
    n_ok = sum(strcmp({results.status}, 'ok'));
    spec_log('done: %d ok, %d error(s)', n_ok, numel(results) - n_ok);
    spec_signal_cache('reset');
end

%% ======================================================================
function regions = expand_regions(ev)
    if strcmp(ev.tipo, 'seizure') && ~isempty(ev.sz)
        regions = {ev.sz.region};
    elseif strcmpi(ev.region, 'both')
        regions = {'HPCr', 'HPCl'};
    else
        regions = {ev.region};
    end
end

function [on, off] = event_times(ev, reg, r)
    if ~isempty(ev.sz)
        k = find(strcmp({ev.sz.region}, reg), 1);
        on = ev.sz(k).start_s;
        off = ev.sz(k).end_s;
        return;
    end
    on = spec_resolve_times('to_t_rel', ev.inicio, ev.time_ref, r);
    off = spec_resolve_times('to_t_rel', ev.fin, ev.time_ref, r);
end

function att = qc_attenuation(root, r)
% Excel Attenuation of this channel, from the subject's 05_summaries/qc_report.csv if present.
    persistent cache
    if isempty(cache)
        cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end
    att = 'NA (no qc_report.csv)';
    key = [root '|' r.subject];
    if ~isKey(cache, key)
        f = dir(fullfile(root, r.subject, '**', 'qc_report.csv'));
        Q = [];
        if isscalar(f)
            try
                Q = read_pipeline_csv(fullfile(f.folder, f.name), 'qc', r.session_start.TimeZone);
            catch
                Q = [];
            end
        end
        cache(key) = Q;
    end
    Q = cache(key);
    if isempty(Q) || ~all(ismember({'source_file', 'region', 'excel_attenuation'}, Q.Properties.VariableNames))
        return;
    end
    m = strcmpi(Q.source_file, r.source_file) & strcmp(Q.region, r.region);
    if any(m)
        att = char(Q.excel_attenuation{find(m, 1)});
        if isempty(att)
            att = '(empty)';
        end
    end
end

function s = apply_overrides(s, ov, tipo)
    if isempty(ov)
        return;
    end
    names = fieldnames(ov);
    for i = 1:numel(names)
        v = ov.(names{i});
        nm = lower(names{i});
        switch nm
            case {'signal_source', 'norm', 'clim', 'cmap', 'freq_scale', 'max_nan_frac'}
                s.(nm) = v;
            case {'win_s', 'step_s'}
                s.seizure.mt.(nm) = v;
            case 'tw'
                s.seizure.mt.TW = v;
            case 'k'
                s.seizure.mt.K = v;
            case {'pad_pre_s', 'pad_post_s', 'postictal_s'}
                s.seizure.(nm) = v;
            case {'fmin', 'fmax'}
                if strcmp(tipo, 'seizure')
                    s.seizure.mt.(nm) = v;
                    s.seizure.morlet.(nm) = v;
                else
                    s.iid.morlet.(nm) = v;
                end
            case 'cycles'
                if strcmp(tipo, 'seizure')
                    s.seizure.morlet.cycles = v;
                else
                    s.iid.morlet.cycles = v;
                end
            case {'baseline', 'baseline_is_t_rel', 'baseline_time_ref'}
                if strcmp(tipo, 'seizure')
                    s.seizure.(nm) = v;
                else
                    s.iid.(nm) = v;
                end
            case {'spike_source', 'polarity', 'align_ms', 'epoch_s'}
                s.iid.(nm) = v;
            case {'threshold_mode', 'threshold_factor', 'threshold_uv', 'upper_uv', 'min_distance_s'}
                fld = strrep(strrep(nm, 'threshold_uv', 'threshold_uV'), 'upper_uv', 'upper_uV');
                s.iid.detect.(fld) = v;
            otherwise
                error('spec_run:BadOverride', 'Unknown per-event override ''%s''.', names{i});
        end
    end
end

function names = preview_names()
    names = {'id', 'tipo', 'subject', 'sesion', 'session_key', 'region', 'status', 'origin', 'txt', ...
        'start_t_rel', 'end_t_rel', 'duration_s', 'start_natus', 'end_natus', 'start_pipeline', 'end_pipeline', ...
        'nan_pct', 'peaks_source', 'seizure_id', 'user_check', 'signal_source'};
end

function row = preview_row(ev, reg, r, ~, status, s)
    sid = NaN; uc = '';
    if ~isempty(ev.sz)
        k = find(strcmp({ev.sz.region}, reg), 1);
        if ~isempty(k)
            sid = ev.sz(k).seizure_id;
            uc = ev.sz(k).user_check;
        end
    end
    if isempty(r)
        row = {ev.id, ev.tipo, ev.subject, ev.sesion, '', reg, status, ev.origin, '', NaN, NaN, NaN, '', '', '', '', ...
            NaN, '', sid, uc, s.signal_source};
        return;
    end
    [n1, p1] = spec_resolve_times('clocks', r.on, r);
    [n2, p2] = spec_resolve_times('clocks', r.off, r);
    peaks = '';
    if strcmp(ev.tipo, 'iid')
        switch s.iid.spike_source
            case 'pipeline'
                if isempty(r.iid_mat)
                    peaks = 'pipeline: MISSING IID_results.mat (event will fail)';
                else
                    peaks = ['pipeline: ' r.iid_mat];
                end
            case 'manual'
                peaks = sprintf('manual (%d)', numel(ev.peaks));
            otherwise
                peaks = sprintf('detect_local (%s)', s.iid.detect.threshold_mode);
        end
    end
    row = {ev.id, ev.tipo, ev.subject, ev.sesion, r.session_key, reg, status, ev.origin, r.txt, r.on, r.off, ...
        r.off - r.on, spec_resolve_times('text', n1), spec_resolve_times('text', n2), ...
        spec_resolve_times('text', p1), spec_resolve_times('text', p2), nan_pct(r, r.on, r.off), peaks, sid, uc, s.signal_source};
end

function p = nan_pct(r, a, b)
% NaN % of [a, b] without reading the signal: 0 when the header has no
% gaps, else from the session's *_gaps.csv; NaN (unknown) if not found.
    p = NaN;
    ng = NaN;
    if isfield(r.hdr.meta, 'n_gaps')
        ng = str2double(r.hdr.meta.n_gaps);
    end
    if ng == 0
        p = 0;
        return;
    end
    g = dir(fullfile(fileparts(r.txt), sprintf('%s_%s_gaps.csv', r.subject, r.session_key)));
    if isempty(g)
        g = dir(fullfile(fileparts(fileparts(r.txt)), '**', sprintf('%s_%s_gaps.csv', r.subject, r.session_key)));
    end
    if isempty(g)
        return;
    end
    try
        G = readtable(fullfile(g(1).folder, g(1).name));
        ov = max(0, min(G.end_s, b) - max(G.start_s, a));
        p = 100 * sum(ov) / (b - a);
    catch
    end
end

function print_preview(P)
    fprintf('\n%-10s %-7s %-7s %-16s %-5s %-24s %-24s %8s  %s\n', 'id', 'tipo', 'subject', 'session', 'reg', ...
        'start (Natus)', 'start (pipeline)', 'dur_s', 'status');
    for i = 1:height(P)
        fprintf('%-10s %-7s %-7s %-16s %-5s %-24s %-24s %8.1f  %s\n', P.id{i}, P.tipo{i}, P.subject{i}, ...
            P.session_key{i}, P.region{i}, P.start_natus{i}, P.start_pipeline{i}, P.duration_s(i), P.status{i});
    end
    fprintf('\n');
end

function lines = cache_lines(res)
    if isfield(res, 'ch')
        srcs = arrayfun(@(c) c.src, res.ch);
    else
        srcs = res.src;
    end
    lines = {};
    for i = 1:numel(srcs)
        if srcs(i).cache_reused
            lines{end+1} = sprintf('signal cache reused for %s', srcs(i).txt); %#ok<AGROW>
        else
            lines{end+1} = sprintf('signal cache built for %s in %.1f s', srcs(i).txt, srcs(i).cache_build_s); %#ok<AGROW>
        end
    end
end

function w = all_warnings(res)
    if isfield(res, 'ch')
        w = unique([res.ch.warnings], 'stable');
    else
        w = res.warnings;
    end
end

function s = where(ME)
    s = '';
    if ~isempty(ME.stack)
        s = sprintf('%s line %d', ME.stack(1).name, ME.stack(1).line);
    end
end

function v = tern(c, a, b)
    if c
        v = a;
    else
        v = b;
    end
end
