%% CONSOLIDATE_RESULTS  One folder with the validated seizure and IID results
% of campaigns C1/C2 (2026) for the 28 validated mice, from EXISTING outputs
% and the user's review. Extraction, cleaning and seizure detection are never
% re-run. Reuses run_campaign_stats.m: its CONFIG block is evaluated as is
% (exclusions, merges, missed seizures, disconnections, IID zones), then the
% same src/stats functions build the seizure set, the observation time and
% the IID re-detection. Decisions of 2026-11-02 that differ from it:
%   - 077-s (C2): labels from its CURRENT seizures_events.csv, not the
%     20261019 copy that run_campaign_stats.m reads (line 54)
%   - 128-s 20260729 reprocessed (A1 / A2 EDFs): included
%   - crisette_focal events: counted as seizures (source pipeline_TP)
%   - IIDs: run_campaign_stats' re-detection; channels without an earlier
%     result with the same zones and code are re-detected (detect_iid only)
% Seizure figures and per-seizure energy metrics: the energy and line-length
% traces are recomputed from 02_clean with the pipeline's own functions
% (seizure_energy_trace, line_length_trace); no detection is re-run.
% Seizure ids (2026-11-03): 1..N per EDF in onset order over the final
% seizures (merges applied); the review numbering is kept as original_id.
%
% Output: <batch 2026>/consolidated_<yyyyMMdd>/  seizures_all.xlsx,
% seizures_summary.xlsx, iid_summary.xlsx, figures/<mouse>/*.fig + .png,
% README.md; _work/ = state, IID runs, worker logs,
% stats/labels_snapshot/<mouse>_seizures_validated.xlsx (the final seizures
% of each mouse) and stats/labels_snapshot_raw/ (the label files as read).
% IID re-detection and figures run in background MATLAB processes (no
% Parallel Toolbox): this same script, with consolidate_role set.
%
% Run:  matlab -batch "run('<repo>/consolidate_results.m')"
% Optional before run(): consolidate_n_workers (per kind, default 3),
% consolidate_date ('yyyyMMdd', default today), consolidate_reuse = true
% (tables and figures again from _work/state.mat: no label reading, and the
% IID tables from _work/iid_summary_source.xlsx -- the iid_summary.xlsx of
% the full run -- so nothing is re-detected; 2026-11-03: the v1 run folder
% that held most reused IID results no longer exists), consolidate_figures =
% false (keep the figures and _work/fig_jobs of an earlier run).
%
% Incremental update (2026-11-04, consolidate_update = 'compute' then
% 'write'): the mice of upd.mice are (re)consolidated into the existing
% folder, every other mouse keeps its rows (upd_hist: one update after another). 'compute': phase A
% over all mice with a fresh label snapshot (_work/update_<tag>/stats), check
% that nothing changed for the other mice (seizures, per-EDF / per-mouse
% rows, IID zones per channel, vs _work/state.mat), IID re-detection of the
% update mice only (with the 055-s AF cases, af_rule_055s) and their new
% seizure figures (into _work/update_<tag>); stops for the user's choice of
% the 055-s case. 'write': with af_rule_055s set, the sheets (rows of the
% other mice checked cell by cell against the files before), the figures, the
% per-mouse validated Excel of the update mice, README.md and _work/state.mat
% (previous files kept under _work/update_<tag>/before). Then run
% consolidate_4OHT.m for its sheets.

if ~exist('consolidate_role', 'var'), consolidate_role = 'main'; end
if ~exist('consolidate_n_workers', 'var'), consolidate_n_workers = 3; end
if ~exist('consolidate_reuse', 'var'), consolidate_reuse = false; end
if ~exist('consolidate_figures', 'var'), consolidate_figures = true; end
if ~exist('consolidate_update', 'var'), consolidate_update = ''; end
repo = fileparts(mfilename('fullpath'));
if ~exist('consolidate_out', 'var')
    if ~exist('consolidate_date', 'var'), consolidate_date = char(datetime('now', 'Format', 'yyyyMMdd')); end
    consolidate_out = fullfile(fileparts(repo), ['consolidated_' consolidate_date]);
end
work = fullfile(consolidate_out, '_work');
state_file = fullfile(work, 'state.mat');
% incremental updates, oldest first: {tag, mice validated since the previous consolidation or whose review
% changed, what changed ('' = new mice)}; consolidate_update runs the last one, against the consolidation as
% the earlier ones left it
upd_hist = {'20261104',  {'016-s', '055-s', '084-s', '113-s'}, ''; ...
            '20261104b', {'025-s'}, ''; ...
            '20261105',  {'055-s', '070-s', '072-s', '076-s', '083-s', '112-s', '113-s', '114-s', '117-s', '124-s'}, ...
                'desconexiones revisadas en video y ND-12 (072-s)'};
upd.tag = upd_hist{end, 1};
upd.mice = upd_hist{end, 2};
upd.desc = upd_hist{end, 3};
upd.dir = fullfile(work, ['update_' upd.tag]);
upd.earlier = fullfile(work, strcat('update_', flip(upd_hist(1:end-1, 1))'));   % newest first
is_upd = any(strcmp(consolidate_update, {'compute', 'write'}));
is_iidh = strcmp(consolidate_update, 'iid_hemisphere');   % one-off: IID per_mouse by hemisphere (2026-11-04)
assert(isempty(consolidate_update) || is_upd || is_iidh, 'consolidate_update: ''compute'', ''write'' or ''iid_hemisphere''');
if is_upd
    state_file = fullfile(upd.dir, 'state.mat');
    wd = fullfile(upd.dir, 'workers');
    fig_out = upd.dir; fig_work = upd.dir;   % new figures and their jobs: copied into the consolidation by 'write'
else
    wd = fullfile(work, 'workers');
    fig_out = consolidate_out; fig_work = work;
end
% 055-s (damaged signal, 2026-11-04): FP whose AF differs between hemispheres. The general rule
% (every other mouse) excludes them in both hemispheres (R4 / R3 of the event's AF union); for
% 055-s: 1 = HPCl's AF for both hemispheres, 2 = each hemisphere its own AF, 3 = every FP as
% "other" (excluded). NaN = not decided yet ('compute' runs the three; 'write' stops).
af_rule_055s = 3;   % the user's choice, 2026-11-04

%% ============================ WORKERS ====================================
if ~strcmp(consolidate_role, 'main')
    restoredefaultpath;
    addpath(genpath(fullfile(repo, 'src')));
    k = consolidate_worker(1); n = consolidate_worker(2);
    msg = 'ok';
    try
        switch consolidate_role
            case 'iid_worker'   % every channel set (the general rule, the 055-s cases) in turn
                q = load(state_file, 'IV');
                for v = q.IV
                    cs_iid_rerun_all(v.cfg, v.Z, k, n);
                end
            case 'fig_worker'
                q = load(state_file, 'cfg', 'SZfig', 'edf_tbl', 'det_rows');
                q.SZ = q.SZfig;
                fig_worker(q, k, n, fig_out, fig_work);
        end
    catch ME
        msg = getReport(ME, 'extended', 'hyperlinks', 'off');
        fprintf(2, '%s\n', msg);
    end
    fid = fopen(fullfile(wd, sprintf('%s_%d_of_%d.done', consolidate_role, k, n)), 'w');
    fprintf(fid, '%s', msg); fclose(fid);
    return;
end

expected.control = {'025-s', '051-s', '056-s', '070-s', '083-s', '109-s'};
controls_added = {'025-s'};   % controls after run_campaign_stats.m's cfg.controls (the user, 2026-11-04)
expected.epileptic = {'001-s', '009-s', '016-s', '020-s', '036-s', '048-s', '055-s', '058-s', '062-s', '072-s', '073-s', ...
    '076-s', '077-s', '084-s', '085-s', '088-s', '091-s', '097-s', '110-s', '112-s', '113-s', '114-s', '117-s', '124-s', '128-s'};
expected.excluded = {'048-s|20260622~ TRAP_fd13b728-0df8-45ad-8f8c-2b79a894d3ef'};   % 084-s 20260706 / 20260715A validated (2026-11-04)
% missed seizures added after run_campaign_stats.m's cfg.nd_source (2026-11-04, ND-12 2026-11-05). Seen on
% video: ND-08/09 during DX-23 (055-s), after the EEG was cut at that disconnection; ND-12 during DX-07
% (072-s, signal = noise). cs_build_seizures expects ND-05/06 as its only video seizures, so these enter as
% missed_eeg and become video below (nd_beyond_eeg: with the EDF of their row; nd_video_in_eeg)
nd_added = {'ND-08', 'missed_eeg'; 'ND-09', 'missed_eeg'; 'ND-10', 'missed_eeg'; 'ND-11', 'missed_eeg'; 'ND-12', 'missed_eeg'};
nd_video_no_eeg = {'ND-08', 'ND-09'};
nd_video_eeg = {'ND-12'};
% disconnections reviewed on video after the end of the recorded EEG whose time counts as observed for
% seizures (the user, 2026-11-05: 055-s on video until 04/06/2026 20:31:39 = DX-23; the rest of that
% disconnection, DX-26, not reviewed). cs_observation_time counts a reviewed DX only inside the analyzed
% EDFs: dx_video_beyond_eeg below. DX-21 (128-s, between the _128s_A1 and _A2 EDFs) is reviewed but does
% not count (the user, 2026-11-05)
dx_video_no_eeg = {'DX-23'};
% notes added to a seizure's comments: {mouse, original_id, text} (2026-11-04)
sz_notes = {'113-s', 'ND-10', ['onset 19 s before the start of the 4-OHT window (09/08/2026 17:27:00 = injection - 2 h), ' ...
    'offset inside it: not in seizures_4-OHT (the onset decides)']};

if is_iidh
    %% ============================ IID per_mouse BY HEMISPHERE ============
    restoredefaultpath;
    addpath(genpath(fullfile(repo, 'src')));
    iid_by_hemisphere(consolidate_out, work, fullfile(work, 'iid_hemisphere_20261104'), af_rule_055s);
    return;
end

if is_upd && strcmp(consolidate_update, 'write')
    %% ============================ UPDATE: STATE OF 'compute' =============
    restoredefaultpath;
    addpath(genpath(fullfile(repo, 'src')));
    q = load(state_file);
    for f = {'cfg', 'S', 'O', 'Z', 'rv', 'edf_tbl', 'ch_tbl', 'det_rows', 'snapinfo', 'IV', 'CV', 'SZfig', 'CHK'}
        assert(isfield(q, f{1}), 'consolidate_results: %s lacks %s -- run consolidate_update = ''compute'' first', state_file, f{1});
    end
    cfg = q.cfg; S = q.S; O = q.O; Z = q.Z; rv = q.rv; edf_tbl = q.edf_tbl; ch_tbl = q.ch_tbl; det_rows = q.det_rows;
    snapinfo = q.snapinfo; IV = q.IV; CV = q.CV; SZfig = q.SZfig; CHK = q.CHK;
    clear q
elseif consolidate_reuse
    %% ============================ REUSE _work ============================
    restoredefaultpath;
    addpath(genpath(fullfile(repo, 'src')));
    q = load(state_file, 'cfg', 'S', 'O', 'Z', 'rv', 'edf_tbl', 'ch_tbl', 'det_rows');
    cfg = q.cfg; S = q.S; O = q.O; Z = q.Z; rv = q.rv; edf_tbl = q.edf_tbl; ch_tbl = q.ch_tbl; det_rows = q.det_rows;
    if ~isempty(whos('-file', state_file, 'iid_C'))   % kept in the state rewritten below (consolidate_4OHT.m)
        q = load(state_file, 'iid_C'); iid_C = q.iid_C;
    end
    clear q
    raw = fullfile(cfg.out_dir, 'labels_snapshot_raw');
    if ~isfolder(raw)   % the per-campaign snapshot this state was built from
        movefile(cfg.labels_snapshot_dir, raw);
    end
    cfg.labels_snapshot_dir = raw;
    M = readtable(fullfile(raw, 'manifest.csv'), 'TextType', 'char', 'Delimiter', ',');
    snapinfo = struct('taken', M.snapshot_taken{1}, 'n', height(M), 'dir', rel_name(raw, consolidate_out));
    iid_src = fullfile(work, 'iid_summary_source.xlsx');   % the IID tables of the full run (2026-11-02)
    if ~isfile(iid_src), copyfile(fullfile(consolidate_out, 'iid_summary.xlsx'), iid_src); end
else
    %% ============================ CONFIG =================================
    % run_campaign_stats.m's CONFIG block, unchanged, with its own overrides
    if ~isfolder(work), mkdir(work); end
    out_dir_override = fullfile(work, 'stats');
    if is_upd, out_dir_override = fullfile(upd.dir, 'stats'); end   % the earlier snapshots and Excel stay as they are
    include_128s_0729_reprocessed_override = true;
    run_threshold_diag_override = false;
    rcs = fileread(fullfile(repo, 'run_campaign_stats.m'));
    k = strfind(rcs, '%% ============================ PHASE A');
    assert(isscalar(k), 'run_campaign_stats.m: PHASE A marker not found once');
    eval(rcs(1:k - 1));
    assert(strcmp(cfg.repo_root, repo), 'cfg.repo_root %s ~= %s', cfg.repo_root, repo);
    cfg.label_file_override = cell(0, 3);   % 077-s: its current seizures_events.csv
    cfg.labels_snapshot_refresh = true;     % one fresh snapshot per run of this script
    for v = {'principal', 'sin_crisette'}   % IID results of the later official runs, reused when zones + code match
        sub = struct('principal', 'iid_rerun_final', 'sin_crisette', 'iid_rerun_final_sin_crisette');
        for d = {'_stats_revision_20261101_final', '_stats_revision_20261030_final'}
            p = fullfile(cfg.batch_root, d{1}, sub.(v{1}));
            if isfolder(p), cfg.iid_rerun_dirs.(v{1}){end+1} = p; end
        end
        if is_upd   % the consolidation's own IID results stay its output root
            cfg.iid_rerun_dirs.(v{1}) = [{fullfile(work, 'stats', sub.(v{1}))}, cfg.iid_rerun_dirs.(v{1})];
        end
    end
    assert(~any(ismember(nd_added(:, 1), cfg.nd_source(:, 1))), 'nd_added: already in run_campaign_stats.m''s cfg.nd_source');
    cfg.nd_source = [cfg.nd_source; nd_added];
    assert(~any(ismember(controls_added, cfg.controls)), 'controls_added: already in run_campaign_stats.m''s cfg.controls');
    cfg.controls = [cfg.controls, controls_added];

    %% ============================ SEIZURES ===============================
    % run_campaign_stats.m phase A, through the IID zones
    cache_src = {cfg.prev_dir, cfg.v1_dir};
    if is_upd, cache_src = [fullfile(upd.earlier, 'stats'), {fullfile(work, 'stats')}, cache_src]; end
    for d = cache_src
        if ~isfolder(fullfile(cfg.out_dir, 'cache')) && isfolder(fullfile(d{1}, 'cache'))
            copyfile(fullfile(d{1}, 'cache'), fullfile(cfg.out_dir, 'cache'));
        end
    end
    [snap, cfg.labels_snapshot_dir] = cs_label_snapshot(cfg);
    raw = fullfile(cfg.out_dir, 'labels_snapshot_raw');   % labels_snapshot/ holds the per-mouse validated Excel
    if isfolder(raw), rmdir(raw, 's'); end
    movefile(cfg.labels_snapshot_dir, raw);
    cfg.labels_snapshot_dir = raw;
    snap = rp128_integrate('snapshot', cfg, snap);
    snapinfo = struct('taken', snap.snapshot_taken{1}, 'n', height(snap), 'dir', rel_name(raw, consolidate_out));
    [edf_tbl, ch_tbl] = cs_load_edf_table(cfg);
    [det_rows, review_rows, det_log] = cs_load_detections(cfg, edf_tbl);
    [edf_tbl, ch_tbl] = rp128_integrate('tables', cfg, edf_tbl, ch_tbl);
    [det_rows, review_rows, rp_log] = rp128_integrate('detections', cfg, edf_tbl, det_rows, review_rows);
    fprintf('%s\n', det_log{:}, rp_log{:});
    events = cs_build_events(det_rows);
    rv = cs_load_review(cfg, edf_tbl);
    [rv.excluded, f_excl] = cs_excluded_edfs(cfg, rv, edf_tbl, events);
    edf_lists = table();
    for c = 1:numel(cfg.campaigns)
        L = readtable(fullfile(cfg.campaigns(c).root, 'edf_list.csv'), 'TextType', 'char', 'Delimiter', ',');
        L.campaign = repmat({cfg.campaigns(c).name}, height(L), 1);
        L.reason = repmat({'error en el analisis (ver log de la campaign)'}, height(L), 1);
        sk = fullfile(cfg.campaigns(c).root, 'skipped_inputs.csv');
        if isfile(sk)
            K = readtable(sk, 'TextType', 'char', 'Delimiter', ',');
            [tf, j] = ismember(strcat(L.animal_id, '|', L.edf_file), strcat(K.animal_id, '|', K.edf_file));
            L.reason(tf) = K.reason(j(tf));
        end
        edf_lists = [edf_lists; table(L.campaign, L.animal_id, L.edf_file, L.reason, ...
            'VariableNames', {'campaign', 'mouse', 'edf', 'reason'})]; %#ok<AGROW>
    end
    edf_lists = rp128_integrate('edf_lists', cfg, edf_lists);
    S = cs_build_seizures(cfg, events, rv, edf_tbl, edf_lists);
    S = nd_beyond_eeg(S, rv, edf_tbl, nd_video_no_eeg);
    S = nd_video_in_eeg(S, rv, nd_video_eeg);
    O = cs_observation_time(cfg, S, rv, edf_tbl);
    O = dx_video_beyond_eeg(O, cfg, rv, edf_tbl, dx_video_no_eeg);
    if ~cfg.git.iid_code_unmodified
        error('consolidate_results:IidCode', 'detect_iid / load_lfp_txt / pipeline_config (or their helpers) differ from HEAD');
    end
    Z = cs_iid_zones(cfg, S, rv, ch_tbl, det_rows);
end

% scope = the user's list (stop otherwise)
inc = S.units(S.units.included, :);
got.control = unique(inc.mouse(strcmp(inc.group, 'control')))';
got.epileptic = unique(inc.mouse(strcmp(inc.group, 'epileptic')))';
got.excluded = sort(strcat(rv.excluded.mouse, '|', rv.excluded.edf))';
for f = {'control', 'epileptic', 'excluded'}
    if ~isequal(sort(got.(f{1})), sort(expected.(f{1})))
        error('consolidate_results:Scope', '%s: got {%s}, expected {%s} -- not consolidated, ask the user', f{1}, ...
            strjoin(setxor(got.(f{1}), expected.(f{1})), ', '), strjoin(expected.(f{1}), ', '));
    end
end
SZ = seizure_table(cfg, S, rv);
for k = 1:size(sz_notes, 1)
    j = find(strcmp(SZ.mouse, sz_notes{k, 1}) & strcmp(SZ.original_id, sz_notes{k, 2}));
    assert(isscalar(j), 'sz_notes: %s %s not found once', sz_notes{k, 1}, sz_notes{k, 2});
    SZ.comments{j} = strjoin([SZ.comments(j), sz_notes(k, 3)], ' | ');
end
% a figure needs EEG at the onset (not the video seizures after the end of 055-s's cut recording)
SZ.has_eeg = arrayfun(@(i) any(strcmp(edf_tbl.mouse, SZ.mouse{i}) & strcmp(edf_tbl.edf, SZ.edf{i}) & ...
    edf_tbl.start_posix <= SZ.start_posix(i) & edf_tbl.end_posix >= SZ.start_posix(i)), (1:height(SZ))');
fprintf('seizures: %d (%s)\n', height(SZ), strjoin(strcat(unique(SZ.source), ':', ...
    arrayfun(@(s) num2str(sum(strcmp(SZ.source, s))), unique(SZ.source), 'UniformOutput', false))', ', '));
if exist('consolidate_stop_after', 'var') && strcmp(consolidate_stop_after, 'seizures')   % check run: state + seizure table only
    save(state_file, 'cfg', 'S', 'O', 'Z', 'rv', 'edf_tbl', 'ch_tbl', 'det_rows', 'SZ', 'snapinfo', '-v7.3');
    return;
end
upd_write = is_upd && strcmp(consolidate_update, 'write');

%% ============================ WORKERS: IID + FIGURES =====================
if ~upd_write
    roles = {};
    SZfig = SZ(SZ.has_eeg, :);
    if is_upd
        % every other mouse as consolidated before; figures only for the update mice's new seizures and the
        % earlier ones renumbered by them (the figure title carries the seizure id)
        CHK = update_checks(load(fullfile(work, 'state.mat'), 'SZ', 'S', 'O', 'Z'), SZ, S, O, Z, upd.mice, cfg.iid_excl_buffer_s);
        SZfig = SZfig(ismember(SZfig.figure, [CHK.renumbered.figure_after; CHK.new_figures]) | ...
            ~cellfun(@(f) isfile(fullfile(consolidate_out, 'figures', [f '.png'])), SZfig.figure), :);
        assert(all(ismember(SZfig.mouse, upd.mice)), 'consolidate_results: a figure to draw outside the update mice');
        roles = {'fig_worker'};
    elseif consolidate_figures
        roles = {'fig_worker'};
    end
    if isempty(SZfig), roles = {}; end
    if ~consolidate_reuse
        [RT, RTcfg] = cs_iid_repro_test(cfg, ch_tbl, cfg.provenance);
        if ~all(RT.reproduce)
            error('consolidate_results:Repro', 'the detect_iid reproducibility test failed: no re-detection');
        end
        if is_upd   % the update mice only; 055-s: its AF case, or the general rule and the three cases if not decided
            rules = af_rule_055s;
            if ~(isscalar(rules) && any(rules == 1:3)), rules = 0:3; end
            IV = iid_sets(cfg, S, rv, ch_tbl, det_rows, Z, upd.mice, fullfile(upd.dir, 'iid_055s'), rules);
        else
            IV = iid_sets(cfg, S, rv, ch_tbl, det_rows, Z, unique(Z.channels.mouse)', fullfile(cfg.out_dir, 'iid_055s'), af_rule_055s);
        end
        roles = [{'iid_worker'}, roles];
    end
    vars = {'cfg', 'S', 'O', 'Z', 'rv', 'edf_tbl', 'ch_tbl', 'det_rows', 'SZ', 'snapinfo', 'SZfig'};
    if ~consolidate_reuse, vars{end+1} = 'IV'; end
    if consolidate_reuse && exist('iid_C', 'var'), vars{end+1} = 'iid_C'; end
    if is_upd, vars{end+1} = 'CHK'; end
    if ~isfolder(fileparts(state_file)), mkdir(fileparts(state_file)); end
    save(state_file, vars{:}, '-v7.3');
    if ~isfolder(wd), mkdir(wd); end
    delete(fullfile(wd, '*.done'));
    if (consolidate_figures || is_upd) && any(strcmp(roles, 'fig_worker'))   % figures are redrawn from scratch (names follow the seizure ids)
        for d = {fullfile(fig_out, 'figures'), fullfile(fig_work, 'fig_jobs')}
            if isfolder(d{1}), rmdir(d{1}, 's'); end
        end
    end
    exe = fullfile(matlabroot, 'bin', 'matlab.exe');
    upd_arg = '';
    if is_upd, upd_arg = sprintf('consolidate_update=''%s''; ', consolidate_update); end
    for r = roles
        for kw = 1:consolidate_n_workers
            cmd = sprintf(['consolidate_role=''%s''; consolidate_worker=[%d %d]; consolidate_out=''%s''; %s' ...
                'run(''%s'')'], r{1}, kw, consolidate_n_workers, consolidate_out, upd_arg, fullfile(repo, 'consolidate_results.m'));
            system(sprintf('start "%s %d" /min "%s" -batch "%s" -logfile "%s"', r{1}, kw, exe, cmd, ...
                fullfile(wd, sprintf('%s_%d.log', r{1}, kw))));
        end
    end
    n_done = 0; t_wait = tic;
    while n_done < numel(roles) * consolidate_n_workers
        pause(60);
        n_done = numel(dir(fullfile(wd, '*.done')));
        fprintf('[%s] workers done %d/%d, figure EDFs done %d\n', char(datetime('now', 'Format', 'HH:mm')), n_done, ...
            numel(roles) * consolidate_n_workers, numel(dir(fullfile(fig_work, 'fig_jobs', '*.mat'))));
        if toc(t_wait) > 24 * 3600, error('consolidate_results:Timeout', 'workers not finished after 24 h (%s)', wd); end
    end
    worker_errors = {};
    for f = dir(fullfile(wd, '*.done'))'
        m = fileread(fullfile(wd, f.name));
        if ~strcmp(m, 'ok'), worker_errors{end+1} = sprintf('%s: %s', f.name, m); end %#ok<SAGROW>
    end
    if ~consolidate_reuse   % collects every set; runs anything a worker left undone
        CV = struct();
        for v = IV
            CV.(v.name) = cs_iid_rerun_all(v.cfg, v.Z, 0);
            c = CV.(v.name)(strcmp(CV.(v.name).variant, 'principal'), :);
            fprintf('IID set %s: %d channel(s) (%d reused from an earlier run or another set, %d run now)\n', v.name, height(c), ...
                sum(strcmp(c.status, 'reutilizado de una corrida anterior')), sum(ismember(c.status, {'reutilizado', 'corrido'})));
        end
    end
    if is_upd   % 'compute' ends here: the user chooses the 055-s case (af_rule_055s), then 'write'
        save(state_file, 'CV', 'worker_errors', '-append');
        fprintf('\n==== update %s computed: %s\n', upd.tag, upd.dir);
        fprintf('  other mice unchanged: %d seizures, %d per-EDF rows, %d IID channels (zones)\n', CHK.n_seizures, CHK.n_edf_rows, CHK.n_iid_channels);
        for i = 1:numel(CHK.changed), fprintf('  earlier seizure of an update mouse changed: %s\n', CHK.changed{i}); end
        for i = 1:height(CHK.renumbered)
            fprintf('  renumbered: %s %s, seizure_id %d -> %d\n', CHK.renumbered.mouse{i}, CHK.renumbered.original_id{i}, ...
                CHK.renumbered.id_before(i), CHK.renumbered.id_after(i));
        end
        for i = 1:height(CHK.new), fprintf('  %s: %d new seizure(s)\n', CHK.new.mouse{i}, CHK.new.n(i)); end
        fprintf('  figures drawn: %d (%d EDF jobs)\n', height(SZfig), numel(dir(fullfile(fig_work, 'fig_jobs', '*.mat'))));
        for i = 1:numel(worker_errors), fprintf(2, '  worker error: %s\n', worker_errors{i}); end
        return;
    end
else
    q = load(state_file, 'worker_errors');
    worker_errors = q.worker_errors;
end
if upd_write   % this run's IIDs for the update mice, every other mouse's rows as consolidated before
    C = iid_final(CV, af_rule_055s);
    I = iid_channel_table(C, S);
    K = cs_iid_channel_intervals(C, S);
    [IE4, IM4] = iid_summary(cfg, I, K, C);
    [IEo, IMo] = iid_from_xlsx(fullfile(consolidate_out, 'iid_summary.xlsx'));
    [IE, IM, n_iid_same] = iid_merge(Z, upd.mice, IE4, IM4, IEo, IMo);
    pr = strcmp(C.variant, 'principal');
    iid_note = sprintf(['IIDs: %s (última actualización) re-detectados con sus zonas actuales (%d canales: %d re-detectados, %d ' ...
        'reutilizados de corridas anteriores con las mismas zonas y código; los %d ya consolidados dan los mismos valores que antes). ' ...
        'El resto de los ratones conserva los valores de la consolidación anterior (%d canales; ver las secciones de actualización).'], ...
        strjoin(upd.mice, ', '), sum(pr), sum(ismember(C.status, {'reutilizado', 'corrido'}) & pr), ...
        sum(strcmp(C.status, 'reutilizado de una corrida anterior') & pr), n_iid_same, height(IE) - sum(pr));
elseif consolidate_reuse   % no IID re-detection at all: the tables of the full run
    [IE, IM] = iid_from_xlsx(iid_src);
    iid_note = sprintf(['IIDs: los mismos valores de la consolidación completa (copia de su iid_summary.xlsx en ' ...
        '`_work/iid_summary_source.xlsx`, %d canales). En esta actualización no se re-detectó nada; los resultados por canal ' ...
        'reutilizados de la corrida v1 (`_stats_revision_20261029`) ya no están en disco.'], height(IE));
else
    C = iid_final(CV, af_rule_055s);
    fprintf('IID re-detection: %d runs (%d reused from an earlier run, %d reused, %d run now)\n', height(C), ...
        sum(strcmp(C.status, 'reutilizado de una corrida anterior')), sum(strcmp(C.status, 'reutilizado')), sum(strcmp(C.status, 'corrido')));
    Corig = rp128_integrate('iid_original', cfg, Z.channels(Z.channels.run, :));
    R = cs_seizure_stats(cfg, S, O);
    I = cs_iid_stats(cfg, S, R, Z, C, Corig);
    K = cs_iid_channel_intervals(C, S);
    [IE, IM] = iid_summary(cfg, I, K, C);
    pr = strcmp(C.variant, 'principal');
    iid_note = sprintf(['Re-detección de IIDs: %d canales (%d reutilizados de corridas anteriores con las mismas zonas y código, ' ...
        '%d re-detectados en esta consolidación).'], sum(pr), sum(strcmp(C.status, 'reutilizado de una corrida anterior') & pr), ...
        sum(ismember(C.status, {'reutilizado', 'corrido'}) & pr));
end
if upd_write   % the new figures and their jobs into the consolidation (renumbered seizures: replaced, old ones kept in before/)
    copy_new_figures(SZfig, upd, consolidate_out, work, SZ, CHK);
end
FJ = fig_results(SZ, work, consolidate_out);
FJ.note(~SZ.has_eeg) = {'no EEG at the onset (video-only seizure after the end of the recorded EEG)'};
mets = {'peak_energy_ratio', 'mean_to_thr', 'auc_to_thr_s'};
for m = 1:3
    for c = 1:2
        SZ.([mets{m} '_' FJ.regions{c}]) = FJ.(mets{m})(:, c);
    end
end

%% ============================ TABLES =====================================
[SE, SM, n_out_obs] = seizure_summary(S, O, SZ);
CMP = compare_20261101(cfg, SM, IM);
out = consolidate_out;
xs = {'seizures_all.xlsx', 'seizures_summary.xlsx', 'iid_summary.xlsx'};
cols = [{'mouse', 'group', 'edf_file', 'seizure_id', 'original_id', 'event_id', 'hemispheres', 'onset', 'offset', 'duration_s'}, ...
    strcat(reshape(repmat(mets, 2, 1), 1, []), '_', repmat(FJ.regions, 1, 3)), {'source', 'light_phase', 'comments'}];
tabs = {{SZ(:, cols)}, {SE, SM}, {IE, IM}};
shs = {{'seizures'}, {'per_EDF', 'per_mouse'}, {'per_mouse_EDF_hemisphere', 'per_mouse'}};
vd = fullfile(work, 'stats', 'labels_snapshot');   % one Excel per mouse = its rows of seizures_all.xlsx
extra = readme_notes(cfg, S, SZ, af_rule_055s, rv);
if upd_write
    % into copies first; the files are replaced only if every other mouse's rows and every other sheet are unchanged
    bk = fullfile(upd.dir, 'before');
    if ~isfolder(bk), mkdir(bk); end
    keep_files = [fullfile(out, [xs, {'README.md'}]), fullfile(work, {'state.mat', 'iid_summary_source.xlsx'}), ...
        fullfile(vd, strcat(upd.mice, '_seizures_validated.xlsx'))];
    for f = keep_files
        if isfile(f{1}) && ~isfile(fullfile(bk, rel_name(f{1}, consolidate_out))), copy_into(f{1}, fullfile(bk, rel_name(f{1}, consolidate_out))); end
    end
    tw = fullfile(upd.dir, 'write');
    if isfolder(tw), rmdir(tw, 's'); end
    mkdir(tw);
    xt = fullfile(tw, xs);
    for f = 1:3
        copyfile(fullfile(out, xs{f}), xt{f});
        for s = 1:numel(shs{f}), writetable(tabs{f}{s}, xt{f}, 'Sheet', shs{f}{s}, 'WriteMode', 'overwritesheet'); end
    end
    n_kept = check_kept(fullfile(out, xs), xt, shs, upd.mice);
    for f = 1:3, copyfile(xt{f}, fullfile(out, xs{f})); end
    for m = upd.mice
        p = fullfile(vd, [m{1} '_seizures_validated.xlsx']);
        if isfile(p), delete(p); end
        writetable(SZ(strcmp(SZ.mouse, m{1}), cols), p, 'Sheet', 'seizures');
    end
    copyfile(fullfile(out, xs{3}), fullfile(work, 'iid_summary_source.xlsx'));   % IID tables for later reuse runs
    keep_iid_results(C, '055-s', fullfile(work, 'stats', 'iid_rerun_final'));
    extra.update = update_notes(upd, CHK, SZfig, n_kept, n_iid_same, af_rule_055s, iid_note, isfield(CV, 'caso1'));
    extra.update_head = sprintf('## Actualización %s: %s', upd.tag, strjoin(upd.mice, ', '));
    if ~isempty(upd.desc), extra.update_head = sprintf('## Actualización %s: %s', upd.tag, upd.desc); end
    % the consolidated channels' complexes (consolidate_4OHT.m): this update's and the earlier updates' others
    iid_C = C;
    sb = fullfile(bk, '_work', 'state.mat');
    if ~isempty(whos('-file', sb, 'iid_C'))
        q = load(sb, 'iid_C');
        iid_C = [q.iid_C(~ismember(q.iid_C.mouse, upd.mice), :); C];
        clear q
    end
    save(fullfile(work, 'state.mat'), 'cfg', 'S', 'O', 'Z', 'rv', 'edf_tbl', 'ch_tbl', 'det_rows', 'SZ', 'snapinfo', 'iid_C', '-v7.3');
else
    for f = xs, if isfile(fullfile(out, f{1})), delete(fullfile(out, f{1})); end, end
    for f = 1:3
        for s = 1:numel(shs{f}), writetable(tabs{f}{s}, fullfile(out, xs{f}), 'Sheet', shs{f}{s}); end
    end
    if ~consolidate_reuse
        copyfile(fullfile(out, xs{3}), fullfile(work, 'iid_summary_source.xlsx'));   % IID tables for later reuse runs
        iid_C = C;
        save(state_file, 'iid_C', '-append');
    end
    if isfolder(vd), rmdir(vd, 's'); end
    mkdir(vd);
    for m = unique(S.units.mouse(S.units.included))'
        writetable(SZ(strcmp(SZ.mouse, m{1}), cols), fullfile(vd, [m{1} '_seizures_validated.xlsx']), 'Sheet', 'seizures');
    end
end
write_readme(fullfile(out, 'README.md'), cfg, S, rv, SZ, FJ, CMP, iid_note, worker_errors, snapinfo, n_out_obs, extra);

%% ============================ REPORT =====================================
fprintf('\n==== consolidated: %s\n', out);
fprintf('  n_TP %d | n_missed %d | n_missed_DX %d | total %d\n', sum(SM.n_TP), sum(SM.n_missed), sum(SM.n_missed_DX), height(SZ));
fprintf('  figures: %d of %d seizures (.fig + .png; %d EDFs with an error)\n', sum(FJ.ok), height(SZ), numel(unique(FJ.job(~FJ.ok))));
fprintf('  recomputed traces vs the pipeline on its accepted rows (%d): max |rel diff| peak_energy_ratio %.2e, ll_ratio %.2e\n', ...
    FJ.check_n, FJ.check_peak, FJ.check_ll);
fprintf('  seizures with onset outside the observed time (not in the rate): %d\n', n_out_obs);
for i = find(~FJ.ok)', fprintf('  figure problem: %s: %s\n', SZ.figure{i}, FJ.note{i}); end
for i = 1:numel(worker_errors), fprintf(2, '  worker error: %s\n', worker_errors{i}); end
disp(CMP(~strcmp(CMP.diferencia, ''), :));

%% ======================================================================
function SZ = seizure_table(cfg, S, rv)
% one row per seizure: run_campaign_stats' final set (TP after merges,
% missed EEG and video seizures) + the crisette_focal events; seizure_id =
% 1..N per EDF in onset order, the review numbering kept as original_id
    ev = S.events;
    tz = cfg.tz;
    rows = {};
    for i = 1:height(S.seizures)
        s = S.seizures(i, :);
        com = {};
        if strcmp(s.source{1}, 'detected_TP')
            ids = strtrim(split(s.origin_ids{1}, '+'));
            mem = [];
            for k = 1:numel(ids)
                j = find(strcmp(ev.event_id, ids{k}));
                assert(isscalar(j), 'consolidate_results: event %s not found once', ids{k});
                mem = [mem; j; find(ev.merged_into == j)]; %#ok<AGROW>
            end
            for k = 1:size(cfg.merge_events, 1)
                c = cfg.merge_events(k, :);
                if any(strcmp(ev.campaign(mem), c{3}) & strcmp(ev.mouse(mem), c{2}) & strcmp(ev.edf(mem), c{4}) & ismember(ev.seizure_id(mem), c{5}))
                    if isscalar(c{5})
                        com{end+1} = sprintf('%s: fragments of seizure %d joined', c{1}, c{5}); %#ok<AGROW>
                    else
                        com{end+1} = sprintf('%s: review seizures %s merged', c{1}, strjoin(arrayfun(@num2str, c{5}, 'UniformOutput', false), '+')); %#ok<AGROW>
                    end
                end
            end
            if numel(ids) > 1, com{end+1} = ['same seizure detected in ' strjoin(ids', ' and ')]; end %#ok<AGROW>
            cm = ev.comments(mem); cm = unique(cm(~cellfun(@isempty, cm)));
            if ~isempty(cm), com{end+1} = ['user comment (AF): ' strjoin(regexprep(cm, '^unrecognized:', '')', '; ')]; end %#ok<AGROW>
            sid = unique(ev.seizure_id(mem))';
            row = {s.mouse{1}, ev.edf{mem(1)}, ev.campaign{mem(1)}, strjoin(arrayfun(@num2str, sid, 'UniformOutput', false), '+'), ...
                hemis(ev.laterality(mem)), 'pipeline_TP'};
        else
            id = s.origin_ids{1};
            k = find(rv.nd.id == string(id), 1);
            com{end+1} = sprintf('%s: %s', id, nd_type_en(char(rv.nd.tipo(k)))); %#ok<AGROW>
            if strcmp(s.source{1}, 'video') && startsWith(lower(rv.nd.analizable(k)), 'no (desconectado')
                com{end+1} = 'video only (no EEG: disconnected, after the end of the recorded EEG)'; %#ok<AGROW>
            elseif strcmp(s.source{1}, 'video')
                com{end+1} = 'video only (signal is noise)'; %#ok<AGROW>
            end
            edfs = split(s.edf{1}, ','); camps = split(s.campaign{1}, ',');
            row = {s.mouse{1}, edfs{1}, camps{1}, id, 'n/a', 'missed_by_pipeline'};
        end
        dx = s.in_dx;
        if iscell(dx), dx = dx{1}; end
        if isempty(dx) || all(ismissing(string(dx))) || strlength(string(dx)) == 0
            dx = '';
        else
            dx = char(string(dx));
        end
        if ~isempty(dx), com{end+1} = sprintf('onset within disconnection %s', dx); end %#ok<AGROW>
        rows(end+1, :) = [row, {dx, s.start_posix, s.end_posix, strjoin(com, ' | ')}]; %#ok<AGROW>
    end
    for i = 1:height(S.crisette)
        c = S.crisette(i, :);
        j = find(strcmp(ev.campaign, c.campaign{1}) & strcmp(ev.mouse, c.mouse{1}) & strcmp(ev.edf, c.edf{1}) & ev.seizure_id == c.seizure_id);
        com = sprintf('%s: focal crisette (not counted by run_campaign_stats)', c.obs_id{1});
        if isempty(ev.label{j}), com = [com '; no AE label, classified in sheet 5_Observaciones']; end
        if ~isempty(ev.comments{j}), com = [com '; user comment (AF): ' regexprep(ev.comments{j}, '^unrecognized:', '')]; end
        rows(end+1, :) = {c.mouse{1}, c.edf{1}, c.campaign{1}, num2str(c.seizure_id), hemis(ev.laterality(j)), 'pipeline_TP', ...
            '', ev.start_posix(j), ev.end_posix(j), com}; %#ok<AGROW>
    end
    SZ = cell2table(rows, 'VariableNames', {'mouse', 'edf', 'fig_campaign', 'original_id', 'hemispheres', 'source', ...
        'in_dx', 'start_posix', 'end_posix', 'comments'});
    SZ = sortrows(SZ, {'mouse', 'edf', 'start_posix'});
    [g, ~] = findgroups(SZ(:, {'mouse', 'edf'}));
    SZ.seizure_id = zeros(height(SZ), 1);
    for k = 1:max(g)
        j = find(g == k);
        SZ.seizure_id(j) = (1:numel(j))';   % sorted by onset within the EDF
    end
    SZ = sortrows(SZ, {'mouse', 'start_posix'});
    SZ.event_id = arrayfun(@(i) sprintf('%s|%s|%d', SZ.mouse{i}, SZ.edf{i}, SZ.seizure_id(i)), (1:height(SZ))', 'UniformOutput', false);
    SZ.group = cellfun(@(m) S.units.group{find(strcmp(S.units.mouse, m), 1)}, SZ.mouse, 'UniformOutput', false);
    SZ.edf_file = strcat(SZ.edf, '.EDF');
    SZ.duration_s = SZ.end_posix - SZ.start_posix;
    fmt = @(p) cellstr(string(datetime(p, 'ConvertFrom', 'posixtime', 'TimeZone', tz), 'dd/MM/yyyy HH:mm:ss.SSS'));
    SZ.onset = fmt(SZ.start_posix);
    SZ.offset = fmt(SZ.end_posix);
    on = hm(cfg.lights_on); off = hm(cfg.lights_off);
    d = datetime(SZ.start_posix, 'ConvertFrom', 'posixtime', 'TimeZone', tz);
    mod_ = hour(d) * 60 + minute(d) + second(d) / 60;
    SZ.light_phase = repmat({'dark'}, height(SZ), 1);
    SZ.light_phase(mod_ >= on & mod_ < off) = {'light'};
    SZ.figure = cell(height(SZ), 1);   % relative path without extension (.fig and .png)
    tags = regexp(SZ.edf, '^[^~]+', 'match', 'once');
    for i = 1:height(SZ)
        tag = tags{i};
        % two EDFs of the mouse from the same day, both with seizures (016-s, 2026-11-04): + the first block of the EDF's id
        if any(strcmp(SZ.mouse, SZ.mouse{i}) & strcmp(tags, tag) & ~strcmp(SZ.edf, SZ.edf{i}))
            id8 = regexp(SZ.edf{i}, '_([0-9a-f]{8})-', 'tokens', 'once');
            tag = [tag '_' id8{1}];
        end
        seg = regexp(SZ.edf{i}, '_128s_(A\d)$', 'tokens', 'once');
        if ~isempty(seg), tag = [tag '_' seg{1}]; end
        kind = 'TP';
        if strcmp(SZ.source{i}, 'missed_by_pipeline'), kind = 'missed'; end
        SZ.figure{i} = fullfile(SZ.mouse{i}, sprintf('%s_%s_sz%02d_%s', SZ.mouse{i}, tag, SZ.seizure_id(i), kind));
    end
    [u, ~, g] = unique(SZ.figure);
    if numel(u) < height(SZ)
        error('consolidate_results:FigureName', 'duplicate figure names: %s', strjoin(u(accumarray(g, 1) > 1)', ', '));
    end
end

function s = nd_type_en(t)
% 'Tipo' of 1_Crisis_no_detectadas (Spanish) in English; stops on an unknown one
    if startsWith(t, 'Crisis no detectada') && contains(t, 'atenuada')
        s = 'seizure not detected by the pipeline (attenuated signal; not even as an FP)';
    elseif startsWith(t, 'Crisis no detectada')
        s = 'seizure not detected by the pipeline (not even as an FP)';
    elseif contains(t, 'electrogr')
        s = 'apparent electrographic seizure, not detected';
    elseif contains(t, 'video') && ~isempty(regexp(t, 'DX-\d+', 'once'))
        s = ['seizure seen on video during disconnection ' regexp(t, 'DX-\d+', 'match', 'once')];
    else
        error('consolidate_results:NdType', 'no English text for the missed-seizure type "%s"', t);
    end
end

function h = hemis(lat)
% union of the events' laterality (bilateral / HPCr_only / HPCl_only)
    r = any(ismember(lat, {'bilateral', 'HPCr_only'}));
    l = any(ismember(lat, {'bilateral', 'HPCl_only'}));
    h = {'n/a', 'HPCl'; 'HPCr', 'HPCr+HPCl'};
    h = h{r + 1, l + 1};
end

function m = hm(s)
    v = sscanf(s, '%d:%d'); m = v(1) * 60 + v(2);
end

function [IE, IM] = iid_from_xlsx(f)
% the IID tables of the full run's iid_summary.xlsx (sheet names of
% 2026-11-02 or the current ones), without the campaign column
    sh = sheetnames(f);
    se = intersect({'per_mouse_EDF_hemisphere', 'por_raton_EDF_hemisferio'}, sh);
    sm = intersect({'per_mouse', 'por_raton'}, sh);
    assert(isscalar(se) && isscalar(sm), 'iid_from_xlsx: sheets of %s not recognized (%s)', f, strjoin(sh, ', '));
    IE = readtable(f, 'Sheet', se{1}, 'TextType', 'char', 'VariableNamingRule', 'preserve');
    if ismember('campaign', IE.Properties.VariableNames), IE.campaign = []; end
    IM = readtable(f, 'Sheet', sm{1}, 'TextType', 'char', 'VariableNamingRule', 'preserve');
end

%% ---------------------------------------------------------------------------
function S = nd_beyond_eeg(S, rv, edf_tbl, ids)
% missed seizures seen on video after the end of the recorded EEG (ND-08/09 of 055-s: the EEG was cut
% at DX-23): cs_build_seizures took them as missed_eeg (nd_added) and found no EDF holding their onset;
% here source = video and the EDF of their row in 1_Crisis_no_detectadas
    for k = 1:numel(ids)
        i = find(strcmp(S.seizures.origin_ids, ids{k}));
        j = find(rv.nd.id == ids{k});
        assert(isscalar(i) && isscalar(j), 'nd_beyond_eeg: %s not found once', ids{k});
        e = find(strcmp(edf_tbl.mouse, S.seizures.mouse{i}) & strcmp(edf_tbl.edf, char(rv.nd.edf(j))));
        assert(isscalar(e), 'nd_beyond_eeg: EDF %s of %s not found once among the analyzed EDFs', rv.nd.edf(j), ids{k});
        assert(isempty(S.seizures.edf{i}) && S.seizures.start_posix(i) > edf_tbl.end_posix(e) && ~isempty(S.seizures.in_dx{i}), ...
            'nd_beyond_eeg: %s is not a seizure after the end of its EDF, inside a disconnection', ids{k});
        S.seizures.source{i} = 'video';
        S.seizures.campaign{i} = edf_tbl.campaign{e};
        S.seizures.edf{i} = edf_tbl.edf{e};
    end
end

function S = nd_video_in_eeg(S, rv, ids)
% missed seizures seen on video during a disconnection inside the recorded EEG (ND-12 of 072-s, DX-07; as
% ND-05/06): cs_build_seizures took them as missed_eeg (nd_added); here source = video, after the check
% cs_build_seizures makes for its video seizures (inside a DX with Video revisado = Sí)
    for k = 1:numel(ids)
        i = find(strcmp(S.seizures.origin_ids, ids{k}));
        j = find(rv.nd.id == ids{k});
        assert(isscalar(i) && isscalar(j) && strcmp(S.seizures.source{i}, 'missed_eeg'), 'nd_video_in_eeg: %s not found once as missed_eeg', ids{k});
        assert(strcmp(S.seizures.edf{i}, char(rv.nd.edf(j))), 'nd_video_in_eeg: %s is not in the EDF of its row', ids{k});
        inv = any(rv.dx.mouse == S.seizures.mouse{i} & rv.dx.video_reviewed & rv.dx.start_posix <= S.seizures.start_posix(i) & ...
            rv.dx.end_posix >= S.seizures.end_posix(i));
        assert(inv, 'nd_video_in_eeg: %s (video) is not inside a DX interval with Video revisado = Sí', ids{k});
        S.seizures.source{i} = 'video';
    end
end

function O = dx_video_beyond_eeg(O, cfg, rv, edf_tbl, ids)
% disconnections reviewed on video after the end of the recorded EEG (055-s, DX-23: the EEG was cut at its
% start; the user, 2026-11-05): cs_observation_time counts a reviewed DX only inside the analyzed EDFs, so
% its part outside them is added here to the mouse's observed seizure time (T_sz, and its clock exposures)
% in its EDF's campaign scope and the merged one, and to that EDF's per-day rows (T_seizure_obs_h)
    tz = cfg.tz;
    for k = 1:numel(ids)
        d = find(rv.dx.id == ids{k});
        assert(isscalar(d) && rv.dx.video_reviewed(d), 'dx_video_beyond_eeg: %s not found once with Video revisado = Sí', ids{k});
        m = char(rv.dx.mouse(d));
        e = find(strcmp(edf_tbl.mouse, m) & strcmp(edf_tbl.edf, char(rv.dx.edf(d))));
        assert(isscalar(e), 'dx_video_beyond_eeg: EDF of %s not found once among the analyzed EDFs', ids{k});
        assert(rv.dx.start_posix(d) >= edf_tbl.end_posix(e) - 1, 'dx_video_beyond_eeg: %s starts before the end of its EDF', ids{k});
        jm = find(strcmp({O.mouse.mouse}, m) & strcmp({O.mouse.scope}, 'merged'));
        jc = find(strcmp({O.mouse.mouse}, m) & strcmp({O.mouse.scope}, edf_tbl.campaign{e}));
        assert(isscalar(jm) && isscalar(jc), 'dx_video_beyond_eeg: observation time of %s not found', m);
        add = cs_iv('subtract', [rv.dx.start_posix(d) rv.dx.end_posix(d)], O.mouse(jm).rec);
        assert(cs_iv('length', add) > 0, 'dx_video_beyond_eeg: %s has no time outside the analyzed EDFs', ids{k});
        for j = [jc jm]
            O.mouse(j).T_sz = cs_iv('union', O.mouse(j).T_sz, add);
            [O.mouse(j).hour_sz_s, ~] = cs_clock_split(O.mouse(j).T_sz, tz, 0:60:1440);
            [O.mouse(j).phase_sz_s, ~] = cs_clock_split(O.mouse(j).T_sz, tz, cfg.phase_edges_min);
        end
        [~, dd] = cs_clock_split(add, tz, [0 1440]);
        for i = 1:height(dd)
            r = find(strcmp(O.edf_day.mouse, m) & strcmp(O.edf_day.edf, edf_tbl.edf{e}) & O.edf_day.date == dd.date(i));
            q = find(strcmp(O.day.mouse, m) & O.day.date == dd.date(i));
            assert(isscalar(r) && isscalar(q), 'dx_video_beyond_eeg: no per-day row of %s on %s', ids{k}, char(dd.date(i)));
            O.edf_day.T_seizure_obs_h(r) = O.edf_day.T_seizure_obs_h(r) + dd.total_s(i) / 3600;
            O.day.T_seizure_obs_h(q) = O.day.T_seizure_obs_h(q) + dd.total_s(i) / 3600;
        end
        O.day.valid_day = O.day.T_seizure_obs_h >= cfg.min_valid_h_day;
        fprintf('%s (%s): +%.2f h observed on video after the end of the EEG\n', ids{k}, m, cs_iv('length', add) / 3600);
    end
end

function IV = iid_sets(cfg, S, rv, ch_tbl, det_rows, Z, mice, root, rules)
% channel sets of the IID re-detection, each {name, cfg, Z} with Z.channels.run = its channels.
% 'principal' = the zones of cs_iid_zones (general rule) on the channels of `mice` but 055-s. 055-s =
% one set per AF case in `rules` (af_rule_055s; 0 = the general rule, for reference): the zones of
% cs_iid_zones with 055-s's FP comments replaced, each case in its own root under `root` (the cases
% reuse each other's channels when the zone intervals are the same)
    m55 = '055-s';
    Zp = Z;
    Zp.channels.run = Z.channels.run & ismember(Z.channels.mouse, setdiff(mice, m55));
    IV = struct('name', {'principal'}, 'cfg', {cfg}, 'Z', {Zp});
    if ~ismember(m55, mice), return; end
    if isempty(rules) || ~all(ismember(rules, 0:3))
        error('consolidate_results:AfRule055s', 'af_rule_055s not decided (1, 2 or 3) -- ask the user');
    end
    % AF of each 055-s FP event per hemisphere: that hemisphere's rows (the event's union if they have none)
    ev = S.events;
    fp = find(strcmp(ev.mouse, m55) & strcmp(ev.label, 'FP'));
    af = struct('HPCr', {ev.comments(fp)}, 'HPCl', {ev.comments(fp)});
    for i = 1:numel(fp)
        r = det_rows(strcmp(det_rows.campaign, ev.campaign{fp(i)}) & strcmp(det_rows.mouse, m55) & ...
            strcmp(det_rows.edf, ev.edf{fp(i)}) & det_rows.seizure_id == ev.seizure_id(fp(i)), :);
        for h = {'HPCr', 'HPCl'}
            c = r.comment(strcmp(r.region, h{1}) & ~cellfun(@isempty, r.comment));
            if ~isempty(c), af.(h{1}){i} = strjoin(unique(c)', ';'); end
        end
    end
    other = cellfun(@(c) strjoin(unique(regexprep(split(c, ';'), '^IIDs$', 'other'))', ';'), ev.comments(fp), 'UniformOutput', false);
    %      rule  set              root             hemispheres        055-s FP comments ([] = as they are)
    def = {0,    'regla_general', 'regla_general', {'HPCr', 'HPCl'}, []; ...
           1,    'caso1',         'caso1',         {'HPCr', 'HPCl'}, af.HPCl; ...
           2,    'caso2_HPCr',    'caso2',         {'HPCr'},         af.HPCr; ...
           2,    'caso2_HPCl',    'caso2',         {'HPCl'},         af.HPCl; ...
           3,    'caso3',         'caso3',         {'HPCr', 'HPCl'}, other};
    sel = find(ismember([def{:, 1}], rules));
    roots = reshape(fullfile(root, unique(def(sel, 3))), 1, []);
    for d = sel
        Zs = Z;
        if ~isempty(def{d, 5})
            Sx = S;
            Sx.events.comments(fp) = def{d, 5};
            Zs = cs_iid_zones(cfg, Sx, rv, ch_tbl, det_rows);
        end
        Zs.channels.run = Zs.channels.run & strcmp(Zs.channels.mouse, m55) & ismember(Zs.channels.region, def{d, 4});
        c = cfg;
        own = fullfile(root, def{d, 3});
        c.iid_rerun_dirs.principal = [{own}, setdiff(roots, {own}, 'stable'), cfg.iid_rerun_dirs.principal];
        c.iid_rerun_dirs.sin_crisette = [{[own '_sin_crisette']}, cfg.iid_rerun_dirs.sin_crisette];
        IV(end+1) = struct('name', def{d, 2}, 'cfg', c, 'Z', Zs); %#ok<AGROW>
    end
end

function C = iid_final(CV, rule)
% the consolidation's IID channels: the general rule for every mouse but 055-s, and 055-s's AF case
    if ~(isscalar(rule) && any(rule == 1:3))
        error('consolidate_results:AfRule055s', 'af_rule_055s not decided (1, 2 or 3) -- ask the user');
    end
    sets = {{'caso1'}, {'caso2_HPCr', 'caso2_HPCl'}, {'caso3'}};
    C = CV.principal;
    for s = sets{rule}
        if isfield(CV, s{1}), C = [C; CV.(s{1})]; end %#ok<AGROW>
    end
end

function I = iid_channel_table(C, S)
% the per-channel columns of cs_iid_stats' I.channels that iid_summary uses, and its pooled per-mouse
% rates (I.mouse: merged; per hemisphere and both), for a subset of channels
    Cm = C(strcmp(C.variant, 'principal'), :);
    grp = cellfun(@(m) S.units.group{find(strcmp(S.units.mouse, m), 1)}, Cm.mouse, 'UniformOutput', false);
    I.channels = table(Cm.campaign, Cm.mouse, grp, Cm.edf, Cm.region, Cm.analyzed_min / 60, Cm.baseline_uV, Cm.lower_uV, ...
        Cm.p97_uV, Cm.n_complexes, Cm.n_single, Cm.n_polyspike, Cm.n_bursts, 'VariableNames', {'campaign', 'mouse', 'group', ...
        'edf', 'region', 'T_iid_h', 'baseline_uV', 'lower_uV', 'p97_real_uV', 'n_complexes', 'n_single', 'n_polyspike', 'n_bursts'});
    [g, M] = findgroups(Cm(:, {'mouse'}));
    M.hemisphere = repmat({'ambos'}, height(M), 1);
    [g2, M2] = findgroups(Cm(:, {'mouse', 'region'}));
    M2 = renamevars(M2, 'region', 'hemisphere');
    M.complexes_per_min = splitapply(@sum, Cm.n_complexes, g) ./ splitapply(@sum, Cm.analyzed_min, g);
    M2.complexes_per_min = splitapply(@sum, Cm.n_complexes, g2) ./ splitapply(@sum, Cm.analyzed_min, g2);
    M = [M; M2];
    M.scope = repmat({'merged'}, height(M), 1);
    I.mouse = M;
end

function [IE, IM, n_same] = iid_merge(Z, mice, IE4, IM4, IEo, IMo)
% IID tables of an update: the rows of `mice` from this run, every other mouse's as consolidated before
% (IEo, IMo); channel rows in the plan's order (as cs_iid_stats gives them). A channel of an update
% mouse that was already consolidated must give the same values (n_same channels checked)
    P = Z.channels(Z.channels.run, :);
    key = @(m, e, h) strcat(m, '|', erase(e, '.EDF'), '|', h);
    kP = key(P.mouse, P.edf, P.region);
    kN = key(IE4.mouse, IE4.edf_file, IE4.hemisphere);
    kO = key(IEo.mouse, IEo.edf_file, IEo.hemisphere);
    new = ismember(P.mouse, mice);
    old_keep = ~ismember(IEo.mouse, mice);
    [tn, jn] = ismember(kP, kN);
    [to, jo] = ismember(kP, kO);
    assert(all(tn(new)) && isequal(sort(jn(new)), (1:height(IE4))'), 'iid_merge: the update''s IID rows are not its plan channels');
    assert(all(to(~new)) && isequal(jo(~new), find(old_keep)), 'iid_merge: the other mice''s IID rows are not their plan channels, in order');
    idx = zeros(height(P), 1);
    idx(new) = jn(new);
    idx(~new) = height(IE4) + jo(~new);
    T = [IE4; IEo];
    IE = T(idx, :);
    ou = find(~old_keep);
    [tf, j] = ismember(kO(ou), kN);
    if ~all(tf)
        error('consolidate_results:IidLost', 'consolidated IID channels no longer in the plan: %s', strjoin(kO(ou(~tf))', ', '));
    end
    vn = IEo.Properties.VariableNames(varfun(@isnumeric, IEo, 'OutputFormat', 'uniform'));
    bad = {};
    for i = 1:numel(ou)
        a = IE4{j(i), vn}; b = IEo{ou(i), vn};
        ok = isequal(isnan(a), isnan(b)) && all(abs(a(~isnan(a)) - b(~isnan(b))) <= 1e-9 * max(1, abs(b(~isnan(b)))));
        if ~ok, bad{end+1} = kO{ou(i)}; end %#ok<AGROW>
    end
    if ~isempty(bad)
        error('consolidate_results:IidDiffers', 'IIDs of already consolidated channels differ from the previous consolidation: %s', strjoin(bad, ', '));
    end
    n_same = numel(ou);
    IM = sortrows([IMo(~ismember(IMo.mouse, mice), :); IM4], {'mouse', 'hemisphere'});   % HPCl before HPCr
end

function iid_by_hemisphere(out, work, bk, af_rule_055s)
% one-off (2026-11-04): the IID per_mouse sheet from one row per mouse (hemispheres summed) to one per
% mouse x hemisphere, from the per-channel results the consolidation holds (nothing re-detected). They
% must give per_mouse_EDF_hemisphere and the pooled per_mouse exactly; then per_mouse is rewritten
% (every other sheet checked unchanged), iid_summary_source.xlsx follows, README.md is regenerated.
% Previous files in bk. per_mouse_4-OHT: consolidate_4OHT.m, run after this.
    q = load(fullfile(work, 'state.mat'), 'cfg', 'S', 'O', 'Z', 'rv', 'SZ', 'snapinfo', 'iid_C');
    cfg = q.cfg; S = q.S; O = q.O; Z = q.Z; rv = q.rv; SZ = q.SZ; snapinfo = q.snapinfo; iid_C = q.iid_C;
    clear q
    xi = fullfile(out, 'iid_summary.xlsx');
    [IEo, IMo] = iid_from_xlsx(xi);
    if ismember('hemisphere', IMo.Properties.VariableNames), error('consolidate_results:IidHemisphere', 'per_mouse is already by hemisphere'); end
    C = iid_stored_channels(cfg, Z, IEo, iid_C, work);
    I = iid_channel_table(C, S);
    K = cs_iid_channel_intervals(C, S);
    [IEc, IMp] = iid_summary(cfg, I, K, C, 'mouse');
    if ~same_values(IEc, IEo) || ~same_values(IMp, IMo)
        error('consolidate_results:IidHemisphere', 'the stored per-channel results do not give the current IID sheets -- nothing written');
    end
    [~, IM] = iid_summary(cfg, I, K, C);
    % into a copy; the file is replaced only if every other sheet is unchanged
    if ~isfolder(bk), mkdir(bk); end
    for f = {fullfile(out, 'iid_summary.xlsx'), fullfile(out, 'README.md'), fullfile(work, 'iid_summary_source.xlsx')}
        [~, nm, ext] = fileparts(f{1});
        if ~isfile(fullfile(bk, [nm ext])), copyfile(f{1}, fullfile(bk, [nm ext])); end
    end
    xt = fullfile(bk, 'iid_summary_new.xlsx');
    copyfile(xi, xt);
    writetable(IM, xt, 'Sheet', 'per_mouse', 'WriteMode', 'overwritesheet');
    sb = sheetnames(xi);
    if ~isequal(sb, sheetnames(xt)), error('consolidate_results:IidHemisphere', 'the sheets would change -- nothing written'); end
    for s = sb'
        if ~strcmp(s{1}, 'per_mouse') && ~isequaln(readcell(xi, 'Sheet', s{1}), readcell(xt, 'Sheet', s{1}))
            error('consolidate_results:IidHemisphere', 'sheet %s would change -- nothing written', s{1});
        end
    end
    movefile(xt, xi, 'f');
    copyfile(xi, fullfile(work, 'iid_summary_source.xlsx'));   % IID tables for later reuse runs
    % README: as the last write, with this change (the IID note of Números kept as it is)
    FJ = fig_results(SZ, work, out);
    FJ.note(~SZ.has_eeg) = {'no EEG at the onset (video-only seizure after the end of the recorded EEG)'};
    [~, SM, n_out_obs] = seizure_summary(S, O, SZ);
    CMP = compare_20261101(cfg, SM, IM);
    t = cellstr(splitlines(string(fileread(fullfile(out, 'README.md')))));
    kn = find(startsWith(t, '## Números'), 1);
    ki = kn + find(startsWith(t(kn + 1:end), '- IIDs: '), 1);
    iid_note = t{ki}(3:end);
    extra = readme_notes(cfg, S, SZ, af_rule_055s, rv);
    extra.update_head = '## Actualización 20261104c: IIDs por hemisferio en `per_mouse`';
    extra.update = {sprintf(['Con `consolidate_update = ''iid_hemisphere''` (`consolidate_results.m`), el %s: la hoja `per_mouse` de ' ...
        '`iid_summary.xlsx` pasa de una fila por ratón (hemisferios sumados) a una por ratón y hemisferio (HPCl, HPCr), con la columna ' ...
        '`hemisphere` y `n_edfs` en lugar de `n_channels`; conteos, tiempos, tasas y amplitud de cada hemisferio por separado. No se ' ...
        're-detectó nada: los resultados por canal guardados reproducen exactamente `per_mouse_EDF_hemisphere` (%d canales) y el ' ...
        '`per_mouse` anterior (%d ratones); las demás hojas no cambian. `per_mouse_4-OHT`: igual, con `consolidate_4OHT.m`.'], ...
        char(datetime('now', 'Format', 'yyyy-MM-dd')), height(IEo), height(IMo)), ...
        sprintf('- Archivos anteriores: `%s/`.', rel_name(bk, out))};
    write_readme(fullfile(out, 'README.md'), cfg, S, rv, SZ, FJ, CMP, iid_note, {}, snapinfo, n_out_obs, extra);
    fprintf('IID per_mouse by hemisphere: %d rows (%d mice) in %s\n', height(IM), numel(unique(IM.mouse)), xi);
end

function C = iid_stored_channels(cfg, Z, IE, iid_C, work)
% one row per row of per_mouse_EDF_hemisphere, with cs_iid_rerun_all's columns that the summaries use,
% from the stored re-detection: the state's iid_C (the updates' mice), the 2026-11-01 run
% (phase3_results.mat) or _work/stats/iid_rerun_final, the first that holds the channel (as
% consolidate_4OHT.m); it must give the row's n_IID and observed_min
    q = load(fullfile(cfg.batch_root, '_stats_revision_20261101_final', 'phase3_results.mat'), 'C');
    C0 = q.C(strcmp(q.C.variant, 'principal'), :);
    U = iid_C(strcmp(iid_C.variant, 'principal'), :);
    P = Z.channels(Z.channels.run, :);
    key = @(m, e, h) strcat(m, '|', erase(e, '.EDF'), '|', h);
    kU = key(U.mouse, U.edf, U.region); k0 = key(C0.mouse, C0.edf, C0.region); kP = key(P.mouse, P.edf, P.region);
    kE = key(IE.mouse, IE.edf_file, IE.hemisphere);
    cols = {'campaign', 'mouse', 'edf', 'region', 'variant', 'start_posix', 'duration_s', 'baseline_uV', 'lower_uV', 'p97_uV', ...
        'analyzed_min', 'n_complexes', 'n_single', 'n_polyspike', 'n_bursts', 'complexes', 'bursts', 'excl'};
    rows = cell(height(IE), 1);
    for i = 1:height(IE)
        jp = find(strcmp(kP, kE{i}));
        ju = find(strcmp(kU, kE{i}));
        j0 = find(strcmp(k0, kE{i}));
        if ~isscalar(jp), error('consolidate_results:IidStored', '%s: not once in the IID plan', kE{i}); end
        if isscalar(ju)
            r = U(ju, cols);
        elseif isscalar(j0)
            r = C0(j0, cols);
        else
            f = dir(fullfile(work, 'stats', 'iid_rerun_final', IE.mouse{i}, erase(IE.edf_file{i}, '.EDF'), ...
                sprintf('*_%s_IID_rerun.mat', IE.hemisphere{i})));
            if ~isscalar(f), error('consolidate_results:IidStored', '%s: no stored re-detection', kE{i}); end
            w = load(fullfile(f.folder, f.name), 'rr');
            s = w.rr.summary; cx = w.rr.spike_complex_table; b = w.rr.burst_table; e = w.rr.exclusion_table;
            r = table(P.campaign(jp), P.mouse(jp), P.edf(jp), P.region(jp), {'principal'}, P.start_posix(jp), P.duration_s(jp), ...
                s.baseline_uV, s.lower_threshold_uV, w.rr.p97_uV, s.analyzed_duration_min, s.total_complexes, s.n_single, s.n_polyspike, ...
                s.n_bursts, {[cx.Start_s cx.End_s cx.N_spikes double(cx.Is_polyspike) cx.Max_amplitude_mV]}, ...
                {[b.Start_s b.End_s b.N_complexes]}, {[e.start_s e.end_s]}, 'VariableNames', cols);
        end
        if r.n_complexes ~= IE.n_IID(i) || abs(r.analyzed_min - IE.observed_min(i)) > 1e-6
            error('consolidate_results:IidStored', '%s: the stored re-detection differs from per_mouse_EDF_hemisphere', kE{i});
        end
        rows{i} = r;
    end
    C = vertcat(rows{:});
end

function ok = same_values(A, B)
% same columns and rows; numbers to 1e-9 relative (NaN == NaN), text exact
    ok = isequal(A.Properties.VariableNames, B.Properties.VariableNames) && height(A) == height(B);
    for v = A.Properties.VariableNames
        if ~ok, break; end
        a = A.(v{1}); b = B.(v{1});
        if isnumeric(a) && isnumeric(b)
            ok = isequal(isnan(a), isnan(b)) && all(abs(a(~isnan(a)) - b(~isnan(b))) <= 1e-9 * max(1, abs(b(~isnan(b)))));
        else
            ok = isequal(a, b);
        end
        if ~ok, fprintf(2, 'same_values: column %s differs\n', v{1}); end
    end
end

function CHK = update_checks(old, SZ, S, O, Z, mice, buf)
% incremental update: every mouse outside `mice` must be as in the previous consolidation (old =
% _work/state.mat): same seizures (every column), per-EDF and per-mouse seizure rows, IID channels
% and zone intervals; an earlier seizure of an update mouse (found by mouse, EDF and original_id) keeps
% its times, and its seizure_id and figure name unless new seizures come before it in its EDF: then its
% id grows by exactly their number (CHK.renumbered); other column changes are listed in CHK.changed
    keep = @(T) ~ismember(T.mouse, mice);
    c = intersect(old.SZ.Properties.VariableNames, SZ.Properties.VariableNames, 'stable');   % not the metrics of a written state
    A = SZ(keep(SZ), c); B = old.SZ(keep(old.SZ), c);
    if ~isequaln(A, B)
        error('consolidate_results:UpdateSeizures', 'seizures of the other mice differ from the previous consolidation (%s)', strjoin(diff_cols(A, B), ', '));
    end
    [SEn, SMn] = seizure_summary(S, O, SZ);
    [SEo, SMo] = seizure_summary(old.S, old.O, old.SZ);
    if ~isequaln(SEn(keep(SEn), :), SEo(keep(SEo), :)) || ~isequaln(SMn(keep(SMn), :), SMo(keep(SMo), :))
        error('consolidate_results:UpdateSummary', 'per_EDF / per_mouse rows of the other mice differ from the previous consolidation');
    end
    key = @(P) strcat(P.campaign, '|', P.mouse, '|', P.edf, '|', P.region);
    Pn = Z.channels(Z.channels.run & keep(Z.channels), :);
    Po = old.Z.channels(old.Z.channels.run & keep(old.Z.channels), :);
    if ~isequal(key(Pn), key(Po))
        error('consolidate_results:UpdateIidPlan', 'IID channels of the other mice differ from the previous consolidation');
    end
    kn = key(Pn);
    for i = 1:height(Pn)
        [sa, ma] = cs_iid_channel_zones(Z.zones, Pn(i, :), buf);
        [sb, mb] = cs_iid_channel_zones(old.Z.zones, Po(i, :), buf);
        if ~isequal(sortrows(sa{:, {'start_s', 'end_s'}}), sortrows(sb{:, {'start_s', 'end_s'}})) || ~isequal(sortrows(ma), sortrows(mb))
            error('consolidate_results:UpdateIidZones', '%s: IID zones differ from the previous consolidation', kn{i});
        end
    end
    Bu = old.SZ(~keep(old.SZ), :);   % earlier seizures of the update mice
    kid = @(T) strcat(T.mouse, '|', T.edf, '|', T.original_id);
    assert(numel(unique(kid(SZ))) == height(SZ), 'update_checks: mouse | EDF | original_id not unique');
    [tf, j] = ismember(kid(Bu), kid(SZ));
    if ~all(tf)
        error('consolidate_results:UpdateLost', 'earlier seizures no longer found: %s', strjoin(Bu.event_id(~tf)', ', '));
    end
    moved = abs(SZ.start_posix(j) - Bu.start_posix) > 1e-6 | abs(SZ.end_posix(j) - Bu.end_posix) > 1e-6;
    if any(moved)
        error('consolidate_results:UpdateMoved', 'earlier seizures with other times: %s', strjoin(Bu.event_id(moved)', ', '));
    end
    nw = SZ(~keep(SZ) & ~ismember(kid(SZ), kid(Bu)), :);
    shift = arrayfun(@(i) sum(strcmp(nw.mouse, Bu.mouse{i}) & strcmp(nw.edf, Bu.edf{i}) & nw.start_posix < Bu.start_posix(i)), (1:height(Bu))');
    bad = SZ.seizure_id(j) ~= Bu.seizure_id + shift | (shift == 0 & ~strcmp(SZ.figure(j), Bu.figure));
    if any(bad)
        error('consolidate_results:UpdateIds', 'earlier seizures with another id or figure than the new seizures explain: %s', strjoin(Bu.event_id(bad)', ', '));
    end
    r = find(shift > 0);
    CHK.renumbered = table(Bu.mouse(r), Bu.edf(r), Bu.original_id(r), Bu.seizure_id(r), SZ.seizure_id(j(r)), Bu.figure(r), SZ.figure(j(r)), ...
        'VariableNames', {'mouse', 'edf', 'original_id', 'id_before', 'id_after', 'figure_before', 'figure_after'});
    CHK.changed = {};
    for v = c
        for i = find(~arrayfun(@(q) isequaln(SZ.(v{1})(j(q)), Bu.(v{1})(q)), (1:height(Bu))'))'
            if shift(i) > 0 && any(strcmp(v{1}, {'seizure_id', 'event_id', 'figure'})), continue; end   % in CHK.renumbered
            CHK.changed{end+1} = sprintf('%s, %s: "%s" -> "%s"', Bu.event_id{i}, v{1}, val2str(Bu.(v{1})(i)), val2str(SZ.(v{1})(j(i))));
        end
    end
    CHK.new = table(unique(nw.mouse), 'VariableNames', {'mouse'});
    CHK.new.n = cellfun(@(m) sum(strcmp(nw.mouse, m)), CHK.new.mouse);
    CHK.new_figures = nw.figure;   % may be a name an earlier seizure had (renumbered)
    CHK.n_seizures = height(A); CHK.n_edf_rows = sum(keep(SEn)); CHK.n_iid_channels = height(Pn); CHK.n_earlier = height(Bu);
end

function d = diff_cols(A, B)
    if height(A) ~= height(B), d = {sprintf('%d vs %d rows', height(A), height(B))}; return; end
    v = A.Properties.VariableNames;
    d = v(~cellfun(@(x) isequaln(A.(x), B.(x)), v));
end

function s = val2str(x)
    if iscell(x), x = x{1}; end
    if isnumeric(x) || islogical(x), s = num2str(x, 12); else, s = char(string(x)); end
end

function n = check_kept(before, after, sheets, mice)
% the update's workbooks (after) against the files they replace (before): the sheets not rewritten
% unchanged; in the rewritten ones the header and every row of the other mice the same cells
    n = 0;
    for f = 1:numel(before)
        sb = sheetnames(before{f}); sa = sheetnames(after{f});
        if ~isequal(sb, sa), error('consolidate_results:KeptSheets', '%s: the sheets would change -- nothing written', before{f}); end
        for s = sb'
            cb = readcell(before{f}, 'Sheet', s{1}); ca = readcell(after{f}, 'Sheet', s{1});
            if ~any(strcmp(sheets{f}, s{1}))
                ok = isequaln(cb, ca);
            else
                km = find(strcmp(cb(1, :), 'mouse'));
                ok = isscalar(km) && isequaln(cb(1, :), ca(1, :));
                if ok
                    ob = cb(2:end, :); oa = ca(2:end, :);
                    ob = ob(~ismember(string(ob(:, km)), mice), :); oa = oa(~ismember(string(oa(:, km)), mice), :);
                    ok = isequaln(ob, oa);
                    n = n + size(ob, 1);
                end
            end
            if ~ok
                error('consolidate_results:KeptRows', '%s, sheet %s: the other mice''s rows (or the sheet) would change -- nothing written', before{f}, s{1});
            end
        end
    end
end

function copy_new_figures(SZfig, upd, out, work, SZ, CHK)
% the update's figures into figures/<mouse>/ and its figure jobs into _work/fig_jobs as
% job_u<tag>_NNN.mat (fig_results reads every job, in name order). An existing file is replaced only if it
% is a name of a renumbered seizure (before or after); a name no seizure uses any more is removed. Every
% replaced or removed file is kept first in _work/update_<tag>/before/figures
    swap = [CHK.renumbered.figure_before; CHK.renumbered.figure_after];
    keep = @(f) copy_once(fullfile(out, 'figures', f), fullfile(upd.dir, 'before', 'figures', f));
    for i = 1:height(SZfig)
        for e = {'.fig', '.png'}
            src = fullfile(upd.dir, 'figures', [SZfig.figure{i} e{1}]);
            dst = fullfile(out, 'figures', [SZfig.figure{i} e{1}]);
            if ~isfile(src), continue; end   % not drawn: reported by fig_results
            if isfile(dst)
                if isequal(file_bytes(src), file_bytes(dst)), continue; end   % a repeated 'write'
                if ~ismember(SZfig.figure{i}, swap)
                    error('consolidate_results:Figure', '%s already in figures/ and different', [SZfig.figure{i} e{1}]);
                end
                keep([SZfig.figure{i} e{1}]);
            end
            copy_into(src, dst);
        end
    end
    for f = setdiff(CHK.renumbered.figure_before, SZ.figure)'
        for e = {'.fig', '.png'}
            p = fullfile(out, 'figures', [f{1} e{1}]);
            if isfile(p), keep([f{1} e{1}]); delete(p); end
        end
    end
    for f = dir(fullfile(upd.dir, 'fig_jobs', 'job_*.mat'))'
        dst = fullfile(work, 'fig_jobs', ['job_u' upd.tag '_' f.name(5:end)]);
        if ~isfile(dst), copyfile(fullfile(f.folder, f.name), dst); end
    end
end

function b = file_bytes(f)
    fid = fopen(f, 'r'); b = fread(fid, inf, '*uint8'); fclose(fid);
end

function keep_iid_results(C, mouse, root)
% the result files of `mouse`'s consolidated channels (055-s: its AF case, computed in the update's
% own roots) also in the consolidation's IID root, like every other channel run there
    for i = find(strcmp(C.mouse, mouse))'
        [~, nm, ext] = fileparts(C.rerun_mat{i});
        dst = fullfile(root, C.mouse{i}, C.edf{i}, [nm ext]);
        if ~isfile(dst), copy_into(C.rerun_mat{i}, dst); end
    end
end

function r = rel_name(p, root)
% p relative to root, with forward slashes (README, backups)
    r = strrep(erase(p, [root filesep]), filesep, '/');
end

function copy_into(src, dst)
    if ~isfolder(fileparts(dst)), mkdir(fileparts(dst)); end
    copyfile(src, dst);
end

function copy_once(src, dst)
% backup: the first copy is kept (a repeated run must not overwrite it with an already replaced file)
    if ~isfile(dst), copy_into(src, dst); end
end

function x = readme_notes(cfg, S, SZ, rule, rv)
% README lines on the mice validated after 2026-11-02 and their special cases (2026-11-04), the
% disconnections reviewed on video and ND-12 (2026-11-05)
    txt = {'el AF de HPCl se aplica a ambos hemisferios', ...
        'cada hemisferio usa su propio AF (los IIDs se cuentan en el hemisferio marcado "IIDs" y se excluyen en el otro)', ...
        'todos sus FP se tratan como "other" (excluidos de los IIDs)'};
    x.af055 = 'sin decidir';
    if isscalar(rule) && any(rule == 1:3), x.af055 = sprintf('caso %d, %s', rule, txt{rule}); end
    x.update = {};
    x.update_head = '';
    x.scope = {};
    inc = S.units.mouse(S.units.included);
    if ismember('025-s', inc)
        x.scope{end+1} = '- 025-s (C2): control, validado e incluido el 04/11/2026; todos sus eventos son FP "other" (excluidos de los IIDs).';
    end
    if ~all(ismember({'016-s', '055-s', '084-s', '113-s'}, inc)), return; end
    nv = SZ(strcmp(SZ.mouse, '055-s') & ~SZ.has_eeg, :);
    d23 = find(rv.dx.id == "DX-23");
    assert(isscalar(d23) && rv.dx.video_reviewed(d23), 'readme_notes: DX-23 not found once with Video revisado = Sí');
    nrv = find(~rv.dx.video_reviewed);
    x.scope = [x.scope, { ...
        ['- Validados después del 02/11 e incluidos el 04/11/2026: 016-s y 055-s (nuevos), 113-s (EDFs 20260807~ y 20260809~, C1) ' ...
        'y 084-s (EDFs 20260706~ y 20260715A~ de C2, antes fuera por etiquetas incompletas, regla A.2).'], ...
        ['- 055-s: la señal está muy dañada en todos sus EDFs (contaminación de línea, atenuación, outliers recortados): sirve para ' ...
        'cuantificar crisis, pero su conteo de IIDs es dudoso; tenerlo en cuenta en la estadística. IIDs: ' x.af055 ' (decisión del ' ...
        '04/11/2026 entre tres casos para sus 21 FP con AF distinto entre hemisferios; en los demás ratones esos FP se excluyen en ambos ' ...
        'hemisferios y los FP "IIDs" se cuentan).'], ...
        sprintf(['- 055-s, EDF 20260604~ (C2): el EEG se cortó en la desconexión del 04/06/2026 15:32:18 (re-corrida del 04/11/2026). El ' ...
        'video se revisó hasta las %s (DX-23; nota del 05/11/2026): ese tramo cuenta como tiempo observado para crisis aunque no haya EEG, y ' ...
        'las %d crisis vistas en él (%s) cuentan en `n_seizures` (como `n_missed_DX`), en la tasa y en la ventana 4-OHT. El resto de la ' ...
        'desconexión, hasta el final del EDF (DX-26), no se revisó y no cuenta. Esas crisis no tienen figura (no hay EEG).'], ...
        char(datetime(rv.dx.end_posix(d23), 'ConvertFrom', 'posixtime', 'TimeZone', cfg.tz, 'Format', 'HH:mm:ss')), height(nv), ...
        strjoin(nv.original_id', ', ')), ...
        sprintf(['- Desconexiones (hoja `2_Desconexiones`, actualizada el 05/11/2026): todas revisadas en video salvo %s. El tiempo de una ' ...
        'desconexión revisada cuenta como observado para crisis (`observed_h`); para los IIDs todas se excluyen. DX-21 (128-s) está revisada ' ...
        'pero cae entre los EDFs `_128s_A1` y `_128s_A2`, fuera del EEG analizado, y no cuenta (decisión del 05/11/2026).'], ...
        strjoin(arrayfun(@(i) sprintf('%s (%s)', char(rv.dx.id(i)), char(rv.dx.mouse(i))), nrv', 'UniformOutput', false), ', ')), ...
        ['- 113-s, EDF 20260807~: DX-24 (desde el 09/08 19:19:59 hasta el fin del registro) y, en 20260809~, DX-25 (desde el ' ...
        'inicio hasta las 19:27:26), fuera de los IIDs (revisadas en video: observadas para crisis). ND-10 y ND-11 = crisis no detectadas ' ...
        'por la atenuación de la señal entre las detecciones 32 y 33; ese tramo atenuado sí cuenta para los IIDs (decisión del 04/11/2026), ' ...
        'aunque puede subestimarlos.']}];
    k12 = find(strcmp(SZ.mouse, '072-s') & strcmp(SZ.original_id, 'ND-12'));
    if isscalar(k12)   % 2026-11-05
        nx = SZ(strcmp(SZ.mouse, '072-s') & strcmp(SZ.edf, SZ.edf{k12}) & SZ.seizure_id > SZ.seizure_id(k12), :);
        x.scope{end+1} = sprintf(['- 072-s, EDF %s: ND-12 = crisis vista en video durante DX-07 (%s–%s; nota del 05/11/2026). Es la crisis %d ' ...
            'de ese EDF: las %d siguientes (%s) subieron un número de `seizure_id` (figuras redibujadas); su `original_id` no cambia.'], ...
            regexp(SZ.edf{k12}, '^[^~]+~', 'match', 'once'), SZ.onset{k12}(1:19), SZ.offset{k12}(12:19), SZ.seizure_id(k12), height(nx), ...
            strjoin(nx.original_id', ', '));
    end
    if ~strcmp(cfg.tz, 'Europe/Paris'), x.scope{end+1} = sprintf('- Zona horaria: %s.', cfg.tz); end
end

function L = update_notes(upd, CHK, SZfig, n_kept, n_iid_same, rule, iid_note, all055)
% README section of the incremental update (the earlier updates' sections are kept by write_readme)
    what = sprintf('se consolidaron %s', strjoin(upd.mice, ', '));
    if ~isempty(upd.desc), what = sprintf('%s; se re-consolidaron %s', upd.desc, strjoin(upd.mice, ', ')); end
    L = {sprintf(['Con `consolidate_update = ''compute''` y luego `''write''` (`consolidate_results.m`), el %s: %s; el resto de ' ...
        'los ratones conserva sus filas de la consolidación anterior.'], char(datetime('now', 'Format', 'yyyy-MM-dd')), what)};
    L{end+1} = sprintf(['- Verificado contra la consolidación anterior (`_work/update_%s/before/_work/state.mat`): las %d crisis, las %d filas ' ...
        'per_EDF y las zonas IID de los %d canales de los demás ratones son idénticas; en las hojas reescritas sus %d filas quedan iguales ' ...
        'celda a celda y las demás hojas no cambian.'], upd.tag, CHK.n_seizures, CHK.n_edf_rows, CHK.n_iid_channels, n_kept);
    nr = height(CHK.renumbered);
    if (CHK.n_earlier > 0 || n_iid_same > 0) && nr == 0
        L{end+1} = sprintf(['- Lo ya consolidado de estos ratones: sus %d crisis con los mismos `seizure_id`, horarios y figuras; sus %d canales ' ...
            'IID con los mismos valores.'], CHK.n_earlier, n_iid_same);
    elseif CHK.n_earlier > 0 || n_iid_same > 0
        L{end+1} = sprintf(['- Lo ya consolidado de estos ratones: sus %d crisis con los mismos horarios y `original_id`, y el mismo `seizure_id` y ' ...
            'figura salvo %d que suben de número por una crisis nueva anterior en su EDF (figuras redibujadas; las anteriores en ' ...
            '`_work/update_%s/before/figures/`); sus %d canales IID con los mismos valores.'], CHK.n_earlier, nr, upd.tag, n_iid_same);
        for i = 1:nr
            L{end+1} = sprintf('  - %s, %s (`original_id` %s): `seizure_id` %d → %d', CHK.renumbered.mouse{i}, CHK.renumbered.edf{i}, ...
                CHK.renumbered.original_id{i}, CHK.renumbered.id_before(i), CHK.renumbered.id_after(i)); %#ok<AGROW>
        end
    end
    for i = 1:numel(CHK.changed), L{end+1} = ['  - cambio: ' CHK.changed{i}]; end %#ok<AGROW>
    nw = 'ninguna';
    if height(CHK.new) > 0
        nw = strjoin(arrayfun(@(i) sprintf('%s %d', CHK.new.mouse{i}, CHK.new.n(i)), (1:height(CHK.new))', 'UniformOutput', false)', ', ');
    end
    L{end+1} = sprintf('- Crisis nuevas: %s; figuras dibujadas: %d.', nw, height(SZfig));
    L{end+1} = ['- ' iid_note];
    if ismember('055-s', upd.mice) && all055
        L{end+1} = sprintf(['- 055-s, IIDs: se usó el caso %d. Los casos 1, 2 y 3 y la regla general se calcularon todos (resultados por ' ...
            'canal en `_work/update_%s/iid_055s/`); los del caso usado, también en `_work/stats/iid_rerun_final/055-s/`.'], rule, upd.tag);
    elseif ismember('055-s', upd.mice)
        L{end+1} = sprintf('- 055-s, IIDs: caso %d (decidido el 04/11/2026).', rule);
    end
    L{end+1} = sprintf('- Archivos anteriores a la actualización: `_work/update_%s/before/`.', upd.tag);
end

%% ---------------------------------------------------------------------------
function fig_worker(q, k, n, out, work)
% worker k of n: one EDF per job (both channels), one .fig + .png per seizure
    jd = fullfile(work, 'fig_jobs');
    if ~isfolder(jd), mkdir(jd); end
    jk = strcat(q.SZ.fig_campaign, '|', q.SZ.mouse, '|', q.SZ.edf);
    jobs = unique(jk);
    for j = find(mod((1:numel(jobs)) - 1, n) + 1 == k)
        jf = fullfile(jd, sprintf('job_%03d.mat', j));
        sel = find(strcmp(jk, jobs{j}));
        t = tic;
        res = render_edf(q, sel, out);
        res.job = jobs{j};
        save(jf, 'res');
        fprintf('[%s] %s: %d figure(s) in %.1f min\n', char(datetime('now', 'Format', 'HH:mm')), jobs{j}, numel(sel), toc(t) / 60);
    end
end

function res = render_edf(q, sel, out)
% traces of both channels of one EDF (pipeline functions, the EDF's own
% config_used.mat); per seizure and channel the energy metrics and ll ratio
% between onset and offset, and the +/- 15 s window for the figure
    win_s = 15;
    sz = q.SZ(sel, :);
    e = q.edf_tbl(strcmp(q.edf_tbl.campaign, sz.fig_campaign{1}) & strcmp(q.edf_tbl.mouse, sz.mouse{1}) & strcmp(q.edf_tbl.edf, sz.edf{1}), :);
    regions = {'HPCl', 'HPCr'};   % left column first, as the pipeline's joint figures
    seg = cell(height(sz), 2);
    info = repmat(struct('label', '', 'note', ''), 1, 2);
    res.check_peak = 0; res.check_ll = 0; res.check_n = 0;
    pcfg = [];
    for c = 1:2
        p = '';
        if height(e) == 1, p = e.(['clean_' regions{c}]){1}; end
        if isempty(p) || ~isfile(p)
            info(c).note = 'no 02_clean file';
            continue;
        end
        try
            if isempty(pcfg)
                w = load(fullfile(fileparts(fileparts(p)), 'config_used.mat'), 'cfg');
                pcfg = w.cfg;
            end
            data = load_lfp_txt(p);
            tr = seizure_energy_trace(data, pcfg);
            if isempty(tr.blocks)
                trimmed = zeros(0, 2);
            else
                trimmed = [vertcat(tr.blocks.trimmed_start), vertcat(tr.blocks.trimmed_end)];
            end
            [ll_full, ll_med] = line_length_trace(tr.bp_full, trimmed, data.fs, pcfg);
        catch ME
            info(c).note = ME.message;
            continue;
        end
        if isfield(data.meta, 'channel_label'), info(c).label = strtrim(data.meta.channel_label); end
        vs = data.signal(data.valid_mask(:));
        span_uV = max(vs) - min(vs);   % seizure_energy_trace normalizes by (max - min) of the valid samples
        t = data.t_rel(:);
        t0 = e.start_posix(1);   % the clean header's session start, as the seizure times (cs_load_edf_table)
        for i = 1:height(sz)
            on_s = sz.start_posix(i) - t0; off_s = sz.end_posix(i) - t0;
            ix = t >= on_s - win_s & t <= off_s + win_s;   % figure window
            iz = t >= on_s & t <= off_s;                   % the seizure: every metric
            es = tr.energy_full(iz);
            m = struct('peak', max(es, [], 'omitnan') / tr.threshold, 'mean', mean(es, 'omitnan') / tr.threshold, ...
                'auc', sum(es, 'omitnan') / data.fs / tr.threshold, 'll', median(ll_full(iz), 'omitnan') / ll_med);
            if all(isnan(es)), m.auc = NaN; end
            seg{i, c} = struct('t', t(ix), 'signal', data.signal(ix), 'bp_uV', tr.bp_full(ix) * span_uV, ...
                'energy', tr.energy_full(ix), 'll', ll_full(ix) / ll_med, 'threshold', tr.threshold, ...
                'll_threshold', pcfg.seizure_robust.ll_threshold, 'on', on_s, 'off', off_s, 'm', m);
        end
        % the recomputed traces must give the pipeline's own peak_energy_ratio and ll_ratio on its accepted rows
        d = q.det_rows(strcmp(q.det_rows.campaign, sz.fig_campaign{1}) & strcmp(q.det_rows.mouse, sz.mouse{1}) & ...
            strcmp(q.det_rows.edf, sz.edf{1}) & strcmp(q.det_rows.region, regions{c}) & strcmp(q.det_rows.detection_status, 'accepted'), :);
        for r = 1:height(d)
            w = t >= d.start_s(r) & t <= d.end_s(r);
            v = max(tr.energy_full(w), [], 'omitnan') / tr.threshold;
            res.check_peak = max(res.check_peak, abs(v - d.peak_energy_ratio(r)) / abs(d.peak_energy_ratio(r)));
            v = median(ll_full(w), 'omitnan') / ll_med;
            res.check_ll = max(res.check_ll, abs(v - d.ll_ratio(r)) / abs(d.ll_ratio(r)));
            res.check_n = res.check_n + 1;
        end
        clear data tr ll_full t vs
    end
    n = height(sz);
    res.figure = sz.figure; res.ok = false(n, 1); res.note = repmat({''}, n, 1);
    res.peak_energy_ratio = NaN(n, 2); res.mean_to_thr = NaN(n, 2); res.auc_to_thr_s = NaN(n, 2); res.ll_ratio = NaN(n, 2);
    for i = 1:n
        for c = 1:2
            if ~isempty(seg{i, c})
                res.peak_energy_ratio(i, c) = seg{i, c}.m.peak; res.mean_to_thr(i, c) = seg{i, c}.m.mean;
                res.auc_to_thr_s(i, c) = seg{i, c}.m.auc; res.ll_ratio(i, c) = seg{i, c}.m.ll;
            end
        end
        try
            draw_seizure(sz(i, :), seg(i, :), info, regions, pcfg, win_s, fullfile(out, 'figures', sz.figure{i}));
            res.ok(i) = true;
            res.note{i} = strjoin({info(~cellfun(@isempty, {info.note})).note}, ' | ');
        catch ME
            res.note{i} = ME.message;
        end
    end
end

function draw_seizure(s, seg, info, regions, pcfg, win_s, base)
% 4 x 2: rows = cleaned LFP, band-passed LFP in uV, energy, line length
% (pipeline row functions plot_trace_voltage / _energy / _segment), columns =
% HPCl, HPCr; x = time from the window start (onset at win_s s), linked;
% onset and offset marked; same y-limits across columns within each row
    fig = figure('Visible', 'off', 'Position', [50, 50, 1600, 1100]);
    ax = gobjects(4, 2); ylo = NaN(4, 2); yhi = NaN(4, 2);
    win = [0, s.duration_s + 2 * win_s];
    for c = 1:2
        for r = 1:4, ax(r, c) = subplot(4, 2, (r - 1) * 2 + c); end
        g = seg{c};
        if isempty(g)
            for r = 1:4, text(ax(r, c), 0.5, 0.5, ['no signal: ' info(c).note], 'HorizontalAlignment', 'center', 'Units', 'normalized'); end
            title(ax(1, c), regions{c}, 'Interpreter', 'none');
            continue;
        end
        t = g.t - (g.on - win_s);
        on = win_s; off = g.off - g.on + win_s;
        y = cell(1, 4);
        y{1} = plot_trace_voltage(ax(1, c), t, g.signal, win);
        y{2} = plot_trace_voltage(ax(2, c), t, g.bp_uV, win);
        ylabel(ax(1, c), ['Voltage (' char(181) 'V)']); ylabel(ax(2, c), ['Voltage (' char(181) 'V)']);
        y{3} = plot_trace_energy(ax(3, c), t, g.energy, win, g.threshold);
        hold(ax(4, c), 'on');
        [a, b] = plot_trace_segment(ax(4, c), t, g.ll, win, [0.6 0.3 0.6]);
        yline(ax(4, c), g.ll_threshold, 'r--', sprintf('Threshold (%.2f)', g.ll_threshold), 'LineWidth', 1.5);
        ylabel(ax(4, c), 'll / ll\_median\_global'); grid(ax(4, c), 'on');
        y{4} = [min(a, g.ll_threshold), max(b, g.ll_threshold)];
        for r = 1:4, ylo(r, c) = y{r}(1); yhi(r, c) = y{r}(2); end
        title(ax(1, c), sprintf('%s (%s)', regions{c}, info(c).label), 'Interpreter', 'none');
        title(ax(2, c), sprintf('Band-passed [%g-%g Hz]', pcfg.seizure.bandpass_band(1), pcfg.seizure.bandpass_band(2)), 'Interpreter', 'none');
        title(ax(3, c), sprintf('Energy Metric: peak/thr=%.1f | mean/thr=%.1f | auc/thr=%.0f s', g.m.peak, g.m.mean, g.m.auc), 'Interpreter', 'none');
        title(ax(4, c), sprintf('Line-Length Ratio: ll=%.2f', g.m.ll), 'Interpreter', 'none');
        for r = 1:4
            xline(ax(r, c), on, '-', 'onset', 'Color', [0.85 0 0], 'LineWidth', 1.5);
            xline(ax(r, c), off, '-', 'offset', 'Color', [0.85 0 0], 'LineWidth', 1.5);
            xlabel(ax(r, c), 'Time (s)');
        end
    end
    linkaxes(ax(:), 'x');
    xlim(ax(1, 1), win);
    for r = 1:4
        lo = min(ylo(r, :), [], 'omitnan'); hi = max(yhi(r, :), [], 'omitnan');
        if isfinite(lo) && isfinite(hi) && hi > lo
            for c = 1:2, if ~isempty(seg{c}), ylim(ax(r, c), [lo hi]); end, end
        end
    end
    sid = sprintf('seizure %d', s.seizure_id);
    if strcmp(s.source{1}, 'missed_by_pipeline'), sid = [sid ' (missed)']; end
    sgtitle(fig, sprintf('%s | %s | %s | %s | onset %s - offset %s (%.1f s)', s.mouse{1}, s.edf_file{1}, sid, ...
        s.light_phase{1}, s.onset{1}, s.offset{1}, s.duration_s), 'Interpreter', 'none', 'FontSize', 11);
    if ~isfolder(fileparts(base)), mkdir(fileparts(base)); end
    exportgraphics(fig, [base '.png'], 'Resolution', 150);
    set(fig, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');   % drawn invisible (batch); opens visible
    savefig(fig, [base '.fig']);
    close(fig);
end

function FJ = fig_results(SZ, work, out)
% per seizure (SZ order): both files written?, note, metrics per channel;
% the trace checks over all jobs. Jobs in name order (job_NNN, then job_u<tag>_NNN by update): a figure
% redrawn by a later update (renumbered seizure) takes that update's result
    n = height(SZ);
    FJ.regions = {'HPCl', 'HPCr'};
    FJ.ok = false(n, 1); FJ.note = repmat({'no worker result'}, n, 1); FJ.job = repmat({''}, n, 1);
    for m = {'peak_energy_ratio', 'mean_to_thr', 'auc_to_thr_s', 'll_ratio'}, FJ.(m{1}) = NaN(n, 2); end
    FJ.check_peak = 0; FJ.check_ll = 0; FJ.check_n = 0;
    fl = dir(fullfile(work, 'fig_jobs', 'job_*.mat'));
    [~, o] = sort({fl.name});
    for f = fl(o)'
        q = load(fullfile(f.folder, f.name), 'res');
        r = q.res;
        [tf, j] = ismember(r.figure, SZ.figure);
        both = cellfun(@(x) isfile(fullfile(out, 'figures', [x '.fig'])) && isfile(fullfile(out, 'figures', [x '.png'])), r.figure(tf));
        FJ.ok(j(tf)) = r.ok(tf) & both;
        FJ.note(j(tf)) = r.note(tf);
        FJ.job(j(tf)) = {r.job};
        for m = {'peak_energy_ratio', 'mean_to_thr', 'auc_to_thr_s', 'll_ratio'}, FJ.(m{1})(j(tf), :) = r.(m{1})(tf, :); end
        FJ.check_peak = max(FJ.check_peak, r.check_peak); FJ.check_ll = max(FJ.check_ll, r.check_ll);
        FJ.check_n = FJ.check_n + r.check_n;
    end
end

%% ---------------------------------------------------------------------------
function [SE, SM, n_out_obs] = seizure_summary(S, O, SZ)
% per EDF and per mouse: observed time = T_seizure_obs (recording minus
% disconnections not reviewed on video, run_campaign_stats), seizures, rate,
% durations, TP / FP / missed (outside and during disconnections)
    om = O.mouse(strcmp({O.mouse.scope}, 'merged'));
    in_obs = false(height(SZ), 1);
    for i = 1:height(SZ)
        T = om(strcmp({om.mouse}, SZ.mouse{i})).T_sz;
        in_obs(i) = any(SZ.start_posix(i) >= T(:, 1) & SZ.start_posix(i) < T(:, 2));
    end
    n_out_obs = sum(~in_obs);
    ed = O.edf_day;
    [g, E] = findgroups(ed(:, {'mouse', 'group', 'campaign', 'edf'}));
    E.observed_h = splitapply(@sum, ed.T_seizure_obs_h, g);
    fp = S.fp_events;
    SE = summary_rows(E, @(r) strcmp(SZ.mouse, r.mouse{1}) & strcmp(SZ.edf, r.edf{1}), ...
        @(r) sum(strcmp(fp.mouse, r.mouse{1}) & strcmp(fp.edf, r.edf{1}) & strcmp(fp.campaign, r.campaign{1})), SZ, in_obs);
    SE.edf_file = strcat(SE.edf, '.EDF');
    SE = SE(:, [{'mouse', 'group', 'edf_file'}, SE.Properties.VariableNames(5:end - 1)]);
    M = unique(ed(:, {'mouse', 'group'}));
    M.n_edfs = arrayfun(@(i) height(unique(ed(strcmp(ed.mouse, M.mouse{i}), {'campaign', 'edf'}))), (1:height(M))');
    M.observed_h = arrayfun(@(i) cs_iv('length', om(strcmp({om.mouse}, M.mouse{i})).T_sz) / 3600, (1:height(M))');
    SM = summary_rows(M, @(r) strcmp(SZ.mouse, r.mouse{1}), @(r) sum(strcmp(fp.mouse, r.mouse{1})), SZ, in_obs);
    for T = {SE, SM}
        bad = T{1}.n_seizures ~= T{1}.n_TP + T{1}.n_missed + T{1}.n_missed_DX;
        assert(~any(bad), 'seizure_summary: n_seizures ~= n_TP + n_missed + n_missed_DX in %d row(s)', sum(bad));
    end
end

function T = summary_rows(T, pick, nfp, SZ, in_obs)
    n = height(T);
    z = NaN(n, 1);
    T.n_seizures = zeros(n, 1); T.seizures_per_24h = z;
    T.dur_mean_s = z; T.dur_median_s = z; T.dur_min_s = z; T.dur_max_s = z;
    T.n_TP = zeros(n, 1); T.n_FP = zeros(n, 1); T.n_missed = zeros(n, 1); T.n_missed_DX = zeros(n, 1);
    tp = strcmp(SZ.source, 'pipeline_TP');
    missed = strcmp(SZ.source, 'missed_by_pipeline');
    in_dx = ~cellfun(@isempty, SZ.in_dx);
    for i = 1:n
        s = pick(T(i, :));
        d = SZ.duration_s(s);
        T.n_seizures(i) = sum(s);
        if T.observed_h(i) > 0, T.seizures_per_24h(i) = sum(s & in_obs) / T.observed_h(i) * 24; end
        if ~isempty(d)
            T.dur_mean_s(i) = mean(d); T.dur_median_s(i) = median(d); T.dur_min_s(i) = min(d); T.dur_max_s(i) = max(d);
        end
        T.n_TP(i) = sum(s & tp);
        T.n_FP(i) = nfp(T(i, :));
        T.n_missed(i) = sum(s & missed & ~in_dx);
        T.n_missed_DX(i) = sum(s & missed & in_dx);
    end
end

%% ---------------------------------------------------------------------------
function [IE, IM] = iid_summary(cfg, I, K, C, by)
% per mouse x EDF x hemisphere (I.channels, re-detection) + light / dark
% from the analyzed intervals and complex onsets (as cs_iid_stats' circadian
% part); per mouse x hemisphere (2026-11-04) = pooled over that hemisphere's
% EDFs (counts / analyzed time), the hemispheres never summed; by = 'mouse'
% pools both hemispheres instead (the per_mouse sheet before, for checks)
    if nargin < 5, by = 'hemisphere'; end
    on = hm(cfg.lights_on); off = hm(cfg.lights_off);
    ch = I.channels;
    key = @(t) strcat(t.campaign, '|', t.mouse, '|', t.edf, '|', t.region);
    kc = key(ch);
    [tf, j] = ismember(kc, key(K));
    assert(all(tf), 'iid_summary: %d channel(s) of I.channels without intervals', sum(~tf));
    Cm = C(strcmp(C.variant, 'principal'), :);
    [~, jc] = ismember(kc, key(Cm));
    n = height(ch);
    L = zeros(n, 4); amp = cell(n, 1);   % light_s, n_light, dark_s, n_dark
    for i = 1:n
        k = j(i);
        p = cs_clock_split(K.inc{k}, cfg.tz, [0 on off 1440]);
        nd = 0; nl = 0;
        if ~isempty(K.on_abs{k})
            dt = datetime(K.on_abs{k}, 'ConvertFrom', 'posixtime', 'TimeZone', cfg.tz);
            m = hour(dt) * 60 + minute(dt) + second(dt) / 60;
            nd = sum(m < on | m >= off); nl = numel(m) - nd;
        end
        L(i, :) = [p(2), nl, p(1) + p(3), nd];
        assert(nl + nd == ch.n_complexes(i), 'iid_summary: complex count mismatch in %s', kc{i});
        amp{i} = Cm.complexes{jc(i)}(:, 5);
    end
    IE = table(ch.mouse, ch.group, strcat(ch.edf, '.EDF'), ch.region, ch.T_iid_h * 60, ch.n_complexes, ...
        'VariableNames', {'mouse', 'group', 'edf_file', 'hemisphere', 'observed_min', 'n_IID'});
    IE = add_rates(IE, L, ch.n_single, ch.n_polyspike, ch.n_bursts, amp);
    IE.baseline_uV = ch.baseline_uV; IE.lower_threshold_uV = ch.lower_uV; IE.p97_real_uV = ch.p97_real_uV;
    if strcmp(by, 'hemisphere')   % rows: mouse, then HPCl before HPCr
        [g, IM] = findgroups(ch(:, {'mouse', 'group', 'region'}));
        IM = renamevars(IM, 'region', 'hemisphere');
        IM.n_edfs = splitapply(@numel, ch.n_complexes, g);
        kim = strcat(IM.mouse, '|', IM.hemisphere);
    else
        [g, IM] = findgroups(ch(:, {'mouse', 'group'}));
        IM.n_channels = splitapply(@numel, ch.n_complexes, g);
        kim = strcat(IM.mouse, '|ambos');
    end
    IM.observed_min = splitapply(@sum, ch.T_iid_h * 60, g);
    IM.n_IID = splitapply(@sum, ch.n_complexes, g);
    Lm = [splitapply(@sum, L(:, 1), g), splitapply(@sum, L(:, 2), g), splitapply(@sum, L(:, 3), g), splitapply(@sum, L(:, 4), g)];
    am = arrayfun(@(k) vertcat(amp{g == k}), (1:max(g))', 'UniformOutput', false);
    IM = add_rates(IM, Lm, splitapply(@sum, ch.n_single, g), splitapply(@sum, ch.n_polyspike, g), splitapply(@sum, ch.n_bursts, g), am);
    % pooled rate = run_campaign_stats' I.mouse (merged; per hemisphere, or both)
    ref = I.mouse(strcmp(I.mouse.scope, 'merged'), :);
    [tf, j] = ismember(kim, strcat(ref.mouse, '|', ref.hemisphere));
    assert(all(tf) && max(abs(IM.IID_per_min - ref.complexes_per_min(j))) < 1e-9, 'iid_summary: per-mouse rate differs from I.mouse');
end

function T = add_rates(T, L, n_single, n_poly, n_bursts, amp)
    T.IID_per_min = T.n_IID ./ T.observed_min;
    T.observed_min_light = L(:, 1) / 60; T.n_IID_light = L(:, 2); T.IID_per_min_light = L(:, 2) ./ (L(:, 1) / 60);
    T.observed_min_dark = L(:, 3) / 60; T.n_IID_dark = L(:, 4); T.IID_per_min_dark = L(:, 4) ./ (L(:, 3) / 60);
    T.n_single = n_single; T.n_polyspike = n_poly; T.pct_polyspike = 100 * n_poly ./ T.n_IID;
    T.n_bursts = n_bursts; T.bursts_per_h = n_bursts ./ (T.observed_min / 60);
    T.median_amplitude_mV = cellfun(@(a) median(a, 'omitnan'), amp);
end

%% ---------------------------------------------------------------------------
function CMP = compare_20261101(cfg, SM, IM)
% per mouse vs the last run_campaign_stats run (2026-11-01); IIDs per hemisphere (IM: one row per
% mouse x hemisphere)
    d = fullfile(cfg.batch_root, '_stats_revision_20261101_final');
    q = load(fullfile(d, 'phase1_results.mat'), 'R');
    w = load(fullfile(d, 'phase3_results.mat'), 'I');
    r = q.R.mouse(strcmp(q.R.mouse.scope, 'merged'), :);
    im = w.I.mouse(strcmp(w.I.mouse.scope, 'merged'), :);
    hs = {'HPCl', 'HPCr'};
    CMP = SM(:, {'mouse', 'n_seizures', 'observed_h', 'seizures_per_24h'});
    CMP.n_seizures_0111 = NaN(height(CMP), 1); CMP.observed_h_0111 = CMP.n_seizures_0111;
    CMP.seizures_per_24h_0111 = CMP.n_seizures_0111;
    for h = hs
        CMP.(['IID_per_min_' h{1} '_0111']) = NaN(height(CMP), 1);
        CMP.(['IID_per_min_' h{1}]) = NaN(height(CMP), 1);
    end
    CMP.diferencia = repmat({''}, height(CMP), 1);
    for i = 1:height(CMP)
        m = CMP.mouse{i};
        for h = hs
            kn = strcmp(IM.mouse, m) & strcmp(IM.hemisphere, h{1});
            if any(kn), CMP.(['IID_per_min_' h{1}])(i) = IM.IID_per_min(kn); end
            ki = strcmp(im.mouse, m) & strcmp(im.hemisphere, h{1});
            if any(ki), CMP.(['IID_per_min_' h{1} '_0111'])(i) = im.complexes_per_min(ki); end
        end
        k = strcmp(r.mouse, m);
        if ~any(k)
            CMP.diferencia{i} = 'no estaba en la corrida del 01/11';
            continue;
        end
        CMP.n_seizures_0111(i) = r.n_total(k) + 0; CMP.observed_h_0111(i) = r.T_seizure_obs_h(k);
        CMP.seizures_per_24h_0111(i) = r.rate_per_24h(k);
        dif = {};
        if CMP.n_seizures(i) ~= CMP.n_seizures_0111(i), dif{end+1} = sprintf('crisis %d -> %d', CMP.n_seizures_0111(i), CMP.n_seizures(i)); end %#ok<AGROW>
        if abs(CMP.observed_h(i) - CMP.observed_h_0111(i)) > 1e-6, dif{end+1} = sprintf('observado %.2f -> %.2f h', CMP.observed_h_0111(i), CMP.observed_h(i)); end %#ok<AGROW>
        for h = hs
            a = CMP.(['IID_per_min_' h{1} '_0111'])(i); b = CMP.(['IID_per_min_' h{1}])(i);
            if abs(b - a) > 1e-9 || isnan(a) ~= isnan(b)
                dif{end+1} = sprintf('IIDs/min %s %.4f -> %.4f', h{1}, a, b); %#ok<AGROW>
            end
        end
        CMP.diferencia{i} = strjoin(dif, '; ');
    end
end

%% ---------------------------------------------------------------------------
function write_readme(file, cfg, S, rv, SZ, FJ, CMP, iid_note, worker_errors, snapinfo, n_out_obs, extra)
    L = {};
    a = @(varargin) sprintf(varargin{:});
    L{end+1} = '# Resultados consolidados: crisis e IIDs, campañas C1 y C2 (2026)';
    L{end+1} = '';
    L{end+1} = a('Generado por `epilepsy_app_natus/consolidate_results.m` el %s (git HEAD %s, árbol con cambios sin commit: %s). No se re-corrió la extracción, la limpieza ni la detección de crisis.', ...
        char(datetime('now', 'TimeZone', cfg.tz, 'Format', 'yyyy-MM-dd HH:mm')), cfg.git.head(1:min(7, end)), mat2str(cfg.git.dirty));
    L{end+1} = '';
    L{end+1} = '## Archivos';
    L{end+1} = ['- `seizures_all.xlsx` (hoja `seizures`): una fila por crisis final. `seizure_id` = número de la crisis dentro de su EDF (1..N por onset, ' ...
        'con las fusiones aplicadas; es el `NN` de las figuras y el último campo de `event_id` = `<ratón>|<EDF>|<seizure_id>`). `original_id` = numeración ' ...
        'del Excel de revisión (`1+2` si se fusionaron) o el ID de la crisis no detectada (`ND-03`). `onset`/`offset` en hora local (Europe/Paris), como en Natus. ' ...
        '`peak_energy_ratio`, `mean_to_thr` y `auc_to_thr_s` por hemisferio: sobre la traza de energía del pipeline entre onset y offset (sin la ventana de ' ...
        '±15 s): máximo / umbral, media / umbral y suma / fs / umbral (s). `source` = `pipeline_TP` (detectadas y etiquetadas TP, más las crisette focal) o ' ...
        '`missed_by_pipeline` (hoja 1_Crisis_no_detectadas). `comments` en inglés.'];
    L{end+1} = ['- `seizures_summary.xlsx`: hojas `per_EDF` y `per_mouse`. `observed_h` = tiempo observado para crisis (registro menos desconexiones; ' ...
        'las desconexiones revisadas en video cuentan como observadas, igual que run_campaign_stats, y en 055-s también el video revisado después ' ...
        'del fin de su EEG, DX-23). `seizures_per_24h` = crisis con onset en ese tiempo / ' ...
        'observed_h × 24. Duraciones sobre todas las crisis (EEG y video). `n_seizures` = `n_TP` (pipeline, con crisette, después de las fusiones) + ' ...
        '`n_missed` (no detectadas fuera de desconexiones) + `n_missed_DX` (no detectadas durante desconexiones, vistas en video), verificado fila a fila. ' ...
        '`n_FP` = detecciones etiquetadas FP.'];
    L{end+1} = ['- `iid_summary.xlsx`: hojas `per_mouse_EDF_hemisphere` y `per_mouse` (desde el 04/11/2026 una fila por ratón y hemisferio, HPCl y ' ...
        'HPCr: los EDFs del hemisferio agrupados, suma de IIDs / suma de su tiempo analizado; los hemisferios no se suman, un IID presente en ambos ' ...
        'cuenta una vez en cada uno; `n_edfs` = EDFs del hemisferio). ' ...
        'IIDs de la re-detección de run_campaign_stats. `observed_min` = tiempo analizado (excluye crisis TP, FP que no son "IIDs", desconexiones y sus ' ...
        'márgenes de ±5 s). Luz = 08:00–19:30. Amplitud: mediana de la amplitud máxima de los complejos; umbrales y P97 por canal. No existen métricas ' ...
        'de intervalo entre IIDs en los outputs.'];
    L{end+1} = ['- `figures/<raton>/<raton>_<fecha-EDF>_sz<NN>_<TP|missed>.fig` y `.png`: una por crisis, `NN` = `seizure_id` (si dos EDFs del ratón ' ...
        'del mismo día tienen crisis, la fecha va seguida del primer bloque del ID del EDF: `016-s_20260327_2884b4a6_sz01_TP`). Columnas HPCl | HPCr; ' ...
        'filas: LFP limpio, LFP filtrado 5–75 Hz en µV (la traza normalizada del pipeline × (máx − mín) del canal), energía y line length. Eje x: ' ...
        'tiempo desde el inicio de la ventana (onset en 15 s, ventana hasta 15 s después del offset), igual en todos los paneles. Títulos: ' ...
        'peak/thr, mean/thr y auc/thr (las métricas de seizures_all) y `ll` = mediana del line length entre onset y offset / mediana global (la fórmula ' ...
        'de `ll_ratio` del pipeline); umbral de line length = `seizure_robust.ll_threshold`.'];
    snapdir = '_work/stats/labels_snapshot_raw';
    if isfield(snapinfo, 'dir'), snapdir = snapinfo.dir; end
    L{end+1} = ['- `_work/stats/labels_snapshot/<raton>_seizures_validated.xlsx`: las crisis finales de cada ratón, mismas filas y columnas que ' ...
        'seizures_all.xlsx. `' snapdir '/`: copia (con SHA-1) de los `seizures_events.csv` de revisión tal como se leyeron.'];
    L{end+1} = '- `_work/`: además, estado intermedio (`state.mat`), resultados de la re-detección de IIDs y logs de los procesos.';
    L{end+1} = '';
    L{end+1} = '## Alcance';
    inc = S.units(S.units.included, :);
    for gname = {'control', 'epileptic'}
        ms = unique(inc.mouse(strcmp(inc.group, gname{1})));
        L{end+1} = a('- %s (%d): %s', gname{1}, numel(ms), strjoin(ms', ', ')); %#ok<AGROW>
    end
    L{end+1} = '- EDFs excluidos:';
    for i = 1:height(rv.excluded)
        L{end+1} = a('  - %s %s %s: %s (%s)', rv.excluded.id{i}, rv.excluded.mouse{i}, rv.excluded.edf{i}, rv.excluded.kind{i}, rv.excluded.reason{i}); %#ok<AGROW>
    end
    L{end+1} = '- 128-s, EDF reprocesado del 29/07 (A1 y A2): incluido. No tiene detecciones; la desconexión DX-21 sigue fuera.';
    L{end+1} = '- 077-s (C2): etiquetas de su `seizures_events.csv` actual (no de la copia `seizures_events_excel_20261019.csv`).';
    L = [L, extra.scope];
    L{end+1} = a('- Etiquetas leídas de una copia tomada el %s (%d archivos).', snapinfo.taken, snapinfo.n);
    L{end+1} = '';
    L{end+1} = '## Criterios';
    L{end+1} = ['- Crisis = eventos con AE = TP (un evento = mismo `seizure_id` del pipeline en HPCr y HPCl; basta TP en uno), más las crisis no detectadas ' ...
        'ND-01..ND-06 y ND-08..ND-12 (ND-02..04 decididas en la fase 1; ND-05, ND-06 y ND-12 vistas solo en video durante DX-06 y DX-07; ND-08 y ND-09 ' ...
        'vistas solo en video durante DX-23, sin EEG; ND-10 y ND-11 no detectadas por atenuación de la señal) y las 3 crisette focal (O-07, O-08, ' ...
        'O-13). ND-07 no tiene horario y queda fuera.'];
    L{end+1} = ['- IIDs: los FP etiquetados "IIDs" cuentan como IID; el resto de FP, las crisis TP y las desconexiones se excluyen (zonas de run_campaign_stats). ' ...
        'Las crisette se excluyen de los IIDs como en run_campaign_stats. Un FP con AF distinto entre hemisferios se excluye en ambos, salvo en 055-s: ' extra.af055 '.'];
    L{end+1} = '- Fusiones (onset = mínimo, offset = máximo; numeración del Excel de revisión):';
    for k = 1:size(cfg.merge_events, 1)
        c = cfg.merge_events(k, :);
        L{end+1} = a('  - %s: %s %s %s, crisis %s', c{1}, c{2}, c{3}, c{4}, strjoin(arrayfun(@num2str, c{5}, 'UniformOutput', false), '+')); %#ok<AGROW>
    end
    L{end+1} = '';
    L{end+1} = '## Números';
    tp = strcmp(SZ.source, 'pipeline_TP'); ms = ~tp; dx = ~cellfun(@isempty, SZ.in_dx);
    L{end+1} = a('- %d crisis: n_TP %d, n_missed %d, n_missed_DX %d.', height(SZ), sum(tp), sum(ms & ~dx), sum(ms & dx));
    L{end+1} = a('- Crisis con onset fuera del tiempo observado (no entran en la tasa): %d.', n_out_obs);
    L{end+1} = a('- Figuras: %d de %d (.fig + .png).', sum(FJ.ok), height(SZ));
    L{end+1} = ['- ' iid_note];
    L{end+1} = a(['- Control de las trazas recalculadas: en las %d filas aceptadas del pipeline reproducen su `peak_energy_ratio` (máx. diferencia relativa ' ...
        '%.1e) y su `ll_ratio` (%.1e).'], FJ.check_n, FJ.check_peak, FJ.check_ll);
    for i = find(~FJ.ok)'
        L{end+1} = a('- Sin figura: %s (%s)', SZ.figure{i}, FJ.note{i}); %#ok<AGROW>
    end
    for i = 1:numel(worker_errors), L{end+1} = a('- Error de proceso: %s', strtok(worker_errors{i}, newline)); end %#ok<AGROW>
    if ~isempty(extra.update)   % the earlier updates' sections of the README as they are, then this one's
        if isfile(file)
            t = cellstr(splitlines(string(fileread(file))))';
            h = find(startsWith(t, '## '));
            for k = find(startsWith(t(h), '## Actualización'))
                e = numel(t);
                if k < numel(h), e = h(k + 1) - 1; end
                sec = t(h(k):e);
                if strcmp(sec{1}, extra.update_head), continue; end   % a repeated 'write' of this update
                L = [L, {''}, sec(1:find(~cellfun(@isempty, strtrim(sec)), 1, 'last'))]; %#ok<AGROW>
            end
        end
        L{end+1} = '';
        L{end+1} = extra.update_head;
        L = [L, extra.update];
    end
    L{end+1} = '';
    L{end+1} = '## Diferencias con run_campaign_stats del 01/11 (`_stats_revision_20261101_final`)';
    L{end+1} = '| ratón | diferencia |';
    L{end+1} = '|---|---|';
    for i = find(~strcmp(CMP.diferencia, ''))'
        L{end+1} = a('| %s | %s |', CMP.mouse{i}, CMP.diferencia{i}); %#ok<AGROW>
    end
    fid = fopen(file, 'w', 'n', 'UTF-8');
    fprintf(fid, '%s\n', L{:});
    fclose(fid);
end
