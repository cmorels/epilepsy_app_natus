function result = run_pipeline_edf(input_folder, cfg)
% RUN_PIPELINE_EDF  Batch orchestrator: EDF -> txt -> quality -> case ->
% precondition -> clean -> seizures -> IID -> consolidated summaries, for
% every .edf in input_folder.
%
%   result = run_pipeline_edf(input_folder, cfg)
%   result = run_pipeline_edf({'path/a.EDF', 'path/b.EDF'}, cfg)
%
% input_folder is a folder (every .edf in it) or a cell array of .edf
% paths (only those; a listed file that does not exist gets a qc_report
% row with the error, the run continues).
%
% Case of each channel: cfg.cases.from_excel = true -> from the recording
% Excel's Attenuation column (utils/resolve_case_excel.m, see
% campaign_config.m); otherwise cases CSV / force / default (resolve_case.m).
%
% Output layout under cfg.paths.output_root:
%   01_txt/        txt per channel + gaps CSV (edf_import.m)
%   02_clean/      *_clean.txt (clean_lfp.m, after precondition_lfp.m's gain + notch)
%   03_seizures/   mat/fig/png per channel (detect_seizures.m)
%   04_iid/        mat/fig/png per channel (detect_iid.m)
%   05_summaries/  the 9 consolidated CSVs + pipeline_summary.xlsx
%   logs/          timestamped run log
%   config_used.mat / config_used.json  (the exact cfg this run used)
%
% Per-channel order (order matters -- see README.md "orden de
% operaciones"): load_lfp_txt -> signal_quality (on the RAW signal,
% immune to 50 Hz) -> resolve_case (CSV / force / default) ->
% precondition_lfp (gain, with dead band) -> clean_lfp (outliers in
% calibrated uV, then notch, writes *_clean.txt) -> detect_seizures ->
% detect_iid.
%
% DEFAULT BEHAVIOR IS UNCHANGED: with pipeline_config() untouched, every
% channel resolves to case 'normal' (gain off, notch off), and
% precondition_lfp.m returns the signal completely untouched (not even a
% x*1 multiplication) -- see KNOWN_ISSUES.md / README.md "sistema de
% casos" for why this is a hard requirement, not an implementation detail.
%
% Robustness: a failing file never aborts the batch. Every stage call is
% individually try/caught; failures are logged (file, region, stage,
% message, line) into qc_report.csv and the run log, and the loop moves on
% with whatever partial results are available for that channel. Each
% file's result rows are also saved to logs/checkpoints/file_NNN.mat as
% soon as the file is done; the final consolidation stamps one TimeZone on
% every datetime column and falls back to text concatenation if a table
% still cannot be stacked (utils/robust_vertcat.m), so a type mismatch can
% no longer throw away the whole run. The checkpoints are deleted once
% 05_summaries/ is written, and kept (with the error in the log) if not.
%
% Resumability (scoped, see README.md):
%   - 02_clean is genuinely skipped (not just not-overwritten) when its
%     output txt already exists, cfg.general.overwrite is false, AND the
%     case resolved for this run matches what's recorded in that file's
%     own header (case_applied, gain_applied, gain_source, notch_applied)
%     -- comparing the header, not re-deriving gain/notch from scratch,
%     is what lets this check stay cheap even for a very long recording.
%     If the case changed, clean_lfp (and only clean_lfp) reruns for that
%     channel; other channels/files are untouched.
%   - 01_txt/03_seizures/04_iid still run every time; they internally
%     refuse to overwrite an existing file, but do not skip the
%     underlying computation.
%
% To combine several separate runs (e.g. one per subject) into one set of
% summaries, see merge_pipeline_runs.m -- this function only consolidates
% files processed within THIS call; a second call with the same
% cfg.paths.output_root overwrites 05_summaries/, it does not append.

    if nargin < 2 || isempty(cfg)
        cfg = pipeline_config();
    end
    [files, missing_inputs, input_desc] = resolve_edf_inputs(input_folder);
    if isempty(files) && isempty(missing_inputs)
        error('run_pipeline_edf:NoFilesFound', 'No .edf files found in: %s', input_desc);
    end

    dirs = make_output_dirs(cfg.paths.output_root);
    log_file = start_log(dirs.logs);
    save_config_used(cfg, cfg.paths.output_root);
    checkpoint_dir = fullfile(dirs.logs, 'checkpoints');
    if isfolder(checkpoint_dir)
        rmdir(checkpoint_dir, 's');  % stale checkpoints of an earlier run in this output_root
    end

    log_line(log_file, sprintf('run_pipeline_edf: %d EDF file(s) from %s', numel(files), input_desc));
    tz = cfg.general.timezone;

    cases_table = table();
    if ~isempty(cfg.cases.file)
        cases_table = load_cases(cfg.cases.file, cfg);
        log_line(log_file, sprintf('Loaded %d case row(s) from %s', height(cases_table), cfg.cases.file));
    end
    row_matched = false(height(cases_table), 1);

    seizure_event_parts = {}; seizure_summary_parts = {};
    iid_event_parts = {}; iid_summary_parts = {}; iid_burst_parts = {};
    gap_parts = {}; qc_parts = {};
    file_subject_map = containers.Map('KeyType', 'char', 'ValueType', 'char');

    for i = 1:numel(missing_inputs)
        [~, mname, mext] = fileparts(missing_inputs{i});
        log_line(log_file, sprintf('  ERROR [%s | | input] listed EDF not found: %s', [mname mext], missing_inputs{i}));
        qc_parts{end+1} = build_qc_row(qc_info_minimal(cfg.edf.subject_id, '', [mname mext], zoned_nat(tz), ...
            {}, {sprintf('input: listed EDF not found (%s)', missing_inputs{i})}, {})); %#ok<AGROW>
    end

    case_tally = containers.Map('KeyType', 'char', 'ValueType', 'double');
    clean_status_tally = containers.Map('KeyType', 'char', 'ValueType', 'double');
    n_unusable = 0;
    n_suggested_mismatch = 0;
    bilateral_on = isfield(cfg, 'bilateral') && ~strcmp(cfg.bilateral.rescue_mode, 'off');
    if bilateral_on
        log_line(log_file, sprintf('bilateral_reconcile: rescue_mode=%s, match_tol_s=%g, exclude_rescued_from_iid=%d', ...
            cfg.bilateral.rescue_mode, cfg.bilateral.match_tol_s, cfg.bilateral.exclude_rescued_from_iid));
    end
    [ll_accept, ll_reject, band_on] = resolve_ll_band(cfg);
    if band_on
        log_line(log_file, sprintf('review band: ll_reject=%.2f, ll_accept=%.2f (robust branch); event categories + figures under 03_seizures/<category>/', ...
            ll_reject, ll_accept));
    end
    if band_on && bilateral_on
        error('run_pipeline_edf:BandAndRescue', ...
            ['cfg.bilateral.rescue_mode=''%s'' together with the review band (ll_accept/ll_reject) is not supported yet. ' ...
             'Set cfg.bilateral.rescue_mode=''off'', or cfg.seizure_robust.ll_accept = ll_reject = [] to switch the band off.'], ...
            cfg.bilateral.rescue_mode);
    end
    % Both need every channel of a recording before any seizure output is written.
    deferred = bilateral_on || band_on;
    any_robust = false;

    for f = 1:numel(files)
        file_path = fullfile(files(f).folder, files(f).name);
        log_line(log_file, sprintf('FILE %d/%d: %s', f, numel(files), files(f).name));

        file_cfg = cfg;
        file_cfg.edf.output_dir = dirs.txt;
        n_before = [numel(seizure_event_parts), numel(seizure_summary_parts), numel(iid_event_parts), ...
            numel(iid_summary_parts), numel(iid_burst_parts), numel(gap_parts), numel(qc_parts)];

        try
            [manifest, gaps_table, ~] = edf_import(file_path, file_cfg);
        catch ME
            log_error(log_file, files(f).name, '', 'edf_import', ME);
            % zoned NaT: an unzoned one cannot be stacked with the other rows
            qc_parts{end+1} = build_qc_row(qc_info_minimal(cfg.edf.subject_id, '', files(f).name, zoned_nat(tz), ...
                {}, {sprintf('edf_import: %s', ME.message)}, {})); %#ok<AGROW>
            write_checkpoint(checkpoint_dir, f, files(f).name, n_before, seizure_event_parts, seizure_summary_parts, ...
                iid_event_parts, iid_summary_parts, iid_burst_parts, gap_parts, qc_parts);
            continue;
        end
        gap_parts{end+1} = gaps_table; %#ok<AGROW>
        if height(manifest) > 0
            file_subject_map(files(f).name) = manifest.subject_id{1};
        end

        channel_ctx = {};
        for c = 1:height(manifest)
            row = manifest(c, :);
            region = row.region{1};
            stages = {}; errors = {}; warnings_list = {};
            raw_data = []; q = []; case_spec = []; clean_data = []; clean_stats = []; notch_info = [];
            seizure_results = []; seizure_mode = '';

            try
                raw_data = load_lfp_txt(row.txt_file{1});
                stages{end+1} = 'import'; %#ok<AGROW>
            catch ME
                log_error(log_file, files(f).name, region, 'load_raw', ME);
                errors{end+1} = sprintf('load_raw: %s', ME.message); %#ok<AGROW>
            end

            if ~isempty(raw_data)
                try
                    q = signal_quality(raw_data, file_cfg);
                    stages{end+1} = 'quality'; %#ok<AGROW>
                    if strcmp(q.quality_class, 'unusable')
                        n_unusable = n_unusable + 1;
                        warnings_list{end+1} = 'quality_class=unusable'; %#ok<AGROW>
                    end
                catch ME
                    log_error(log_file, files(f).name, region, 'signal_quality', ME);
                    errors{end+1} = sprintf('signal_quality: %s', ME.message); %#ok<AGROW>
                end
            end

            if ~isempty(q)
                try
                    if isfield(cfg.cases, 'from_excel') && cfg.cases.from_excel
                        [case_spec, case_note] = resolve_case_excel(manifest_attenuation(row), cfg);
                        if ~isempty(case_note)
                            warnings_list{end+1} = case_note; %#ok<AGROW>
                            log_line(log_file, sprintf('  [%s] WARNING %s', region, case_note));
                        end
                        matched_idx = NaN;
                    else
                        [~, txt_name, txt_ext] = fileparts(row.txt_file{1});
                        candidate_names = {files(f).name, [txt_name txt_ext]};
                        [case_spec, matched_idx] = resolve_case(candidate_names, region, cases_table, cfg);
                    end
                    if ~isnan(matched_idx)
                        row_matched(matched_idx) = true;
                    end
                    case_tally(case_spec.case_applied) = get_or_zero(case_tally, case_spec.case_applied) + 1;
                    if ~isempty(case_spec.suggested_case) && ~strcmp(case_spec.suggested_case, case_spec.case_applied)
                        n_suggested_mismatch = n_suggested_mismatch + 1;
                    end
                    stages{end+1} = 'case'; %#ok<AGROW>
                    fprintf('  [%s] case=%s (source=%s)%s\n', region, case_spec.case_applied, case_spec.case_source, ...
                        suggestion_note(case_spec));
                catch ME
                    log_error(log_file, files(f).name, region, 'resolve_case', ME);
                    errors{end+1} = sprintf('resolve_case: %s', ME.message); %#ok<AGROW>
                end
            end

            if ~isempty(case_spec)
                seizure_mode = case_spec.seizure_mode;  % known even if clean fails (bilateral_reconcile reports such a channel)
                try
                    [clean_data, clean_stats, clean_status, notch_info] = run_clean_stage( ...
                        raw_data, q, case_spec, file_cfg, dirs.clean);
                    stages{end+1} = 'clean'; %#ok<AGROW>
                    clean_status_tally(clean_status) = get_or_zero(clean_status_tally, clean_status) + 1;
                    log_line(log_file, sprintf('  [%s] clean_lfp: %s (case=%s)', region, clean_status, case_spec.case_applied));
                    if clean_stats.auto_switch
                        warnings_list{end+1} = clean_stats.switched_reason; %#ok<AGROW>
                    end
                    gain_note = gain_warning(case_spec, clean_data);
                    if ~isempty(gain_note)
                        warnings_list{end+1} = gain_note; %#ok<AGROW>
                        log_line(log_file, sprintf('  [%s] WARNING %s', region, gain_note));
                    end
                catch ME
                    log_error(log_file, files(f).name, region, 'clean_lfp', ME);
                    errors{end+1} = sprintf('clean_lfp: %s', ME.message); %#ok<AGROW>
                end
            end

            if ~isempty(clean_data)
                seizure_mode = case_spec.seizure_mode;
                reconciled_channel = bilateral_on && strcmp(seizure_mode, 'robust');
                if reconciled_channel
                    % The robust detector's own figures are numbered per channel;
                    % with reconciliation on they go to a throwaway folder and
                    % save_bilateral_seizure_figures.m writes the real ones.
                    % Legacy channels are never reconciled: their figures stay.
                    detector_dir = tempname();
                    file_cfg.seizure.output_dir = detector_dir;
                else
                    file_cfg.seizure.output_dir = dirs.seizures;
                end
                try
                    if strcmp(seizure_mode, 'robust')
                        seizure_results = detect_seizures_robust(clean_data, file_cfg);
                    else
                        seizure_results = detect_seizures(clean_data, file_cfg);
                    end
                    stages{end+1} = 'seizures'; %#ok<AGROW>
                    if seizure_results.qc.n_blocks_rejected_short > 0
                        warnings_list{end+1} = sprintf('%d block(s) rejected as too short for seizure detection', ...
                            seizure_results.qc.n_blocks_rejected_short); %#ok<AGROW>
                    end
                catch ME
                    log_error(log_file, files(f).name, region, 'detect_seizures', ME);
                    errors{end+1} = sprintf('detect_seizures(%s): %s', seizure_mode, ME.message); %#ok<AGROW>
                end
                if reconciled_channel && isfolder(detector_dir)
                    rmdir(detector_dir, 's');
                end
            end

            ctx = struct('row', row, 'region', region, 'source_name', files(f).name, ...
                'stages', {stages}, 'errors', {errors}, 'warnings_list', {warnings_list}, ...
                'q', q, 'case_spec', case_spec, 'clean_data', clean_data, 'clean_stats', clean_stats, ...
                'notch_info', notch_info, 'seizure_results', seizure_results, 'seizure_mode', seizure_mode, ...
                'session_start', pick_session_start(raw_data, clean_data, cfg.general.timezone));
            raw_data = []; %#ok<NASGU>
            any_robust = any_robust || strcmp(seizure_mode, 'robust');

            if deferred
                % IID and every seizure output wait until all channels of this
                % recording are detected (bilateral_reconcile.m / bilateral_events.m).
                channel_ctx{end+1} = ctx; %#ok<AGROW>
            else
                out = finalize_channel(ctx, file_cfg, dirs, log_file, false, false);
                [seizure_event_parts, seizure_summary_parts, iid_event_parts, iid_summary_parts, iid_burst_parts, qc_parts] = ...
                    append_parts(out, seizure_event_parts, seizure_summary_parts, iid_event_parts, iid_summary_parts, iid_burst_parts, qc_parts);
            end
        end

        if deferred && ~isempty(channel_ctx)
            if bilateral_on
                channel_ctx = reconcile_file_channels(channel_ctx, file_cfg, dirs, log_file);
            end
            outs = cell(1, numel(channel_ctx));
            for c = 1:numel(channel_ctx)
                outs{c} = finalize_channel(channel_ctx{c}, file_cfg, dirs, log_file, bilateral_on, band_on);
            end
            if band_on
                outs = finalize_file_events(outs, channel_ctx, file_cfg, dirs, log_file);
            end
            for c = 1:numel(outs)
                [seizure_event_parts, seizure_summary_parts, iid_event_parts, iid_summary_parts, iid_burst_parts, qc_parts] = ...
                    append_parts(outs{c}, seizure_event_parts, seizure_summary_parts, iid_event_parts, iid_summary_parts, iid_burst_parts, qc_parts);
            end
        end
        clear channel_ctx outs
        write_checkpoint(checkpoint_dir, f, files(f).name, n_before, seizure_event_parts, seizure_summary_parts, ...
            iid_event_parts, iid_summary_parts, iid_burst_parts, gap_parts, qc_parts);
    end

    for i = 1:numel(row_matched)
        if ~row_matched(i)
            warning('run_pipeline_edf:UnmatchedCaseRow', ...
                'cases.file row %d (source_file="%s") never matched a processed file/channel.', ...
                i, cases_table.source_file{i});
        end
    end

    % Consolidation: robust_vertcat never throws (see header); the review
    % sheet and the writer are guarded so whatever can be written is.
    events_on = band_on && any_robust;  % event columns only if the run had robust channels (legacy-only: just ll_status)
    if bilateral_on
        seizures_events = robust_vertcat(seizure_event_parts, @() empty_seizure_events_table_bilateral(tz), tz, 'seizures_events');
        seizures_summary = robust_vertcat(seizure_summary_parts, @() empty_seizure_summary_table_bilateral(tz), tz, 'seizures_summary');
    elseif band_on
        ev_kinds = {'ll_status'};
        if events_on
            ev_kinds{end+1} = 'seizures_events';
            seizure_summary_parts = cellfun(@(t) add_event_columns(t, 'seizures_summary'), seizure_summary_parts, 'UniformOutput', false);
        end
        for k = ev_kinds
            seizure_event_parts = cellfun(@(t) add_event_columns(t, k{1}), seizure_event_parts, 'UniformOutput', false);
        end
        seizures_events = robust_vertcat(seizure_event_parts, @() empty_seizure_events_table_events(tz, events_on), tz, 'seizures_events');
        seizures_summary = robust_vertcat(seizure_summary_parts, @() empty_seizure_summary_table_events(tz, events_on), tz, 'seizures_summary');
    else
        seizures_events = robust_vertcat(seizure_event_parts, @() empty_seizure_events_table(tz), tz, 'seizures_events');
        seizures_summary = robust_vertcat(seizure_summary_parts, @() empty_seizure_summary_table(tz), tz, 'seizures_summary');
    end
    iid_events = robust_vertcat(iid_event_parts, @() empty_iid_events_table(tz), tz, 'iid_events');
    iid_summary = robust_vertcat(iid_summary_parts, @() empty_iid_summary_table(tz), tz, 'iid_summary');
    iid_bursts = robust_vertcat(iid_burst_parts, @() empty_iid_bursts_table(tz), tz, 'iid_bursts');
    gaps_summary = robust_vertcat(gap_parts, @() empty_gaps_table(tz), tz, 'gaps_summary');
    qc_report = robust_vertcat(qc_parts, @() empty_qc_table(tz), tz, 'qc_report');

    try
        if bilateral_on
            natus_review_sheet = build_natus_review_sheet_bilateral(seizures_events, iid_bursts, gaps_summary, file_subject_map, tz);
        elseif events_on
            natus_review_sheet = build_natus_review_sheet_events(seizures_events, iid_bursts, gaps_summary, file_subject_map, tz);
        else
            natus_review_sheet = build_natus_review_sheet(seizures_events, iid_bursts, gaps_summary, file_subject_map, tz);
        end
    catch ME
        log_line(log_file, sprintf('ERROR building natus_review_sheet (written empty, every other summary unaffected): %s', ME.message));
        natus_review_sheet = table(cell(0, 1), 'VariableNames', {'natus_review_sheet_failed'});
    end

    try
        paths = write_all_summaries(dirs.summaries, seizures_events, seizures_summary, ...
            iid_events, iid_summary, iid_bursts, gaps_summary, qc_report, natus_review_sheet);
        if isfolder(checkpoint_dir)
            rmdir(checkpoint_dir, 's');
        end
    catch ME
        log_line(log_file, sprintf('ERROR writing 05_summaries: %s -- per-file results kept in %s', ME.message, checkpoint_dir));
        paths = struct();
    end

    log_line(log_file, sprintf('DONE: %d seizure(s), %d IID complex(es), %d burst(s), %d error(s) across %d file(s)', ...
        height(seizures_events), height(iid_events), height(iid_bursts), count_errors(qc_report), numel(files)));
    print_case_summary(case_tally, clean_status_tally, n_unusable, n_suggested_mismatch);
    fclose(log_file);

    result = struct();
    result.dirs = dirs;
    result.seizures_events = seizures_events;
    result.seizures_summary = seizures_summary;
    result.iid_events = iid_events;
    result.iid_summary = iid_summary;
    result.iid_bursts = iid_bursts;
    result.gaps_summary = gaps_summary;
    result.qc_report = qc_report;
    result.natus_review_sheet = natus_review_sheet;
    result.paths = paths;
    result.bilateral_on = bilateral_on;
end

%% ======================================================================
% Per-channel tail of the pipeline (IID + every CSV row of that channel),
% shared by the unreconciled path (called right after each channel's
% detection, exactly the order the pipeline always had) and the
% reconciled one (called after bilateral_reconcile.m, once per channel).
function out = finalize_channel(ctx, file_cfg, dirs, log_file, bilateral_on, band_on)
    row = ctx.row;
    region = ctx.region;
    stages = ctx.stages; errors = ctx.errors; warnings_list = ctx.warnings_list;
    clean_data = ctx.clean_data;
    seizure_results = ctx.seizure_results;
    iid_results = [];

    if ~isempty(clean_data)
        file_cfg.iid.output_dir = dirs.iid;
        try
            seizures_for_iid = seizures_for_iid_of(seizure_results, file_cfg, bilateral_on, band_on);
            iid_results = detect_iid(clean_data, seizures_for_iid, file_cfg);
            stages{end+1} = 'iid';
        catch ME
            log_error(log_file, ctx.source_name, region, 'detect_iid', ME);
            errors{end+1} = sprintf('detect_iid: %s', ME.message);
        end
    end

    session_start = ctx.session_start;
    out = struct('seizure_events', {{}}, 'seizure_summary', {{}}, 'iid_events', {{}}, ...
        'iid_summary', {{}}, 'iid_bursts', {{}}, 'qc', {{}});

    if ~isempty(seizure_results)
        if bilateral_on
            if isempty(clean_data)  % failed robust channel reported by bilateral_reconcile (no signal to check gaps on)
                valid_mask = []; t_rel = [];
            else
                valid_mask = clean_data.valid_mask; t_rel = clean_data.t_rel;
            end
            out.seizure_events = {build_seizure_event_rows_bilateral( ...
                seizure_results.seizures, region, row.subject_id{1}, session_start, row.source_file{1}, ...
                valid_mask, t_rel, file_cfg.seizure.edge_trim_s, ctx.seizure_mode)};
            out.seizure_summary = {build_seizure_summary_row_bilateral(row, session_start, seizure_results, file_cfg)};
        else
            T = build_seizure_event_rows( ...
                seizure_results.seizures, region, row.subject_id{1}, session_start, row.source_file{1}, ...
                clean_data.valid_mask, clean_data.t_rel, file_cfg.seizure.edge_trim_s, ctx.seizure_mode);
            summary_sr = seizure_results;
            if band_on
                % ll_status per row (legacy rows: always 'accepted'); every
                % PRE-EXISTING summary column keeps counting accepted rows only,
                % in_band rows are counted apart (finalize_file_events).
                S = seizure_results.seizures;
                if ismember('ll_status', S.Properties.VariableNames)
                    T.ll_status = S.ll_status;
                    summary_sr.seizures = S(strcmp(S.ll_status, 'accepted'), :);
                else
                    T.ll_status = repmat({'accepted'}, height(T), 1);
                end
            end
            out.seizure_events = {T};
            out.seizure_summary = {build_seizure_summary_row(row, session_start, summary_sr, file_cfg)};
        end
    end
    if ~isempty(iid_results)
        out.iid_events = {build_iid_event_rows( ...
            iid_results.spike_complex_table, iid_results.burst_table, region, row.subject_id{1}, session_start, row.source_file{1})};
        out.iid_summary = {build_iid_summary_row(row, session_start, iid_results)};
        out.iid_bursts = {build_iid_burst_rows(iid_results.burst_table, region, row.subject_id{1}, row.source_file{1}, session_start.TimeZone)};
    end

    info = qc_info_full(row, session_start, stages, errors, warnings_list, ...
        clean_data, ctx.clean_stats, seizure_results, ctx.q, ctx.case_spec, ctx.notch_info, file_cfg);
    out.qc = {build_qc_row(info)};
end

function outs = finalize_file_events(outs, ctxs, cfg, dirs, log_file)
% Review band on: group this recording's robust detections into events
% (bilateral_events.m), draw the per-event individual + joint figures into
% 03_seizures/<category>/ (save_event_figures.m), and append the per-channel
% event counts / rates to each channel's summary row. A recording with no
% robust channel is left with ll_status only (legacy rows are unchanged).
    robust = cellfun(@(c) strcmp(c.seizure_mode, 'robust'), ctxs);
    if ~any(robust)
        return;
    end
    has_rows = cellfun(@(o) ~isempty(o.seizure_events), outs);
    parts = cellfun(@(o) o.seizure_events{1}, outs(has_rows), 'UniformOutput', false);
    heights = cellfun(@height, parts);
    if isempty(parts)
        E = add_event_columns(empty_seizure_events_table(cfg.general.timezone), 'll_status');
    else
        E = vertcat(parts{:});
    end
    [E, notes] = bilateral_events(E, cfg);

    % one column per robust channel with a usable clean signal, detected or not
    channels = struct('region', {}, 'subject_id', {}, 'session_start', {}, 'label', {}, 'fs', {}, ...
        't_rel', {}, 'signal', {}, 'trace', {});
    for i = find(robust)
        c = ctxs{i};
        if isempty(c.clean_data)
            continue;
        end
        try
            tr = seizure_energy_trace(c.clean_data, cfg);
            if isempty(tr.blocks)
                trimmed = zeros(0, 2);
            else
                trimmed = [vertcat(tr.blocks.trimmed_start), vertcat(tr.blocks.trimmed_end)];
            end
            [ll_full, ll_med] = line_length_trace(tr.bp_full, trimmed, c.clean_data.fs, cfg);
        catch ME
            notes{end+1} = sprintf('[%s] no traces for the event figures (%s): column omitted', c.region, ME.message); %#ok<AGROW>
            continue;
        end
        label = '';
        if isfield(c.clean_data.meta, 'channel_label')
            label = strtrim(c.clean_data.meta.channel_label);
        end
        channels(end+1) = struct('region', c.region, 'subject_id', c.row.subject_id{1}, ...
            'session_start', c.session_start, 'label', label, 'fs', c.clean_data.fs, ...
            't_rel', c.clean_data.t_rel, 'signal', c.clean_data.signal, ...
            'trace', struct('energy_full', tr.energy_full, 'bp_full', tr.bp_full, 'll_full', ll_full, ...
                'll_median_global', ll_med, 'threshold', tr.threshold)); %#ok<AGROW>
    end
    try
        [E, fnotes] = save_event_figures(E, channels, dirs.seizures, cfg);
        notes = [notes, fnotes];
    catch ME
        log_error(log_file, ctxs{1}.source_name, '', 'save_event_figures', ME);
    end
    for i = 1:numel(notes)
        log_line(log_file, ['  [events] ' notes{i}]);
    end
    ev = E(strcmp(E.seizure_mode, 'robust') & ~isnan(E.event_id), :);
    n_cat = @(cat) numel(unique(ev.event_id(strcmp(ev.category, cat))));
    log_line(log_file, sprintf('  [events] %d event(s): %d Crisis, %d Candidates, %d Candidates_in_band', ...
        numel(unique(ev.event_id)), n_cat('Crisis'), n_cat('Candidates'), n_cat('Candidates_in_band')));

    % split E back into the per-channel parts, in their original order
    k = 0;
    idx_rows = find(has_rows);
    for j = 1:numel(idx_rows)
        outs{idx_rows(j)}.seizure_events = {E(k + (1:heights(j)), :)};
        k = k + heights(j);
    end

    % per-channel summary columns (event_columns('seizures_summary'))
    for i = 1:numel(outs)
        if isempty(outs{i}.seizure_summary)
            continue;
        end
        S = outs{i}.seizure_summary{1};
        mine = E(strcmp(E.region, ctxs{i}.region) & strcmp(E.seizure_mode, 'robust'), :);
        total_s = 60 * S.total_duration_min;
        t_crisis = sum(mine.duration_s(strcmp(mine.category, 'Crisis')));
        t_cc = sum(mine.duration_s(ismember(mine.category, {'Crisis', 'Candidates'})));
        S.n_events_crisis = n_cat('Crisis');
        S.n_events_candidates = n_cat('Candidates');
        S.n_events_candidates_in_band = n_cat('Candidates_in_band');
        S.n_accepted_this_channel = nnz(strcmp(mine.ll_status, 'accepted'));
        S.n_in_band_this_channel = nnz(strcmp(mine.ll_status, 'in_band'));
        S.total_seizure_time_s_crisis = t_crisis;
        S.pct_time_in_seizure_crisis = 100 * t_crisis / total_s;
        S.total_seizure_time_s_crisis_candidates = t_cc;
        S.pct_time_in_seizure_crisis_candidates = 100 * t_cc / total_s;
        outs{i}.seizure_summary = {S};
    end
end

function [se, ss, ie, is, ib, qc] = append_parts(out, se, ss, ie, is, ib, qc)
    se = [se, out.seizure_events];
    ss = [ss, out.seizure_summary];
    ie = [ie, out.iid_events];
    is = [is, out.iid_summary];
    ib = [ib, out.iid_bursts];
    qc = [qc, out.qc];
end

function T = seizures_for_iid_of(seizure_results, cfg, bilateral_on, band_on)
% Unreconciled: the detector's table, as always -- except that with the
% review band on, 'in_band' rows are dropped unless
% cfg.bilateral.exclude_in_band_from_iid (default false: unconfirmed events
% must not remove time from the spike analysis). Reconciled: every
% reported row (accepted + rescued + imputed) when
% cfg.bilateral.exclude_rescued_from_iid, else only the accepted ones;
% the exclusion zone's source_id is then the SHARED seizure_id.
    if isempty(seizure_results)
        T = [];
        return;
    end
    if ~bilateral_on
        T = seizure_results.seizures;
        if band_on && ismember('ll_status', T.Properties.VariableNames) && ~cfg.bilateral.exclude_in_band_from_iid
            T = T(strcmp(T.ll_status, 'accepted'), :);
        end
        return;
    end
    S = seizure_results.seizures;
    if ~cfg.bilateral.exclude_rescued_from_iid
        S = S(ismember(S.detection_status, {'accepted', 'not_reconciled'}), :);
    end
    T = table(S.seizure_id, S.start_s, S.end_s, S.duration_s, S.start_abs, S.end_abs, S.block_id, ...
        'VariableNames', {'id', 'start_s', 'end_s', 'duration_s', 'start_abs', 'end_abs', 'block_id'});
end

%% ======================================================================
function ctxs = reconcile_file_channels(ctxs, file_cfg, dirs, log_file)
% Runs bilateral_reconcile.m over the robust-branch channels of ONE
% recording, then writes each reconciled channel's figures + .mat into
% 03_seizures/. Legacy channels are left exactly as detect_seizures.m
% produced them (own figures, own per-channel ids) and only tagged
% 'not_reconciled' so the CSV schema stays uniform.
    n = numel(ctxs);
    pc = struct('region', cell(1, n), 'subject_id', [], 'source_file', [], 'session_start', [], ...
        'seizure_mode', [], 'seizure_results', [], 'data', []);
    for i = 1:n
        c = ctxs{i};
        pc(i).region = c.region;
        pc(i).subject_id = c.row.subject_id{1};
        pc(i).source_file = c.row.source_file{1};
        pc(i).session_start = c.session_start;
        pc(i).seizure_mode = c.seizure_mode;
        pc(i).seizure_results = c.seizure_results;
        pc(i).data = c.clean_data;
    end

    try
        rec = bilateral_reconcile(pc, file_cfg);
    catch ME
        for i = 1:n
            log_error(log_file, ctxs{i}.source_name, ctxs{i}.region, 'bilateral_reconcile', ME);
            ctxs{i}.errors{end+1} = sprintf('bilateral_reconcile: %s', ME.message);
            if ~isempty(ctxs{i}.seizure_results)
                ctxs{i}.seizure_results = mark_unreconciled(ctxs{i}.seizure_results, 'failed');
            end
        end
        return;
    end

    for i = 1:n
        if isempty(rec(i).trace)
            sr = rec(i).seizure_results;
            if ~isempty(sr) && isfield(sr, 'detection_failed')
                % failed robust channel: one imputed row per event, no metrics, no figures
                ctxs{i}.seizure_results = sr;
                ctxs{i}.stages{end+1} = 'bilateral';
                ctxs{i}.warnings_list{end+1} = sprintf('seizure detection failed: %d event(s) reported as imputed without metrics', ...
                    sr.bilateral.n_imputed);
                log_line(log_file, sprintf('  [%s] bilateral_reconcile: detection FAILED on this channel -> %d imputed row(s) without metrics; %d event(s) across %d channel(s)', ...
                    ctxs{i}.region, sr.bilateral.n_imputed, sr.bilateral.n_events_session, sr.bilateral.n_channels_total));
            elseif ~isempty(ctxs{i}.seizure_results)
                % legacy channel: left exactly as detect_seizures.m produced it
                ctxs{i}.seizure_results = mark_unreconciled(ctxs{i}.seizure_results, ctxs{i}.seizure_mode);
            end
            continue;
        end
        sr = rec(i).seizure_results;
        [~, fname] = fileparts(ctxs{i}.clean_data.file);
        base = strrep(fname, '_clean', '');
        sr.mat_file = ''; sr.fig_file = ''; sr.png_file = '';
        try
            fl = save_bilateral_seizure_figures(dirs.seizures, base, file_cfg, ctxs{i}.clean_data, rec(i).trace, ...
                sr.seizures, ctxs{i}.seizure_mode, ctxs{i}.region);
            sr.fig_file = fl.panorama_fig;
            sr.png_file = fl.panorama_png;
        catch ME
            warning('run_pipeline_edf:BilateralFigures', 'Could not save reconciled seizure figures for %s: %s', base, ME.message);
            ctxs{i}.warnings_list{end+1} = sprintf('bilateral figures failed: %s', ME.message);
        end
        if height(sr.seizures) > 0
            sr.mat_file = fullfile(dirs.seizures, [base '_seizures.mat']);
            seizure_results = sr;
            save(sr.mat_file, 'seizure_results');
        end
        b = sr.bilateral;
        log_line(log_file, sprintf('  [%s] bilateral_reconcile: %d accepted, %d rescued, %d imputed -> %d row(s); %d event(s) across %d channel(s)', ...
            ctxs{i}.region, b.n_accepted, b.n_rescued, b.n_imputed, b.n_reported, b.n_events_session, b.n_channels_total));
        ctxs{i}.stages{end+1} = 'bilateral';
        ctxs{i}.seizure_results = sr;
    end
end

function sr = mark_unreconciled(sr, why)
% Channel not reconciled (legacy branch, or reconciliation failed): keep
% the detector's own rows, under their per-channel id, tagged
% 'not_reconciled' (never 'accepted'), so the output schema stays uniform
% and nothing claims a comparison across channels that did not happen.
% Idempotent: a table already carrying detection_status is left alone.
    if isfield(sr, 'bilateral')
        return;
    end
    S = sr.seizures;
    n = height(S);
    for col = {'over_max_duration', 'll_ratio', 'peak_energy_ratio', 'hf_ratio_db', 'envelope_cv'}
        if ~ismember(col{1}, S.Properties.VariableNames)
            if strcmp(col{1}, 'over_max_duration')
                S.(col{1}) = false(n, 1);
            else
                S.(col{1}) = nan(n, 1);
            end
        end
    end
    S.seizure_id = S.id;
    S.detection_status = repmat({'not_reconciled'}, n, 1);
    S.accepted_in_n_channels = nan(n, 1);
    S.accepted_in_regions = repmat({''}, n, 1);
    S.n_channels_total = nan(n, 1);
    S.is_bilateral_accepted = false(n, 1);
    S.fragmented = false(n, 1);
    S.ref_start_s = S.start_s;
    S.ref_end_s = S.end_s;
    S.rejected_by = repmat({''}, n, 1);
    sr.seizures_accepted = sr.seizures;
    sr.seizures = S;
    sr.bilateral = struct('rescue_mode', ['not_reconciled:' why], 'n_accepted', n, 'n_reported', n, 'n_events', n, ...
        'n_rescued', 0, 'n_imputed', 0, 'n_channels_total', NaN, 'n_events_session', NaN);
end

function T = build_seizure_event_rows_bilateral(seizures, region, subject_id, session_start, source_file, valid_mask, t_rel, edge_trim_s, seizure_mode)
% Reconciled counterpart of build_seizure_event_rows: same 18 columns in the
% same order, where seizure_id is now the SHARED event id, plus
% bilateral_columns('seizures_events') appended (channel_seizure_id holds
% the detector's own per-channel index, NaN on rescued/imputed rows).
    n = height(seizures);
    if n == 0
        T = empty_seizure_events_table_bilateral(session_start.TimeZone);
        return;
    end
    T = seizures;
    T.channel_seizure_id = T.id;
    T.id = [];
    T.subject_id = repmat({subject_id}, n, 1);
    T.region = repmat({region}, n, 1);
    T.session_start = repmat(session_start, n, 1);
    T.source_file = repmat({source_file}, n, 1);
    T.adjacent_to_gap = compute_adjacent_to_gap(seizures.start_s, seizures.end_s, valid_mask, t_rel, edge_trim_s);
    T.seizure_mode = repmat({seizure_mode}, n, 1);
    base_names = empty_seizure_events_table(session_start.TimeZone).Properties.VariableNames;
    T = T(:, [base_names, bilateral_columns('seizures_events')]);
end

function T = build_seizure_summary_row_bilateral(row, session_start, seizure_results, cfg)
% Every pre-existing column keeps its meaning: it describes what THIS
% channel's detector accepted (seizures_accepted). The appended columns
% give the same counts/time/percentage twice, *_accepted and *_reported,
% so an arrastrada crisis never silently inflates a per-channel rate.
    base = seizure_results;
    base.seizures = seizure_results.seizures_accepted;
    T = build_seizure_summary_row(row, session_start, base, cfg);

    S = seizure_results.seizures;
    acc = ismember(S.detection_status, {'accepted', 'not_reconciled'});
    total_s = row.total_duration_s;
    t_acc = sum(S.duration_s(acc));
    t_rep = sum(S.duration_s);
    T.n_seizures_accepted = nnz(acc);
    T.n_seizures_reported = height(S);
    T.n_events_reported = numel(unique(S.seizure_id));
    T.n_rescued = nnz(strcmp(S.detection_status, 'rescued'));
    T.n_imputed = nnz(strcmp(S.detection_status, 'imputed'));
    T.total_seizure_time_s_accepted = t_acc;
    T.pct_time_in_seizure_accepted = 100 * t_acc / total_s;
    T.total_seizure_time_s_reported = t_rep;
    T.pct_time_in_seizure_reported = 100 * t_rep / total_s;
end

%% ======================================================================
function dirs = make_output_dirs(output_root)
    names = {'01_txt', '02_clean', '03_seizures', '04_iid', '05_summaries', 'logs'};
    for i = 1:numel(names)
        d = fullfile(output_root, names{i});
        if ~isfolder(d)
            mkdir(d);
        end
    end
    dirs = struct('root', output_root, 'txt', fullfile(output_root, '01_txt'), ...
        'clean', fullfile(output_root, '02_clean'), 'seizures', fullfile(output_root, '03_seizures'), ...
        'iid', fullfile(output_root, '04_iid'), 'summaries', fullfile(output_root, '05_summaries'), ...
        'logs', fullfile(output_root, 'logs'));
end

function fid = start_log(logs_dir)
    log_path = fullfile(logs_dir, sprintf('run_%s.log', datestr(now, 'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
    fid = fopen(log_path, 'w');
    log_line(fid, sprintf('run_pipeline_edf started %s', datestr(now))); %#ok<TNOW1,DATST>
end

function log_line(fid, msg)
    line = sprintf('[%s] %s', datestr(now, 'HH:MM:SS'), msg); %#ok<TNOW1,DATST>
    fprintf('%s\n', line);
    fprintf(fid, '%s\n', line);
end

function log_error(fid, source_file, region, stage, ME)
    if isempty(ME.stack)
        loc = 'n/a';
    else
        loc = sprintf('%s:%d', ME.stack(1).name, ME.stack(1).line);
    end
    log_line(fid, sprintf('  ERROR [%s | %s | %s] %s (at %s)', source_file, region, stage, ME.message, loc));
end

function save_config_used(cfg, output_root)
    save(fullfile(output_root, 'config_used.mat'), 'cfg');
    try
        json_txt = jsonencode(jsonify(cfg), 'PrettyPrint', true);
        fid = fopen(fullfile(output_root, 'config_used.json'), 'w');
        fprintf(fid, '%s', json_txt);
        fclose(fid);
    catch ME
        warning('run_pipeline_edf:ConfigJson', 'Could not write config_used.json: %s', ME.message);
    end
end

function v = jsonify(v)
% Recursively replace containers.Map (not supported by jsonencode) with a
% plain key/value struct array, so config_used.json can be written.
    if isa(v, 'containers.Map')
        k = keys(v); val = values(v);
        if isempty(k)
            v = struct('key', {}, 'value', {});
        else
            v = struct('key', k(:), 'value', val(:));
        end
    elseif isstruct(v)
        for fn = fieldnames(v)'
            for i = 1:numel(v)
                v(i).(fn{1}) = jsonify(v(i).(fn{1}));
            end
        end
    end
end

function files_out = discover_edf_files(folder)
    files = [dir(fullfile(folder, '*.edf')); dir(fullfile(folder, '*.EDF'))];
    if isempty(files)
        files_out = files;
        return;
    end
    full_paths = fullfile({files.folder}, {files.name});
    [~, ia] = unique(full_paths, 'stable');
    files_out = files(ia);
end

function [files, missing, desc] = resolve_edf_inputs(input)
% A folder -> every .edf in it. A cell array (or a single path ending in
% .edf) -> exactly those files; the ones that do not exist are returned in
% `missing` so the caller reports them instead of stopping.
    missing = {};
    if (ischar(input) || (isstring(input) && isscalar(input))) && isfolder(input)
        desc = char(input);
        files = discover_edf_files(desc);
        return;
    end
    if ischar(input) || isstring(input)
        input = cellstr(input);
        [~, ~, ext] = fileparts(input{1});
        if isscalar(input) && ~strcmpi(ext, '.edf')
            error('run_pipeline_edf:BadInput', 'Not a folder: %s', input{1});
        end
    end
    if ~iscell(input)
        error('run_pipeline_edf:BadInput', 'input must be a folder or a cell array of .edf paths.');
    end
    desc = sprintf('a list of %d EDF file(s)', numel(input));
    files = [];
    for i = 1:numel(input)
        d = dir(char(input{i}));
        if isempty(d) || d(1).isdir
            missing{end+1} = char(input{i}); %#ok<AGROW>
        else
            files = [files; d(1)]; %#ok<AGROW>
        end
    end
end

function t = zoned_nat(tz)
    t = NaT;
    t.TimeZone = tz;
end

function v = manifest_attenuation(row)
% The Excel's Attenuation cell edf_import.m put in the manifest ('' when
% the channel mode is not 'log').
    if ismember('attenuation', row.Properties.VariableNames)
        v = char(row.attenuation{1});
    else
        v = '';
    end
end

function note = gain_warning(case_spec, clean_data)
% qc_report warning when a gain case could not apply its automatic gain,
% or used the all-regions fallback reference (precondition_lfp.m).
    note = '';
    if ~strcmp(case_spec.gain_mode, 'auto') || ~isfield(clean_data, 'meta')
        return;
    end
    m = clean_data.meta;
    if strcmp(meta_char(m, 'reference_source'), 'none') || isnan(meta_num(m, 'gain_estimate_raw'))
        note = 'gain NOT applied (no reference amplitude, or no measurable amplitude in this channel)';
    elseif strcmp(meta_char(m, 'reference_source'), 'fallback')
        note = 'gain reference taken from all regions pooled (no non-attenuated channel of this region in the run)';
    end
end

function write_checkpoint(checkpoint_dir, f, source_name, n_before, varargin)
% Saves the result rows THIS file added (one .mat per file, so the cost
% does not grow with the run). Recovery only: deleted once 05_summaries/
% is written.
    names = {'seizure_events', 'seizure_summary', 'iid_events', 'iid_summary', 'iid_bursts', 'gaps', 'qc'};
    cp = struct('file_index', f, 'source_file', source_name);
    for k = 1:numel(names)
        parts = varargin{k};
        cp.(names{k}) = parts(n_before(k) + 1:end);
    end
    try
        if ~isfolder(checkpoint_dir)
            mkdir(checkpoint_dir);
        end
        save(fullfile(checkpoint_dir, sprintf('file_%03d.mat', f)), '-struct', 'cp');
    catch ME
        warning('run_pipeline_edf:CheckpointFailed', 'Could not save the checkpoint of %s: %s', source_name, ME.message);
    end
end

function n = count_errors(qc_report)
    try
        v = qc_report.n_errors;
        if iscell(v)
            v = str2double(v);
        end
        n = sum(v(~isnan(v)));
    catch
        n = NaN;
    end
end

%% ======================================================================
function [clean_data, clean_stats, status, notch_info] = run_clean_stage(raw_data, q, case_spec, cfg, clean_dir)
% Case-aware selective reprocessing: precondition_lfp.m is cheap (gain is
% just a scalar multiply or a no-op), so it always runs to find out what
% gain_applied WOULD be; only the potentially-expensive clean_lfp.m
% (outlier interpolation + notch filtfilt) is skipped, and only when the
% existing file's own header (peeked without reading its, possibly huge,
% data section) already matches on every field that would otherwise
% change: case_applied, gain_applied, gain_source, notch_applied.
    [~, name, ~] = fileparts(raw_data.file);
    expected_clean = fullfile(clean_dir, [name '_clean.txt']);

    cond = precondition_lfp(raw_data, q, case_spec, cfg);
    predicted_notch_active = resolve_notch_active(case_spec.notch_mode, q.quality_class);

    file_exists = exist(expected_clean, 'file') == 2;
    if file_exists && ~cfg.general.overwrite
        existing_meta = peek_header(expected_clean);
        if header_matches_case(existing_meta, case_spec, cond, predicted_notch_active)
            clean_data = load_lfp_txt(expected_clean);
            clean_stats = struct('auto_switch', false, 'switched_reason', '', 'outlier_pct', NaN);
            if isfield(clean_data.meta, 'outlier_pct')
                clean_stats.outlier_pct = str2double(clean_data.meta.outlier_pct);
            end
            status = 'skipped';
            notch_info = struct('applied', predicted_notch_active, ...
                'blocks_skipped', field_num(clean_data.meta, 'notch_blocks_skipped', 0));
            return;
        end
        status = 'reprocessed';
    elseif file_exists
        status = 'reprocessed';
    else
        status = 'new';
    end

    run_cfg = cfg;
    run_cfg.clean.output_dir = clean_dir;
    run_cfg.general.overwrite = true;  % already decided to (re)write above
    clean_result = clean_lfp(cond, run_cfg, case_spec);
    clean_data = load_lfp_txt(clean_result.txt_file);
    clean_stats = clean_result.stats;
    notch_info = clean_result.notch;
end

function meta = peek_header(path)
    fid = fopen(path, 'rt');
    if fid == -1
        meta = struct();
        return;
    end
    hdr = parse_header(fid);
    fclose(fid);
    meta = hdr.fields;
end

function tf = header_matches_case(meta, case_spec, cond, predicted_notch_active)
    if ~isfield(meta, 'case_applied')
        tf = false;
        return;
    end
    existing_gain_applied = str2double(field_char(meta, 'gain_applied', 'NaN'));
    existing_notch_applied = strcmpi(field_char(meta, 'notch_applied', 'false'), 'true');
    tf = strcmp(strtrim(meta.case_applied), case_spec.case_applied) && ...
         strcmp(field_char(meta, 'gain_source', 'off'), cond.gain_source) && ...
         gain_equal(existing_gain_applied, cond.gain_applied) && ...
         (existing_notch_applied == predicted_notch_active);
end

function tf = gain_equal(a, b)
    if isnan(a) && isnan(b)
        tf = true;
    else
        tf = abs(a - b) < 1e-9 * max(1, abs(b));
    end
end

function v = field_char(s, name, default)
    if isfield(s, name)
        v = strtrim(s.(name));
    else
        v = default;
    end
end

function v = field_num(s, name, default)
    if isfield(s, name)
        v = str2double(s.(name));
        if isnan(v)
            v = default;
        end
    else
        v = default;
    end
end

function session_start = pick_session_start(raw_data, clean_data, tz)
% Always returns a TimeZone-bearing datetime, even when the value itself
% is NaT (unknown session start) -- otherwise vertcat-ing datetime columns
% across files where some have a real session_start and others don't
% fails ("cannot concatenate a zoned datetime with an unzoned one").
    if ~isempty(clean_data)
        session_start = clean_data.session_start;
    elseif ~isempty(raw_data)
        session_start = raw_data.session_start;
    else
        session_start = NaT;
    end
    if isempty(session_start.TimeZone)
        session_start.TimeZone = tz;
    end
end

function v = get_or_zero(map, key)
    if isKey(map, key)
        v = map(key);
    else
        v = 0;
    end
end

function s = suggestion_note(case_spec)
    if isempty(case_spec.suggested_case) || strcmp(case_spec.suggested_case, case_spec.case_applied)
        s = '';
    else
        s = sprintf(' [suggested_case=%s]', case_spec.suggested_case);
    end
end

function print_case_summary(case_tally, clean_status_tally, n_unusable, n_suggested_mismatch)
    fprintf('\n--- case summary ---\n');
    ck = keys(case_tally);
    for i = 1:numel(ck)
        fprintf('  case %-12s : %d channel(s)\n', ck{i}, case_tally(ck{i}));
    end
    sk = keys(clean_status_tally);
    for i = 1:numel(sk)
        fprintf('  clean_lfp %-11s: %d channel(s)\n', sk{i}, clean_status_tally(sk{i}));
    end
    fprintf('  unusable (quality_class)      : %d channel(s)\n', n_unusable);
    fprintf('  suggested_case != case_applied: %d channel(s)\n', n_suggested_mismatch);
end

%% ======================================================================
function tf = compute_adjacent_to_gap(start_s, end_s, valid_mask, t_rel, tol_s)
    gap_idx = mask_to_segments(~valid_mask);
    n = numel(start_s);
    tf = false(n, 1);
    if isempty(gap_idx)
        return;
    end
    gap_starts_t = t_rel(gap_idx(:, 1));
    gap_ends_t = t_rel(gap_idx(:, 2));
    for i = 1:n
        near = any(abs(start_s(i) - gap_ends_t) <= tol_s) || any(abs(start_s(i) - gap_starts_t) <= tol_s) || ...
               any(abs(end_s(i) - gap_ends_t) <= tol_s) || any(abs(end_s(i) - gap_starts_t) <= tol_s);
        tf(i) = near;
    end
end

function T = build_seizure_event_rows(seizures, region, subject_id, session_start, source_file, valid_mask, t_rel, edge_trim_s, seizure_mode)
% seizure_mode ('legacy'|'robust') tags every row for traceability. The
% robust branch's .seizures carries 5 extra confidence columns (see
% detect_seizures_robust.m); the legacy branch's does not, so they are
% backfilled here (over_max_duration=false, the rest NaN) -- this keeps
% seizures_events.csv's original 12 columns and their VALUES exactly as
% they were for case='normal' (seizure_mode='legacy' throughout), while
% giving every row (legacy included) the same seizure_mode + confidence
% column set. See README.md "dos ramas" for why the file grows a column
% rather than legacy rows simply omitting it (mixed-mode files must vertcat).
    n = height(seizures);
    if n == 0
        T = empty_seizure_events_table(session_start.TimeZone);
        return;
    end
    T = seizures;
    T.Properties.VariableNames{'id'} = 'seizure_id';
    T.subject_id = repmat({subject_id}, n, 1);
    T.region = repmat({region}, n, 1);
    T.session_start = repmat(session_start, n, 1);
    T.source_file = repmat({source_file}, n, 1);
    T.adjacent_to_gap = compute_adjacent_to_gap(seizures.start_s, seizures.end_s, valid_mask, t_rel, edge_trim_s);
    T.seizure_mode = repmat({seizure_mode}, n, 1);

    confidence_cols = {'over_max_duration', 'll_ratio', 'peak_energy_ratio', 'hf_ratio_db', 'envelope_cv'};
    for i = 1:numel(confidence_cols)
        col = confidence_cols{i};
        if ~ismember(col, T.Properties.VariableNames)
            if strcmp(col, 'over_max_duration')
                T.(col) = false(n, 1);
            else
                T.(col) = nan(n, 1);
            end
        end
    end

    T = T(:, [{'subject_id', 'region', 'session_start', 'source_file', 'seizure_id', 'start_s', 'end_s', ...
        'duration_s', 'start_abs', 'end_abs', 'block_id', 'adjacent_to_gap', 'seizure_mode'}, confidence_cols]);
end

function T = build_seizure_summary_row(row, session_start, seizure_results, cfg)
    seizures = seizure_results.seizures;
    n_seizures = height(seizures);
    if n_seizures > 0
        total_seizure_time_s = sum(seizures.duration_s);
        mean_duration_s = mean(seizures.duration_s);
        min_duration_s = min(seizures.duration_s);
        max_duration_s = max(seizures.duration_s);
    else
        total_seizure_time_s = 0; mean_duration_s = NaN; min_duration_s = NaN; max_duration_s = NaN;
    end
    total_duration_min = row.total_duration_s / 60;
    valid_duration_min = row.valid_duration_s / 60;

    T = table(row.subject_id, row.region, session_start, row.source_file, ...
        total_duration_min, valid_duration_min, row.n_gaps, total_duration_min - valid_duration_min, ...
        n_seizures, total_seizure_time_s, 100 * total_seizure_time_s / (60 * total_duration_min), ...
        mean_duration_s, min_duration_s, max_duration_s, ...
        seizure_results.metrics.median_energy, seizure_results.metrics.threshold, seizure_results.metrics.pct_above, ...
        seizure_results.metrics.n_segments, seizure_results.metrics.n_rejected, ...
        cfg.seizure.bandpass_band(1), cfg.seizure.bandpass_band(2), cfg.seizure.power_exponent, ...
        cfg.seizure.window_sec, cfg.seizure.median_factor, cfg.seizure.min_seizure_duration, ...
        'VariableNames', {'subject_id', 'region', 'session_start', 'source_file', 'total_duration_min', 'valid_duration_min', ...
        'n_gaps', 'gap_duration_min', 'n_seizures', 'total_seizure_time_s', 'pct_time_in_seizure', ...
        'mean_duration_s', 'min_duration_s', 'max_duration_s', 'median_energy', 'threshold_value', 'pct_above_thr', ...
        'n_segments', 'n_rejected', 'bandpass_low', 'bandpass_high', 'power_exponent', 'window_s', 'median_factor', 'min_seizure_duration'});
end

function [in_burst, burst_id] = find_burst_membership(complex_start_s, burst_table)
    n = numel(complex_start_s);
    in_burst = false(n, 1);
    burst_id = nan(n, 1);
    for b = 1:height(burst_table)
        mask = complex_start_s >= burst_table.Start_s(b) & complex_start_s <= burst_table.End_s(b);
        in_burst(mask) = true;
        burst_id(mask) = b;
    end
end

function T = build_iid_event_rows(spike_complex_table, burst_table, region, subject_id, session_start, source_file)
    n = height(spike_complex_table);
    if n == 0
        T = empty_iid_events_table(session_start.TimeZone);
        return;
    end
    T = spike_complex_table;
    classification = repmat({'single'}, n, 1);
    classification(T.Is_polyspike) = {'polyspike'};
    T.classification = classification;

    [in_burst, burst_id] = find_burst_membership(T.Start_s, burst_table);
    T.in_burst = in_burst;
    T.burst_id = burst_id;

    T.subject_id = repmat({subject_id}, n, 1);
    T.region = repmat({region}, n, 1);
    T.session_start = repmat(session_start, n, 1);
    T.source_file = repmat({source_file}, n, 1);

    T.Properties.VariableNames{'Complex_ID'} = 'complex_id';
    T.Properties.VariableNames{'Start_s'} = 'start_s';
    T.Properties.VariableNames{'End_s'} = 'end_s';
    T.Properties.VariableNames{'Duration_ms'} = 'duration_ms';
    T.Properties.VariableNames{'N_spikes'} = 'n_spikes';
    T.Properties.VariableNames{'Max_amplitude_mV'} = 'max_amplitude_mV';
    T.Properties.VariableNames{'Mean_amplitude_mV'} = 'mean_amplitude_mV';

    T = T(:, {'subject_id', 'region', 'session_start', 'source_file', 'complex_id', 'start_s', 'end_s', ...
        'start_abs', 'end_abs', 'duration_ms', 'n_spikes', 'classification', 'max_amplitude_mV', 'mean_amplitude_mV', ...
        'in_burst', 'burst_id'});
end

function T = build_iid_summary_row(row, session_start, iid_results)
    s = iid_results.summary;
    T = table(row.subject_id, row.region, session_start, row.source_file, ...
        s.total_duration_min, s.analyzed_duration_min, s.excluded_duration_min, s.n_exclusion_zones, ...
        s.baseline_uV, s.lower_threshold_uV, s.upper_threshold_uV, ...
        s.total_peaks, s.total_complexes, s.n_single, s.n_polyspike, s.pct_polyspike, ...
        s.complexes_per_min, s.single_per_min, s.polyspikes_per_min, s.n_bursts, s.bursts_per_hour, s.mean_spikes_per_polyspike, ...
        'VariableNames', {'subject_id', 'region', 'session_start', 'source_file', 'total_duration_min', 'analyzed_duration_min', ...
        'excluded_duration_min', 'n_exclusion_zones', 'baseline_uV', 'lower_threshold_uV', 'upper_threshold_uV', ...
        'total_peaks', 'total_complexes', 'n_single', 'n_polyspike', 'pct_polyspike', 'complexes_per_min', 'single_per_min', ...
        'polyspikes_per_min', 'n_bursts', 'bursts_per_hour', 'mean_spikes_per_polyspike'});
end

function T = build_iid_burst_rows(burst_table, region, subject_id, source_file, tz)
    n = height(burst_table);
    if n == 0
        T = empty_iid_bursts_table(tz);
        return;
    end
    T = burst_table;
    T.subject_id = repmat({subject_id}, n, 1);
    T.region = repmat({region}, n, 1);
    T.Properties.VariableNames{'Start_s'} = 'start_s';
    T.Properties.VariableNames{'End_s'} = 'end_s';
    T.Properties.VariableNames{'N_complexes'} = 'n_complexes';
    T.Properties.VariableNames{'Duration_s'} = 'duration_s';
    T = T(:, {'subject_id', 'region', 'start_s', 'end_s', 'start_abs', 'end_abs', 'duration_s', 'n_complexes'});
end

%% ======================================================================
% qc_report row assembly. info is a struct with one field per qc_report
% column (see build_qc_row's VariableNames list) -- built either minimally
% (a whole EDF failed to import, nothing else ran) or fully (whatever
% stages of one channel ran, each field NaN/''/false if that stage never
% reached).
function info = qc_info_minimal(subject_id, region, source_file, session_start, stages, errors, warnings_list)
    info = qc_info_blank();
    info.subject_id = subject_id;
    info.region = region;
    info.source_file = source_file;
    info.session_start = session_start;
    info.stages_completed = strjoin(stages, ',');
    info.n_errors = numel(errors);
    info.error_messages = strjoin(errors, '; ');
    info.n_warnings = numel(warnings_list);
    info.warning_messages = strjoin(warnings_list, '; ');
end

function info = qc_info_blank()
    info = struct( ...
        'subject_id', '', 'region', '', 'source_file', '', 'session_start', NaT, ...
        'stages_completed', '', 'n_errors', 0, 'error_messages', '', 'n_blocks_rejected_short', NaN, ...
        'outlier_pct', NaN, 'nan_pct', NaN, 'n_warnings', 0, 'warning_messages', '', ...
        'case_applied', '', 'case_source', '', 'suggested_case', '', 'quality_class', '', ...
        'sigma_band_uV', NaN, 'reference_used', NaN, 'reference_source', '', 'gain_estimate_raw', NaN, ...
        'gain_applied', NaN, 'gain_source', '', 'sensitivity_equivalent_uV_per_mm', NaN, ...
        'line_ratio_db', NaN, 'line_ratio_p95_db', NaN, 'line_ratio_max_db', NaN, 'pct_time_line_high', NaN, ...
        'notch_applied', false, 'quantization_step_uV', NaN, 'snr_quantization_db', NaN, 'adc_codes_span', NaN, ...
        'pct_clipped', NaN, 'flat_fraction', NaN, 'notch_blocks_skipped', NaN, ...
        'seizure_threshold_mode', '', 'iid_threshold_mode', '', 'excel_attenuation', '');
end

function info = qc_info_full(row, session_start, stages, errors, warnings_list, clean_data, clean_stats, seizure_results, q, case_spec, notch_info, cfg)
    info = qc_info_blank();
    info.subject_id = row.subject_id{1};
    info.region = row.region{1};
    info.source_file = row.source_file{1};
    info.session_start = session_start;
    info.stages_completed = strjoin(stages, ',');
    info.n_errors = numel(errors);
    info.error_messages = strjoin(errors, '; ');
    info.n_warnings = numel(warnings_list);
    info.warning_messages = strjoin(warnings_list, '; ');
    info.nan_pct = 100 * (1 - row.n_valid_samples / row.n_samples);
    info.excel_attenuation = manifest_attenuation(row);

    if ~isempty(clean_stats)
        info.outlier_pct = clean_stats.outlier_pct;
    end
    if ~isempty(seizure_results)
        info.n_blocks_rejected_short = seizure_results.qc.n_blocks_rejected_short;
    end
    info.seizure_threshold_mode = cfg.seizure.threshold_mode;
    info.iid_threshold_mode = cfg.iid.threshold_mode;

    if ~isempty(q)
        info.quality_class = q.quality_class;
        info.sigma_band_uV = q.sigma_band_uV;
        info.line_ratio_db = q.line_ratio_db;
        info.line_ratio_p95_db = q.line_ratio_p95_db;
        info.line_ratio_max_db = q.line_ratio_max_db;
        info.pct_time_line_high = q.pct_time_line_high;
        info.quantization_step_uV = q.quantization_step_uV;
        info.snr_quantization_db = q.snr_quantization_db;
        info.adc_codes_span = q.adc_codes_span;
        info.pct_clipped = q.pct_clipped;
        info.flat_fraction = q.flat_fraction;
    end

    if ~isempty(case_spec)
        info.case_applied = case_spec.case_applied;
        info.case_source = case_spec.case_source;
        info.suggested_case = case_spec.suggested_case;
    end

    % gain/reference fields: clean_data is load_lfp_txt's own return value
    % (skipped or freshly reprocessed, either way re-read from the written
    % header), so .meta carries these as text regardless of which path was
    % taken -- no need to reach into precondition_lfp's intermediate struct.
    if ~isempty(clean_data) && isfield(clean_data, 'meta')
        m = clean_data.meta;
        info.gain_estimate_raw = meta_num(m, 'gain_estimate_raw');
        info.gain_applied = meta_num(m, 'gain_applied');
        info.gain_source = meta_char(m, 'gain_source');
        info.reference_used = meta_num(m, 'reference_used');
        info.reference_source = meta_char(m, 'reference_source');
        info.sensitivity_equivalent_uV_per_mm = meta_num(m, 'sensitivity_equivalent_uV_per_mm');
    end

    if ~isempty(notch_info)
        info.notch_applied = notch_info.applied;
        info.notch_blocks_skipped = notch_info.blocks_skipped;
    end
end

function v = meta_num(m, name)
    if isfield(m, name)
        v = str2double(m.(name));
    else
        v = NaN;
    end
end

function v = meta_char(m, name)
    if isfield(m, name)
        v = strtrim(m.(name));
    else
        v = '';
    end
end

function T = build_qc_row(info)
    T = table({info.subject_id}, {info.region}, {info.source_file}, info.session_start, ...
        {info.stages_completed}, info.n_errors, {info.error_messages}, info.n_blocks_rejected_short, ...
        info.outlier_pct, info.nan_pct, info.n_warnings, {info.warning_messages}, ...
        {info.case_applied}, {info.case_source}, {info.suggested_case}, {info.quality_class}, ...
        info.sigma_band_uV, info.reference_used, {info.reference_source}, info.gain_estimate_raw, ...
        info.gain_applied, {info.gain_source}, info.sensitivity_equivalent_uV_per_mm, ...
        info.line_ratio_db, info.line_ratio_p95_db, info.line_ratio_max_db, info.pct_time_line_high, ...
        info.notch_applied, info.quantization_step_uV, info.snr_quantization_db, info.adc_codes_span, ...
        info.pct_clipped, info.flat_fraction, info.notch_blocks_skipped, ...
        {info.seizure_threshold_mode}, {info.iid_threshold_mode}, {info.excel_attenuation}, ...
        'VariableNames', qc_column_names());
end

% qc_column_names, vertcat_or_empty, and every empty_*_table function
% used to be duplicated here as local functions, which silently SHADOWED
% the shared src/utils/ versions for every call made from within this
% file (MATLAB resolves a call to a file's own local function before its
% path) -- see write_all_summaries's history for how that class of bug
% was first found (a fix to the shared copy alone had no effect on this
% file's own output). Removed; this file now uses the shared
% src/utils/qc_column_names.m, vertcat_or_empty.m and empty_*_table.m
% (single source of truth, also used by read_pipeline_csv.m).

%% ======================================================================
function T = build_natus_review_sheet(seizures_events, iid_bursts, gaps_summary, file_subject_map, tz)
    rows_abs = [seizures_events.start_abs; iid_bursts.start_abs; gaps_summary.start_abs];
    rows_dur = [seizures_events.duration_s; iid_bursts.duration_s; gaps_summary.duration_s];
    rows_type = [repmat({'seizure'}, height(seizures_events), 1); ...
                 repmat({'iid_burst'}, height(iid_bursts), 1); ...
                 repmat({'gap'}, height(gaps_summary), 1)];
    rows_region = [seizures_events.region; iid_bursts.region; repmat({''}, height(gaps_summary), 1)];
    gap_subjects = cellfun(@(s) lookup_subject(s, file_subject_map), gaps_summary.source_file, 'UniformOutput', false);
    rows_subject = [seizures_events.subject_id; iid_bursts.subject_id; gap_subjects];

    n = numel(rows_abs);
    if n == 0
        T = table('Size', [0 7], ...
            'VariableTypes', {'datetime', 'cell', 'cell', 'double', 'cell', 'cell', 'cell'}, ...
            'VariableNames', {'abs_time', 'clock_time', 'event_type', 'duration_s', 'region', 'subject_id', 'natus_confirmed'});
        T.abs_time.TimeZone = tz;
        return;
    end

    [abs_sorted, order] = sort(rows_abs);
    clock_time = cellstr(string(abs_sorted, 'HH:mm:ss'));

    T = table(abs_sorted, clock_time, rows_type(order), rows_dur(order), rows_region(order), rows_subject(order), ...
        repmat({''}, n, 1), ...
        'VariableNames', {'abs_time', 'clock_time', 'event_type', 'duration_s', 'region', 'subject_id', 'natus_confirmed'});
end

function subj = lookup_subject(source_file, file_subject_map)
    if isKey(file_subject_map, source_file)
        subj = file_subject_map(source_file);
    else
        subj = '';
    end
end

