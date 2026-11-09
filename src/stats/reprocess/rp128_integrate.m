function varargout = rp128_integrate(step, cfg, varargin)
% RP128_INTEGRATE  Hooks of run_campaign_stats.m for
% cfg.include_128s_0729_reprocessed = true: the reprocessed A1 / A2 EDFs of
% 128-s 20260729 (batch 2026\128-s\reprocess_20260729_A1A2, campaign C2)
% replace the original EDF (exclusion ED-17); DX-21 stays (it covers the
% gap between the two EDFs). Nothing here runs with the flag false.
%
%   cfg                       = rp128_integrate('config', cfg)
%   snap                      = rp128_integrate('snapshot', cfg, snap)
%   [edf_tbl, ch_tbl]         = rp128_integrate('tables', cfg, edf_tbl, ch_tbl)
%   [det_rows, review, log]   = rp128_integrate('detections', cfg, edf_tbl, det_rows, review)
%   edf_lists                 = rp128_integrate('edf_lists', cfg, edf_lists)
%   RD                        = rp128_integrate('iid_plan', cfg, RD)
%   [Corig, flags]            = rp128_integrate('iid_original', cfg, P)
%   [edf_tbl_diag, flags]     = rp128_integrate('threshold_diag_edfs', cfg, edf_tbl)
%
% The reprocessed labels are read like any other unit's: the pipeline's
% own seizures_events_pipeline.csv (full precision) joined row by row with
% the labeled seizures_events.csv (AE user_check / AF comments), both
% copied into this run's labels_snapshot. cs_load_detections and
% cs_iid_original are reused unchanged on a small campaign-shaped view
% (<out_dir>/rp128_view) holding copies of those files.

    switch step
        case 'config'
            o = rp128_paths();
            cfg.rp128 = struct('dir', o.out, 'campaign', 'C2', 'mouse', o.subject, 'replaced_edf', o.stem, ...
                'edfs', {{o.seg.stem}}, 'ed_id', 'ED-17');
            cfg.edf_excluded(strcmp(cfg.edf_excluded(:, 1), cfg.rp128.ed_id), :) = [];
            varargout = {cfg};
        case 'snapshot'
            varargout = {snapshot(cfg, varargin{1})};
        case 'tables'
            [varargout{1}, varargout{2}] = tables(cfg, varargin{:});
        case 'detections'
            [varargout{1}, varargout{2}, varargout{3}] = detections(cfg, varargin{:});
        case 'edf_lists'
            L = varargin{1};
            r = cfg.rp128;
            k = find(strcmp(L.campaign, r.campaign) & strcmp(L.mouse, r.mouse) & strcmp(L.edf, r.replaced_edf));
            new = repmat(L(k(1), :), numel(r.edfs), 1);
            new.edf = r.edfs(:);
            L(k, :) = [];
            varargout = {[L; new]};
        case 'iid_plan'
            RD = varargin{1};
            k = is_rp(cfg, RD) & strcmp(RD.accion, 're-detectar');
            RD.motivo(k) = {'canal nuevo: 128-s 20260729 reprocesado (A1/A2), reemplaza a ED-17'};
            RD.explicado(k) = true;
            varargout = {RD};
        case 'iid_original'
            [varargout{1}, varargout{2}] = iid_original(cfg, varargin{1});
        case 'threshold_diag_edfs'
            E = varargin{1};
            k = is_rp(cfg, E);
            fl = cell(0, 3);
            if any(k)
                fl(end+1, :) = {'info', '6.3 sin EDF reprocesados', ['128-s 20260729 A1/A2 (reprocesados) no entran en el ' ...
                    'diagnostico informativo de umbrales 6.3: sus *_seizures.mat no estan en la carpeta de la C2']};
            end
            varargout = {E(~k, :), fl};
        otherwise
            error('rp128_integrate:Step', 'unknown step %s', step);
    end
end

function k = is_rp(cfg, T)
    r = cfg.rp128;
    k = strcmp(T.campaign, r.campaign) & strcmp(T.mouse, r.mouse) & ismember(T.edf, r.edfs);
end

% ---- label snapshot ---------------------------------------------------------
function snap = snapshot(cfg, snap)
    r = cfg.rp128;
    dst_dir = fullfile(cfg.labels_snapshot_dir, r.campaign, r.mouse);
    files = {'seizures_events.csv', 'rp128_seizures_events.csv'; 'seizures_events_pipeline.csv', 'rp128_seizures_events_pipeline.csv'};
    man = fullfile(cfg.labels_snapshot_dir, 'manifest.csv');
    for i = 1:size(files, 1)
        if any(strcmp(snap.campaign, r.campaign) & strcmp(snap.mouse, r.mouse) & strcmp(snap.file, files{i, 2}))
            continue;   % already in this run's snapshot (reused)
        end
        src = fullfile(r.dir, '05_summaries', files{i, 1});
        if ~isfile(src), error('rp128_integrate:NoLabels', 'missing %s', src); end
        if ~isfolder(dst_dir), mkdir(dst_dir); end
        dst = fullfile(dst_dir, files{i, 2});
        copyfile(src, dst);
        d = dir(src);
        row = {r.campaign, r.mouse, files{i, 2}, src, dst, sha1_file(dst), d.bytes, ...
            char(datetime(d.datenum, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
            char(datetime('now', 'TimeZone', cfg.tz, 'Format', 'yyyy-MM-dd HH:mm:ss Z'))};
        M = readtable(man, 'TextType', 'char', 'Delimiter', ',', 'DatetimeType', 'text');
        M = [M; cell2table(row, 'VariableNames', M.Properties.VariableNames)]; %#ok<AGROW>
        writetable(M, man);
        snap = [snap; cell2table([row, {false}], 'VariableNames', snap.Properties.VariableNames)]; %#ok<AGROW>
    end
end

function h = sha1_file(f)
    fid = fopen(f, 'r');
    b = fread(fid, Inf, '*uint8');
    fclose(fid);
    md = java.security.MessageDigest.getInstance('SHA-1');
    h = lower(reshape(dec2hex(typecast(md.digest(b), 'uint8'))', 1, []));
end

% ---- 02_clean -> ch_tbl / edf_tbl (same columns as cs_load_edf_table) ---------
function [edf_tbl, ch_tbl] = tables(cfg, edf_tbl, ch_tbl)
    r = cfg.rp128;
    drop = @(T) strcmp(T.campaign, r.campaign) & strcmp(T.mouse, r.mouse) & strcmp(T.edf, r.replaced_edf);
    ch_tbl(drop(ch_tbl), :) = [];
    edf_tbl(drop(edf_tbl), :) = [];
    rows = {};
    ff = dir(fullfile(r.dir, '02_clean', '*_clean.txt'));
    for k = 1:numel(ff)
        f = fullfile(ff(k).folder, ff(k).name);
        h = cs_read_clean_header(f);
        edf = regexprep(h.source_file, '\.EDF$', '', 'ignorecase');
        if ~ismember(edf, r.edfs), continue; end
        fs = str2double(h.fs);
        ns = str2double(h.n_samples);
        gaps = zeros(0, 2);
        ngaps = str2double(h.n_gaps);
        if ngaps > 0
            g = readtable(fullfile(r.dir, '01_txt', regexprep(ff(k).name, '_HPC[rl]_clean\.txt$', '_gaps.csv')));
            gaps = [g.start_s g.end_s];
        end
        rows(end+1, :) = {r.campaign, r.mouse, edf, h.region, f, str2double(h.session_start_unix), ...
            h.session_start_datetime, fs, ns, ns / fs, ngaps, {gaps}, getf(h, 'case_applied'), ...
            str2double(getf(h, 'gain_applied')), getf(h, 'quality_class')}; %#ok<AGROW>
    end
    C = cell2table(rows, 'VariableNames', ch_tbl.Properties.VariableNames);
    if height(C) ~= 2 * numel(r.edfs)
        error('rp128_integrate:Clean', '%d clean files for %d reprocessed EDFs in %s', height(C), numel(r.edfs), r.dir);
    end
    ch_tbl = [ch_tbl; C];
    for e = 1:numel(r.edfs)
        x = C(strcmp(C.edf, r.edfs{e}), :);
        s = min(x.start_posix);
        t = max(x.start_posix + x.duration_s);
        new = {r.campaign, r.mouse, r.edfs{e}, s, t, t - s, height(x), char(x.file(strcmp(x.region, 'HPCr'))), ...
            char(x.file(strcmp(x.region, 'HPCl'))), {vertcat(x.gaps{:})}, max(abs(x.start_posix - s)), true};
        edf_tbl = [edf_tbl; cell2table(new, 'VariableNames', edf_tbl.Properties.VariableNames)]; %#ok<AGROW>
    end
end

function v = getf(h, k)
    if isfield(h, k), v = h.(k); else, v = ''; end
end

% ---- detections: cs_load_detections on the view --------------------------------
function [det_rows, review, log] = detections(cfg, edf_tbl, det_rows, review)
    r = cfg.rp128;
    det_rows(strcmp(det_rows.campaign, r.campaign) & strcmp(det_rows.mouse, r.mouse) & ...
        strcmp(det_rows.edf, r.replaced_edf), :) = [];
    snapdir = fullfile(cfg.labels_snapshot_dir, r.campaign, r.mouse);
    pri = fullfile(snapdir, 'rp128_seizures_events_pipeline.csv');
    if numel(readlines(pri, 'EmptyLineRule', 'skip')) < 2
        % no detection at all (the 2026-10-30 run): nothing to join or label;
        % cs_load_detections cannot build its table from zero rows
        new = det_rows([], :);
        rev = review(1, :);
        rev.campaign = {r.campaign}; rev.mouse = {r.mouse};
        rev.n_rows = 0; rev.n_rows_labeled = 0;
        rev.ae_column = {'31:user_check'}; rev.af_column = {'32:comments'};
        rev.row_status = {'no_detections'};
        log = {};
    else
        view = make_view(cfg, {pri, fullfile('merged', '05_summaries', 'seizures_events.csv'); ...
            fullfile(snapdir, 'rp128_seizures_events.csv'), fullfile(r.mouse, '05_summaries', 'seizures_events.csv')});
        cv = cfg;
        cv.campaigns = struct('name', r.campaign, 'root', view);
        cv.label_file_override = cell(0, 3);
        cv.labels_snapshot_dir = '';
        [new, rev, log] = cs_load_detections(cv, edf_tbl);
    end
    if any(~ismember(new.edf, r.edfs))
        error('rp128_integrate:Edf', 'reprocessed labels hold rows of other EDFs');
    end
    det_rows = [det_rows; new];
    rev.label_file = {'rp128_seizures_events.csv (reprocess_20260729_A1A2)'};
    review = [review; rev];
    log = strcat({'rp128: '}, log);
    log{end+1} = sprintf('rp128: %s %s %s replaced by %s (%d rows, %d labeled)', r.campaign, r.mouse, r.replaced_edf, ...
        strjoin(r.edfs, ' + '), height(new), sum(~cellfun(@isempty, new.label)));
end

function view = make_view(cfg, pairs)
% copies of small files in a campaign-shaped folder: {source, path in the view}
    view = fullfile(cfg.out_dir, 'rp128_view');
    for i = 1:size(pairs, 1)
        dst = fullfile(view, pairs{i, 2});
        if ~isfolder(fileparts(dst)), mkdir(fileparts(dst)); end
        copyfile(pairs{i, 1}, dst);
    end
end

% ---- the campaigns' own IID results: the reprocessed channels from their run -------
function [Corig, flags] = iid_original(cfg, P)
    r = cfg.rp128;
    k = is_rp(cfg, P);
    [C1, flags] = cs_iid_original(cfg, P(~k, :));
    if ~any(k)
        Corig = C1;
        return;
    end
    pairs = {fullfile(r.dir, '05_summaries', 'iid_summary.csv'), fullfile('merged', '05_summaries', 'iid_summary.csv')};
    for e = 1:numel(r.edfs)
        mf = dir(fullfile(r.dir, '04_iid', r.edfs{e}, '*_IID_results.mat'));
        for i = 1:numel(mf)
            pairs(end+1, :) = {fullfile(mf(i).folder, mf(i).name), fullfile(r.mouse, '04_iid', r.edfs{e}, mf(i).name)}; %#ok<AGROW>
        end
    end
    view = make_view(cfg, pairs);
    cv = cfg;
    cv.campaigns = struct('name', r.campaign, 'root', view);
    [C2, f2] = cs_iid_original(cv, P(k, :));
    Corig = [C1; C2];
    Corig([find(~k); find(k)], :) = Corig;   % back to P's row order
    flags = [flags; f2];
end
