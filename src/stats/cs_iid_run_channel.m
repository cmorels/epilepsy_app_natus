function rr = cs_iid_run_channel(clean, seizure_zones, manual_zones, out_dir, opts)
% CS_IID_RUN_CHANNEL  Wrapper around src/detect_iid.m (unmodified) for the
% IID re-detection of section 5.3: runs detect_iid on an EXISTING
% *_clean.txt with an explicit list of exclusion zones instead of the
% seizure detector's list. Baseline, thresholds and every rate are
% recomputed by detect_iid on the samples these zones leave; nothing saved
% by the campaign (thresholds, masks) is reused.
%
%   rr = cs_iid_run_channel(clean, seizure_zones, manual_zones, out_dir, opts)
%
%   clean         : path of the *_clean.txt, or the struct load_lfp_txt.m
%                   returned for it
%   seizure_zones : table (id, start_s, end_s) in the file's t_rel seconds,
%                   WITHOUT buffer -> detect_iid's seizures argument (it adds
%                   cfg.iid.exclusion_buffer_s on both sides)
%   manual_zones  : Nx2 [start_s end_s] in t_rel seconds, buffer ALREADY
%                   added -> cfg.iid.exclusion_zones_manual (no buffer there)
%   out_dir       : folder for detect_iid's own outputs (created; existing
%                   files of the same channel are overwritten)
%   opts.cfg      : pipeline cfg (pipeline_config()); only cfg.iid is used
%   opts.figures  : false -> skip detect_iid's panorama / polyspike figures
%   opts.provenance : struct stored with the result (git hash, file SHA-1s)
%   opts.zone_info  : table describing the zones (stored with the result)
%
% Outputs in out_dir: detect_iid's <base>_IID_results.mat (written by
% detect_iid only when there is >= 1 complex) and, always,
% <base>_IID_rerun.mat with the struct rr returned here (summary, zones,
% complexes, bursts, real P97 of the included samples (cs_iid_p97.m),
% cfg.iid, provenance, timings).
%
% opts.figures = false: a temporary decimate_minmax.m that throws is put
% first on the path for the duration of the call. In detect_iid,
% decimate_minmax is called only as the first statement of its figure
% helper, inside a try/catch that runs AFTER every result is computed and
% the .mat saved; detect_iid catches the error (its warning
% detect_iid:FigureSaveFailed is silenced here) and returns the same
% results. Saves ~4 min per channel (the panorama loops over every
% polyspike against the whole signal).

    cfg = opts.cfg;
    cfg.iid.output_dir = out_dir;
    cfg.iid.exclusion_zones_manual = reshape(manual_zones, [], 2);
    cfg.general.overwrite = true;
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    t_load = 0;
    if ischar(clean) || isstring(clean)
        tl = tic;
        data = load_lfp_txt(clean);
        t_load = toc(tl);
    else
        data = clean;
    end

    if ~opts.figures
        tmp = tempname;
        mkdir(tmp);
        fid = fopen(fullfile(tmp, 'decimate_minmax.m'), 'w');
        fprintf(fid, 'function varargout = decimate_minmax(varargin) %%#ok<STOUT>\n');
        fprintf(fid, 'error(''cs_iid_run_channel:FiguresOff'', ''IID figures disabled by cs_iid_run_channel'');\nend\n');
        fclose(fid);
        addpath(tmp, '-begin');
        ws = warning('off', 'detect_iid:FigureSaveFailed');
        cleaner = onCleanup(@() restore_path(tmp, ws));
    end

    td = tic;
    r = detect_iid(data, seizure_zones, cfg);
    t_detect = toc(td);
    if ~opts.figures
        clear cleaner
    end
    % real P97 behind the quantized (and capped) baseline; the rebuilt mask
    % must give detect_iid's own analyzed duration
    q = cs_iid_p97(data, r.exclusion_table, cfg.iid.baseline.max_uV);
    if abs(q.analyzed_min - r.summary.analyzed_duration_min) > 1e-9
        error('cs_iid_run_channel:Mask', '%s: rebuilt inclusion mask %.6f min ~= detect_iid %.6f min', ...
            data.file, q.analyzed_min, r.summary.analyzed_duration_min);
    end

    [~, fname] = fileparts(data.file);
    base = strrep(fname, '_clean', '');
    rr = struct();
    rr.clean_file = data.file;
    rr.base = base;
    rr.summary = r.summary;
    rr.p97_uV = q.p97_uV;
    rr.frac_below_grid_top = q.frac_below_cap;
    rr.exclusion_table = r.exclusion_table;
    rr.spike_complex_table = r.spike_complex_table;
    rr.burst_table = r.burst_table;
    rr.seizure_zones = seizure_zones;
    rr.manual_zones = cfg.iid.exclusion_zones_manual;
    rr.zone_info = [];
    if isfield(opts, 'zone_info')
        rr.zone_info = opts.zone_info;
    end
    rr.cfg_iid = cfg.iid;
    rr.provenance = [];
    if isfield(opts, 'provenance')
        rr.provenance = opts.provenance;
    end
    rr.figures = opts.figures;
    rr.iid_results_mat = r.mat_file;
    rr.load_s = t_load;
    rr.detect_s = t_detect;
    rr.finished = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    rr.rerun_mat = fullfile(out_dir, [base '_IID_rerun.mat']);
    save(rr.rerun_mat, 'rr');
end

function restore_path(tmp, ws)
    rmpath(tmp);
    warning(ws);
    rmdir(tmp, 's');
end
