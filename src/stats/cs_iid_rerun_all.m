function C = cs_iid_rerun_all(cfg, Z, worker, n_workers)
% CS_IID_RERUN_ALL  IID re-detection (section 5.3) on every channel of the
% plan (Z.channels.run), through cs_iid_run_channel.m, with the zones of
% cs_iid_zones.m. Outputs in cfg.iid_rerun_dirs.principal{1}/<mouse>/<EDF>/.
%
% Crisette: every channel a crisette_focal zone touches is run a second
% time WITHOUT the crisette zones (cfg.iid_rerun_dirs.sin_crisette{1}),
% only to count the IIDs that would fall inside them; the main run keeps
% them excluded.
%
% Incremental / resumable (cs_iid_previous.m): a channel already computed
% with exactly the same zone intervals and the same code -- in this run's
% folder or in an earlier run's (the other roots of cfg.iid_rerun_dirs,
% e.g. v1) -- is read back instead of re-run.
%
%   C = cs_iid_rerun_all(cfg, Z, 0)        run what is missing, then collect
%   cs_iid_rerun_all(cfg, Z, k, n)         worker k of n (k = 1..n): runs
%                                          only job j with mod(j-1, n)+1 == k,
%                                          collects nothing (C = [])
%
% C: one row per (channel, variant) with the summary, P97, complexes
% ([Start_s End_s N_spikes Is_polyspike Max_amplitude_mV], file-relative s),
% bursts ([Start_s End_s N_complexes]) and the exclusion intervals actually
% applied ([start_s end_s], from detect_iid's exclusion_table).

    if nargin < 3
        worker = 0;
        n_workers = 1;
    end
    pc = pipeline_config();
    if pc.iid.exclusion_buffer_s ~= cfg.iid_excl_buffer_s
        error('cs_iid_rerun_all:Buffer', 'pipeline_config exclusion_buffer_s ~= cfg.iid_excl_buffer_s');
    end
    if ~cfg.provenance.iid_code_unmodified
        error('cs_iid_rerun_all:IidCode', 'detect_iid or its dependencies differ from HEAD');
    end
    buf = cfg.iid_excl_buffer_s;
    opts = struct('cfg', pc, 'figures', cfg.iid_rerun_figures, 'provenance', cfg.provenance);

    P = Z.channels(Z.channels.run, :);
    is_cr = strcmp(Z.zones.rule, 'crisette');
    jobs = cell(0, 2);
    for i = 1:height(P)
        jobs(end+1, :) = {i, 'principal'}; %#ok<AGROW>
        if any(ismember(P.zone_ids{i}, Z.zones.zone_id(is_cr)))
            jobs(end+1, :) = {i, 'sin_crisette'}; %#ok<AGROW>
        end
    end

    rows = cell(size(jobs, 1), 1);
    for j = 1:size(jobs, 1)
        mine = worker == 0 || mod(j - 1, n_workers) + 1 == worker;
        if ~mine
            continue;
        end
        ch = P(jobs{j, 1}, :);
        variant = jobs{j, 2};
        zones = Z.zones;
        if strcmp(variant, 'sin_crisette')
            zones = zones(~is_cr, :);
        end
        [sz, man, ids] = cs_iid_channel_zones(zones, ch, buf);
        dirs = cfg.iid_rerun_dirs.(variant);   % {this run's root, earlier runs' roots...}
        out = fullfile(dirs{1}, ch.mouse{1}, ch.edf{1});

        [rr, src] = cs_iid_previous(dirs, ch, sz, man, cfg.provenance);
        status = 'reutilizado';
        if ~isempty(rr) && ~startsWith(src, dirs{1})
            status = 'reutilizado de una corrida anterior';
        end
        if isempty(rr)
            o = opts;
            o.zone_info = zones(ismember(zones.zone_id, ids), {'zone_id', 'rule', 'origin', 'start_local', 'end_local'});
            rr = cs_iid_run_channel(ch.file{1}, sz, man, out, o);
            status = 'corrido';
            fprintf('[%s] job %d/%d %s %s %s %s (%s): %d complejos, analizado %.1f min, P97 %.0f uV, carga %.0f s + deteccion %.0f s\n', ...
                char(datetime('now', 'Format', 'HH:mm:ss')), j, size(jobs, 1), ch.campaign{1}, ch.mouse{1}, ch.edf{1}(1:min(10, end)), ...
                ch.region{1}, variant, rr.summary.total_complexes, rr.summary.analyzed_duration_min, rr.p97_uV, rr.load_s, rr.detect_s);
        end
        if worker == 0
            rows{j} = summarize(rr, ch, variant, status, ids);
        end
    end
    if worker ~= 0
        C = [];
        return;
    end
    C = vertcat(rows{:});
end

function T = summarize(rr, ch, variant, status, ids)
    s = rr.summary;
    cx = rr.spike_complex_table;
    b = rr.burst_table;
    e = rr.exclusion_table;
    T = table(ch.campaign, ch.mouse, ch.edf, ch.region, {variant}, {status}, ch.start_posix, ch.duration_s, ...
        s.baseline_uV, s.lower_threshold_uV, rr.p97_uV, rr.frac_below_grid_top, s.total_duration_min, s.analyzed_duration_min, ...
        s.total_complexes, s.n_single, s.n_polyspike, s.n_bursts, ...
        {[cx.Start_s cx.End_s cx.N_spikes double(cx.Is_polyspike) cx.Max_amplitude_mV]}, ...
        {[b.Start_s b.End_s b.N_complexes]}, {[e.start_s e.end_s]}, {ids}, rr.load_s, rr.detect_s, {rr.rerun_mat}, ...
        'VariableNames', {'campaign', 'mouse', 'edf', 'region', 'variant', 'status', 'start_posix', 'duration_s', ...
        'baseline_uV', 'lower_uV', 'p97_uV', 'frac_below_130uV', 'total_min', 'analyzed_min', ...
        'n_complexes', 'n_single', 'n_polyspike', 'n_bursts', 'complexes', 'bursts', 'excl', 'zone_ids', ...
        'load_s', 'detect_s', 'rerun_mat'});
end
