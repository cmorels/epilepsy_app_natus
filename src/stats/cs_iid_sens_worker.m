function cs_iid_sens_worker(out_dir, k, n)
% CS_IID_SENS_WORKER  One background MATLAB process of phase C (v2): worker
% k of n over the IID plan of phase B (out_dir/phase2_results.mat: cfg, Z;
% out_dir/phase3_results.mat: C). Results are cached per channel in
% out_dir/cache/sens/; run_campaign_stats.m then collects them.
%
%   matlab -batch "restoredefaultpath; addpath(genpath('<repo>/src')); cs_iid_sens_worker('<out_dir>', 1, 4)"

    q = load(fullfile(out_dir, 'phase2_results.mat'), 'cfg', 'Z');
    c = load(fullfile(out_dir, 'phase3_results.mat'), 'C');
    cfg = q.cfg;
    cfg.run_started = char(datetime('now', 'TimeZone', cfg.tz, 'Format', 'yyyy-MM-dd HH:mm:ss Z'));
    [~, cfg.provenance] = cs_provenance(cfg);
    fprintf('cs_iid_sens_worker %d/%d: started %s\n', k, n, cfg.run_started);
    cs_iid_sens_all(cfg, q.Z, c.C, k, n);
    fprintf('cs_iid_sens_worker %d/%d: done %s\n', k, n, char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
end
