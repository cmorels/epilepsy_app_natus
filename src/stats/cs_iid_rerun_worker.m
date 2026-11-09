function cs_iid_rerun_worker(results_mat, k, n)
% CS_IID_RERUN_WORKER  Entry point of one background MATLAB process of the
% IID re-detection: worker k of n over the plan saved by phase 2
% (phase2_results.mat: cfg, Z). Provenance is recomputed here, at the time
% the worker runs. run_campaign_stats.m (phase 3) then reuses every
% finished channel and runs whatever is still missing.
%
%   matlab -batch "restoredefaultpath; addpath(genpath('<repo>/src')); cs_iid_rerun_worker('<out>/phase2_results.mat', 1, 3)"

    q = load(results_mat, 'cfg', 'Z');
    cfg = q.cfg;
    cfg.run_started = char(datetime('now', 'TimeZone', cfg.tz, 'Format', 'yyyy-MM-dd HH:mm:ss Z'));
    [~, cfg.provenance] = cs_provenance(cfg);
    fprintf('cs_iid_rerun_worker %d/%d: HEAD %s, started %s\n', k, n, cfg.provenance.git_head, cfg.run_started);
    cs_iid_rerun_all(cfg, q.Z, k, n);
    fprintf('cs_iid_rerun_worker %d/%d: done %s\n', k, n, char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
end
