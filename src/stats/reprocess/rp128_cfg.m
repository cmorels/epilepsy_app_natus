function cfg = rp128_cfg(o, log_file)
% RP128_CFG  Configuration of the C2 for 128-s, built exactly as
% resume_campaign_20261014.m / rerun_077s_campaign_20261014.m did
% (campaign_config + recording log + the gain reference SAVED by the C2,
% injected explicitly; compute_run_reference is NOT called), with this
% run's recording log and output folder.

    cfg = campaign_config();
    cfg.edf.channels.log_file = log_file;
    R = readtable(o.c2_ref_csv, 'TextType', 'char');
    is_fb = strcmp(R.region, '(all regions, fallback)');
    ref_map = containers.Map(R.region(~is_fb), num2cell(R.reference_uV(~is_fb)));
    if ref_map.Count > 0
        cfg.quality.reference_sigma_uV = ref_map;
    end
    cfg.quality.reference_sigma_fallback_uV = R.reference_uV(is_fb);
    cfg.edf.subject_id = o.subject;
    cfg.paths.output_root = o.out;
end
