function [baseline_uV, lower_uV, upper_uV] = iids_compute_baseline(abs_mV_included, cfg)
% IIDS_COMPUTE_BASELINE  Verbatim copy of the local function compute_baseline of src/detect_iid.m
% (git HEAD 3672166), renamed only so it does not shadow anything on the path.
% Used by cs_iid_sens_channel.m (threshold sensitivity, v2 section 5.2);
% the lower threshold is already a parameter where it applies.

    if isempty(abs_mV_included)
        error('detect_iid:NoIncludedSamples', 'No included (non-excluded, non-gap) samples remain for baseline calculation.');
    end

    sv = cfg.iid.baseline.start_uV;
    Percent = 0;
    while Percent < cfg.iid.baseline.target_fraction && sv <= cfg.iid.baseline.max_uV
        thr_mV = sv / 1000;
        Percent = mean(abs_mV_included < thr_mV);
        if Percent < cfg.iid.baseline.target_fraction
            sv = sv + cfg.iid.baseline.step_uV;
        end
    end

    if sv > cfg.iid.baseline.max_uV
        baseline_uV = cfg.iid.baseline.max_uV - cfg.iid.baseline.step_uV;
    else
        baseline_uV = sv - cfg.iid.baseline.step_uV;
    end

    lower_uV = cfg.iid.lower_threshold_factor * baseline_uV;
    upper_uV = cfg.iid.upper_threshold_uV;
end
