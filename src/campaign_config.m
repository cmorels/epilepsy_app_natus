function cfg = campaign_config()
% CAMPAIGN_CONFIG  Configuration of the EDF campaign (decided 2026-10-03).
% Starts from pipeline_config() -- whose defaults stay untouched, so any
% other run is reproducible as before -- and sets only what the campaign
% changes. run_all_subjects.m uses it; config_used.json of every run
% records the full result.
%
%   1. Seizures: robust branch (detect_seizures_robust.m) on EVERY channel;
%      the legacy detector is not used.
%   2. Bilateral reconciliation ON (rescue_and_impute), review band OFF
%      (single ll_ratio cut at cfg.seizure_robust.ll_threshold = 1.75).
%      The two cannot be combined (see README.md).
%   3. IID exclusion zones cover reconciled seizures in BOTH channels
%      (accepted, rescued and imputed rows).
%   4. 50 Hz notch on EVERY channel, whatever its measured line noise,
%      stopband 49-51 Hz. It is applied in clean_lfp.m, so seizure AND IID
%      detection both work on the notched signal.
%   5. Gain from the Excel (EEG_recording_log.xlsx), column Attenuation:
%      severe / mild / "unsure, perhaps mild" (or any unknown value) ->
%      attenuated -> notch + automatic gain; no / none / empty -> notch only.
%      Applies to both channels of the animal in that EDF.
%   6. Automatic gain reference: median amplitude (sigma_band_uV) of the
%      NON-attenuated channels of the same run, per region, computed by
%      compute_run_reference.m before processing (run_all_subjects.m sets
%      cfg.quality.reference_sigma_uV / reference_sigma_fallback_uV). If a
%      region has none, all regions pooled; if there are none at all, no
%      gain (warning in qc_report.csv). The run never stops for this.
%   7. A complete import of an EDF is reused instead of re-read (the
%      reference pre-pass imports every EDF once).

    cfg = pipeline_config();

    % Channels per (animal, EDF) from the Excel.
    cfg.edf.channels.mode = 'log';
    cfg.edf.reuse_import = true;

    % 1 + 4 + 5: one profile per Excel answer, both robust + notch.
    cfg.cases.profiles.excel_clean = struct('gain_mode', 'off', 'notch_mode', 'on', 'seizure_mode', 'robust');
    cfg.cases.profiles.excel_attenuated = struct('gain_mode', 'auto', 'notch_mode', 'on', 'seizure_mode', 'robust');
    cfg.cases.from_excel = true;
    cfg.cases.excel.case_clean = 'excel_clean';
    cfg.cases.excel.case_attenuated = 'excel_attenuated';
    cfg.cases.excel.attenuation_no = {'', 'no', 'none'};
    cfg.cases.excel.attenuation_yes = {'severe', 'mild', 'unsure, perhaps mild'};

    % 4: notch 50 +/- 1 Hz, no harmonics.
    cfg.precondition.notch_freqs = 50;
    cfg.precondition.notch_harmonics = [];
    cfg.precondition.notch_halfwidth_hz = 1;

    % 2: reconciliation on, review band off.
    cfg.seizure_robust.ll_accept = [];
    cfg.seizure_robust.ll_reject = [];
    cfg.bilateral.rescue_mode = 'rescue_and_impute';

    % 3: IID excludes reconciled seizures in both channels.
    cfg.bilateral.exclude_rescued_from_iid = true;
end
