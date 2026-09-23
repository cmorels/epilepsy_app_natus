function result = clean_lfp(data, cfg, case_spec)
% CLEAN_LFP  Outlier removal, THEN mains notch (see below), THEN write.
% Outlier removal ports the amplitude-threshold logic from
% complete_pipeline_seizures.m (STEP 2) verbatim: same fixed thresholds,
% same k_factor/MAD auto thresholds, same >0.1% auto-switch rule, same
% linear interpolation of flagged samples.
%
% Adaptation for gaps (NaN): every statistic (median, MAD, std, range) is
% computed over valid (non-NaN) samples only, gap samples are never
% flagged as outliers, and outlier interpolation is done independently
% within each contiguous valid block so it can never reach across a gap
% (see interpolate_block_outliers below). Gap NaNs are copied through to
% the output unchanged.
%
%   result = clean_lfp(data, cfg)
%   result = clean_lfp(data, cfg, case_spec)
%
%   data      : struct from load_lfp_txt.m (loaded from an edf_import.m
%               txt, or any txt following the same header convention).
%               When it has already been through precondition_lfp.m, its
%               gain_applied/quality/etc. fields are carried into this
%               clean file's header; when called standalone without one
%               (the 2-arg form), the case system is treated as inactive
%               ('normal': gain off, notch off) -- unchanged from before
%               the case system existed.
%   cfg       : struct from pipeline_config.m
%   case_spec : struct from resolve_case.m (gain_mode, notch_mode,
%               case_applied, case_source, suggested_case -- only
%               notch_mode and the reporting fields are used here; gain
%               was already applied upstream by precondition_lfp.m).
%               Defaults to the inactive 'normal' case if omitted.
%
% Outliers are removed BEFORE the notch, not after: large transients
% (movement artifacts, the outliers themselves) can ring a narrow IIR
% bandstop filter, so removing them first keeps that ringing from ever
% being introduced.
%
% OUTPUT (struct result):
%   .signal, .valid_mask, .fs, .t_rel, .session_start, .meta, .file
%   .txt_file   : path written (or that would be written, if skipped)
%   .stats      : baseline_uV, noise_uV, signal_std, threshold_type,
%                 lower_threshold_uV, upper_threshold_uV, k_factor,
%                 auto_switch, switched_reason, n_outliers, outlier_pct,
%                 signal_range_original, signal_range_clean
%   .outlier_info : [global_idx, original_value] for every flagged sample
%   .notch      : applied, freqs_used, blocks_skipped

    if nargin < 3 || isempty(case_spec)
        case_spec = struct('gain_mode', 'off', 'gain', NaN, 'notch_mode', 'off', ...
            'case_applied', 'normal', 'case_source', 'default', 'suggested_case', '');
    end

    signal = data.signal(:);
    fs = data.fs;
    valid_mask = data.valid_mask(:);
    valid_signal = signal(valid_mask);
    n_valid = numel(valid_signal);

    if n_valid < 2
        error('clean_lfp:NotEnoughValidSamples', ...
            '%s has fewer than 2 valid (non-gap) samples; cannot clean.', data.file);
    end

    %% ---- signal statistics (valid samples only) --------------------------
    baseline_uV = median(valid_signal);
    noise_uV = mad(valid_signal, 1);
    if noise_uV == 0
        warning('clean_lfp:ZeroMAD', 'MAD = 0 for %s; using minimum noise floor (1 uV).', data.file);
        noise_uV = 1;
    end
    signal_std = std(valid_signal);
    obs_min = min(valid_signal);
    obs_max = max(valid_signal);

    min_auto_uV = baseline_uV - cfg.clean.k_factor * noise_uV;
    max_auto_uV = baseline_uV + cfg.clean.k_factor * noise_uV;

    %% ---- decide thresholds (auto-switch logic, ported as-is) -------------
    trigger_fixed = cfg.clean.use_fixed_thresholds;
    switched_reason = '';

    if ~cfg.clean.use_fixed_thresholds
        [~, test_info] = interpolate_block_outliers(signal, valid_mask, min_auto_uV, max_auto_uV);
        outlier_pct_test = 100 * size(test_info, 1) / n_valid;
        if outlier_pct_test > cfg.clean.auto_switch_pct
            trigger_fixed = true;
            switched_reason = sprintf('Auto-switched (%.2f%% > %.2f%%)', outlier_pct_test, cfg.clean.auto_switch_pct);
        end
    end

    if trigger_fixed
        min_allowed = cfg.clean.fixed_lower_uV;
        max_allowed = cfg.clean.fixed_upper_uV;
        if isempty(switched_reason)
            threshold_type = 'FIXED (user config)';
        else
            threshold_type = 'FIXED';
        end
    else
        min_allowed = min_auto_uV;
        max_allowed = max_auto_uV;
        threshold_type = sprintf('AUTOMATIC (k=%.1f)', cfg.clean.k_factor);
    end

    %% ---- apply final thresholds, block-respecting interpolation ----------
    [clean_signal, outlier_info] = interpolate_block_outliers(signal, valid_mask, min_allowed, max_allowed);
    n_outliers = size(outlier_info, 1);
    outlier_pct = 100 * n_outliers / n_valid;

    fprintf('clean_lfp: %s -> %d/%d valid samples flagged (%.4f%%), thresholds=[%.1f, %.1f] uV (%s)\n', ...
        data.file, n_outliers, n_valid, outlier_pct, min_allowed, max_allowed, threshold_type);

    %% ---- notch (after outliers -- see header comment for why) ------------
    if strcmp(case_spec.notch_mode, 'auto') && ~(isfield(data, 'quality') && isfield(data.quality, 'quality_class'))
        warning('clean_lfp:AutoNotchNoQuality', ...
            '%s: notch_mode=''auto'' but no quality info is available (run precondition_lfp.m first); leaving notch off.', ...
            data.file);
        notch_active = false;
    else
        notch_active = resolve_notch_active(case_spec.notch_mode, quality_or_defaults(data).quality_class);
    end

    if notch_active
        [final_signal, notch_blocks_skipped, notch_freqs_used] = apply_notch_blocks(clean_signal, valid_mask, fs, cfg);
        fprintf('clean_lfp: %s -> notch applied at [%s] Hz (%d block(s) skipped, too short)\n', ...
            data.file, num2str(notch_freqs_used), notch_blocks_skipped);
    else
        final_signal = clean_signal;
        notch_blocks_skipped = 0;
        notch_freqs_used = [];
    end

    clean_valid = final_signal(valid_mask);
    clean_min = min(clean_valid);
    clean_max = max(clean_valid);

    %% ---- header: carry input header, append cleaning + case/quality/notch fields
    q = quality_or_defaults(data);
    clean_pairs = { ...
        'clean_source_file',      data.file; ...
        'threshold_type',         threshold_type; ...
        'k_factor',                cfg.clean.k_factor; ...
        'auto_switch',             ~isempty(switched_reason); ...
        'auto_switch_reason',      pick(switched_reason, 'none'); ...
        'lower_threshold_uV',      min_allowed; ...
        'upper_threshold_uV',      max_allowed; ...
        'n_outliers',              n_outliers; ...
        'outlier_pct',             outlier_pct; ...
        'signal_range_original_uV', sprintf('[%.6g, %.6g]', obs_min, obs_max); ...
        'signal_range_clean_uV',   sprintf('[%.6g, %.6g]', clean_min, clean_max); ...
        'case_applied',            case_spec.case_applied; ...
        'case_source',             case_spec.case_source; ...
        'suggested_case',          case_spec.suggested_case; ...
        'gain_mode',               case_spec.gain_mode; ...
        'notch_mode',              case_spec.notch_mode; ...
        'gain_applied',            field_or(data, 'gain_applied', 1); ...
        'gain_estimate_raw',       field_or(data, 'gain_estimate_raw', NaN); ...
        'gain_source',             field_or(data, 'gain_source', 'off'); ...
        'reference_used',          field_or(data, 'reference_used', NaN); ...
        'reference_source',        field_or(data, 'reference_source', 'off'); ...
        'sensitivity_equivalent_uV_per_mm', field_or(data, 'sensitivity_equivalent_uV_per_mm', 100); ...
        'sigma_band_uV_original',  q.sigma_band_uV; ...
        'mad_uV_original',         q.mad_uV; ...
        'line_ratio_db',           q.line_ratio_db; ...
        'line_ratio_p95_db',       q.line_ratio_p95_db; ...
        'pct_time_line_high',      q.pct_time_line_high; ...
        'quality_class',           q.quality_class; ...
        'quantization_step_uV',    q.quantization_step_uV; ...
        'snr_quantization_db',     q.snr_quantization_db; ...
        'adc_codes_span',          q.adc_codes_span; ...
        'pct_clipped',             q.pct_clipped; ...
        'flat_fraction',           q.flat_fraction; ...
        'notch_applied',           notch_active; ...
        'notch_freqs',             num2str(notch_freqs_used); ...
        'notch_halfwidth_hz',      cfg.precondition.notch_halfwidth_hz; ...
        'notch_order',             cfg.precondition.notch_order; ...
        'notch_blocks_skipped',    notch_blocks_skipped; ...
    };
    header_pairs = [data.header_pairs; clean_pairs];

    %% ---- resolve output path, copy sibling gaps CSV if unambiguous -------
    out_dir = cfg.clean.output_dir;
    if isempty(out_dir)
        out_dir = pwd;
    end
    if ~isfolder(out_dir)
        mkdir(out_dir);
    end

    [~, name, ~] = fileparts(data.file);
    out_name = [name '_clean.txt'];
    txt_path = fullfile(out_dir, out_name);

    if cfg.general.overwrite || exist(txt_path, 'file') ~= 2
        write_lfp_txt(txt_path, header_pairs, final_signal);
    end
    copy_sibling_gaps_csv(data.folder, out_dir, cfg.general.overwrite);

    %% ---- assemble result ---------------------------------------------
    result = struct();
    result.signal = final_signal;
    result.valid_mask = valid_mask;
    result.fs = fs;
    result.t_rel = data.t_rel;
    result.session_start = data.session_start;
    result.meta = data.meta;
    result.file = data.file;
    result.txt_file = txt_path;
    result.outlier_info = outlier_info;
    result.stats = struct( ...
        'baseline_uV', baseline_uV, ...
        'noise_uV', noise_uV, ...
        'signal_std', signal_std, ...
        'threshold_type', threshold_type, ...
        'lower_threshold_uV', min_allowed, ...
        'upper_threshold_uV', max_allowed, ...
        'k_factor', cfg.clean.k_factor, ...
        'auto_switch', ~isempty(switched_reason), ...
        'switched_reason', switched_reason, ...
        'n_outliers', n_outliers, ...
        'outlier_pct', outlier_pct, ...
        'signal_range_original_uV', [obs_min, obs_max], ...
        'signal_range_clean_uV', [clean_min, clean_max]);
    result.case_spec = case_spec;
    result.notch = struct('applied', notch_active, 'freqs_used', notch_freqs_used, 'blocks_skipped', notch_blocks_skipped);
end

%% ======================================================================
function [clean_signal, outlier_info] = interpolate_block_outliers(signal, valid_mask, min_allowed, max_allowed)
% Same amplitude-threshold + linear-interpolation rule as the original
% remove_outliers_amplitude, applied independently within each contiguous
% valid (non-gap) block so interpolation never reaches across a NaN gap.

    clean_signal = signal;
    outlier_info = zeros(0, 2);

    blocks = mask_to_segments(valid_mask);
    for i = 1:size(blocks, 1)
        s = blocks(i, 1);
        e = blocks(i, 2);
        blk = signal(s:e);

        is_out = blk < min_allowed | blk > max_allowed;
        out_idx = find(is_out);
        if isempty(out_idx)
            continue;
        end

        good_idx = find(~is_out);
        if numel(good_idx) < 2
            warning('clean_lfp:TooFewValidPoints', ...
                'Block [%d %d] has too few in-range points to interpolate; left unchanged.', s, e);
            continue;
        end

        interp_vals = interp1(good_idx, blk(good_idx), out_idx, 'linear', 'extrap');
        clean_signal(s - 1 + out_idx) = interp_vals;
        outlier_info = [outlier_info; s - 1 + out_idx, blk(out_idx)]; %#ok<AGROW>
    end
end

function copy_sibling_gaps_csv(src_folder, dst_folder, overwrite)
    matches = dir(fullfile(src_folder, '*_gaps.csv'));
    if numel(matches) ~= 1
        return;  % none, or ambiguous -- not fatal, just skip (best-effort convenience copy)
    end
    dst_path = fullfile(dst_folder, matches(1).name);
    if overwrite || exist(dst_path, 'file') ~= 2
        copyfile(fullfile(matches(1).folder, matches(1).name), dst_path);
    end
end

function v = pick(value, default_if_empty)
    if isempty(value)
        v = default_if_empty;
    else
        v = value;
    end
end

function v = field_or(s, name, default)
    if isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end

function q = quality_or_defaults(data)
% precondition_lfp.m attaches data.quality (a signal_quality.m struct);
% when clean_lfp.m is called standalone without it (the 2-arg form), the
% quality header fields are written as NaN/'' rather than left out, so
% every clean.txt has the same header shape regardless of how it was made.
    if isfield(data, 'quality')
        q = data.quality;
        return;
    end
    q = struct('sigma_band_uV', NaN, 'mad_uV', NaN, 'line_ratio_db', NaN, ...
        'line_ratio_p95_db', NaN, 'pct_time_line_high', NaN, 'quality_class', '', ...
        'quantization_step_uV', NaN, 'snr_quantization_db', NaN, 'adc_codes_span', NaN, ...
        'pct_clipped', NaN, 'flat_fraction', NaN);
end
