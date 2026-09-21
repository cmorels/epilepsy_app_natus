function result = clean_lfp(data, cfg)
% CLEAN_LFP  Outlier removal stage. Ports the amplitude-threshold outlier
% logic from complete_pipeline_seizures.m (STEP 2) verbatim: same fixed
% thresholds, same k_factor/MAD auto thresholds, same >0.1% auto-switch
% rule, same linear interpolation of flagged samples.
%
% Adaptation for gaps (NaN): every statistic (median, MAD, std, range) is
% computed over valid (non-NaN) samples only, gap samples are never
% flagged as outliers, and outlier interpolation is done independently
% within each contiguous valid block so it can never reach across a gap
% (see interpolate_block_outliers below). Gap NaNs are copied through to
% the output unchanged.
%
%   result = clean_lfp(data, cfg)
%
%   data : struct from load_lfp_txt.m (loaded from an edf_import.m txt,
%          or any txt following the same header convention)
%   cfg  : struct from pipeline_config.m
%
% OUTPUT (struct result):
%   .signal, .valid_mask, .fs, .t_rel, .session_start, .meta, .file
%   .txt_file   : path written (or that would be written, if skipped)
%   .stats      : baseline_uV, noise_uV, signal_std, threshold_type,
%                 lower_threshold_uV, upper_threshold_uV, k_factor,
%                 auto_switch, switched_reason, n_outliers, outlier_pct,
%                 signal_range_original, signal_range_clean
%   .outlier_info : [global_idx, original_value] for every flagged sample

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

    clean_valid = clean_signal(valid_mask);
    clean_min = min(clean_valid);
    clean_max = max(clean_valid);

    fprintf('clean_lfp: %s -> %d/%d valid samples flagged (%.4f%%), thresholds=[%.1f, %.1f] uV (%s)\n', ...
        data.file, n_outliers, n_valid, outlier_pct, min_allowed, max_allowed, threshold_type);

    %% ---- header: carry input header, append cleaning fields --------------
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
        write_lfp_txt(txt_path, header_pairs, clean_signal);
    end
    copy_sibling_gaps_csv(data.folder, out_dir, cfg.general.overwrite);

    %% ---- assemble result ---------------------------------------------
    result = struct();
    result.signal = clean_signal;
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
