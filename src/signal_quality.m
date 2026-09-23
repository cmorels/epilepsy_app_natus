function q = signal_quality(data, cfg)
% SIGNAL_QUALITY  Measures signal quality on a raw (pre-gain) channel.
% Never modifies the signal -- measuring is free and always runs, on
% every channel, regardless of what (if anything) the case system later
% decides to do about what it measures (see README.md "sistema de casos").
%
%   q = signal_quality(data, cfg)
%
%   data : struct from load_lfp_txt.m (the RAW, not gain-adjusted, signal)
%   cfg  : struct from pipeline_config.m (uses cfg.quality.*)
%
% OUTPUT (struct q): sigma_band_uV, mad_uV, line_ratio_db,
% line_ratio_p95_db, line_ratio_max_db, pct_time_line_high,
% line_ratio_series (struct with .t_s, .ratio_db), quantization_step_uV,
% quantization_step_source ('header'|'estimated'), snr_quantization_db,
% adc_codes_span, n_distinct_values, pct_clipped, pct_nan, flat_fraction,
% quality_class.
%
% sigma_band_uV design (see README.md for the full rationale):
%   - restricted to cfg.quality.scale_band ([15 70] Hz default): measures
%     what the pipeline actually analyzes, not the whole spectrum;
%   - excludes cfg.quality.line_exclusion_band ([48 52] Hz default): makes
%     the metric immune to mains contamination, so attenuation and line
%     contamination can be diagnosed independently even when both are
%     present at once;
%   - median ACROSS windows, not a single estimate over the whole
%     recording: a seizure or a handful of spikes can't inflate it, so it
%     represents the background, not the events.

    signal = data.signal(:);
    fs = data.fs;
    valid_mask = data.valid_mask(:);

    window_samples = round(cfg.quality.welch_window_s * fs);
    noverlap_samples = round(window_samples / 2);
    win = hamming(window_samples);

    line_band = cfg.quality.line_freq_hz + [-1 1];
    lower_neighbor = cfg.quality.line_freq_hz + [-5 -2];
    upper_neighbor = cfg.quality.line_freq_hz + [2 5];

    raw_blocks = mask_to_segments(valid_mask);
    all_p_scale = zeros(1, 0);
    all_ratio_db = zeros(1, 0);
    all_t = zeros(1, 0);

    for b = 1:size(raw_blocks, 1)
        s = raw_blocks(b, 1);
        e = raw_blocks(b, 2);
        if (e - s + 1) < window_samples
            continue;
        end

        seg = signal(s:e);
        [~, f, t_block, P] = spectrogram(seg, win, noverlap_samples, [], fs, 'psd');

        p_scale = band_power_from_psd(P, f, cfg.quality.scale_band) - band_power_from_psd(P, f, cfg.quality.line_exclusion_band);
        p_scale = max(p_scale, 0);
        p_line = band_power_from_psd(P, f, line_band);
        p_lo = band_power_from_psd(P, f, lower_neighbor);
        p_hi = band_power_from_psd(P, f, upper_neighbor);
        ratio_db = 10 * log10(p_line ./ median([p_lo; p_hi], 1));

        all_p_scale = [all_p_scale, p_scale]; %#ok<AGROW>
        all_ratio_db = [all_ratio_db, ratio_db]; %#ok<AGROW>
        all_t = [all_t, data.t_rel(s) + t_block(:)']; %#ok<AGROW>
    end

    if isempty(all_p_scale)
        sigma_band_uV = NaN;
        line_ratio_db = NaN;
        line_ratio_p95_db = NaN;
        line_ratio_max_db = NaN;
        pct_time_line_high = NaN;
    else
        sigma_band_uV = sqrt(median(all_p_scale));
        line_ratio_db = median(all_ratio_db);
        line_ratio_p95_db = simple_prctile(all_ratio_db, 95);
        line_ratio_max_db = max(all_ratio_db);
        pct_time_line_high = 100 * mean(all_ratio_db > cfg.quality.line_ratio_thr_db);
    end

    valid_signal = signal(valid_mask);
    mad_uV = 1.4826 * mad(valid_signal, 1);

    [quantization_step_uV, quantization_step_source] = resolve_quantization_step(data.meta, valid_signal);
    snr_quantization_db = 20 * log10(sigma_band_uV / quantization_step_uV);
    adc_codes_span = (simple_prctile(valid_signal, 75) - simple_prctile(valid_signal, 25)) / quantization_step_uV;

    n_distinct_values = numel(unique(valid_signal));
    pct_nan = 100 * (1 - sum(valid_mask) / numel(valid_mask));
    pct_clipped = compute_pct_clipped(valid_signal, data.meta, quantization_step_uV);
    flat_fraction = compute_flat_fraction(signal, raw_blocks, fs);

    quality_class = classify_quality(snr_quantization_db, line_ratio_db, pct_time_line_high, ...
        flat_fraction, pct_clipped, cfg);

    q = struct();
    q.sigma_band_uV = sigma_band_uV;
    q.mad_uV = mad_uV;
    q.line_ratio_db = line_ratio_db;
    q.line_ratio_p95_db = line_ratio_p95_db;
    q.line_ratio_max_db = line_ratio_max_db;
    q.pct_time_line_high = pct_time_line_high;
    q.line_ratio_series = struct('t_s', all_t(:), 'ratio_db', all_ratio_db(:));
    q.quantization_step_uV = quantization_step_uV;
    q.quantization_step_source = quantization_step_source;
    q.snr_quantization_db = snr_quantization_db;
    q.adc_codes_span = adc_codes_span;
    q.n_distinct_values = n_distinct_values;
    q.pct_clipped = pct_clipped;
    q.pct_nan = pct_nan;
    q.flat_fraction = flat_fraction;
    q.quality_class = quality_class;
end

%% ======================================================================
function bp = band_power_from_psd(P, f, band)
    mask = f >= band(1) & f <= band(2);
    bp = trapz(f(mask), P(mask, :), 1);
end

function [step_uV, source] = resolve_quantization_step(meta, valid_signal)
    if isfield(meta, 'quantization_step_uV')
        v = str2double(meta.quantization_step_uV);
        if ~isnan(v) && v > 0
            step_uV = v;
            source = 'header';
            return;
        end
    end

    u = unique(valid_signal);
    if numel(u) < 2
        step_uV = NaN;
    else
        d = diff(sort(u));
        d = d(d > 0);
        if isempty(d)
            step_uV = NaN;
        else
            step_uV = min(d);
        end
    end
    source = 'estimated';
end

function pct_clipped = compute_pct_clipped(valid_signal, meta, quantization_step_uV)
    if ~isfield(meta, 'physical_min_uV') || ~isfield(meta, 'physical_max_uV')
        pct_clipped = NaN;
        return;
    end
    phys_min = str2double(meta.physical_min_uV);
    phys_max = str2double(meta.physical_max_uV);
    if isnan(phys_min) || isnan(phys_max) || isnan(quantization_step_uV)
        pct_clipped = NaN;
        return;
    end
    tol = 0.5 * quantization_step_uV;
    is_clipped = abs(valid_signal - phys_min) <= tol | abs(valid_signal - phys_max) <= tol;
    pct_clipped = 100 * sum(is_clipped) / numel(valid_signal);
end

function flat_fraction = compute_flat_fraction(signal, raw_blocks, fs)
    min_run = round(0.5 * fs);
    n_flat = 0;
    n_valid = 0;
    for b = 1:size(raw_blocks, 1)
        s = raw_blocks(b, 1);
        e = raw_blocks(b, 2);
        seg = signal(s:e);
        n_valid = n_valid + numel(seg);
        if numel(seg) < min_run
            continue;
        end
        is_same_as_prev = [false; diff(seg) == 0];
        run_id = cumsum(~is_same_as_prev);
        run_lengths = accumarray(run_id, 1);
        n_flat = n_flat + sum(run_lengths(run_lengths >= min_run));
    end
    if n_valid == 0
        flat_fraction = NaN;
    else
        flat_fraction = n_flat / n_valid;
    end
end

function quality_class = classify_quality(snr_quantization_db, line_ratio_db, pct_time_line_high, ...
        flat_fraction, pct_clipped, cfg)

    is_unusable = isnan(snr_quantization_db) || snr_quantization_db < cfg.quality.min_snr_quant_db || ...
        (~isnan(flat_fraction) && flat_fraction > 0.5) || (~isnan(pct_clipped) && pct_clipped > 5);
    if is_unusable
        quality_class = 'unusable';
        return;
    end

    is_low = snr_quantization_db < cfg.quality.low_amp_snr_db;
    is_line = (line_ratio_db > cfg.quality.line_ratio_thr_db) || (pct_time_line_high > cfg.quality.line_intermittent_pct);

    if is_low && is_line
        quality_class = 'low_amplitude+line';
    elseif is_low
        quality_class = 'low_amplitude';
    elseif is_line
        quality_class = 'line_contaminated';
    else
        quality_class = 'ok';
    end
end
