function T = make_cases_template(input, out_csv, cfg)
% MAKE_CASES_TEMPLATE  Generate a pre-filled cases CSV by running
% signal_quality.m over every channel found, across an entire cohort.
%
%   T = make_cases_template(input, out_csv, cfg)
%
%   input   : a folder; a folder searched recursively if
%             cfg.quality.recursive_scan is true; a cell array of
%             folders; or a cell array of specific .txt paths. Writes ONE
%             csv spanning every animal scanned (subject_id is a column),
%             so one cases file can serve a whole cohort even though EDF
%             exports are organized one folder per animal.
%   out_csv : path to write the CSV to
%   cfg     : struct from pipeline_config.m
%
% CRITICAL DESIGN DECISION: the 'case' column is filled with 'normal' on
% EVERY row, always. The algorithm's recommendation goes in a SEPARATE
% 'suggested_case' column that the pipeline never reads. A freshly
% generated CSV is therefore, as-is, equivalent to today's behavior; the
% user promotes only the rows they decide to, row by row, comparing the
% suggestion against what they actually observed reviewing in Natus. That
% turns agreement between suggested_sensitivity_uV_per_mm and the
% sensitivity actually used in Natus into an independent check on the
% calibration.
%
% Reference resolution for the suggestions, per channel, in this order:
%   1. cfg.quality.reference_sigma_uV (containers.Map by region, or a
%      scalar) if it covers this channel's region -> reference_source='config'
%   2. else the MEDIAN sigma_band_uV of same-region channels in THIS scan
%      -> reference_source='cohort_median'. Self-calibrating: most animals
%      in a cohort are fine, so the median is robust and the problematic
%      ones fall below it. Warns if that region has fewer than
%      cfg.quality.min_reference_channels channels in the scan.
%   3. if a region has exactly one channel in the whole scan, a
%      "cohort median" of one channel is just that channel itself, which
%      is meaningless as a reference -- suggested_gain is left empty and
%      suggested_case falls back to the ABSOLUTE snr_quantization_db
%      criterion alone.
% Scanning the full cohort in one call (see cfg.quality.recursive_scan)
% is what makes this self-calibration reliable; scanning a single animal
% would make its own median trivially itself.

    files = resolve_txt_input_list(input, cfg.quality.recursive_scan);
    if isempty(files)
        error('make_cases_template:NoFiles', 'No .txt files found for the given input.');
    end

    rows = struct('subject_id', {}, 'source_file', {}, 'region', {}, 'sigma_band_uV', {}, ...
        'line_ratio_db', {}, 'line_ratio_p95_db', {}, 'pct_time_line_high', {}, ...
        'snr_quantization_db', {}, 'adc_codes_span', {}, 'quality_class', {});

    for i = 1:numel(files)
        try
            data = load_lfp_txt(files{i});
            q = signal_quality(data, cfg);
            [~, name, ext] = fileparts(files{i});
            rows(end+1) = struct( ... %#ok<AGROW>
                'subject_id', field_or(data.meta, 'subject_id', ''), ...
                'source_file', [name ext], ...
                'region', field_or(data.meta, 'region', ''), ...
                'sigma_band_uV', q.sigma_band_uV, ...
                'line_ratio_db', q.line_ratio_db, ...
                'line_ratio_p95_db', q.line_ratio_p95_db, ...
                'pct_time_line_high', q.pct_time_line_high, ...
                'snr_quantization_db', q.snr_quantization_db, ...
                'adc_codes_span', q.adc_codes_span, ...
                'quality_class', q.quality_class);
        catch ME
            warning('make_cases_template:FileFailed', 'Skipping %s: %s', files{i}, ME.message);
        end
    end

    if isempty(rows)
        error('make_cases_template:NoUsableFiles', 'No file could be read successfully.');
    end

    regions = unique({rows.region});
    region_n = containers.Map('KeyType', 'char', 'ValueType', 'double');
    region_median = containers.Map('KeyType', 'char', 'ValueType', 'double');
    for r = 1:numel(regions)
        region = regions{r};
        vals = [rows(strcmp({rows.region}, region)).sigma_band_uV];
        vals = vals(~isnan(vals));
        region_n(region) = numel(vals);
        if ~isempty(vals)
            region_median(region) = median(vals);
        else
            region_median(region) = NaN;
        end
        if numel(vals) > 0 && numel(vals) < cfg.quality.min_reference_channels
            warning('make_cases_template:FewChannels', ...
                'Region "%s" has only %d channel(s) in this scan (< cfg.quality.min_reference_channels = %d); its cohort-median reference may be unreliable.', ...
                region, numel(vals), cfg.quality.min_reference_channels);
        end
    end

    n = numel(rows);
    case_col = repmat({'normal'}, n, 1);
    suggested_case = cell(n, 1);
    reference_used = nan(n, 1);
    reference_source = cell(n, 1);
    suggested_gain = nan(n, 1);
    suggested_sensitivity = nan(n, 1);

    for i = 1:n
        region = rows(i).region;
        [ref_val, ref_source, usable] = resolve_reference(region, region_n, region_median, cfg.quality.reference_sigma_uV);
        reference_used(i) = ref_val;
        reference_source{i} = ref_source;

        if usable
            suggested_gain(i) = ref_val / rows(i).sigma_band_uV;
            suggested_sensitivity(i) = 100 / suggested_gain(i);
            is_attenuated = (rows(i).sigma_band_uV < ref_val / cfg.quality.attenuation_ratio_thr) || ...
                (rows(i).snr_quantization_db < cfg.quality.low_amp_snr_db);
        else
            is_attenuated = rows(i).snr_quantization_db < cfg.quality.low_amp_snr_db;
        end

        is_line = (rows(i).line_ratio_db > cfg.quality.line_ratio_thr_db) || ...
            (rows(i).pct_time_line_high > cfg.quality.line_intermittent_pct);

        if strcmp(rows(i).quality_class, 'unusable')
            suggested_case{i} = 'unusable';
        elseif is_attenuated && is_line
            suggested_case{i} = 'both';
        elseif is_attenuated
            suggested_case{i} = 'attenuated';
        elseif is_line
            suggested_case{i} = 'line';
        else
            suggested_case{i} = 'normal';
        end
    end

    T = table({rows.subject_id}', {rows.source_file}', {rows.region}', case_col, suggested_case, ...
        [rows.sigma_band_uV]', reference_used, reference_source, [rows.line_ratio_db]', [rows.line_ratio_p95_db]', ...
        [rows.pct_time_line_high]', [rows.snr_quantization_db]', [rows.adc_codes_span]', {rows.quality_class}', ...
        suggested_gain, suggested_sensitivity, nan(n, 1), repmat({''}, n, 1), repmat({''}, n, 1), ...
        'VariableNames', {'subject_id', 'source_file', 'region', 'case', 'suggested_case', ...
        'sigma_band_uV', 'reference_used', 'reference_source', 'line_ratio_db', 'line_ratio_p95_db', ...
        'pct_time_line_high', 'snr_quantization_db', 'adc_codes_span', 'quality_class', ...
        'suggested_gain', 'suggested_sensitivity_uV_per_mm', 'gain', 'notch', 'notes'});

    writetable(T, out_csv);

    n_suggested_nonnormal = sum(~strcmp(suggested_case, 'normal'));
    fprintf('make_cases_template: %d channel(s) scanned, %d region(s), %d row(s) with a non-normal suggested_case.\n', ...
        n, numel(regions), n_suggested_nonnormal);
    fprintf('  Written to: %s\n', out_csv);
    fprintf('  Every row''s "case" column is ''normal'' -- edit it by hand before using cfg.cases.file.\n');
end

function v = field_or(s, name, default)
    if isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end

function [ref_val, ref_source, usable] = resolve_reference(region, region_n, region_median, cfg_ref)
    if isa(cfg_ref, 'containers.Map')
        if isKey(cfg_ref, region)
            ref_val = cfg_ref(region);
            ref_source = 'config';
            usable = true;
            return;
        end
    elseif isnumeric(cfg_ref) && isscalar(cfg_ref) && ~isnan(cfg_ref)
        ref_val = cfg_ref;
        ref_source = 'config';
        usable = true;
        return;
    end

    ref_source = 'cohort_median';
    if isKey(region_n, region) && region_n(region) >= 2
        ref_val = region_median(region);
        usable = true;
    else
        ref_val = NaN;
        usable = false;
    end
end
