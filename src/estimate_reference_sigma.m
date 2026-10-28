function ref = estimate_reference_sigma(input, cfg)
% ESTIMATE_REFERENCE_SIGMA  Calibrates the amplitude reference used by
% precondition_lfp.m / make_cases_template.m from the recordings the user
% considers GOOD, instead of hardcoding a number.
%
%   ref = estimate_reference_sigma(input, cfg)
%
%   input : a folder; a folder searched recursively if
%           cfg.quality.recursive_scan is true; a cell array of folders;
%           or a cell array of specific .txt paths. EDF exports are
%           organized one folder per animal, so scanning several folders
%           (or one recursively) at once is the normal case.
%   cfg   : struct from pipeline_config.m (uses cfg.quality.*)
%
% Groups by REGION across every animal scanned -- "what does a healthy
% HPCleft look like in this montage", not "what does this one animal look
% like", since different regions can have genuinely different baseline
% amplitude. Never writes config; only prints ready-to-paste lines and
% returns the numbers.
%
% OUTPUT (struct ref): .per_region (table: region, n, median, iqr, p25,
% p75, min, max), .per_channel (table: subject_id, region, source_file,
% sigma_band_uV, snr_quantization_db, deviation_multiple -- sorted by how
% far each channel deviates from its region's median, most first).

    files = resolve_txt_input_list(input, cfg.quality.recursive_scan);
    if isempty(files)
        error('estimate_reference_sigma:NoFiles', 'No .txt files found for the given input.');
    end

    rows = cell(0, 5);  % subject_id, region, source_file, sigma_band_uV, snr_quantization_db
    for i = 1:numel(files)
        try
            data = load_lfp_txt(files{i});
            q = signal_quality(data, cfg);
            subject_id = field_or(data.meta, 'subject_id', '');
            region = field_or(data.meta, 'region', '');
            rows(end+1, :) = {subject_id, region, data.file, q.sigma_band_uV, q.snr_quantization_db}; %#ok<AGROW>
        catch ME
            warning('estimate_reference_sigma:FileFailed', 'Skipping %s: %s', files{i}, ME.message);
        end
    end

    if isempty(rows)
        error('estimate_reference_sigma:NoUsableFiles', 'No file could be read successfully.');
    end

    per_channel = cell2table(rows, 'VariableNames', {'subject_id', 'region', 'source_file', 'sigma_band_uV', 'snr_quantization_db'});

    regions = unique(per_channel.region, 'stable');
    region_rows = cell(0, 8);
    per_channel.deviation_multiple = nan(height(per_channel), 1);

    fprintf('\n=== estimate_reference_sigma: %d channel(s), %d region(s) ===\n', height(per_channel), numel(regions));

    for r = 1:numel(regions)
        region = regions{r};
        mask = strcmp(per_channel.region, region);
        vals = per_channel.sigma_band_uV(mask);
        vals = vals(~isnan(vals));
        n = numel(vals);

        if n == 0
            warning('estimate_reference_sigma:NoValidChannels', 'Region "%s" has no channels with a usable sigma_band_uV.', region);
            continue;
        end
        if n < cfg.quality.min_reference_channels
            warning('estimate_reference_sigma:FewChannels', ...
                'Region "%s" has only %d channel(s) (< cfg.quality.min_reference_channels = %d); its reference may be unreliable.', ...
                region, n, cfg.quality.min_reference_channels);
        end

        med = median(vals);
        p25 = simple_prctile(vals, 25);
        p75 = simple_prctile(vals, 75);
        region_rows(end+1, :) = {region, n, med, p75 - p25, p25, p75, min(vals), max(vals)}; %#ok<AGROW>
        per_channel.deviation_multiple(mask) = per_channel.sigma_band_uV(mask) / med;

        fprintf('  %-12s n=%-3d  median=%7.2f uV  IQR=%6.2f  [%6.2f, %6.2f]\n', region, n, med, p75 - p25, min(vals), max(vals));

        fig = figure('Name', sprintf('sigma_band_uV -- %s', region), 'Visible', 'on');
        histogram(vals, 'FaceColor', [0.3 0.5 0.7]); %#ok<NASGU>
        xline(med, 'r--', 'median', 'LineWidth', 1.5);
        title(sprintf('%s: sigma\\_band\\_uV (n=%d)', region, n), 'Interpreter', 'tex');
        xlabel('sigma\_band\_uV', 'Interpreter', 'tex'); ylabel('channels'); grid on;
    end

    per_region = cell2table(region_rows, 'VariableNames', {'region', 'n', 'median', 'iqr', 'p25', 'p75', 'min', 'max'});

    sort_key = abs(log2(per_channel.deviation_multiple));
    [~, order] = sort(sort_key, 'descend', 'MissingPlacement', 'last');
    per_channel = per_channel(order, :);

    fprintf('\n--- channels sorted by deviation from their region''s median (most first) ---\n');
    disp(per_channel);

    fprintf('\n--- config lines ready to paste ---\n');
    all_vals = per_channel.sigma_band_uV(~isnan(per_channel.sigma_band_uV));
    fprintf('%% scalar (overall median across all regions pooled together -- use only if you have one region, or as a rough fallback):\n');
    fprintf('cfg.quality.reference_sigma_uV = %.2f;\n', median(all_vals));
    fprintf('%% per-region (recommended):\n');
    fprintf('cfg.quality.reference_sigma_uV = containers.Map(%s, %s);\n', ...
        cellstr_literal(per_region.region), num_literal(per_region.median));

    ref = struct();
    ref.per_region = per_region;
    ref.per_channel = per_channel;
end

function v = field_or(s, name, default)
    if isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end

function s = cellstr_literal(c)
    parts = cellfun(@(x) ['''' x ''''], c, 'UniformOutput', false);
    s = ['{' strjoin(parts, ', ') '}'];
end

function s = num_literal(v)
    parts = arrayfun(@(x) sprintf('%.2f', x), v, 'UniformOutput', false);
    s = ['{' strjoin(parts, ', ') '}'];
end
