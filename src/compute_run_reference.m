function ref = compute_run_reference(jobs, cfg, out_dir)
% COMPUTE_RUN_REFERENCE  Pre-pass of a campaign run: imports every EDF of
% the run once and measures the reference amplitude for the automatic
% gain from the NON-attenuated channels (Excel Attenuation column).
%
%   ref = compute_run_reference(jobs, cfg, out_dir)
%
%   jobs    : struct array, one per EDF: .subject_id (exactly as in the
%             Excel, e.g. '005-s'), .edf_path, .output_root (that animal's
%             cfg.paths.output_root; the txt go to <output_root>/01_txt and
%             are reused by run_pipeline_edf.m through cfg.edf.reuse_import)
%   cfg     : campaign_config() (+ cfg.edf.channels.log_file)
%   out_dir : where reference_channels.csv / reference_per_region.csv are
%             written ('' = not written)
%
% Reference per region = median sigma_band_uV (signal_quality.m) of that
% region's non-attenuated channels in THIS run (any number >= 1).
% ref.fallback = median over every non-attenuated channel of the run, for
% a region with none; NaN if the run has none at all (precondition_lfp.m
% then applies no gain, with a warning). An EDF or channel that fails is
% recorded with its error and skipped -- this function never stops.
%
% OUTPUT (struct ref): .map (containers.Map region -> reference uV, only
% regions with a value), .fallback (uV or NaN), .channels (table: one row
% per channel), .per_region (table: region, n_channels, reference_uV).

    rows = struct('subject_id', {}, 'source_file', {}, 'region', {}, 'excel_attenuation', {}, ...
        'attenuated', {}, 'sigma_band_uV', {}, 'status', {});

    for j = 1:numel(jobs)
        job = jobs(j);
        [~, stem, ext] = fileparts(job.edf_path);
        fprintf('[reference %d/%d] %s | %s\n', j, numel(jobs), job.subject_id, [stem ext]);
        jcfg = cfg;
        jcfg.edf.subject_id = job.subject_id;
        jcfg.edf.output_dir = fullfile(job.output_root, '01_txt');
        try
            manifest = edf_import(job.edf_path, jcfg);
        catch ME
            rows(end+1) = make_row(job.subject_id, [stem ext], '', '', NaN, NaN, ['edf_import failed: ' ME.message]); %#ok<AGROW>
            fprintf(2, '  edf_import failed: %s\n', ME.message);
            continue;
        end
        for c = 1:height(manifest)
            att_text = '';
            if ismember('attenuation', manifest.Properties.VariableNames)
                att_text = char(manifest.attenuation{c});
            end
            is_att = excel_attenuation(att_text, cfg);
            region = manifest.region{c};
            if is_att
                rows(end+1) = make_row(job.subject_id, manifest.source_file{c}, region, att_text, 1, NaN, 'attenuated: not used'); %#ok<AGROW>
                continue;
            end
            try
                data = load_lfp_txt(manifest.txt_file{c});
                q = signal_quality(data, cfg);
                status = 'used';
                if isnan(q.sigma_band_uV) || q.sigma_band_uV <= 0
                    status = 'no measurable amplitude: not used';
                end
                rows(end+1) = make_row(job.subject_id, manifest.source_file{c}, region, att_text, 0, q.sigma_band_uV, status); %#ok<AGROW>
                clear data
            catch ME
                rows(end+1) = make_row(job.subject_id, manifest.source_file{c}, region, att_text, 0, NaN, ['measure failed: ' ME.message]); %#ok<AGROW>
                fprintf(2, '  [%s] measure failed: %s\n', region, ME.message);
            end
        end
    end

    names = {'subject_id', 'source_file', 'region', 'excel_attenuation', 'attenuated', 'sigma_band_uV', 'status'};
    if isempty(rows)
        channels = table(cell(0, 1), cell(0, 1), cell(0, 1), cell(0, 1), zeros(0, 1), zeros(0, 1), cell(0, 1), ...
            'VariableNames', names);
    else
        % attenuated: 1 / 0, NaN when the EDF could not be imported
        channels = table({rows.subject_id}', {rows.source_file}', {rows.region}', {rows.excel_attenuation}', ...
            [rows.attenuated]', [rows.sigma_band_uV]', {rows.status}', 'VariableNames', names);
    end

    used = strcmp(channels.status, 'used');
    ref_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
    regions = unique(channels.region(used));
    per_region = table(cell(0, 1), zeros(0, 1), zeros(0, 1), 'VariableNames', {'region', 'n_channels', 'reference_uV'});
    for r = 1:numel(regions)
        vals = channels.sigma_band_uV(used & strcmp(channels.region, regions{r}));
        ref_map(regions{r}) = median(vals);
        per_region(end+1, :) = {regions{r}, numel(vals), median(vals)}; %#ok<AGROW>
    end
    if any(used)
        fallback = median(channels.sigma_band_uV(used));
    else
        fallback = NaN;
    end
    per_region(end+1, :) = {'(all regions, fallback)', nnz(used), fallback};

    fprintf('\nReference amplitude (sigma_band_uV of non-attenuated channels):\n');
    for r = 1:height(per_region)
        fprintf('  %-26s n=%-3d %.3f uV\n', per_region.region{r}, per_region.n_channels(r), per_region.reference_uV(r));
    end
    if ~any(used)
        fprintf(2, '  No non-attenuated channel in this run: attenuated channels will get NO gain.\n');
    end

    if nargin >= 3 && ~isempty(out_dir)
        try
            if ~isfolder(out_dir)
                mkdir(out_dir);
            end
            writetable(channels, fullfile(out_dir, 'reference_channels.csv'));
            writetable(per_region, fullfile(out_dir, 'reference_per_region.csv'));
        catch ME
            warning('compute_run_reference:WriteFailed', 'Could not write the reference CSVs: %s', ME.message);
        end
    end

    ref = struct('map', ref_map, 'fallback', fallback, 'channels', channels, 'per_region', per_region);
end

function r = make_row(subject_id, source_file, region, att_text, attenuated, sigma, status)
    r = struct('subject_id', subject_id, 'source_file', source_file, 'region', region, ...
        'excel_attenuation', att_text, 'attenuated', attenuated, 'sigma_band_uV', sigma, 'status', status);
end
