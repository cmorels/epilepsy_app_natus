function [T, cfg_cmp] = cs_iid_repro_test(cfg, ch_tbl, provenance)
% CS_IID_REPRO_TEST  Section 5.3 reproducibility test, before any
% re-detection: on cfg.iid_repro_channels, cs_iid_run_channel.m fed with
% the seizure detector's ORIGINAL list must reproduce that campaign's
% 04_iid exactly.
%
% Original list = what run_pipeline_edf.m's seizures_for_iid_of passed to
% detect_iid in the campaign (bilateral branch with
% cfg.bilateral.exclude_rescued_from_iid = true, asserted from the mouse's
% config_used.mat): every reconciled row (accepted, rescued, imputed) of
% the channel's 03_seizures/<EDF>/<base>_seizures.mat (no file = no rows).
%
% Two variants per channel, one load of the *_clean.txt:
%   A  list through detect_iid's seizures argument (as the pipeline did)
%   B  the same list as manual zones with the buffer added by the caller
%      (the path the re-detection uses for every non-seizure zone)
% Compared with == (no tolerance): baseline, lower/upper threshold,
% analyzed/total duration, n peaks/complexes/single/polyspike/bursts,
% the start/end/n_spikes/max amplitude of every complex, and the start/end
% of every exclusion zone.
%
% Results are cached per channel in cfg.out_dir/cache (key: channel + SHA-1
% of detect_iid.m and of the wrapper/test code); outputs of the runs go to
% cfg.out_dir/iid_repro_test. Also reported: the real P97 of |signal| over
% the included samples (what the 100:5:130 uV baseline grid estimates).
%
% T       : one row per channel x variant
% cfg_cmp : cfg.iid of each campaign's config_used.mat vs pipeline_config()

    pc = pipeline_config();
    sc = provenance.stats_code;
    sha = [provenance.iid_code{strcmp(provenance.iid_code(:, 1), 'detect_iid.m'), 2}, ...
        sc{strcmp(sc(:, 1), 'cs_iid_run_channel.m'), 2}, sc{strcmp(sc(:, 1), 'cs_iid_repro_test.m'), 2}, ...
        sc{strcmp(sc(:, 1), 'cs_iid_p97.m'), 2}];
    buf = pc.iid.exclusion_buffer_s;
    if buf ~= cfg.iid_excl_buffer_s
        error('cs_iid_repro_test:Buffer', 'pipeline_config exclusion_buffer_s (%g) ~= cfg.iid_excl_buffer_s (%g)', buf, cfg.iid_excl_buffer_s);
    end
    rows = {};
    cfg_cmp = {};
    for k = 1:size(cfg.iid_repro_channels, 1)
        c = cfg.iid_repro_channels(k, :);
        j = find(strcmp(ch_tbl.campaign, c{1}) & strcmp(ch_tbl.mouse, c{2}) & startsWith(ch_tbl.edf, c{3}) & strcmp(ch_tbl.region, c{4}));
        if numel(j) ~= 1
            error('cs_iid_repro_test:Channel', '%s %s %s %s: %d channels match', c{:}, numel(j));
        end
        ch = ch_tbl(j, :);
        root = cfg.campaigns(strcmp({cfg.campaigns.name}, c{1})).root;
        [~, nm] = fileparts(ch.file{1});
        base = strrep(nm, '_clean', '');
        key = sprintf('%s_%s_%s_%s', c{1}, c{2}, c{3}, c{4});
        cache = fullfile(cfg.out_dir, 'cache', ['repro_' key '.mat']);

        used = load(fullfile(root, c{2}, 'config_used.mat'), 'cfg');
        used = used.cfg;
        cfg_cmp(end+1, :) = {c{1}, c{2}, compare_iid_cfg(used.iid, pc.iid)}; %#ok<AGROW>

        if isfile(cache)
            q = load(cache, 'res');
            if strcmp(q.res.sha, sha)
                rows = [rows; q.res.rows]; %#ok<AGROW>
                continue;
            end
        end

        if strcmp(used.bilateral.rescue_mode, 'off') || ~used.bilateral.exclude_rescued_from_iid
            error('cs_iid_repro_test:Branch', '%s %s: campaign not run with bilateral reconciliation + exclude_rescued_from_iid', c{1}, c{2});
        end
        sz_file = fullfile(root, c{2}, '03_seizures', ch.edf{1}, [base '_seizures.mat']);
        if isfile(sz_file)
            q = load(sz_file, 'seizure_results');
            X = q.seizure_results.seizures;
            sz = table(X.seizure_id, X.start_s, X.end_s, 'VariableNames', {'id', 'start_s', 'end_s'});
            statuses = strjoin(cellfun(@(s) sprintf('%s=%d', s, sum(strcmp(X.detection_status, s))), ...
                unique(X.detection_status)', 'UniformOutput', false), ', ');
        else
            sz = table(zeros(0, 1), zeros(0, 1), zeros(0, 1), 'VariableNames', {'id', 'start_s', 'end_s'});
            statuses = 'sin *_seizures.mat (0 filas)';
        end
        o = load(fullfile(root, c{2}, '04_iid', ch.edf{1}, [base '_IID_results.mat']), 'iid_results');
        o = o.iid_results;

        tl = tic;
        data = load_lfp_txt(ch.file{1});
        t_load = toc(tl);
        opts = struct('cfg', pc, 'figures', false, 'provenance', provenance);
        res_rows = {};
        for v = 'AB'
            out = fullfile(cfg.out_dir, 'iid_repro_test', v, c{2}, ch.edf{1});
            if v == 'A'
                n = cs_iid_run_channel(data, sz, zeros(0, 2), out, opts);
            else
                n = cs_iid_run_channel(data, sz([], :), [sz.start_s - buf, sz.end_s + buf], out, opts);
            end
            res_rows(end+1, :) = compare_runs(c, v, height(sz), statuses, o, n, t_load); %#ok<AGROW>
        end
        clear data
        res = struct('sha', sha, 'rows', {res_rows});
        save(cache, 'res');
        rows = [rows; res_rows]; %#ok<AGROW>
    end
    T = cell2table(rows, 'VariableNames', {'campaign', 'mouse', 'edf', 'region', 'variante', 'n_crisis_detector', ...
        'filas_por_estado', 'baseline_uV_orig', 'baseline_uV_wrapper', 'lower_uV_orig', 'lower_uV_wrapper', ...
        'n_complejos_orig', 'n_complejos_wrapper', 'n_polyspike_orig', 'n_polyspike_wrapper', 'n_bursts_orig', ...
        'n_bursts_wrapper', 'analyzed_min_orig', 'analyzed_min_wrapper', 'n_zonas_orig', 'n_zonas_wrapper', ...
        'complejos_identicos', 'zonas_identicas', 'campos_distintos', 'reproduce', 'p97_real_uV', ...
        'frac_bajo_130uV', 'carga_s', 'deteccion_s'});
    cfg_cmp = cell2table(cfg_cmp, 'VariableNames', {'campaign', 'mouse', 'diferencias_cfg_iid'});
end

function row = compare_runs(c, v, n_sz, statuses, o, n, t_load)
    so = o.summary;
    sn = n.summary;
    f = {'baseline_uV', 'lower_threshold_uV', 'upper_threshold_uV', 'total_duration_min', 'analyzed_duration_min', ...
        'excluded_duration_min', 'n_exclusion_zones', 'total_peaks', 'total_complexes', 'n_single', 'n_polyspike', ...
        'n_bursts', 'complexes_per_min', 'single_per_min', 'polyspikes_per_min', 'bursts_per_hour'};
    bad = f(~cellfun(@(x) isequaln(so.(x), sn.(x)), f));
    cols = {'Start_s', 'End_s', 'N_spikes', 'Max_amplitude_mV', 'Is_polyspike'};
    cx_same = isequal(o.spike_complex_table{:, cols}, n.spike_complex_table{:, cols});
    zo = o.exclusion_table;
    zn = n.exclusion_table;
    z_same = isequal(zo{:, {'start_s', 'end_s'}}, zn{:, {'start_s', 'end_s'}});
    if v == 'A'
        z_same = z_same && isequal(zo.reason, zn.reason);
    end
    ok = isempty(bad) && cx_same && z_same;
    row = {c{1}, c{2}, c{3}, c{4}, v, n_sz, statuses, so.baseline_uV, sn.baseline_uV, so.lower_threshold_uV, ...
        sn.lower_threshold_uV, so.total_complexes, sn.total_complexes, so.n_polyspike, sn.n_polyspike, so.n_bursts, ...
        sn.n_bursts, so.analyzed_duration_min, sn.analyzed_duration_min, so.n_exclusion_zones, sn.n_exclusion_zones, ...
        cx_same, z_same, strjoin(bad, ','), ok, n.p97_uV, n.frac_below_grid_top, t_load, n.detect_s};
end

function s = compare_iid_cfg(a, b)
% Differences between two cfg.iid structs, ignoring output_dir (set per run).
    fa = setdiff(fieldnames(a), {'output_dir'});
    fb = setdiff(fieldnames(b), {'output_dir'});
    d = [cellfun(@(x) ['solo en campaign: ' x], setdiff(fa, fb), 'UniformOutput', false); ...
         cellfun(@(x) ['solo en actual: ' x], setdiff(fb, fa), 'UniformOutput', false)];
    d = d(:);
    for x = intersect(fa, fb)'
        if ~isequaln(a.(x{1}), b.(x{1}))
            d{end+1, 1} = sprintf('%s: campaign=%s actual=%s', x{1}, val(a.(x{1})), val(b.(x{1}))); %#ok<AGROW>
        end
    end
    if isempty(d)
        s = 'identico';
    else
        s = strjoin(d', '; ');
    end
end

function s = val(x)
    if isstruct(x)
        s = jsonencode(x);
    elseif ischar(x)
        s = x;
    else
        s = mat2str(x);
    end
end
