function [T, V, F] = cs_iid_sens_all(cfg, Z, C, worker, n_workers, only)
% CS_IID_SENS_ALL  Phase C of v2 over every channel of the IID plan:
% cs_iid_sens_channel.m per channel (one read + filter), cached in
% cfg.out_dir/cache/sens/ (reused when zones, sensitivity settings and code
% are the same). Worker mode as cs_iid_rerun_all.m.
%
%   [T, V, F] = cs_iid_sens_all(cfg, Z, C, 0)         run what is missing, collect
%   cs_iid_sens_all(cfg, Z, C, k, n)                  worker k of n, no collection
%   [...] = cs_iid_sens_all(cfg, Z, C, 0, 1, idx)     only plan channels idx (test)
%
% T: one row per channel (percentiles, baseline, whether the 'actual'
%    variant reproduces the main run C exactly); V: channel x variant
%    counts; F: flat segments (t_rel s + absolute POSIX).

    if nargin < 4, worker = 0; n_workers = 1; end
    P = Z.channels(Z.channels.run, :);
    if nargin >= 6 && ~isempty(only), P = P(only, :); end
    Cm = C(strcmp(C.variant, 'principal'), :);
    keyC = strcat(Cm.campaign, '|', Cm.mouse, '|', Cm.edf, '|', Cm.region);
    cache_dir = fullfile(cfg.out_dir, 'cache', 'sens');
    if ~isfolder(cache_dir), mkdir(cache_dir); end
    code = sens_code_key(cfg);
    rows = cell(height(P), 1); vrows = cell(height(P), 1); frows = cell(height(P), 1);
    for i = 1:height(P)
        if worker ~= 0 && mod(i - 1, n_workers) + 1 ~= worker
            continue;
        end
        ch = P(i, :);
        [sz, man] = cs_iid_channel_zones(Z.zones, ch, cfg.iid_excl_buffer_s);
        f = fullfile(cache_dir, regexprep(sprintf('%s_%s_%s_%s', ch.campaign{1}, ch.mouse{1}, ch.edf{1}, ch.region{1}), '[^A-Za-z0-9_-]', '_'));
        f = [f '.mat'];
        r = [];
        if isfile(f)
            q = load(f, 'r', 'code');
            if strcmp(q.code, code) && isequal(sortrows(reshape(q.r.sz{:, {'start_s', 'end_s'}}, [], 2)), ...
                    sortrows(reshape(sz{:, {'start_s', 'end_s'}}, [], 2))) && isequal(sortrows(q.r.man), sortrows(reshape(man, [], 2)))
                r = q.r;
            end
        end
        if isempty(r)
            r = cs_iid_sens_channel(cfg, ch, sz, man);
            save(f, 'r', 'code');
            fprintf('[%s] sens %d/%d %s %s %s %s: %.0f s (carga %.0f s), actual %d complejos, P99 |x_bp| %.0f uV, 312.5 uV = P%.2f, %d tramos planos\n', ...
                char(datetime('now', 'Format', 'HH:mm:ss')), i, height(P), ch.campaign{1}, ch.mouse{1}, ch.edf{1}(1:min(10, end)), ...
                ch.region{1}, r.total_s, r.load_s, r.variants.n_complexes(1), r.pct_bp_uV(1), r.pctile_of_312_5, size(r.flat, 1));
        end
        if worker ~= 0, continue; end
        j = find(strcmp(keyC, sprintf('%s|%s|%s|%s', ch.campaign{1}, ch.mouse{1}, ch.edf{1}, ch.region{1})));
        a = r.variants(1, :);
        repro = numel(j) == 1 && a.n_complexes == Cm.n_complexes(j) && a.n_single == Cm.n_single(j) && ...
            a.n_polyspike == Cm.n_polyspike(j) && a.n_bursts == Cm.n_bursts(j) && abs(r.analyzed_min - Cm.analyzed_min(j)) < 1e-9 ...
            && r.baseline_uV == Cm.baseline_uV(j);
        rows{i} = table(ch.campaign, ch.mouse, ch.edf, ch.region, r.analyzed_min / 60, r.baseline_uV, r.lower_uV, r.p97_x_uV, ...
            r.pct_bp_uV(1), r.pct_bp_uV(2), r.pct_bp_uV(3), r.pctile_of_312_5, size(r.flat, 1), repro, r.total_s, ...
            'VariableNames', {'campaign', 'mouse', 'edf', 'region', 'T_iid_h', 'baseline_uV', 'umbral_actual_uV', 'p97_x_uV', ...
            'p99_xbp_uV', 'p99_5_xbp_uV', 'p99_9_xbp_uV', 'percentil_xbp_de_312_5uV', 'n_tramos_planos', 'actual_reproduce', 'tiempo_s'});
        v = r.variants;
        v = [repmat(table(ch.campaign, ch.mouse, ch.edf, ch.region, r.analyzed_min, 'VariableNames', ...
            {'campaign', 'mouse', 'edf', 'region', 'analyzed_min'}), height(v), 1), v];
        vrows{i} = v;
        nf = size(r.flat, 1);
        frows{i} = table(repmat(ch.campaign, nf, 1), repmat(ch.mouse, nf, 1), repmat(ch.edf, nf, 1), repmat(ch.region, nf, 1), ...
            r.flat(:, 1), r.flat(:, 2), r.flat(:, 1) + ch.start_posix, r.flat(:, 2) + ch.start_posix, ...
            'VariableNames', {'campaign', 'mouse', 'edf', 'region', 'start_s', 'end_s', 'start_posix', 'end_posix'});
    end
    if worker ~= 0
        T = []; V = []; F = [];
        return;
    end
    T = vertcat(rows{:});
    V = vertcat(vrows{:});
    F = vertcat(frows{:});
end

function k = sens_code_key(cfg)
% SHA-1 of the files that compute a channel's sensitivity + its settings.
    sc = cfg.provenance.stats_code;
    f = sc(startsWith(sc(:, 1), 'iids_') | strcmp(sc(:, 1), 'cs_iid_sens_channel.m'), :);
    k = [strjoin(sort(strcat(f(:, 1), ':', f(:, 2)))', ';'), '|', ...
        cfg.provenance.iid_code{strcmp(cfg.provenance.iid_code(:, 1), 'detect_iid.m'), 2}, '|', ...
        mat2str([cfg.flat_thr_uV cfg.flat_min_s cfg.sens_fixed_uV cfg.sens_pct_bp])];
end
