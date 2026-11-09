function [Corig, flags] = cs_iid_original(cfg, P)
% CS_IID_ORIGINAL  The campaigns' own IID results (04_iid, never modified)
% for the channels of the re-detection plan P (Z.channels rows that run):
%   - summary: the channel's row of <campaign>/merged/05_summaries/iid_summary.csv
%     (every analyzed channel, also those with 0 complexes)
%   - complexes and exclusion zones: 04_iid/<EDF>/<base>_IID_results.mat,
%     which detect_iid writes only with >= 1 complex
% A .mat whose summary disagrees with the CSV row (a file left over by an
% earlier run and never overwritten) is not used -> flag; its complexes
% are then unknown (NaN in the before/after decomposition).
%
% Corig: one row per channel of P, columns aligned with cs_iid_rerun_all.m.

    flags = cell(0, 3);
    sumcsv = cell(numel(cfg.campaigns), 1);
    for c = 1:numel(cfg.campaigns)
        f = fullfile(cfg.campaigns(c).root, 'merged', '05_summaries', 'iid_summary.csv');
        T = readtable(f, 'TextType', 'char', 'Delimiter', ',', 'VariableNamingRule', 'preserve');
        T.edf = regexprep(T.source_file, '\.EDF$', '', 'ignorecase');
        sumcsv{c} = T;
    end
    n = height(P);
    rows = cell(n, 1);
    for i = 1:n
        ch = P(i, :);
        ci = find(strcmp({cfg.campaigns.name}, ch.campaign{1}));
        T = sumcsv{ci};
        k = find(strcmp(T.subject_id, ch.mouse{1}) & strcmp(T.region, ch.region{1}) & strcmp(T.edf, ch.edf{1}));
        if isempty(k)
            error('cs_iid_original:NoSummary', '%s %s %s %s: no row in iid_summary.csv', ch.campaign{1}, ch.mouse{1}, ch.edf{1}, ch.region{1});
        end
        if numel(k) > 1
            flags(end+1, :) = {'info', 'iid_summary.csv con filas repetidas', sprintf('%s %s %s %s: %d filas, se usa la ultima', ...
                ch.campaign{1}, ch.mouse{1}, ch.edf{1}, ch.region{1}, numel(k))}; %#ok<AGROW>
            k = k(end);
        end
        s = T(k, :);
        [~, nm] = fileparts(ch.file{1});
        mf = fullfile(cfg.campaigns(ci).root, ch.mouse{1}, '04_iid', ch.edf{1}, [strrep(nm, '_clean', '') '_IID_results.mat']);
        cx = zeros(0, 5);
        ex = zeros(0, 2);
        src = 'csv (0 complejos)';
        known = true;
        if isfile(mf)
            q = load(mf, 'iid_results');
            r = q.iid_results;
            same = r.summary.total_complexes == s.total_complexes && abs(r.summary.analyzed_duration_min - s.analyzed_duration_min) < 1e-6;
            if same
                t = r.spike_complex_table;
                cx = [t.Start_s t.End_s t.N_spikes double(t.Is_polyspike) t.Max_amplitude_mV];
                e = r.exclusion_table;
                ex = [e.start_s e.end_s];
                src = 'mat';
            else
                known = false;
                src = '.mat descartado (no coincide con iid_summary.csv)';
                flags(end+1, :) = {'revisar', '04_iid .mat no coincide con iid_summary.csv', sprintf(['%s %s %s %s: .mat %d complejos / %.2f min, ' ...
                    'CSV %d / %.2f min; se usa el CSV'], ch.campaign{1}, ch.mouse{1}, ch.edf{1}, ch.region{1}, r.summary.total_complexes, ...
                    r.summary.analyzed_duration_min, s.total_complexes, s.analyzed_duration_min)}; %#ok<AGROW>
            end
        elseif s.total_complexes > 0
            known = false;
            src = 'sin .mat';
            flags(end+1, :) = {'revisar', '04_iid sin .mat', sprintf('%s %s %s %s: CSV %d complejos pero no hay .mat', ...
                ch.campaign{1}, ch.mouse{1}, ch.edf{1}, ch.region{1}, s.total_complexes)}; %#ok<AGROW>
        end
        rows{i} = table(ch.campaign, ch.mouse, ch.edf, ch.region, s.baseline_uV, s.lower_threshold_uV, s.total_duration_min, ...
            s.analyzed_duration_min, s.n_exclusion_zones, s.total_complexes, s.n_single, s.n_polyspike, s.n_bursts, {cx}, {ex}, ...
            known, {src}, 'VariableNames', {'campaign', 'mouse', 'edf', 'region', 'baseline_uV', 'lower_uV', 'total_min', ...
            'analyzed_min', 'n_zones', 'n_complexes', 'n_single', 'n_polyspike', 'n_bursts', 'complexes', 'excl', ...
            'complexes_known', 'source'});
    end
    Corig = vertcat(rows{:});
end
