function [edf_tbl, ch_tbl] = cs_load_edf_table(cfg)
% CS_LOAD_EDF_TABLE  Every analyzed (campaign, mouse, EDF) from the
% *_clean.txt headers in 02_clean (header lines only).
%
%   ch_tbl : one row per channel file (campaign, mouse, edf, region, file,
%            start_posix, fs, n_samples, duration_s, n_gaps, gain/case info)
%   edf_tbl: one row per (campaign, mouse, edf) with both regions:
%            start_posix, end_posix (= start + n_samples/fs), duration_s,
%            n_regions, clean_HPCr/clean_HPCl paths, gaps (Nx2 file-relative s)

    rows = {};
    for c = 1:numel(cfg.campaigns)
        camp = cfg.campaigns(c);
        md = dir(camp.root);
        md = md([md.isdir] & ~cellfun(@isempty, regexp({md.name}, '^\d{3}-s$', 'once')));
        for m = 1:numel(md)
            mouse = md(m).name;
            cdir = fullfile(camp.root, mouse, '02_clean');
            ff = dir(fullfile(cdir, '*_clean.txt'));
            for k = 1:numel(ff)
                f = fullfile(cdir, ff(k).name);
                h = cs_read_clean_header(f);
                edf = regexprep(h.source_file, '\.EDF$', '', 'ignorecase');
                fs = str2double(h.fs);
                ns = str2double(h.n_samples);
                gaps = zeros(0, 2);
                ngaps = str2double(h.n_gaps);
                if ngaps > 0
                    gf = fullfile(camp.root, mouse, '01_txt', regexprep(ff(k).name, '_HPC[rl]_clean\.txt$', '_gaps.csv'));
                    g = readtable(gf);
                    gaps = [g.start_s g.end_s];
                end
                rows(end+1, :) = {camp.name, mouse, edf, h.region, f, str2double(h.session_start_unix), ...
                    h.session_start_datetime, fs, ns, ns / fs, ngaps, {gaps}, getf(h, 'case_applied'), ...
                    str2double(getf(h, 'gain_applied')), getf(h, 'quality_class')}; %#ok<AGROW>
            end
        end
    end
    ch_tbl = cell2table(rows, 'VariableNames', {'campaign', 'mouse', 'edf', 'region', 'file', 'start_posix', ...
        'session_start_txt', 'fs', 'n_samples', 'duration_s', 'n_gaps', 'gaps', 'case_applied', 'gain_applied', 'quality_class'});

    [g, keys] = findgroups(ch_tbl(:, {'campaign', 'mouse', 'edf'}));
    n = height(keys);
    start_posix = splitapply(@min, ch_tbl.start_posix, g);
    end_posix = splitapply(@max, ch_tbl.start_posix + ch_tbl.duration_s, g);
    n_regions = splitapply(@numel, ch_tbl.region, g);
    clean_HPCr = repmat({''}, n, 1);
    clean_HPCl = repmat({''}, n, 1);
    gaps = repmat({zeros(0, 2)}, n, 1);
    start_mismatch_s = zeros(n, 1);
    for i = 1:height(ch_tbl)
        j = g(i);
        if strcmp(ch_tbl.region{i}, 'HPCr')
            clean_HPCr{j} = ch_tbl.file{i};
        elseif strcmp(ch_tbl.region{i}, 'HPCl')
            clean_HPCl{j} = ch_tbl.file{i};
        end
        gaps{j} = [gaps{j}; ch_tbl.gaps{i}];
        start_mismatch_s(j) = max(start_mismatch_s(j), abs(ch_tbl.start_posix(i) - start_posix(j)));
    end
    edf_tbl = [keys, table(start_posix, end_posix, end_posix - start_posix, n_regions, clean_HPCr, clean_HPCl, ...
        gaps, start_mismatch_s, 'VariableNames', {'start_posix', 'end_posix', 'duration_s', 'n_regions', ...
        'clean_HPCr', 'clean_HPCl', 'gaps', 'start_mismatch_s'})];
    edf_tbl.analyzed = ~cellfun(@isempty, edf_tbl.clean_HPCr) & ~cellfun(@isempty, edf_tbl.clean_HPCl);
end

function v = getf(h, k)
    if isfield(h, k)
        v = h.(k);
    else
        v = '';
    end
end
