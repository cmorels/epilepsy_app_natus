function T = cs_label_changes(cfg, snap, det_rows)
% CS_LABEL_CHANGES  Label files of the snapshot vs the labels v1 used
% (det_rows of cfg.v1_dir/phase1_results.mat): per file, whether any row's
% AE/AF (raw text) differs, how many rows, and in which EDFs. Rows are
% matched by (campaign, mouse, EDF, region, seizure_id, start_s).
%
%   T: the snapshot manifest + changed_vs_v1, n_rows_changed, edfs_changed,
%      modified_after_v1 (file time later than v1's last run start)

    q = load(fullfile(cfg.v1_dir, 'phase1_results.mat'), 'det_rows', 'cfg');
    v1 = q.det_rows;
    v1_start = datetime(q.cfg.run_started, 'InputFormat', 'yyyy-MM-dd HH:mm:ss Z', 'TimeZone', cfg.tz);
    key = @(R) strcat(R.campaign, '|', R.mouse, '|', R.edf, '|', R.region, '|', ...
        arrayfun(@(x) sprintf('%d', x), R.seizure_id, 'UniformOutput', false), '|', ...
        arrayfun(@(x) sprintf('%.3f', x), R.start_s, 'UniformOutput', false));
    k_new = key(det_rows);
    k_old = key(v1);
    [tf, j] = ismember(k_new, k_old);
    same = false(height(det_rows), 1);
    same(tf) = strcmp(det_rows.ae_raw(tf), v1.ae_raw(j(tf))) & strcmp(det_rows.af_raw(tf), v1.af_raw(j(tf)));
    T = snap;
    n = height(T);
    T.changed_vs_v1 = false(n, 1);
    T.n_rows_changed = zeros(n, 1);
    T.edfs_changed = repmat({''}, n, 1);
    T.modified_after_v1 = false(n, 1);
    T.read_by_analysis = true(n, 1);   % false: a second file of the unit, not the one cfg.label_file_override names
    for i = 1:n
        ov = cfg.label_file_override;
        for k = 1:size(ov, 1)
            if strcmp(ov{k, 1}, T.campaign{i}) && strcmp(ov{k, 2}, T.mouse{i})
                T.read_by_analysis(i) = strcmp(T.file{i}, ov{k, 3});
            end
        end
        if ~strcmp(T.file{i}, 'seizures_events.csv') && ~any(strcmp(ov(:, 3), T.file{i}))
            T.read_by_analysis(i) = false;
        end
        if ~T.read_by_analysis(i)
            continue;
        end
        s = strcmp(det_rows.campaign, T.campaign{i}) & strcmp(det_rows.mouse, T.mouse{i});
        ch = s & ~same;
        T.n_rows_changed(i) = sum(ch);
        T.changed_vs_v1(i) = any(ch) || sum(strcmp(v1.campaign, T.campaign{i}) & strcmp(v1.mouse, T.mouse{i})) ~= sum(s);
        T.edfs_changed{i} = strjoin(unique(det_rows.edf(ch))', ' | ');
        T.modified_after_v1(i) = datetime(T.modified{i}, 'InputFormat', 'yyyy-MM-dd HH:mm:ss', 'TimeZone', cfg.tz) > v1_start;
    end
end
