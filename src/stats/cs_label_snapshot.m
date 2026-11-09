function [snap, snap_dir] = cs_label_snapshot(cfg)
% CS_LABEL_SNAPSHOT  One snapshot of every label file (each mouse's
% 05_summaries/seizures_events*.csv, both campaigns), copied to
% cfg.out_dir/labels_snapshot/<campaign>/<mouse>/<file> with its SHA-1,
% size and modification time. The whole run reads labels from here, never
% from the campaign folders (which are edited while the analysis runs).
%
% An existing snapshot is reused as is unless cfg.labels_snapshot_refresh
% (so separate MATLAB processes of the same run read the same labels);
% its manifest is labels_snapshot/manifest.csv.

    snap_dir = fullfile(cfg.out_dir, 'labels_snapshot');
    man = fullfile(snap_dir, 'manifest.csv');
    if isfile(man) && ~cfg.labels_snapshot_refresh
        snap = readtable(man, 'TextType', 'char', 'Delimiter', ',', 'DatetimeType', 'text');
        snap.reused = true(height(snap), 1);
        return;
    end
    if isfolder(snap_dir)
        rmdir(snap_dir, 's');
    end
    rows = cell(0, 9);
    taken = char(datetime('now', 'TimeZone', cfg.tz, 'Format', 'yyyy-MM-dd HH:mm:ss Z'));
    for c = 1:numel(cfg.campaigns)
        camp = cfg.campaigns(c);
        lf = dir(fullfile(camp.root, '*-s', '05_summaries', 'seizures_events*.csv'));
        for k = 1:numel(lf)
            [~, mouse] = fileparts(fileparts(lf(k).folder));
            src = fullfile(lf(k).folder, lf(k).name);
            dst_dir = fullfile(snap_dir, camp.name, mouse);
            if ~isfolder(dst_dir), mkdir(dst_dir); end
            dst = fullfile(dst_dir, lf(k).name);
            copyfile(src, dst);
            rows(end+1, :) = {camp.name, mouse, lf(k).name, src, dst, sha1_file(dst), lf(k).bytes, ...
                char(datetime(lf(k).datenum, 'ConvertFrom', 'datenum', 'Format', 'yyyy-MM-dd HH:mm:ss')), taken}; %#ok<AGROW>
        end
    end
    snap = cell2table(rows, 'VariableNames', {'campaign', 'mouse', 'file', 'source_path', 'snapshot_path', 'sha1', ...
        'bytes', 'modified', 'snapshot_taken'});
    writetable(snap, man);
    snap.reused = false(height(snap), 1);
end

function h = sha1_file(f)
    fid = fopen(f, 'r');
    b = fread(fid, Inf, '*uint8');
    fclose(fid);
    md = java.security.MessageDigest.getInstance('SHA-1');
    h = lower(reshape(dec2hex(typecast(md.digest(b), 'uint8'))', 1, []));
end
