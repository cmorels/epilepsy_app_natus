function [rows, review, log] = cs_load_detections(cfg, edf_tbl)
% CS_LOAD_DETECTIONS  Every detection row of both campaigns with the
% user's labels attached.
%
% Times and features come from <campaign>/merged/05_summaries/
% seizures_events.csv, written by the pipeline BEFORE the review (full
% precision). Labels (AE user_check / AF comments) come from each mouse's
% reviewed 05_summaries/seizures_events.csv, which Excel re-saved (seconds
% and hours lost in the datetime text columns, start_s rounded to ~5e-5 s).
% The two are joined row by row in file order; every row's (source_file,
% region, seizure_id) must match and |start_s difference| < 1e-3 s, else
% this errors. Label columns are found by header NAME (112-s/114-s carry a
% duplicate user_check in AF and the comments in AG; some files say
% 'comment'; 124-s has an unnamed comments column at position 32).
%
% Absolute time = session start of the EDF (*_clean.txt header,
% session_start_unix) + start_s, i.e. elapsed time, checked against the
% pipeline's own start_abs text (max difference reported in log).

    log = {};
    all_rows = {};
    rev = {};
    for c = 1:numel(cfg.campaigns)
        camp = cfg.campaigns(c);
        mf = fullfile(camp.root, 'merged', '05_summaries', 'seizures_events.csv');
        opts = detectImportOptions(mf, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
        opts = setvartype(opts, opts.VariableNames, 'char');
        M = readtable(mf, opts);
        mice = unique(M.subject_id, 'stable');
        md = dir(camp.root);
        md = md([md.isdir] & ~cellfun(@isempty, regexp({md.name}, '^\d{3}-s$', 'once')));
        mice = unique([mice; {md.name}'], 'stable');
        for m = 1:numel(mice)
            mouse = mice{m};
            Mm = M(strcmp(M.subject_id, mouse), :);
            lf = 'seizures_events.csv';
            ov = cfg.label_file_override;
            for k = 1:size(ov, 1)
                if strcmp(ov{k, 1}, camp.name) && strcmp(ov{k, 2}, mouse)
                    lf = ov{k, 3};
                end
            end
            lfile = fullfile(camp.root, mouse, '05_summaries', lf);
            if isfield(cfg, 'labels_snapshot_dir') && ~isempty(cfg.labels_snapshot_dir)
                lfile = fullfile(cfg.labels_snapshot_dir, camp.name, mouse, lf);   % cs_label_snapshot.m
            end
            [lab, ae_col, af_col] = read_labels(lfile, height(Mm));
            % ---- join by order, verify keys
            if height(lab.keys) ~= height(Mm)
                error('cs_load_detections:RowCount', '%s %s: %d label rows vs %d merged rows', camp.name, mouse, height(lab.keys), height(Mm));
            end
            for i = 1:height(Mm)
                ok = strcmp(lab.keys.source_file{i}, Mm.source_file{i}) && strcmp(lab.keys.region{i}, Mm.region{i}) ...
                    && str2double(lab.keys.seizure_id{i}) == str2double(Mm.seizure_id{i}) ...
                    && abs(str2double(lab.keys.start_s{i}) - str2double(Mm.start_s{i})) < 1e-3;
                if ~ok
                    error('cs_load_detections:KeyMismatch', '%s %s row %d does not match the merged row', camp.name, mouse, i);
                end
            end
            n_lab = 0;
            for i = 1:height(Mm)
                edf = regexprep(Mm.source_file{i}, '\.EDF$', '', 'ignorecase');
                j = find(strcmp(edf_tbl.campaign, camp.name) & strcmp(edf_tbl.mouse, mouse) & strcmp(edf_tbl.edf, edf));
                if numel(j) ~= 1
                    error('cs_load_detections:NoEdf', '%s %s %s: EDF not in 02_clean', camp.name, mouse, edf);
                end
                s0 = edf_tbl.start_posix(j);
                st = str2double(Mm.start_s{i});
                en = str2double(Mm.end_s{i});
                [label, ae_flag] = norm_label(lab.ae{i});
                [comment, mapped] = norm_comment(lab.af{i}, cfg.comment_map);
                n_lab = n_lab + ~isempty(label);
                all_rows(end+1, :) = {camp.name, mouse, edf, Mm.region{i}, str2double(Mm.seizure_id{i}), ...
                    Mm.channel_seizure_id{i}, Mm.detection_status{i}, st, en, s0 + st, s0 + en, ...
                    str2double(Mm.ll_ratio{i}), str2double(Mm.peak_energy_ratio{i}), Mm.accepted_in_regions{i}, ...
                    lab.ae{i}, lab.af{i}, label, comment, ae_flag, mapped, Mm.start_abs{i}}; %#ok<AGROW>
            end
            if height(Mm) == 0
                status = 'no_detections';
            elseif isempty(ae_col)
                status = 'not_reviewed';
            elseif n_lab == height(Mm)
                status = 'reviewed';
            elseif n_lab == 0
                status = 'not_reviewed';
            else
                status = 'partial';
            end
            rev(end+1, :) = {camp.name, mouse, lf, height(Mm), n_lab, ae_col, af_col, status}; %#ok<AGROW>
        end
    end
    rows = cell2table(all_rows, 'VariableNames', {'campaign', 'mouse', 'edf', 'region', 'seizure_id', ...
        'channel_seizure_id', 'detection_status', 'start_s', 'end_s', 'start_posix', 'end_posix', 'll_ratio', ...
        'peak_energy_ratio', 'accepted_in_regions', 'ae_raw', 'af_raw', 'label', 'comment', 'ae_unrecognized', ...
        'comment_mapped', 'start_abs_pipeline'});
    review = cell2table(rev, 'VariableNames', {'campaign', 'mouse', 'label_file', 'n_rows', 'n_rows_labeled', ...
        'ae_column', 'af_column', 'row_status'});

    % ---- absolute-time check against the pipeline's own text
    t = datetime(rows.start_abs_pipeline, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss.SSS', 'TimeZone', cfg.tz, 'Locale', 'en_US');
    d = abs(posixtime(t) - rows.start_posix);
    log{end+1} = sprintf('abs-time check: %d rows, max |session_start+start_s - pipeline start_abs| = %.4f s', height(rows), max([d; 0]));
    if any(d > 0.002)
        error('cs_load_detections:TimeCheck', '%d rows differ from the pipeline start_abs by > 2 ms', sum(d > 0.002));
    end
end

function [lab, ae_col, af_col] = read_labels(file, n_expected)
    [hdr, body] = cs_read_csv_raw(file);
    hn = lower(strtrim(hdr));
    keep = ~cellfun(@(s) isempty(strtrim(s)), body(:, 1));  % drop trailing empty rows (110-s)
    body = body(keep, :);
    col = @(name) find(strcmp(hn, name), 1);
    lab.keys = table(body(:, col('source_file')), body(:, col('region')), body(:, col('seizure_id')), ...
        body(:, col('start_s')), 'VariableNames', {'source_file', 'region', 'seizure_id', 'start_s'});
    ia = col('user_check');
    ic = find(strcmp(hn, 'comments') | strcmp(hn, 'comment'), 1);
    if isempty(ic) && ~isempty(ia) && numel(hdr) >= 32 && isempty(hn{32}) && ia == 31
        ic = 32;  % unnamed comments column (124-s)
    end
    ae_col = ''; af_col = '';
    lab.ae = repmat({''}, size(body, 1), 1);
    lab.af = repmat({''}, size(body, 1), 1);
    if ~isempty(ia)
        ae_col = sprintf('%d:%s', ia, hdr{ia});
        lab.ae = body(:, ia);
    end
    if ~isempty(ic)
        af_col = sprintf('%d:%s', ic, hdr{ic});
        lab.af = body(:, ic);
    end
    if size(body, 1) ~= n_expected && n_expected > 0 && size(body, 1) == 0
        error('cs_load_detections:Empty', '%s has no rows', file);
    end
end

function [label, unrec] = norm_label(s)
    s = upper(strtrim(s));
    unrec = false;
    if any(strcmp(s, {'TP', 'FP'}))
        label = s;
    else
        label = '';
        unrec = ~isempty(s);
    end
end

function [c, mapped] = norm_comment(s, cmap)
    s0 = lower(strtrim(s));
    mapped = '';
    switch s0
        case ''
            c = '';
        case {'iid', 'iids'}
            c = 'IIDs';
        case {'grooming', 'scratching', 'grooming/scratching', 'grooming / scratching', 'grooming-scratching', 'scratching/grooming'}
            c = 'grooming/scratching';
        case {'disconnect', 'disconnected', 'disconnection'}
            c = 'disconnect';
        case {'other', 'others'}
            c = 'other';
        otherwise
            k = find(strcmp(cmap(:, 1), s0), 1);
            if ~isempty(k)
                c = cmap{k, 2};
                mapped = s0;
            else
                c = ['unrecognized:' strtrim(s)];
            end
    end
end
