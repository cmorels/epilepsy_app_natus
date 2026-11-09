function r = spec_resolve_files(root, subject, sesion, region, signal_source, tz)
% SPEC_RESOLVE_FILES  (subject, session, region) -> the pipeline files of
% that channel, found RECURSIVELY under <root>\<subject>\ (no fixed
% layout: 01_txt\, 02_clean\, one folder per EDF, ... all work).
%
%   r = spec_resolve_files(root, subject, sesion, region, signal_source, tz)
%   spec_resolve_files('reset')            % forget the cached file index
%
%   subject : exact folder / file prefix ('005-s' ~= '005')
%   sesion  : 'yyyyMMdd_HHmmss' (the session part of the file names) or
%             the EDF file name (with or without extension), matched
%             against the "# source_file" header line of the txt files
%   region  : 'HPCr' | 'HPCl' | ... (one region; 'both' is expanded by the caller)
%   signal_source : 'clean' | 'raw' -- which txt r.txt points to
%   tz      : TimeZone used only if a header has no "timezone" line
%
% r: subject, region, session_key, clean, raw, txt, iid_mat ('' when
%    absent), source_file, session_start (zoned), wall_start (unzoned:
%    the clock time Natus shows at t_rel = 0), fs, n_samples, hdr (header
%    of r.txt), clean_hdr (header of the clean txt, [] if none), notes.
%
% Errors (never guesses): subject folder missing, session not found,
% more than one file for the same (subject, session, region, kind),
% missing txt for the requested signal_source.

    persistent index
    if ischar(root) && strcmp(root, 'reset')
        index = [];
        peek('reset');
        r = [];
        return;
    end
    if isempty(index)
        index = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end

    subject = strtrim(char(subject));
    sesion = strtrim(char(sesion));
    region = strtrim(char(region));
    sub_dir = fullfile(root, subject);
    if ~isfolder(sub_dir)
        error('spec_resolve_files:NoSubject', 'Subject folder not found: %s', sub_dir);
    end

    key = [root '|' subject];
    if ~isKey(index, key)
        index(key) = build_index(sub_dir, subject);
    end
    T = index(key);

    % ---- session: key or EDF name ------------------------------------
    notes = {};
    if ~isempty(regexp(sesion, '^\d{8}_\d{6}$', 'once'))
        session_key = sesion;
        if ~any(strcmp(T.session_key, session_key))
            error('spec_resolve_files:NoSession', 'No txt for subject %s, session %s under %s', ...
                subject, session_key, sub_dir);
        end
    else
        session_key = session_from_edf(T, sesion, subject, sub_dir);
    end

    rows = T(strcmp(T.session_key, session_key) & strcmp(T.region, region), :);
    clean = pick_one(rows(strcmp(rows.kind, 'clean'), :), 'clean txt', subject, session_key, region);
    raw = pick_one(rows(strcmp(rows.kind, 'raw'), :), 'raw txt (01_txt)', subject, session_key, region);
    iid = pick_one(rows(strcmp(rows.kind, 'iid'), :), 'IID_results.mat', subject, session_key, region);

    switch signal_source
        case 'clean'
            if isempty(clean)
                error('spec_resolve_files:NoClean', ['No %s_%s_%s_clean.txt under %s. Run the pipeline''s ' ...
                    'import + cleaning first (run_pipeline_edf / run_all_subjects).'], subject, session_key, region, sub_dir);
            end
            txt = clean;
        case 'raw'
            if isempty(raw)
                error('spec_resolve_files:NoRaw', ['No %s_%s_%s.txt (01_txt) under %s. Run the pipeline''s ' ...
                    'import first (edf_import via run_pipeline_edf / run_all_subjects).'], subject, session_key, region, sub_dir);
            end
            if isempty(clean)
                error('spec_resolve_files:NoCleanForRaw', ['signal_source = ''raw'' needs the clean txt too ' ...
                    '(its header gives gain_applied / notch_applied): no %s_%s_%s_clean.txt under %s.'], ...
                    subject, session_key, region, sub_dir);
            end
            txt = raw;
        otherwise
            error('spec_resolve_files:BadSource', 'signal_source must be ''clean'' or ''raw'' (got ''%s'').', signal_source);
    end

    hdr = peek(txt);
    clean_hdr = [];
    if ~isempty(clean)
        clean_hdr = peek(clean);
    end

    ss = hdr.session_start;
    if ~isdatetime(ss) || isnat(ss)
        error('spec_resolve_files:NoSessionStart', '%s has no usable session_start in its header.', txt);
    end
    if ~isfield(hdr.meta, 'timezone') || isempty(strtrim(hdr.meta.timezone))
        ss.TimeZone = tz;
        notes{end+1} = sprintf('header has no timezone; session_start interpreted in %s', tz);
    end

    r = struct();
    r.subject = subject;
    r.region = region;
    r.session_key = session_key;
    r.clean = clean;
    r.raw = raw;
    r.txt = txt;
    r.iid_mat = iid;
    r.source_file = field_text(hdr.meta, 'source_file');
    r.session_start = ss;
    r.wall_start = datetime(ss.Year, ss.Month, ss.Day, ss.Hour, ss.Minute, ss.Second);
    r.wall_start.Format = 'dd/MM/yyyy HH:mm:ss.SSS';
    r.fs = hdr.fs;
    r.n_samples = str2double(field_text(hdr.meta, 'n_samples'));
    r.hdr = hdr;
    r.clean_hdr = clean_hdr;
    r.notes = notes;
end

%% ======================================================================
function h = peek(p)
% Header of p via spec_peek_header, cached for the run.
    persistent peeks
    if strcmp(p, 'reset') || isempty(peeks)
        peeks = containers.Map('KeyType', 'char', 'ValueType', 'any');
        if strcmp(p, 'reset')
            h = [];
            return;
        end
    end
    if ~isKey(peeks, p)
        peeks(p) = spec_peek_header(p);
    end
    h = peeks(p);
end

function k = session_from_edf(T, edf, subject, sub_dir)
    [~, stem, ext] = fileparts(edf);
    if ~any(strcmpi(ext, {'.edf', ''}))
        stem = edf;   % a name with dots that is not an extension
    end
    cand = T(ismember(T.kind, {'clean', 'raw'}), :);
    hit = false(height(cand), 1);
    for i = 1:height(cand)
        h = peek(cand.path{i});
        sf = field_text(h.meta, 'source_file');
        [~, sf_stem] = fileparts(sf);
        hit(i) = strcmpi(strtrim(sf_stem), strtrim(stem));
    end
    keys_hit = unique(cand.session_key(hit));
    if isempty(keys_hit)
        error('spec_resolve_files:NoSession', ['No txt of subject %s under %s comes from EDF "%s" ' ...
            '(checked the "# source_file" header of every txt). Run the pipeline''s import + cleaning for it first.'], ...
            subject, sub_dir, edf);
    end
    if numel(keys_hit) > 1
        error('spec_resolve_files:AmbiguousSession', 'EDF "%s" maps to several sessions of %s: %s', ...
            edf, subject, strjoin(keys_hit, ', '));
    end
    k = keys_hit{1};
end

%% ======================================================================
function T = build_index(sub_dir, subject)
    esc = regexptranslate('escape', subject);
    txt = dir(fullfile(sub_dir, '**', [subject '_*.txt']));
    mat = dir(fullfile(sub_dir, '**', [subject '_*_IID_results.mat']));
    files = [txt; mat];
    n = numel(files);
    kind = cell(n, 1); session_key = cell(n, 1); region = cell(n, 1); path = cell(n, 1);
    keep = false(n, 1);
    for i = 1:n
        f = files(i);
        if contains(f.folder, '06_spectrograms')
            continue;
        end
        tk = regexp(f.name, ['^' esc '_(\d{8}_\d{6})_(.+?)(_clean\.txt|\.txt|_IID_results\.mat)$'], 'tokens', 'once');
        if isempty(tk)
            continue;
        end
        switch tk{3}
            case '_clean.txt'
                kind{i} = 'clean';
            case '.txt'
                kind{i} = 'raw';
            otherwise
                kind{i} = 'iid';
        end
        session_key{i} = tk{1};
        region{i} = tk{2};
        path{i} = fullfile(f.folder, f.name);
        keep(i) = true;
    end
    T = table(kind(keep), session_key(keep), region(keep), path(keep), ...
        'VariableNames', {'kind', 'session_key', 'region', 'path'});
end

function p = pick_one(rows, what, subject, session_key, region)
    if height(rows) == 0
        p = '';
    elseif height(rows) == 1
        p = rows.path{1};
    else
        error('spec_resolve_files:Ambiguous', 'More than one %s for %s / %s / %s:\n  %s', ...
            what, subject, session_key, region, strjoin(rows.path, sprintf('\n  ')));
    end
end

function s = field_text(meta, name)
    if isfield(meta, name)
        s = strtrim(char(meta.(name)));
    else
        s = '';
    end
end
