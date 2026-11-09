function T = reorganize_by_edf(campaign_root, apply)
% REORGANIZE_BY_EDF  Move the files of each animal's 03_seizures/ and 04_iid/
% of a campaign output (run_all_subjects.m) into one subfolder per EDF,
% named after the EDF file name without extension:
%
%   <campaign>/<animal>/03_seizures/<EDF>/...
%   <campaign>/<animal>/04_iid/<EDF>/...
%   e.g. 128-s/03_seizures/20260801~ TRAP_39376c78-f43d-49fd-97f7-c768d2a533db/128-s_20260801_185338_HPCl_seizures.mat
%
%   T = reorganize_by_edf(campaign_root)          % DRY RUN: prints what it would do, moves/writes nothing
%   T = reorganize_by_edf(campaign_root, true)    % move + log CSV
%
% campaign_root : a campaign output, one subfolder per animal. Subfolders
%                 without 03_seizures/ or 04_iid/ (merged/, ...) are ignored.
%
% A file belongs to an EDF through its name prefix <animal>_<yyyyMMdd>_<HHmmss>_
% (the session start, edf_import.m) and the "# source_file = ..." header
% line of the 01_txt/<animal>_<yyyyMMdd>_<HHmmss>_<region>.txt edf_import.m
% wrote for that session. Nothing is guessed: a file without that prefix
% (e.g. Thumbs.db) stays where it is ('raiz'), and so does a file whose
% session has no txt or whose txt files name different EDFs ('ambiguo').
% Never overwrites: an existing destination is reported ('conflicto') and
% that file is not moved. Nothing is deleted. Safe to run twice: files
% already in their EDF subfolder are not looked at again.
%
% With apply = true every file looked at is logged, as it is processed, to
% <campaign_root>/reorganize_by_edf_<yyyyMMdd_HHmmss>.csv (raton, etapa,
% origen, destino, estado); undo_reorganize_by_edf.m moves them back from it.
%
% OUTPUT T: table raton, etapa, origen, destino, estado ('a mover' (dry run)
% | 'movido' | 'raiz' | 'ambiguo' | 'conflicto' | 'error: <message>').

    if nargin < 2
        apply = false;
    end
    if ~isfolder(campaign_root)
        error('reorganize_by_edf:NoFolder', 'Campaign folder not found: %s', campaign_root);
    end
    stages = {'03_seizures', '04_iid'};
    log_file = fullfile(campaign_root, ['reorganize_by_edf_' datestr(now, 'yyyymmdd_HHMMSS') '.csv']); %#ok<TNOW1,DATST>

    fid = -1;
    if apply
        fid = fopen(log_file, 'w', 'n', 'UTF-8');
        if fid < 0
            error('reorganize_by_edf:LogOpen', 'Could not create the log %s', log_file);
        end
        closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
        fprintf(fid, 'raton,etapa,origen,destino,estado\n');
    end

    rows = cell(0, 5);
    visited = cell(0, 2);
    animals = animal_folders(campaign_root, stages);
    for a = 1:numel(animals)
        id = animals{a};
        adir = fullfile(campaign_root, id);
        edf_of = session_map(adir, id);
        edf_stems = unique(values(edf_of));
        rx = ['^' regexptranslate('escape', id) '_(\d{8}_\d{6})_'];
        for s = 1:numel(stages)
            sdir = fullfile(adir, stages{s});
            if ~isfolder(sdir)
                continue;
            end
            visited(end+1, :) = {id, stages{s}}; %#ok<AGROW>
            d = dir(sdir);
            d = d(~ismember({d.name}, {'.', '..'}));
            for k = 1:numel(d)
                src = fullfile(sdir, d(k).name);
                if d(k).isdir
                    if ismember(d(k).name, edf_stems)
                        continue;  % an EDF subfolder (earlier run)
                    end
                    row = {id, stages{s}, src, '', 'raiz'};
                else
                    tok = regexp(d(k).name, rx, 'tokens', 'once');
                    if isempty(tok)
                        row = {id, stages{s}, src, '', 'raiz'};
                    elseif ~isKey(edf_of, tok{1})
                        row = {id, stages{s}, src, '', 'ambiguo'};
                    else
                        dst_dir = fullfile(sdir, edf_of(tok{1}));
                        dst = fullfile(dst_dir, d(k).name);
                        if exist(dst, 'file') ~= 0  % a file or a folder
                            st = 'conflicto';
                        elseif ~apply
                            st = 'a mover';
                        else
                            st = move_one(src, dst_dir, dst);
                        end
                        row = {id, stages{s}, src, dst, st};
                    end
                end
                rows(end+1, :) = row; %#ok<AGROW>
                if apply
                    write_row(fid, row);
                end
            end
        end
    end

    T = cell2table(rows, 'VariableNames', {'raton', 'etapa', 'origen', 'destino', 'estado'});
    print_summary(T, visited, apply, campaign_root);
    if apply
        fprintf('Log: %s\n', log_file);
    else
        fprintf('DRY RUN: nada movido ni escrito. Para aplicar: reorganize_by_edf(campaign_root, true)\n');
    end
end

%% ======================================================================
function names = animal_folders(campaign_root, stages)
    d = dir(campaign_root);
    d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
    names = {};
    for k = 1:numel(d)
        p = fullfile(campaign_root, d(k).name);
        if any(cellfun(@(s) isfolder(fullfile(p, s)), stages))
            names{end+1} = d(k).name; %#ok<AGROW>
        end
    end
end

function edf_of = session_map(adir, id)
% '<yyyyMMdd>_<HHmmss>' -> EDF name without extension, from the source_file
% header line of every 01_txt/<id>_<yyyyMMdd>_<HHmmss>_<region>.txt. A
% session whose txt files give no source_file or different EDFs is left
% out, so its files are reported 'ambiguo'.
    edf_of = containers.Map('KeyType', 'char', 'ValueType', 'char');
    bad = {};
    rx = ['^' regexptranslate('escape', id) '_(\d{8}_\d{6})_'];
    t = dir(fullfile(adir, '01_txt', '*.txt'));
    for k = 1:numel(t)
        tok = regexp(t(k).name, rx, 'tokens', 'once');
        if isempty(tok)
            continue;
        end
        [~, stem] = fileparts(header_source_file(fullfile(t(k).folder, t(k).name)));
        key = tok{1};
        if isempty(stem) || (isKey(edf_of, key) && ~strcmp(edf_of(key), stem))
            bad{end+1} = key; %#ok<AGROW>
        elseif ~isKey(edf_of, key)
            edf_of(key) = stem;
        end
    end
    bad = unique(bad);
    for i = 1:numel(bad)
        if isKey(edf_of, bad{i})
            remove(edf_of, bad{i});
        end
    end
end

function src = header_source_file(txt_path)
    src = '';
    fid = fopen(txt_path, 'r');
    if fid < 0
        return;
    end
    for i = 1:30
        ln = fgetl(fid);
        if ~ischar(ln) || ~startsWith(ln, '#')
            break;
        end
        if startsWith(ln, '# source_file = ')
            src = strtrim(ln(numel('# source_file = ') + 1:end));
            break;
        end
    end
    fclose(fid);
end

function st = move_one(src, dst_dir, dst)
    if ~isfolder(dst_dir)
        [ok, msg] = mkdir(dst_dir);
        if ~ok
            st = ['error: ' msg];
            return;
        end
    end
    [ok, msg] = movefile(src, dst);
    if ok
        st = 'movido';
    else
        st = ['error: ' msg];
    end
end

function write_row(fid, row)
    q = cellfun(@(v) ['"' strrep(v, '"', '""') '"'], row, 'UniformOutput', false);
    fprintf(fid, '%s\n', strjoin(q, ','));
end

function print_summary(T, visited, apply, campaign_root)
    move_state = ternary(apply, 'movido', 'a mover');
    fprintf('\nreorganize_by_edf: %s (%s)\n', campaign_root, ternary(apply, 'APLICADO', 'DRY RUN'));
    fprintf('%-8s %-12s %8s %6s %6s %9s %10s %6s\n', 'raton', 'etapa', move_state, 'n_edf', 'raiz', 'ambiguos', 'conflictos', 'error');
    for i = 1:size(visited, 1)
        r = strcmp(T.raton, visited{i, 1}) & strcmp(T.etapa, visited{i, 2});
        mv = r & strcmp(T.estado, move_state);
        n_edf = numel(unique(cellfun(@fileparts, T.destino(mv), 'UniformOutput', false)));
        fprintf('%-8s %-12s %8d %6d %6d %9d %10d %6d\n', visited{i, 1}, visited{i, 2}, nnz(mv), n_edf, ...
            nnz(r & strcmp(T.estado, 'raiz')), nnz(r & strcmp(T.estado, 'ambiguo')), ...
            nnz(r & strcmp(T.estado, 'conflicto')), nnz(r & startsWith(T.estado, 'error')));
    end
    fprintf('TOTAL: %d %s, %d raiz, %d ambiguos, %d conflictos, %d errores (%d raton(es))\n', ...
        nnz(strcmp(T.estado, move_state)), move_state, nnz(strcmp(T.estado, 'raiz')), nnz(strcmp(T.estado, 'ambiguo')), ...
        nnz(strcmp(T.estado, 'conflicto')), nnz(startsWith(T.estado, 'error')), numel(unique(visited(:, 1))));
    other = find(~strcmp(T.estado, move_state));
    if ~isempty(other)
        fprintf('Se quedan donde estan:\n');
        for i = other'
            fprintf('  %-10s %s\n', T.estado{i}, T.origen{i}(numel(campaign_root) + 2:end));
        end
    end
end

function v = ternary(c, a, b)
    if c
        v = a;
    else
        v = b;
    end
end
