function T = undo_reorganize_by_edf(log_csv, apply)
% UNDO_REORGANIZE_BY_EDF  Move back what reorganize_by_edf.m moved, from
% the log CSV it wrote (<campaign>/reorganize_by_edf_<yyyyMMdd_HHmmss>.csv).
%
%   T = undo_reorganize_by_edf(log_csv)          % DRY RUN: prints what it would do, moves/writes nothing
%   T = undo_reorganize_by_edf(log_csv, true)    % move back + log CSV
%
% Only the rows with estado = 'movido' are considered. Each file goes from
% destino back to origen only when destino exists and origen does not:
% nothing is overwritten and nothing is deleted (the EDF subfolders stay,
% empty). Safe to run twice: a file already back is reported 'ya devuelto'.
%
% With apply = true every row is logged, as it is processed, to
% <log_csv folder>/<log_csv name>_undo_<yyyyMMdd_HHmmss>.csv (same columns).
%
% OUTPUT T: table raton, etapa, origen, destino, estado ('a devolver' (dry
% run) | 'devuelto' | 'ya devuelto' | 'conflicto' (both exist) |
% 'no encontrado' (neither exists) | 'error: <message>').

    if nargin < 2
        apply = false;
    end
    if exist(log_csv, 'file') ~= 2
        error('undo_reorganize_by_edf:NoLog', 'Log not found: %s', log_csv);
    end
    cols = {'raton', 'etapa', 'origen', 'destino', 'estado'};
    opts = detectImportOptions(log_csv, 'Delimiter', ',', 'VariableNamingRule', 'preserve');
    missing = setdiff(cols, opts.VariableNames);
    if ~isempty(missing)
        error('undo_reorganize_by_edf:BadLog', 'Log %s is missing column(s): %s', log_csv, strjoin(missing, ', '));
    end
    opts.SelectedVariableNames = cols;
    opts = setvartype(opts, cols, 'char');
    L = readtable(log_csv, opts);
    L = L(strcmp(L.estado, 'movido'), :);

    [folder, name] = fileparts(log_csv);
    out_file = fullfile(folder, [name '_undo_' datestr(now, 'yyyymmdd_HHMMSS') '.csv']); %#ok<TNOW1,DATST>
    fid = -1;
    if apply
        fid = fopen(out_file, 'w', 'n', 'UTF-8');
        if fid < 0
            error('undo_reorganize_by_edf:LogOpen', 'Could not create the log %s', out_file);
        end
        closer = onCleanup(@() fclose(fid)); %#ok<NASGU>
        fprintf(fid, 'raton,etapa,origen,destino,estado\n');
    end

    rows = cell(height(L), 5);
    for i = 1:height(L)
        src = L.origen{i};
        dst = L.destino{i};
        at_dst = exist(dst, 'file') == 2;
        at_src = exist(src, 'file') ~= 0;
        if at_dst && at_src
            st = 'conflicto';
        elseif at_dst && ~apply
            st = 'a devolver';
        elseif at_dst
            [ok, msg] = movefile(dst, src);
            if ok
                st = 'devuelto';
            else
                st = ['error: ' msg];
            end
        elseif at_src
            st = 'ya devuelto';
        else
            st = 'no encontrado';
        end
        rows(i, :) = {L.raton{i}, L.etapa{i}, src, dst, st};
        if apply
            write_row(fid, rows(i, :));
        end
    end

    T = cell2table(rows, 'VariableNames', cols);
    [sts, ~, ic] = unique(T.estado);
    fprintf('\nundo_reorganize_by_edf: %s (%s), %d fila(s) movido\n', log_csv, ternary(apply, 'APLICADO', 'DRY RUN'), height(T));
    for k = 1:numel(sts)
        fprintf('  %-14s %d\n', sts{k}, nnz(ic == k));
    end
    if apply
        fprintf('Log: %s\n', out_file);
    else
        fprintf('DRY RUN: nada movido ni escrito. Para aplicar: undo_reorganize_by_edf(log_csv, true)\n');
    end
end

function write_row(fid, row)
    q = cellfun(@(v) ['"' strrep(v, '"', '""') '"'], row, 'UniformOutput', false);
    fprintf(fid, '%s\n', strjoin(q, ','));
end

function v = ternary(c, a, b)
    if c
        v = a;
    else
        v = b;
    end
end
