function T = rename_seizure_figures_by_status(output_root, apply)
% RENAME_SEIZURE_FIGURES_BY_STATUS  Put each reconciled seizure zoom
% figure's status in its file name, for outputs written before
% save_bilateral_seizure_figures.m did it itself (2026-10-04):
%
%   {base}_seizure03.png  ->  {base}_imputed_seizure03.png      (and .fig)
%   e.g. 005-s_20260329_171443_HPCl_seizure03.png -> 005-s_20260329_171443_HPCl_imputed_seizure03.png
%
%   T = rename_seizure_figures_by_status(output_root)          % preview only, nothing renamed
%   T = rename_seizure_figures_by_status(output_root, true)    % rename
%
% output_root : a campaign output (one subfolder per animal, each with
%               03_seizures/) or a single run_pipeline_edf.m output_root.
%
% The status of each zoom comes from {base}_seizures.mat in the same
% folder (seizure_results.seizures: seizure_id, detection_status = accepted
% / rescued / imputed), never from the file name or dates. One zoom is one
% (seizure_id, channel); if a channel has several rows for the same id
% (fragmented) they are all 'accepted'. Panoramas ({base}_seizures.*) are
% not touched: they hold every status.
%
% Safe to run twice: a zoom that already carries its status is skipped. A
% zoom listed in the .mat but missing on disk (e.g. a figure that could not
% be saved) is reported, and an old-style zoom with no .mat entry is
% reported as 'no status found' and left alone. With apply = true the full
% list is written to <output_root>/renamed_files.csv (old and new name of
% every file), which is what you need to undo it.
%
% OUTPUT T: table folder, old_name, new_name, status, action
% ('renamed' | 'to rename' | 'already renamed' | 'missing on disk' |
% 'no status found').

    if nargin < 2
        apply = false;
    end
    statuses = {'accepted', 'rescued', 'imputed'};

    dirs = seizure_dirs(output_root);
    rows = cell(0, 5);
    for i = 1:numel(dirs)
        d = dirs{i};
        handled = {};
        mats = dir(fullfile(d, '*_seizures.mat'));
        for k = 1:numel(mats)
            base = regexprep(mats(k).name, '_seizures\.mat$', '');
            try
                S = load(fullfile(d, mats(k).name), 'seizure_results');
                sz = S.seizure_results.seizures;
            catch ME
                warning('rename_seizure_figures_by_status:BadMat', 'Skipping %s: %s', fullfile(d, mats(k).name), ME.message);
                continue;
            end
            if ~all(ismember({'seizure_id', 'detection_status'}, sz.Properties.VariableNames))
                continue;  % not a reconciled run: its figures carry no status
            end
            for sid = unique(sz.seizure_id)'
                status = sz.detection_status{find(sz.seizure_id == sid, 1)};
                tag = regexprep(status, '[^\w\-]', '_');
                for ext = {'png', 'fig'}
                    old_name = sprintf('%s_seizure%02d.%s', base, sid, ext{1});
                    new_name = sprintf('%s_%s_seizure%02d.%s', base, tag, sid, ext{1});
                    handled{end+1} = old_name; %#ok<AGROW>
                    if exist(fullfile(d, new_name), 'file') == 2
                        action = 'already renamed';
                    elseif exist(fullfile(d, old_name), 'file') == 2
                        if apply
                            [ok, msg] = movefile(fullfile(d, old_name), fullfile(d, new_name));
                            if ok
                                action = 'renamed';
                            else
                                action = ['rename failed: ' msg];
                            end
                        else
                            action = 'to rename';
                        end
                    else
                        action = 'missing on disk';
                    end
                    rows(end+1, :) = {d, old_name, new_name, status, action}; %#ok<AGROW>
                end
            end
        end

        % old-style zooms that no .mat explains
        status_re = strjoin(statuses, '|');
        files = dir(fullfile(d, '*_seizure*.*'));
        for k = 1:numel(files)
            n = files(k).name;
            is_zoom = ~isempty(regexp(n, '_seizure\d+\.(png|fig)$', 'once'));
            has_status = ~isempty(regexp(n, ['_(' status_re ')_seizure\d+\.(png|fig)$'], 'once'));
            if is_zoom && ~has_status && ~ismember(n, handled)
                rows(end+1, :) = {d, n, '', '', 'no status found'}; %#ok<AGROW>
            end
        end
    end

    T = cell2table(rows, 'VariableNames', {'folder', 'old_name', 'new_name', 'status', 'action'});

    [acts, ~, ia] = unique(T.action);
    fprintf('rename_seizure_figures_by_status: %d folder(s) with 03_seizures, %d file(s) (%s)\n', ...
        numel(dirs), height(T), ternary(apply, 'APPLIED', 'preview only, nothing renamed'));
    for a = 1:numel(acts)
        fprintf('  %-16s %d\n', acts{a}, nnz(ia == a));
    end
    known = T.status(~cellfun(@isempty, T.status));
    if ~isempty(known)
        [sts, ~, is] = unique(known);
        c = [sts'; num2cell(accumarray(is, 1)')];
        fprintf('  by status (files):'); fprintf(' %s=%d', c{:}); fprintf('\n');
    end

    if apply && height(T) > 0
        out_csv = fullfile(output_root, 'renamed_files.csv');
        writetable(T, out_csv);
        fprintf('  list written to %s\n', out_csv);
    end
end

function dirs = seizure_dirs(output_root)
    if isfolder(fullfile(output_root, '03_seizures'))
        dirs = {fullfile(output_root, '03_seizures')};
        return;
    end
    d = dir(output_root);
    d = d([d.isdir] & ~ismember({d.name}, {'.', '..'}));
    dirs = {};
    for k = 1:numel(d)
        p = fullfile(output_root, d(k).name, '03_seizures');
        if isfolder(p)
            dirs{end+1} = p; %#ok<AGROW>
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
