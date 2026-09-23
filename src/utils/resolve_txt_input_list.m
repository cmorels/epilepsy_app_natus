function files = resolve_txt_input_list(input, recursive)
% RESOLVE_TXT_INPUT_LIST  Flexible input resolution shared by
% estimate_reference_sigma.m and make_cases_template.m: a folder, a folder
% searched recursively, a cell array of folders, or a cell array of
% specific .txt paths -- because EDF/txt exports are typically organized
% one folder per animal, and a cohort-wide tool needs to see all of them
% at once.
%
%   files = resolve_txt_input_list(input, recursive)
%
% Returns a column cellstr of absolute .txt paths, deduplicated.

    if iscell(input)
        files = cell(0, 1);
        for i = 1:numel(input)
            files = [files; resolve_one(input{i}, recursive)]; %#ok<AGROW>
        end
    else
        files = resolve_one(input, recursive);
    end

    files = unique(files, 'stable');
end

function files = resolve_one(item, recursive)
    if ~(ischar(item) || isstring(item))
        error('resolve_txt_input_list:BadInput', ...
            'Each input must be a folder path or a .txt file path (char/string).');
    end
    item = char(item);

    if isfolder(item)
        if recursive
            d = dir(fullfile(item, '**', '*.txt'));
        else
            d = dir(fullfile(item, '*.txt'));
        end
        d = d(~[d.isdir]);
        files = fullfile({d.folder}, {d.name})';
    elseif exist(item, 'file') == 2
        files = {item};
    else
        warning('resolve_txt_input_list:NotFound', 'Input not found, skipping: %s', item);
        files = cell(0, 1);
    end
end
