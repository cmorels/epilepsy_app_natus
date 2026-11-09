function R = spec_read_text_csv(file)
% SPEC_READ_TEXT_CSV  A comma-separated file read ENTIRELY as text, with
% line 1 as the header -- never type-guessed, never header-guessed (a
% reviewed seizures_events.csv may have a trailing comma, rows shorter
% than the header, duplicate or empty column names).
%
%   R = spec_read_text_csv(file)
%
% Column names are the header cells as written (trimmed); an empty name
% becomes Var<k>, a repeated one gets a suffix (_1, _2, ...). Every column
% is a cellstr; missing trailing cells are ''.

    fid = fopen(file, 'r', 'n', 'UTF-8');
    if fid < 0
        error('spec_read_text_csv:FileOpenError', 'Could not open %s', file);
    end
    first = fgetl(fid);
    fclose(fid);
    if ~ischar(first)
        error('spec_read_text_csv:Empty', '%s is empty.', file);
    end
    first = regexprep(first, '^\xEF\xBB\xBF|^\x{FEFF}', '');   % UTF-8 BOM
    names = strtrim(strsplit(first, ','));
    for k = 1:numel(names)
        names{k} = regexprep(names{k}, '^"(.*)"$', '$1');
        if isempty(names{k})
            names{k} = sprintf('Var%d', k);
        end
    end
    names = matlab.lang.makeUniqueStrings(names);
    n = numel(names);
    opts = delimitedTextImportOptions('NumVariables', n, 'VariableNames', names, 'Delimiter', ',', ...
        'DataLines', [2 Inf], 'VariableTypes', repmat({'char'}, 1, n), 'ExtraColumnsRule', 'ignore', ...
        'EmptyLineRule', 'skip', 'Encoding', 'UTF-8');
    opts.VariableNamingRule = 'preserve';
    opts = setvaropts(opts, names, 'WhitespaceRule', 'trim', 'FillValue', '');
    R = readtable(file, opts);
end
