function [header, rows] = cs_read_csv_raw(file)
% CS_READ_CSV_RAW  Header + body of a comma-separated file as raw text
% cells (no type guessing, duplicate/blank header names kept as-is). Used
% for the seizures_events.csv files the user re-saved with Excel, whose
% headers can repeat (user_check twice) or be blank.
%
%   header : 1 x n cell of char
%   rows   : m x n cell of char (short rows padded with '')

    txt = fileread(file);
    if ~isempty(txt) && txt(1) == char(65279)
        txt = txt(2:end);
    end
    lines = regexp(txt, '\r\n|\n|\r', 'split');
    while ~isempty(lines) && isempty(strtrim(lines{end}))
        lines(end) = [];
    end
    if isempty(lines)
        header = {};
        rows = cell(0, 0);
        return;
    end
    parsed = cellfun(@split_line, lines, 'UniformOutput', false);
    header = parsed{1};
    n = max(cellfun(@numel, parsed));
    header(end+1:n) = {''};
    rows = repmat({''}, numel(parsed) - 1, n);
    for i = 2:numel(parsed)
        f = parsed{i};
        rows(i-1, 1:numel(f)) = f;
    end
end

function f = split_line(s)
    if ~contains(s, '"')
        f = strsplit(s, ',', 'CollapseDelimiters', false);
        return;
    end
    f = {};
    cur = '';
    inq = false;
    k = 1;
    while k <= numel(s)
        c = s(k);
        if inq
            if c == '"'
                if k < numel(s) && s(k+1) == '"'
                    cur(end+1) = '"'; %#ok<AGROW>
                    k = k + 1;
                else
                    inq = false;
                end
            else
                cur(end+1) = c; %#ok<AGROW>
            end
        elseif c == '"'
            inq = true;
        elseif c == ','
            f{end+1} = cur; %#ok<AGROW>
            cur = '';
        else
            cur(end+1) = c; %#ok<AGROW>
        end
        k = k + 1;
    end
    f{end+1} = cur;
end
