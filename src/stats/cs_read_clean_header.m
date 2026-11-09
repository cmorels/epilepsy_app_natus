function h = cs_read_clean_header(file)
% CS_READ_CLEAN_HEADER  The '# key = value' header of a *_clean.txt as a
% struct of raw strings (stops at the first non-'#' line; the signal is
% never read). Keys keep their original names (valid field names here).

    fid = fopen(file, 'r');
    if fid < 0
        error('cs_read_clean_header:Open', 'Cannot open %s', file);
    end
    c = onCleanup(@() fclose(fid));
    h = struct();
    while true
        ln = fgetl(fid);
        if ~ischar(ln) || isempty(ln) || ln(1) ~= '#'
            break;
        end
        tok = regexp(ln, '^#\s*([^=]+?)\s*=\s*(.*)$', 'tokens', 'once');
        if isempty(tok)
            continue;
        end
        key = matlab.lang.makeValidName(tok{1});
        if ~isfield(h, key)  % first occurrence wins (quantization_step_uV is written twice)
            h.(key) = strtrim(tok{2});
        end
    end
end
