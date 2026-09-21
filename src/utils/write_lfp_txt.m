function write_lfp_txt(file_path, header_pairs, signal)
% WRITE_LFP_TXT  Writer counterpart of load_lfp_txt.m / parse_header.m.
% Emits the shared "# key = value" header block followed by one sample
% per line ('%.6f', literal 'NaN' preserved at gaps), matching the format
% the original Intan-export scripts used.
%
% INPUTS
%   file_path    : output .txt path
%   header_pairs : Nx2 cell array {key, value}; value is char/string
%                  (printed as-is) or a numeric/logical scalar (printed
%                  with '%.10g' / 'true'|'false').
%   signal       : numeric vector, one sample per line.

    if ~(iscell(header_pairs) && size(header_pairs, 2) == 2)
        error('write_lfp_txt:BadHeader', 'header_pairs must be an Nx2 cell array of {key, value}.');
    end
    if ~isnumeric(signal)
        error('write_lfp_txt:BadSignal', 'signal must be numeric.');
    end

    fid = fopen(file_path, 'w');
    if fid == -1
        error('write_lfp_txt:FileOpenError', 'Could not open file for writing: %s', file_path);
    end

    cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    for i = 1:size(header_pairs, 1)
        key = header_pairs{i, 1};
        fprintf(fid, '# %s = %s\n', key, format_value(header_pairs{i, 2}, key));
    end

    fprintf(fid, '%.6f\n', signal(:));
end

function s = format_value(v, key)
    if ischar(v) || isstring(v)
        s = char(v);
    elseif islogical(v) && isscalar(v)
        if v
            s = 'true';
        else
            s = 'false';
        end
    elseif isnumeric(v) && isscalar(v)
        s = sprintf('%.17g', v);   % 17 significant digits round-trips any double exactly (needed for session_start_unix)
    else
        error('write_lfp_txt:BadValue', 'Unsupported header value type for key "%s".', key);
    end
end
