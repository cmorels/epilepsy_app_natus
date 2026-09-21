function hdr = parse_header(fid)
% PARSE_HEADER  Generic reader for the "# key = value" header block used by
% every txt file in this pipeline (EDF-origin and legacy Intan-origin
% alike). Consumes consecutive leading lines starting with '#' from the
% already-open file id FID, and leaves FID positioned at the first
% non-header line (ready for textscan of the numeric data).
%
% OUTPUT (struct hdr):
%   .fields  : struct with one dynamic field per header key. Every value
%              is kept as trimmed char, deliberately NOT auto-typed: an
%              identifier that happens to look numeric (e.g. a zero-padded
%              subject_id "097") must not silently become the number 97.
%              Callers that need a number (fs, n_samples, ...) convert
%              explicitly with str2double and validate the result -- see
%              load_lfp_txt.m for the pattern.
%   .pairs   : Nx2 cell array {original_key, raw_value_string}, in file
%              order, for lossless passthrough when a later stage needs to
%              re-emit the same header (e.g. clean_lfp.m).
%   .n_lines : number of header lines consumed.

    fields = struct();
    pairs = cell(0, 2);
    n_lines = 0;

    pos = ftell(fid);
    line = fgetl(fid);

    while ischar(line) && startsWith(strtrim(line), '#')
        n_lines = n_lines + 1;

        tokens = regexp(strtrim(line), '^#\s*([^=]+?)\s*=\s*(.*)$', 'tokens', 'once');
        if ~isempty(tokens)
            key_raw = tokens{1};
            value_raw = strtrim(tokens{2});

            pairs(end+1, :) = {key_raw, value_raw}; %#ok<AGROW>
            fields.(sanitize_field_name(key_raw)) = value_raw;
        end

        pos = ftell(fid);
        line = fgetl(fid);
    end

    % Rewind to just after the last consumed header line (fgetl already
    % advanced past a non-header line, so undo that one step).
    fseek(fid, pos, 'bof');

    hdr = struct('fields', fields, 'pairs', {pairs}, 'n_lines', n_lines);
end

function name = sanitize_field_name(key_raw)
    name = regexprep(strtrim(key_raw), '[^a-zA-Z0-9_]', '_');
    if isempty(name) || ~isstrprop(name(1), 'alpha')
        name = ['f_' name];
    end
end

