function h = spec_peek_header(txt_path)
% SPEC_PEEK_HEADER  Header of a pipeline txt WITHOUT reading its data
% section: the header lines plus one dummy sample are copied to a temp
% file and parsed by load_lfp_txt.m, so fs / session_start / meta are
% exactly what the pipeline itself would get.
%
%   h = spec_peek_header(txt_path)
%   h: fs, session_start (zoned datetime or NaT), meta (struct, text values),
%      header_pairs, file, folder, n_lines

    txt_path = char(txt_path);
    fid = fopen(txt_path, 'rt');
    if fid < 0
        error('spec_peek_header:FileOpenError', 'Could not open %s', txt_path);
    end
    hdr = parse_header(fid);
    fseek(fid, 0, 'bof');
    lines = cell(hdr.n_lines, 1);
    for k = 1:hdr.n_lines
        lines{k} = fgetl(fid);
    end
    fclose(fid);

    tmp = [tempname '.txt'];
    ft = fopen(tmp, 'w');
    fprintf(ft, '%s\n', lines{:});
    fprintf(ft, '0\n');
    fclose(ft);
    cleaner = onCleanup(@() delete(tmp));
    d = load_lfp_txt(tmp);

    [folder, name, ext] = fileparts(txt_path);
    h = struct('fs', d.fs, 'session_start', d.session_start, 'meta', d.meta, ...
        'header_pairs', {d.header_pairs}, 'file', [name ext], 'folder', folder, 'n_lines', hdr.n_lines);
end
