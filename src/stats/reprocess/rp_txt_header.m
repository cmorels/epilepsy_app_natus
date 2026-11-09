function h = rp_txt_header(p)
% RP_TXT_HEADER  '# key = value' header of a pipeline txt (01_txt / 02_clean) as a struct of strings.
    fid = fopen(p, 'r');
    if fid < 0, error('rp_txt_header:Open', 'Cannot open %s', p); end
    cl = onCleanup(@() fclose(fid));
    h = struct();
    while true
        l = fgetl(fid);
        if ~ischar(l) || ~startsWith(l, '#'), break; end
        kv = regexp(l, '^#\s*(\w+)\s*=\s*(.*)$', 'tokens', 'once');
        if ~isempty(kv), h.(kv{1}) = strtrim(kv{2}); end
    end
end
