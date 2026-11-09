function g = cs_git_info(repo_root, files)
% CS_GIT_INFO  Code provenance for every output: HEAD hash, whether the
% working tree differs from it (porcelain status), and the SHA-1 of each
% given file as it was when the run started (uncommitted code is thus
% still identifiable).

    [~, h] = system(sprintf('git -C "%s" rev-parse HEAD', repo_root));
    [~, st] = system(sprintf('git -C "%s" status --porcelain', repo_root));
    g.head = strtrim(h);
    g.dirty = ~isempty(strtrim(st));
    g.status = strtrim(st);
    g.files = cell(numel(files), 2);
    md = java.security.MessageDigest.getInstance('SHA-1');
    for i = 1:numel(files)
        fid = fopen(files{i}, 'r');
        b = fread(fid, Inf, '*uint8');
        fclose(fid);
        md.reset();
        d = typecast(md.digest(b), 'uint8');
        [~, nm, ex] = fileparts(files{i});
        g.files(i, :) = {[nm ex], lower(reshape(dec2hex(d)', 1, []))};
    end
end
