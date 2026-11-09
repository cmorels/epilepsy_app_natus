function g = spec_git_info()
% SPEC_GIT_INFO  Code provenance stored in every output: HEAD hash of the
% repo this module lives in, and whether the working tree has
% uncommitted changes ('' / false when git is not available).

    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    g = struct('repo', repo, 'head', '', 'dirty', false);
    try
        [st1, h] = system(sprintf('git -C "%s" rev-parse HEAD', repo));
        [st2, d] = system(sprintf('git -C "%s" status --porcelain', repo));
        if st1 == 0
            g.head = strtrim(h);
        end
        if st2 == 0
            g.dirty = ~isempty(strtrim(d));
        end
    catch
    end
end
