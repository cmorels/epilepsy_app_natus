function cs = cs_code_state(cfg)
% CS_CODE_STATE  Records the exact code of an uncommitted run in
% cfg.out_dir: git_head.txt (git rev-parse HEAD) and code_state.patch =
% `git diff HEAD` (tracked changes) followed by every untracked file of
% src/stats/ and run_campaign_stats.m as a new-file diff (git diff
% --no-index /dev/null <file>). Applying the patch on HEAD rebuilds the code.

    r = cfg.repo_root;
    [~, head] = system(sprintf('git -C "%s" rev-parse HEAD', r));
    head = strtrim(head);
    [~, tracked] = system(sprintf('git -C "%s" diff HEAD 2>NUL', r));
    [~, un] = system(sprintf('git -C "%s" ls-files --others --exclude-standard -- src/stats run_campaign_stats.m', r));
    un = strtrim(splitlines(strtrim(un)));
    un = un(~cellfun(@isempty, un));
    parts = {tracked};
    for i = 1:numel(un)
        [~, d] = system(sprintf('git -C "%s" diff --no-index -- /dev/null "%s" 2>NUL', r, un{i}));
        d = regexprep(d, '^warning:[^\n]*\n', '', 'lineanchors');
        parts{end+1} = d; %#ok<AGROW>
    end
    cs.head = head;
    cs.untracked = un;
    cs.patch_file = fullfile(cfg.out_dir, 'code_state.patch');
    fid = fopen(cs.patch_file, 'w');
    fprintf(fid, '%s', strjoin(parts, ''));
    fclose(fid);
    fid = fopen(fullfile(cfg.out_dir, 'git_head.txt'), 'w');
    fprintf(fid, '%s\n', head);
    fclose(fid);
end
