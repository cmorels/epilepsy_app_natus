function [git, prov] = cs_provenance(cfg)
% CS_PROVENANCE  Code provenance recorded with every output: git HEAD and
% dirty flag, SHA-1 of every stats file (run_campaign_stats.m + src/stats,
% so uncommitted code is still identifiable) and of the IID code the
% re-detection runs (detect_iid.m and its dependencies), plus whether that
% IID code is unmodified with respect to HEAD.
%
%   git  : cs_git_info.m struct + iid_code (Nx2 {file, sha1}), iid_code_unmodified
%   prov : the compact struct stored in every re-detection output

    stats_files = [{fullfile(cfg.repo_root, 'run_campaign_stats.m')}; ...
        arrayfun(@(d) fullfile(d.folder, d.name), dir(fullfile(cfg.repo_root, 'src', 'stats', '**', '*.m')), 'UniformOutput', false)];   % incl. iid_sens/, plots/
    git = cs_git_info(cfg.repo_root, stats_files);
    iid_code = {'detect_iid.m', 'load_lfp_txt.m', 'pipeline_config.m', 'utils/parse_header.m', 'utils/mask_to_segments.m', ...
        'utils/filter_by_blocks.m', 'utils/rel_to_abs_time.m', 'utils/ensure_findpeaks_signal_toolbox.m', 'utils/decimate_minmax.m'};
    gi = cs_git_info(cfg.repo_root, fullfile(cfg.repo_root, 'src', iid_code));
    [~, st] = system(sprintf('git -C "%s" status --porcelain -- %s', cfg.repo_root, strjoin(strcat('"src/', iid_code, '"'), ' ')));
    git.iid_code = gi.files;
    git.iid_code_unmodified = isempty(strtrim(st));
    prov = struct('git_head', git.head, 'git_dirty', git.dirty, 'iid_code', {git.iid_code}, ...
        'iid_code_unmodified', git.iid_code_unmodified, 'stats_code', {git.files}, 'run_started', cfg.run_started);
end
