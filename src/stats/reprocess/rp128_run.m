function result = rp128_run()
% RP128_RUN  Phase 3: recording log override (only the 2 new EDFs of
% 128-s, same Port/channel mapping and Attenuation as the master row of
% A1 -> same case as the C2), provenance (git HEAD, diff, status, a copy of
% src/stats/reprocess), then run_pipeline_edf on the A1 and A2 EDFs with
% the C2 configuration and the C2 gain reference injected explicitly.
% The master EEG_recording_log.xlsx is only read.

    o = rp128_paths();
    addpath(fullfile(o.repo, 'src'), fullfile(o.repo, 'src', 'utils'));
    ts = char(datetime('now'), 'yyyyMMdd_HHmmss');
    diary(fullfile(o.out, sprintf('run_console_%s.log', ts)));
    cl = onCleanup(@() diary('off'));

    % ---- recording log override ----------------------------------------
    M = load_recording_log(o.master_log);
    m = M(strcmp(M.filename, o.stem) & strcmp(M.animal_id, o.subject), :);
    T = cell2table([{o.seg.stem}', {o.subject; o.subject}, {o.seg.port}', repmat(m.hpcr, 2, 1), repmat(m.hpcl, 2, 1), ...
        repmat(m.attenuation, 2, 1)], 'VariableNames', {'Filename', 'Animal ID', 'Port', 'HPCr_channel', 'HPCl_channel', 'Attenuation'});
    if isfile(o.log_override), delete(o.log_override); end
    writetable(T, o.log_override);
    L = load_recording_log(o.log_override);
    cfg_c2 = load(o.c2_cfg_mat, 'cfg');
    want = resolve_case_excel(m.attenuation{1}, cfg_c2.cfg);
    for i = 1:height(L)
        got = resolve_case_excel(L.attenuation{i}, cfg_c2.cfg);
        if ~strcmp(got.case_applied, want.case_applied)
            error('rp128_run:Case', '%s: case %s, the C2 used %s', L.filename{i}, got.case_applied, want.case_applied);
        end
    end
    fprintf('Override log %s:\n', o.log_override); disp(L);
    fprintf('Case of both EDFs: %s (as 128-s 20260729 in the C2)\n', want.case_applied);

    % ---- configuration ---------------------------------------------------
    cfg = rp128_cfg(o, o.log_override);
    d = rp_cfg_diff(cfg, cfg_c2.cfg, {'paths.output_root', 'edf.channels.log_file'});
    if ~isempty(d), error('rp128_run:Config', 'configuration differs from the C2 in: %s', strjoin(d, ', ')); end
    fprintf('Configuration identical to the C2''s config_used.mat of 128-s (except output_root and log_file)\n');
    fprintf('Gain reference injected (from %s): HPCr %.10g, HPCl %.10g, fallback %.10g uV\n', o.c2_ref_csv, ...
        cfg.quality.reference_sigma_uV('HPCr'), cfg.quality.reference_sigma_uV('HPCl'), cfg.quality.reference_sigma_fallback_uV);

    % ---- provenance ------------------------------------------------------
    git = @(a) evalc(sprintf('system(''git -C "%s" %s'');', o.repo, a));
    write_txt(fullfile(o.out, 'git_head.txt'), git('rev-parse HEAD'));
    write_txt(fullfile(o.out, 'git_diff.patch'), git('diff'));
    write_txt(fullfile(o.out, 'git_status.txt'), git('status --porcelain'));
    cd_ = fullfile(o.out, 'code_reprocess');
    if ~isfolder(cd_), mkdir(cd_); end
    copyfile(fullfile(o.repo, 'src', 'stats', 'reprocess', '*.m'), cd_);

    % ---- run ---------------------------------------------------------------
    fprintf('run_pipeline_edf: %s\n', strjoin({o.seg.edf}, ' | '));
    result = run_pipeline_edf({o.seg.edf}, cfg);
end

function write_txt(p, s)
    fid = fopen(p, 'w'); fprintf(fid, '%s', s); fclose(fid);
end
