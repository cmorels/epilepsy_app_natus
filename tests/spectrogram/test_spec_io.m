function tests = test_spec_io
% TEST_SPEC_IO  Files, selection, dry-run, duplicates, outputs, on a
% synthetic folder tree (run: runtests('tests/spectrogram')).
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(repo, 'src'), fullfile(repo, 'src', 'utils'), fullfile(repo, 'src', 'spectrogram'), ...
        fileparts(mfilename('fullpath')));
    rng(2);
    root = tempname;
    mkdir(root);
    fs = 512;
    n = 300 * fs;
    t = (0:n - 1)' / fs;
    tz = 'Europe/Paris';
    ssA = datetime(2026, 7, 1, 12, 0, 0, 'TimeZone', tz);
    ssB = datetime(2026, 7, 2, 9, 30, 15.5, 'TimeZone', tz);
    sz_win = t >= 120 & t < 150;                         % one "seizure" at 120-150 s in session A
    for reg = {'HPCr', 'HPCl'}
        x = 50 * randn(n, 1) + 80 * sin(2 * pi * 6 * t);
        x(sz_win) = x(sz_win) + 400 * sin(2 * pi * 12 * t(sz_win));
        spikes = 40:3:100;
        for p = spikes
            x = x - 500 * exp(-(t - p) .^ 2 / (2 * 0.012 ^ 2));
        end
        % session A: one folder per EDF + 01_txt raw (gain 2)
        make_spec_session(fullfile(root, '097', 'EDFA'), '097', ssA, reg{1}, x, fs, 'source_file', 'EDFA.EDF', ...
            'raw_folder', fullfile(root, '097', '01_txt'), 'gain', 2);
        make_spec_session(fullfile(root, '097', 'EDFB'), '097', ssB, reg{1}, x(1:60 * fs), fs, 'source_file', 'EDFB.EDF');
        % a different animal whose id contains '097': never matched for '097'
        make_spec_session(fullfile(root, '097-s', 'EDFA'), '097-s', ssA, reg{1}, x(1:10 * fs), fs, 'source_file', 'EDFA.EDF');
    end
    % minimal IID_results.mat for session A, HPCr
    iid_results = struct();
    iid_results.individual_peaks_table = table((40:3:100)', 0.5 * ones(21, 1), 10 * ones(21, 1), 0.4 * ones(21, 1), ...
        ones(21, 1), 'VariableNames', {'Time_s', 'Amplitude_mV', 'Width_samples', 'Prominence_mV', 'block_id'});
    iid_results.exclusion_table = table(1, "seizure", 1, 115, 155, 40, ...
        'VariableNames', {'zone_id', 'reason', 'source_id', 'start_s', 'end_s', 'duration_s'});
    mkdir(fullfile(root, '097', '04_iid', 'EDFA'));
    save(fullfile(root, '097', '04_iid', 'EDFA', '097_20260701_120000_HPCr_IID_results.mat'), 'iid_results');

    % seizures_events: (1) pipeline-shaped CSV, (2) reviewed copy (Excel-broken datetimes + labels)
    hdr = 'subject_id,region,session_start,source_file,seizure_id,start_datetime,end_datetime,start_s,end_s,duration_s,start_abs,end_abs,block_id,adjacent_to_gap,seizure_mode,over_max_duration,ll_ratio,peak_energy_ratio,hf_ratio_db,envelope_cv';
    fmt = 'dd-MMM-yyyy HH:mm:ss.SSS';
    rowsP = {};
    rowsR = {};
    for reg = {'HPCr', 'HPCl'}
        sa = ssA + seconds(120); sa.Format = fmt;
        ea = ssA + seconds(150); ea.Format = fmt;
        ssf = ssA; ssf.Format = fmt;
        rowsP{end+1} = sprintf('097,%s,%s,EDFA.EDF,1,01/07/2026 12:02:00,01/07/2026 12:02:30,120,150,30,%s,%s,1,0,robust,0,3,100,-18,1', ...
            reg{1}, char(ssf), char(sa), char(ea)); %#ok<AGROW>
        rowsR{end+1} = sprintf('097,%s,00:00.0,EDFA.EDF,1,01/07/2026 12:02:00,01/07/2026 12:02:30,120,150,30,02:00.0,02:30.0,1,0,robust,0,3,100,-18,1,TP,ok', reg{1}); %#ok<AGROW>
        rowsR{end+1} = sprintf('097,%s,00:00.0,EDFA.EDF,2,01/07/2026 12:04:00,01/07/2026 12:04:10,240,250,10,04:00.0,04:10.0,1,0,robust,0,3,100,-18,1,FP,grooming', reg{1}); %#ok<AGROW>
    end
    write_lines(fullfile(root, 'seizures_pipeline.csv'), [{hdr}, rowsP]);
    write_lines(fullfile(root, 'seizures_reviewed.csv'), [{[hdr ',user_check,comments,']}, rowsR]);

    tc.TestData.root = root;
    tc.TestData.fs = fs;
    tc.TestData.cfg = base_cfg(root);
end

function teardownOnce(tc)
    spec_resolve_files('reset');
    spec_signal_cache('reset');
    close all force
    rmdir(tc.TestData.root, 's');
end

function cfg = base_cfg(root)
    cfg = spectrogram_config();
    cfg.spec.pipeline_root = root;
    cfg.spec.output_root = fullfile(root, 'out');
    cfg.spec.seizure.pad_pre_s = 30;
    cfg.spec.seizure.pad_post_s = 30;
end

function write_lines(f, lines)
    fid = fopen(f, 'w');
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);
end

%% ---- file resolution -------------------------------------------------------
function testResolveByKeyAndByEdfName(tc)
    spec_resolve_files('reset');
    a = spec_resolve_files(tc.TestData.root, '097', '20260701_120000', 'HPCr', 'clean', 'Europe/Paris');
    b = spec_resolve_files(tc.TestData.root, '097', 'EDFA.EDF', 'HPCr', 'clean', 'Europe/Paris');
    c = spec_resolve_files(tc.TestData.root, '097', 'EDFA', 'HPCr', 'clean', 'Europe/Paris');
    tc.verifyEqual(a.clean, b.clean);
    tc.verifyEqual(a.clean, c.clean);
    tc.verifyTrue(contains(a.clean, [filesep '097' filesep 'EDFA' filesep]));
    tc.verifyNotEmpty(a.raw);
    tc.verifyNotEmpty(a.iid_mat);
    tc.verifyEmpty(spec_resolve_files(tc.TestData.root, '097', '20260701_120000', 'HPCl', 'clean', 'Europe/Paris').iid_mat);
    tc.verifyEqual(a.session_key, '20260701_120000');
end

function testMissingSessionAndSubjectErrors(tc)
    spec_resolve_files('reset');
    tc.verifyError(@() spec_resolve_files(tc.TestData.root, '097', '20260709_000000', 'HPCr', 'clean', 'Europe/Paris'), ...
        'spec_resolve_files:NoSession');
    tc.verifyError(@() spec_resolve_files(tc.TestData.root, '097', 'EDFZ.EDF', 'HPCr', 'clean', 'Europe/Paris'), ...
        'spec_resolve_files:NoSession');
    tc.verifyError(@() spec_resolve_files(tc.TestData.root, '97', '20260701_120000', 'HPCr', 'clean', 'Europe/Paris'), ...
        'spec_resolve_files:NoSubject');
    tc.verifyError(@() spec_resolve_files(tc.TestData.root, '097', '20260702_093015', 'HPCr', 'raw', 'Europe/Paris'), ...
        'spec_resolve_files:NoRaw');
end

function testAmbiguousFileErrors(tc)
    root = tc.TestData.root;
    src = fullfile(root, '097', 'EDFB', '097_20260702_093015_HPCr_clean.txt');
    dupdir = fullfile(root, '097', 'copy');
    mkdir(dupdir);
    copyfile(src, dupdir);
    cleanup = onCleanup(@() rmdir(dupdir, 's'));
    spec_resolve_files('reset');
    tc.verifyError(@() spec_resolve_files(root, '097', '20260702_093015', 'HPCr', 'clean', 'Europe/Paris'), ...
        'spec_resolve_files:Ambiguous');
    clear cleanup
    spec_resolve_files('reset');
end

function testCacheEqualsLoadLfpTxt(tc)
    spec_resolve_files('reset');
    spec_signal_cache('reset');
    r = spec_resolve_files(tc.TestData.root, '097', '20260701_120000', 'HPCr', 'clean', 'Europe/Paris');
    cdir = fullfile(tc.TestData.root, 'cache_test');
    info = spec_signal_cache('open', r.clean, cdir, 10000);   % small chunks: several textscan passes
    d = load_lfp_txt(r.clean);
    tc.verifyEqual(info.n, numel(d.signal));
    tc.verifyTrue(isequaln(spec_signal_cache('read', info, 1, info.n), d.signal));
    tc.verifyTrue(isequaln(spec_signal_cache('read', info, 1000, 2000), d.signal(1000:2000)));
    x = spec_signal_cache('read', info, -9, 10);
    tc.verifyTrue(all(isnan(x(1:10))) && isequal(x(11:20), d.signal(1:10)));
    tc.verifyEqual(info.session_start, d.session_start);
    spec_signal_cache('reset');
    info2 = spec_signal_cache('open', r.clean, cdir, 10000);
    tc.verifyTrue(info2.reused, 'second open reuses the binary cache');
    spec_signal_cache('reset');
end

function testRawSourceAppliesGain(tc)
    spec_resolve_files('reset');
    spec_signal_cache('reset');
    cfg = tc.TestData.cfg;
    rc = spec_resolve_files(tc.TestData.root, '097', '20260701_120000', 'HPCr', 'clean', 'Europe/Paris');
    rr = spec_resolve_files(tc.TestData.root, '097', '20260701_120000', 'HPCr', 'raw', 'Europe/Paris');
    cd_ = fullfile(tc.TestData.root, 'cache_raw');
    cfg.spec.signal_source = 'clean';
    a = spec_load_segment(rc, 10, 20, cfg, cd_);
    cfg.spec.signal_source = 'raw';
    b = spec_load_segment(rr, 10, 20, cfg, cd_);
    tc.verifyEqual(b.gain_applied, 2);
    tc.verifyEqual(b.x, a.x, 'AbsTol', 2e-6, 'raw x gain == clean when nothing was interpolated / notched');
    tc.verifyEqual(a.n_outliers_segment, 0);
    spec_signal_cache('reset');
end

%% ---- seizures_events -----------------------------------------------------------
function testReadSeizuresTables(tc)
    [P, ~] = spec_read_seizures_table(fullfile(tc.TestData.root, 'seizures_pipeline.csv'), 'Europe/Paris');
    tc.verifyEqual(P.subject_id{1}, '097', 'subject_id stays text');
    tc.verifyFalse(any(isnat(P.start_abs)), 'pipeline CSV: typed datetimes');
    tc.verifyTrue(all(cellfun(@isempty, P.user_check)));
    [R, notes] = spec_read_seizures_table(fullfile(tc.TestData.root, 'seizures_reviewed.csv'), 'Europe/Paris');
    tc.verifyEqual(R.subject_id{1}, '097');
    tc.verifyEqual(R.user_check, {'TP'; 'FP'; 'TP'; 'FP'});
    tc.verifyEqual(R.comments{2}, 'grooming');
    tc.verifyTrue(all(isnat(R.start_abs)), 'Excel-broken start_abs -> NaT');
    tc.verifyTrue(any(contains(notes, 'read_pipeline_csv could not read it')));
end

function testSelectByIdAndFilter(tc)
    cfg = tc.TestData.cfg;
    spec_resolve_files('reset');
    sel = struct('seizures_file', fullfile(tc.TestData.root, 'seizures_reviewed.csv'), 'region', 'both');
    sel.ids = {'097', 'EDFA.EDF', [1 2]};
    ev = spec_select_seizures(sel, cfg);
    tc.verifyEqual(numel(ev), 2);
    tc.verifyEqual(ev(1).region, 'both');
    tc.verifyEqual({ev(1).sz.region}, {'HPCr', 'HPCl'});
    tc.verifyEqual([ev(1).inicio ev(1).fin], [120 150]);
    sel.ids = {'097', '20260701_120000', 1};
    tc.verifyEqual(numel(spec_select_seizures(sel, cfg)), 1);
    sel = rmfield(sel, 'ids');
    sel.filter = struct('user_check', 'tp');
    ev = spec_select_seizures(sel, cfg);
    tc.verifyEqual(numel(ev), 1);
    tc.verifyEqual(ev(1).sz(1).user_check, 'TP');
    sel.region = 'HPCl';
    ev = spec_select_seizures(sel, cfg);
    tc.verifyEqual({ev.region}, {'HPCl'});
    sel.region = 'both';
    sel.filter = struct('duration_s', [0 20]);
    ev = spec_select_seizures(sel, cfg);
    tc.verifyEqual(ev(1).id, 'sz002');
end

%% ---- dry-run, duplicates, full run, redraw ----------------------------------------
function testDryRunAndDuplicates(tc)
    cfg = tc.TestData.cfg;
    cfg.spec.dry_run = true;
    cfg.spec.output_root = fullfile(tc.TestData.root, 'out_dry');
    manual = {
        'a1', 'seizure', '097', '20260701_120000', 'both', '01/07/2026 12:02:00.000', '01/07/2026 12:02:30.000', 'natus'
        'a2', 'seizure', '097', 'EDFA.EDF', 'HPCr', 120, 150, 't_rel'
        'a3', 'iid', '097', '20260709_000000', 'HPCr', 1, 2, 't_rel'
        };
    out = spec_run(cfg, spec_read_event_list(manual));
    P = out.preview;
    tc.verifyEqual(height(P), 4);
    tc.verifyEqual(P.status(1:2), {'ok'; 'ok'});
    tc.verifyTrue(startsWith(P.status{3}, 'duplicate'), 'a2 HPCr = a1 HPCr');
    tc.verifyTrue(startsWith(P.status{4}, 'error'));
    tc.verifyEqual(P.start_t_rel(1), 120, 'AbsTol', 1e-6);
    tc.verifyEqual(P.nan_pct(1), 0);
    tc.verifyTrue(isfile(fullfile(cfg.spec.output_root, '06_spectrograms_summary', '_selection_preview.csv')));
    tc.verifyFalse(isfolder(fullfile(cfg.spec.output_root, '097')), 'dry run computes / writes no event output');
end

function testFullRunAndRedraw(tc)
    cfg = tc.TestData.cfg;
    cfg.spec.dry_run = false;
    cfg.spec.output_root = fullfile(tc.TestData.root, 'out_full');
    cfg.spec.seizures_file = fullfile(tc.TestData.root, 'seizures_reviewed.csv');
    sel = struct('seizures_file', cfg.spec.seizures_file, 'region', 'both');
    sel.ids = {'097', 'EDFA.EDF', 1};
    ev = spec_select_seizures(sel, cfg);
    manual = {'iid01', 'iid', '097', '20260701_120000', 'HPCr', 35, 105, 't_rel'};
    ev = [ev, spec_read_event_list(manual)];
    out = spec_run(cfg, ev);
    tc.verifyEqual({out.results.status}, {'ok', 'ok'});
    m1 = out.results(1).paths.mat;
    tc.verifyTrue(isfile(m1) && isfile(out.results(1).paths.png) && isfile(out.results(1).paths.fig));
    tc.verifyTrue(contains(m1, fullfile('097', '06_spectrograms', '097_20260701_120000', 'sz001_both.mat')));
    S = load(m1);
    tc.verifyEqual(numel(S.res.ch), 2);
    tc.verifyEqual(S.res.ch(1).sz.user_check, 'TP');
    f = S.res.features;
    tc.verifyEqual(f.peak_freq_rel_hz(1), 12, 'AbsTol', 1, 'ictal 12 Hz rhythm re baseline');
    tc.verifyEqual(f.peak_freq_osc_hz(1), 12, 'AbsTol', 1, 'ictal 12 Hz rhythm above the 1/f background');
    tc.verifyEqual(f.domfreq_mean_hz(1), 12, 'AbsTol', 1, 'dominant frequency (aperiodic) follows the 12 Hz rhythm');
    tc.verifyGreaterThan(f.rel_power_theta_db(1), 6);
    fig = spec_plot_seizure(S.res);                      % redraw from the .mat only
    tc.verifyTrue(isgraphics(fig));
    close(fig);
    S2 = load(out.results(2).paths.mat);
    tc.verifyEqual(S2.res.features.n_detected, 21);
    tc.verifyEqual(S2.res.features.n_used, 21, 'no peak in the 115-155 s exclusion zone');
    tc.verifyEqual(S2.res.peaks.polarity, 'neg');
    fig = spec_plot_iid(S2.res);
    tc.verifyTrue(isgraphics(fig));
    close(fig);
    sumdir = fullfile(cfg.spec.output_root, '06_spectrograms_summary');
    tc.verifyTrue(isfile(fullfile(sumdir, 'spec_features_seizures.csv')) && isfile(fullfile(sumdir, 'spec_features_iid.csv')));
    tc.verifyTrue(isfile(fullfile(sumdir, 'spec_run_log.txt')));
    % nothing written outside 06_*
    d = dir(fullfile(cfg.spec.output_root, '097'));
    tc.verifyEqual(sort(setdiff({d.name}, {'.', '..'})), {'06_spectrograms'});
end
