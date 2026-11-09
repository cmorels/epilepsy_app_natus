function tests = test_spec_time
% TEST_SPEC_TIME  Natus clock vs pipeline clock across the 29/03/2026
% daylight-saving change (run: runtests('tests/spectrogram')).
    tests = functiontests(localfunctions);
end

function setupOnce(tc)
    repo = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(repo, 'src'), fullfile(repo, 'src', 'utils'), fullfile(repo, 'src', 'spectrogram'), ...
        fileparts(mfilename('fullpath')));
    tc.TestData.root = tempname;
    mkdir(tc.TestData.root);
    % session starting 27/03/2026 10:00:00 (CET, +01:00); the header claims 4 days
    % of data, the file holds 10 samples (resolution only reads headers)
    ss = datetime(2026, 3, 27, 10, 0, 0, 'TimeZone', 'Europe/Paris');
    make_spec_session(fullfile(tc.TestData.root, '016-s', 'EDF1'), '016-s', ss, 'HPCr', zeros(10, 1), 512, ...
        'source_file', 'EDF1.EDF', 'n_samples_header', 4 * 86400 * 512);
    spec_resolve_files('reset');
    tc.TestData.r = spec_resolve_files(tc.TestData.root, '016-s', '20260327_100000', 'HPCr', 'clean', 'Europe/Paris');
end

function teardownOnce(tc)
    spec_resolve_files('reset');
    rmdir(tc.TestData.root, 's');
end

function testWallStart(tc)
    r = tc.TestData.r;
    tc.verifyEqual(spec_resolve_times('text', r.wall_start), '27/03/2026 10:00:00.000');
end

function testBeforeChangeClocksAgree(tc)
    r = tc.TestData.r;
    [n, p] = spec_resolve_times('clocks', 3600.25, r);
    tc.verifyEqual(spec_resolve_times('text', n), '27/03/2026 11:00:00.250');
    tc.verifyEqual(spec_resolve_times('text', p), '27/03/2026 11:00:00.250');
end

function testAfterChangePipelineIsOneHourAhead(tc)
    r = tc.TestData.r;
    t = 3 * 86400 + 600.5;                    % 30/03/2026, after the spring change
    [n, p] = spec_resolve_times('clocks', t, r);
    tc.verifyEqual(spec_resolve_times('text', n), '30/03/2026 10:10:00.500');
    tc.verifyEqual(spec_resolve_times('text', p), '30/03/2026 11:10:00.500');
    wall_p = datetime(p.Year, p.Month, p.Day, p.Hour, p.Minute, p.Second);
    tc.verifyEqual(wall_p - n, hours(1));
end

function testNatusTextToTRel(tc)
    r = tc.TestData.r;
    t = spec_resolve_times('to_t_rel', '30/03/2026 10:10:00.500', 'natus', r);
    tc.verifyEqual(t, 3 * 86400 + 600.5, 'AbsTol', 1e-6);
end

function testPipelineTextToTRel(tc)
    r = tc.TestData.r;
    t = spec_resolve_times('to_t_rel', '30/03/2026 11:10:00.500', 'pipeline', r);
    tc.verifyEqual(t, 3 * 86400 + 600.5, 'AbsTol', 1e-6);
end

function testNumberIsTRel(tc)
    r = tc.TestData.r;
    tc.verifyEqual(spec_resolve_times('to_t_rel', 1234.5, 'natus', r), 1234.5);
    tc.verifyEqual(spec_resolve_times('to_t_rel', '1234.5', 'natus', r), 1234.5);
end

function testNonexistentSpringTimeErrors(tc)
    r = tc.TestData.r;
    tc.verifyError(@() spec_resolve_times('to_t_rel', '29/03/2026 02:30:00.000', 'pipeline', r), ...
        'spec_resolve_times:NonexistentTime');
end

function testAmbiguousAutumnTimeErrors(tc)
    r = tc.TestData.r;
    tc.verifyError(@() spec_resolve_times('to_t_rel', '25/10/2026 02:30:00.000', 'pipeline', r), ...
        'spec_resolve_times:AmbiguousTime');
end

function testSummerSessionMatchesSeizuresEvents(tc)
    % 073-s, seizure 1 HPCr: session 06/07/2026 13:00:48.599 + start_s 59657.41797
    % -> start_datetime 07/07/2026 05:35:06 in the reviewed seizures_events
    ss = datetime(1783335648.599, 'ConvertFrom', 'posixtime', 'TimeZone', 'Europe/Paris');
    r = struct('session_start', ss, 'wall_start', datetime(ss.Year, ss.Month, ss.Day, ss.Hour, ss.Minute, ss.Second));
    [n, p] = spec_resolve_times('clocks', 59657.41797, r);
    tc.verifyLessThan(abs(seconds(n - datetime(2026, 7, 7, 5, 35, 6.017))), 2e-3);
    wall_p = datetime(p.Year, p.Month, p.Day, p.Hour, p.Minute, p.Second);
    tc.verifyLessThan(abs(seconds(wall_p - n)), 1e-6);   % no DST change in July: same clock
end
