function D = rp128_compat_diag()
% RP128_COMPAT_DIAG  Read-only follow-up of Phase 1's A1 vs A2 STOP: is
% A2's larger amplitude a port difference or the hour after handling?
%   1. hourly profile: broadband MAD and the pipeline's sigma_band_uV
%      (15-70 Hz without 48-52, signal_quality.m) of A1 (start -> 17:18:39)
%      and A2 (17:23:55 -> end);
%   2. same clock hours on 30/07 in 20260730 (128-s on A2 there too);
%   3. sigma_band_uV of every 128-s EDF in the C1 / C2 qc_report.csv.

    o = rp128_paths();
    addpath(fullfile(o.repo, 'src'), fullfile(o.repo, 'src', 'utils'));
    cfg = rp128_cfg(o, o.master_log);
    lab = {'EEG A1C2', 'EEG A1C4', 'EEG A2C2', 'EEG A2C4'};
    E = rp_edf_read(o.src_edf, lab);
    fs = E.sig(1).fs;
    t0 = rp_start(E, o.tz);
    k1 = rp_first_sample(t0, o.cut1, fs);
    k2 = rp_first_sample(t0, o.cut2, fs);
    h = 3600 * fs;

    % 1. hourly profile, hours counted back from 17:18:39 (A1) and forward from 17:23:55 (A2)
    rows = {};
    for r = 1:2
        for j = floor((k1 - 1) / h):-1:1
            idx = k1 - j * h:k1 - (j - 1) * h - 1;
            rows(end+1, :) = [{o.regions{r}, 'A1', t0 + seconds((idx(1) - 1) / fs)}, m(E, r, idx, fs, cfg)]; %#ok<AGROW>
        end
        for j = 0:floor((size(E.x, 1) - k2 + 1) / h) - 1
            idx = k2 + j * h:k2 + (j + 1) * h - 1;
            rows(end+1, :) = [{o.regions{r}, 'A2', t0 + seconds((idx(1) - 1) / fs)}, m(E, r + 2, idx, fs, cfg)]; %#ok<AGROW>
        end
    end
    D.hourly = cell2table(rows, 'VariableNames', {'region', 'port', 'hour_start', 'MAD_uV', 'P97_abs_uV', 'sigma_band_uV'});
    D.hourly.hour_start.Format = 'dd/MM HH:mm';

    % 2. 20260730 (A2) at the clock hours of the two Phase 1 windows
    N = rp_edf_read(o.next_edf, {'EEG A2C2', 'EEG A2C4'});
    n0 = rp_start(N, o.tz);
    rows = {};
    for w = {o.cut1 - hours(1), o.cut2}
        a = rp_first_sample(n0, w{1} + days(1), fs);
        idx = a:a + h - 1;
        for r = 1:2
            rows(end+1, :) = [{o.regions{r}, N.sig(r).label, w{1} + days(1)}, m(N, r, idx, fs, cfg)]; %#ok<AGROW>
        end
    end
    D.next_day = cell2table(rows, 'VariableNames', {'region', 'channel', 'hour_start', 'MAD_uV', 'P97_abs_uV', 'sigma_band_uV'});
    D.next_day.hour_start.Format = 'dd/MM HH:mm';

    % 3. qc_report of the campaigns
    rows = {};
    for c = {'pipeline_output_campaign_20261006', 'pipeline_output_campaign_20261014'}
        q = readtable(fullfile(o.batch, c{1}, o.subject, '05_summaries', 'qc_report.csv'), 'TextType', 'char', ...
            'VariableNamingRule', 'preserve');
        for i = 1:height(q)
            rows(end+1, :) = {c{1}(end-7:end), q.source_file{i}, q.region{i}, q.sigma_band_uV(i)}; %#ok<AGROW>
        end
    end
    D.qc = cell2table(rows, 'VariableNames', {'campaign', 'edf', 'region', 'sigma_band_uV'});

    disp(D.hourly); disp(D.next_day); disp(D.qc);
end

function v = m(E, c, idx, fs, cfg)
    s = E.sig(c);
    x = rp_edf_phys(E.x(idx, c), s);
    step = abs(s.pmax - s.pmin) / (s.dmax - s.dmin);
    data = struct('signal', x, 'fs', fs, 't_rel', (0:numel(x) - 1)' / fs, 'valid_mask', true(numel(x), 1), ...
        'meta', struct('quantization_step_uV', num2str(step, 17), 'physical_min_uV', num2str(s.pmin), ...
        'physical_max_uV', num2str(s.pmax)));
    q = signal_quality(data, cfg);
    v = {mad(x, 1), prctile(abs(x), 97), q.sigma_band_uV};
end
