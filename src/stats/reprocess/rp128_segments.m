function S = rp128_segments()
% RP128_SEGMENTS  Phase 2 (option 1): write the A1-only and the A2-only EDF
% of 128-s 20260729 into the output folder, copying the source's DIGITAL
% samples untouched, with the original's sample times:
%   A1: first sample of the recording (10:07:51.457) up to 17:18:39 (exclusive)
%   A2: 17:23:55 up to the original's last sample (13:02:11.457, exclusive)
% 1 s records, so each segment ends < 1 s before its cut (whole records
% only); the source's annotations inside each segment are carried over.
% Verifies: every sample identical to the source, first-sample times,
% MATLAB's edfread reads them as edf_import will, sample counts.

    o = rp128_paths();
    if ~isfolder(o.out), mkdir(o.out); end
    lab = [o.seg.labels];
    E = rp_edf_read(o.src_edf, lab);
    fs = E.sig(1).fs;
    t0 = rp_start(E, o.tz);
    h0 = datetime([E.hdr.startdate ' ' E.hdr.starttime], 'InputFormat', 'dd.MM.yy HH.mm.ss', 'TimeZone', o.tz);
    tal0_ns = int64(round(E.tal0_s * 1e9));
    ns_per_sample = int64(1e9 / fs);                      % 1953125 ns at 512 Hz, exact
    k1 = rp_first_sample(t0, o.cut1, fs);
    k2 = rp_first_sample(t0, o.cut2, fs);
    c2h = rp_txt_header(o.c2_txt.HPCr);
    n_orig = str2double(c2h.n_samples);
    ann = source_annotations(E);

    first = [1, k2];
    avail = [k1 - 1, n_orig - k2 + 1];          % samples of each segment up to its cut
    S = struct();
    for i = 1:2
        g = o.seg(i);
        n = o.spr * floor(avail(i) / o.spr);
        idx = first(i):first(i) + n - 1;
        ch = find(ismember(lab, g.labels));
        tot_ns = tal0_ns + int64(first(i) - 1) * ns_per_sample;      % first sample, after the source header second
        shift_s = double(idivide(tot_ns, int64(1e9), 'floor'));
        frac_ns = double(tot_ns - int64(shift_s) * int64(1e9));
        seg_t = [double(tot_ns) / 1e9, double(tot_ns) / 1e9 + n / fs];   % s after the source header start
        a = ann(ann.onset_s >= seg_t(1) & ann.onset_s < seg_t(2), :);
        tals = [num2cell(a.onset_s - shift_s), a.bytes];

        W = struct('patient', E.hdr.patient, 'recording', E.hdr.recording, 'start', h0 + seconds(shift_s), ...
            'frac_ns', frac_ns, 'spr', o.spr, 'rec_dur_s', o.spr / fs, 'sig', E.sig(ismember({E.sig.label}, g.labels)), ...
            'x', E.x(idx, ch), 'ann_spr', o.ann_spr);
        W.tals = tals;
        rp_edf_write(g.edf, W);

        % ---- verification ----
        R = rp_edf_read(g.edf, g.labels);
        v.identical = isequal(R.x, E.x(idx, ch));
        v.start = rp_start(R, o.tz);
        v.start_expected = t0 + seconds((first(i) - 1) / fs);
        v.start_err_us = 1e6 * seconds(v.start - v.start_expected);
        v.n = n;
        v.end = v.start + seconds(n / fs);
        tt = edfread(g.edf, 'SelectedSignals', g.labels(1), 'SelectedDataRecords', 1:2, 'TimeOutputType', 'duration');
        v.edfread_rt0_s = seconds(tt.("Record Time")(1));
        y = tt{1, 1}{1};
        v.edfread_max_diff_uV = max(abs(y(:) - rp_edf_phys(E.x(idx(1:numel(y)), ch(1)), E.sig(ch(1)))));
        if i == 1, tp = datetime(2026, 7, 29, 10, 8, 0, 'TimeZone', o.tz);
        else, tp = datetime(2026, 7, 29, 17, 30, 0, 'TimeZone', o.tz); end
        ks = rp_first_sample(t0, tp, fs);
        v.probe = struct('time', t0 + seconds((ks - 1) / fs), 'source', E.x(ks, ch), 'new', R.x(ks - first(i) + 1, :), ...
            'uV', rp_edf_phys(R.x(ks - first(i) + 1, :), E.sig(ch(1))));
        d = find(any(diff(R.x, 1, 1) ~= 0, 2), 1);   % first change of value (flat start after the reconnection?)
        v.first_change = v.start + seconds((d - 1) / fs);
        v.n_annotations = height(a);
        v.annotations = a(:, {'onset_s', 'time', 'text'});
        S.(g.name) = v;
        fprintf('%s: %s\n  %s -> %s, %d samples (%.3f s), identical to source: %d, start error %.3f us\n', g.name, g.edf, ...
            char(v.start, 'yyyy-MM-dd HH:mm:ss.SSSSSS'), char(v.end, 'yyyy-MM-dd HH:mm:ss.SSSSSS'), n, n / fs, ...
            v.identical, v.start_err_us);
        fprintf('  edfread: first Record Time %.9f s, max |diff| vs source %.2e uV; probe %s: source %s / new %s (%s uV)\n', ...
            v.edfread_rt0_s, v.edfread_max_diff_uV, char(v.probe.time, 'HH:mm:ss.SSSSSS'), mat2str(v.probe.source), ...
            mat2str(v.probe.new), mat2str(v.probe.uV, 6));
        fprintf('  first sample whose value changes: %s; %d annotation(s) carried over\n', ...
            char(v.first_change, 'HH:mm:ss.SSS'), v.n_annotations);
    end
    S.counts = struct('n_orig', n_orig, 'A1', S.A1.n, 'A1_trim', avail(1) - S.A1.n, 'gap_DX21', k2 - k1, ...
        'A2', S.A2.n, 'A2_trim', avail(2) - S.A2.n);
    c = S.counts;
    S.counts.sum_ok = c.A1 + c.A1_trim + c.gap_DX21 + c.A2 + c.A2_trim == n_orig;
    fprintf(['Samples: A1 %d + trimmed %d (%.3f s) + DX-21 %d (%.3f s) + A2 %d + trimmed %d (%.3f s) = %d ' ...
        '(original %d, %.0f s): %d\n'], c.A1, c.A1_trim, c.A1_trim / fs, c.gap_DX21, c.gap_DX21 / fs, c.A2, c.A2_trim, ...
        c.A2_trim / fs, c.A1 + c.A1_trim + c.gap_DX21 + c.A2 + c.A2_trim, n_orig, n_orig / fs, S.counts.sum_ok);
    ann.time.Format = 'yyyy-MM-dd HH:mm:ss.SSS';
    writetable(ann(:, {'record', 'onset_s', 'time', 'text'}), fullfile(o.out, 'annotations_source_edf.csv'));
    S.ok = S.A1.identical && S.A2.identical && S.counts.sum_ok && abs(S.A1.start_err_us) < 1 && abs(S.A2.start_err_us) < 1 ...
        && S.A1.edfread_max_diff_uV < 1e-9 && S.A2.edfread_max_diff_uV < 1e-9;
end

function A = source_annotations(E)
% every TAL of the source except the time keeping ones; onsets unwrapped
% (the export's time keeping restarts at 86400 s)
    rows = {};
    t0 = rp_start(E, 'Europe/Paris');
    for r = 1:size(E.ann, 2)
        b = E.ann(:, r)';
        z = [0, find(b == 0)];
        for t = 1:numel(z) - 1
            tal = b(z(t) + 1:z(t + 1) - 1);
            if isempty(tal), continue; end
            p = find(tal == 20, 1);
            onset = char(tal(1:p - 1));
            rest = tal(p + 1:end);
            if all(rest == 20), continue; end          % time keeping TAL
            if any(onset == char(21)), error('rp128_segments:Duration', 'TAL with a duration in record %d', r); end
            on = str2double(onset);
            rec_t = (r - 1) * E.hdr.rec_dur_s;
            on = on + 86400 * round((rec_t - on) / 86400);
            txt = native2unicode(rest(rest ~= 20 | [diff(rest) ~= 0, true]), 'UTF-8');
            rows(end+1, :) = {r, on, rest, strtrim(regexprep(strrep(txt, char(20), ' | '), '[\r\n]+', ' '))}; %#ok<AGROW>
        end
    end
    A = cell2table(rows, 'VariableNames', {'record', 'onset_s', 'bytes', 'text'});
    A.time = t0 - seconds(E.tal0_s) + seconds(A.onset_s);
end
