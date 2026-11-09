function [DW, M] = cs_day_windows(cfg, S, O, rv, K, edf_tbl, edf_lists)
% CS_DAY_WINDOWS  "Days" of v2 (section 6.1), per in-scope mouse:
%   t0 = first sample of the mouse's first EDF -- any campaign, any state:
%        analyzed (02_clean header), excluded, unit out of scope, or never
%        analyzed (EDF header, cs_edf_header_span.m)
%   day k = [t0 + 24(k-1) h, t0 + 24k h), independent of the calendar date.
%
% DW: one row per mouse x day with any recording (analyzed or not):
%   rec_h (any EDF), excluded_h, pending_h, T_eeg_h, T_seizure_obs_h, T_iid_h
%   (analyzed IID time, mean of an EDF's channels, summed over EDFs), seizures
%   by source (onsets inside T_seizure_obs, as v1), crisettes, IIDs (both
%   hemispheres) and rates, plus enters_assoc (T_iid and T_seizure_obs >=
%   cfg.min_valid_h_day) and enters_violin (T_seizure_obs >= cfg.violin_min_h).
% M: per mouse, the interval sets the raster draws (POSIX s).

    tz = cfg.tz;
    om = O.mouse(strcmp({O.mouse.scope}, 'merged'));
    rows = {};
    M = struct('mouse', {}, 'group', {}, 't0', {}, 'n_days', {}, 'rec_all', {}, 'T_eeg', {}, 'T_sz', {}, 'dx', {}, ...
        'excluded', {}, 'pending', {}, 'seizures', {}, 'crisette', {});
    for m = 1:numel(om)
        o = om(m);
        mouse = o.mouse;
        E = edf_tbl(strcmp(edf_tbl.mouse, mouse), :);
        spans = [E.start_posix E.end_posix];
        % never-analyzed EDFs of the mouse (listed, no 02_clean): header span
        L = edf_lists(strcmp(edf_lists.mouse, mouse) & ~ismember(edf_lists.edf, E.edf), :);
        na = zeros(0, 2);
        for i = 1:height(L)
            f = fullfile(cfg.edf_root, mouse, [L.edf{i} '.EDF']);
            if ~isfile(f), f = fullfile(cfg.edf_root, mouse, [L.edf{i} '.edf']); end
            [a, b] = cs_edf_header_span(f, tz);
            if isfinite(a), na(end+1, :) = [a b]; end %#ok<AGROW>
        end
        rec_all = cs_iv('union', [spans; na]);
        t0 = min(rec_all(:, 1));
        nd = ceil((max(rec_all(:, 2)) - t0) / 86400);
        % excluded / pending parts of the recording (not in the analyzed record)
        pend = zeros(0, 2);
        x = rv.excluded(strcmp(rv.excluded.mouse, mouse) & strcmp(rv.excluded.kind, 'pendiente de reprocesamiento'), :);
        if ~isempty(x), pend = cs_iv('union', [x.start_posix x.end_posix]); end
        excl = cs_iv('subtract', cs_iv('subtract', rec_all, o.rec), pend);
        sz = S.seizures(strcmp(S.seizures.mouse, mouse), :);
        cr = S.crisette(strcmp(S.crisette.mouse, mouse), :);
        kc = find(strcmp(K.mouse, mouse));
        ekey = strcat(K.campaign(kc), '|', K.edf(kc));
        [ue, ~, ge] = unique(ekey);
        M(end+1) = struct('mouse', mouse, 'group', o.group, 't0', t0, 'n_days', nd, 'rec_all', rec_all, 'T_eeg', o.T_eeg, ...
            'T_sz', o.T_sz, 'dx', o.dx, 'excluded', excl, 'pending', pend, 'seizures', sz, 'crisette', cr); %#ok<AGROW>
        for k = 1:nd
            w = t0 + [(k - 1) k] * 86400;
            rec_h = cs_iv('length', cs_iv('intersect', rec_all, w)) / 3600;
            if rec_h <= 0, continue; end
            T_sz = cs_iv('length', cs_iv('intersect', o.T_sz, w)) / 3600;
            T_eeg = cs_iv('length', cs_iv('intersect', o.T_eeg, w)) / 3600;
            ch_h = 0; n_iid = 0; T_iid = 0;
            for e = 1:numel(ue)
                cc = kc(ge == e);
                he = arrayfun(@(c) cs_iv('length', cs_iv('intersect', K.inc{c}, w)) / 3600, cc);
                T_iid = T_iid + mean(he);
                ch_h = ch_h + sum(he);
                n_iid = n_iid + sum(arrayfun(@(c) sum(K.on_abs{c} >= w(1) & K.on_abs{c} < w(2)), cc));
            end
            s = sz(sz.start_posix >= w(1) & sz.start_posix < w(2), :);
            s = s(cs_iv('overlap', [s.start_posix s.start_posix + 1e-6], o.T_sz) > 0, :);
            ncr = sum(cr.start_posix >= w(1) & cr.start_posix < w(2));
            rows(end+1, :) = {mouse, o.group, k, char(datetime(w(1), 'ConvertFrom', 'posixtime', 'TimeZone', tz, ...
                'Format', 'yyyy-MM-dd HH:mm')), rec_h, cs_iv('length', cs_iv('intersect', excl, w)) / 3600, ...
                cs_iv('length', cs_iv('intersect', pend, w)) / 3600, T_eeg, T_sz, T_iid, height(s), ...
                sum(strcmp(s.source, 'detected_TP')), sum(strcmp(s.source, 'missed_eeg')), sum(strcmp(s.source, 'video')), ncr, ...
                n_iid, ch_h, n_iid / ch_h, height(s) / T_sz * 24}; %#ok<AGROW>
        end
    end
    DW = cell2table(rows, 'VariableNames', {'mouse', 'group', 'day', 'window_start_local', 'rec_h', 'excluded_h', 'pending_h', ...
        'T_eeg_h', 'T_seizure_obs_h', 'T_iid_h', 'n_seizures', 'n_detected_TP', 'n_missed_eeg', 'n_video', 'n_crisette', ...
        'n_iid', 'iid_channel_h', 'iid_per_h', 'seizures_per_24h'});
    DW.enters_assoc = DW.T_iid_h >= cfg.min_valid_h_day & DW.T_seizure_obs_h >= cfg.min_valid_h_day;
    DW.enters_violin = DW.T_seizure_obs_h >= cfg.violin_min_h;
end
