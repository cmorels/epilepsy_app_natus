function [A, D] = cs_iid_daily_assoc(cfg, S, K, DW, M)
% CS_IID_DAILY_ASSOC  Section 6.2 of v2: within-animal association of daily
% seizures with the daily IID rate (epileptic mice, days of cs_day_windows.m),
% in three versions:
%   (i)   no post-ictal exclusion
%   (ii)  the cfg.postictal_excl_min(1) min after the END of every final
%         seizure (any source; also when the window leaves a disconnection)
%         removed from the IIDs and from T_iid
%   (iii) the same with cfg.postictal_excl_min(2) min
% For each: rmcorr and Poisson GLMM (cs_within_assoc.m), same day and with
% lag (IIDs of day k vs seizures of day k+1, consecutive valid days only).
% A day is valid when T_iid (of that version) and T_seizure_obs are both
% >= cfg.min_valid_h_day.
%
% A: comparison table; D: the mouse-day table of every version.

    versions = [0, cfg.postictal_excl_min(:)'];
    vname = {'(i) sin exclusion post-ictal', sprintf('(ii) excluye %d min post-ictales', cfg.postictal_excl_min(1)), ...
        sprintf('(iii) excluye %d min post-ictales', cfg.postictal_excl_min(2))};
    rows = {};
    D = table();
    for vi = 1:numel(versions)
        dd = {};
        for m = 1:numel(M)
            if ~strcmp(M(m).group, 'epileptic'), continue; end
            mouse = M(m).mouse;
            sz = S.seizures(strcmp(S.seizures.mouse, mouse), :);
            post = cs_iv('union', [sz.end_posix, sz.end_posix + versions(vi) * 60]);
            kc = find(strcmp(K.mouse, mouse));
            ekey = strcat(K.campaign(kc), '|', K.edf(kc));
            [ue, ~, ge] = unique(ekey);
            dw = DW(strcmp(DW.mouse, mouse), :);
            for r = 1:height(dw)
                w = M(m).t0 + [(dw.day(r) - 1) dw.day(r)] * 86400;
                T_iid = 0; ch_h = 0; n = 0;
                for e = 1:numel(ue)
                    cc = kc(ge == e);
                    he = zeros(numel(cc), 1);
                    for j = 1:numel(cc)
                        x = cs_iv('subtract', cs_iv('intersect', K.inc{cc(j)}, w), post);
                        he(j) = cs_iv('length', x) / 3600;
                        t = K.on_abs{cc(j)};
                        t = t(t >= w(1) & t < w(2));
                        n = n + sum(~cs_in_iv(t, post));
                    end
                    T_iid = T_iid + mean(he);
                    ch_h = ch_h + sum(he);
                end
                dd(end+1, :) = {vname{vi}, mouse, dw.day(r), T_iid, n, n / ch_h, dw.n_seizures(r), dw.T_seizure_obs_h(r)}; %#ok<AGROW>
            end
        end
        Dv = cell2table(dd, 'VariableNames', {'version', 'mouse', 'day', 'T_iid_h', 'n_iid', 'iid_per_h', 'n_seizures', 'T_seizure_obs_h'});
        Dv.valid = Dv.T_iid_h >= cfg.min_valid_h_day & Dv.T_seizure_obs_h >= cfg.min_valid_h_day;
        D = [D; Dv]; %#ok<AGROW>
        same = Dv(Dv.valid, :);
        r1 = cs_within_assoc(same, 'mismo dia');
        % lag: IIDs of day k (valid) vs seizures of day k+1 (valid)
        [tf, j] = ismember(strcat(Dv.mouse, '|', string(Dv.day + 1)), strcat(Dv.mouse, '|', string(Dv.day)));
        ok = Dv.valid & tf;
        ok(ok) = Dv.valid(j(ok));
        lag = table(Dv.mouse(ok), Dv.iid_per_h(ok), Dv.n_seizures(j(ok)), Dv.T_seizure_obs_h(j(ok)), ...
            'VariableNames', {'mouse', 'iid_per_h', 'n_seizures', 'T_seizure_obs_h'});
        r2 = cs_within_assoc(lag, 'IID dia k vs crisis dia k+1');
        rr = [r1; r2];
        rows = [rows; [repmat(vname(vi), size(rr, 1), 1), rr]]; %#ok<AGROW>
    end
    A = cell2table(rows, 'VariableNames', {'version', 'analisis', 'metodo', 'datos', 'conjunto', 'n_raton_dias', 'n_ratones', ...
        'estimacion', 'ci95_lo', 'ci95_hi', 'tipo_ic', 'p', 'nota'});
end
