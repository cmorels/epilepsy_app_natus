function O = cs_observation_time(cfg, S, rv, edf_tbl)
% CS_OBSERVATION_TIME  Denominators (section 5.2), as interval sets in
% POSIX seconds, per in-scope mouse and per scope (C1, C2, merged):
%
%   rec   = union of the mouse's analyzed EDF spans in that scope, minus
%           excluded EDFs (rv.excluded, cs_excluded_edfs.m: DX "EDF completo",
%           cfg.edf_excluded, EDFs without complete AE)
%   T_eeg = rec - NaN gaps - union of the mouse's DX intervals
%   T_sz  = T_eeg + (DX intervals with Video revisado = Si, within rec)
%
% Overlapping intervals are unioned before subtracting, so nothing is
% subtracted twice. Also returns the mouse x EDF x local-day table, the
% mouse x day table and clock-bin exposures (hours, light/dark phases).

    tz = cfg.tz;
    units = S.units(S.units.included, :);
    mice = unique(units.mouse);
    scopes = {'C1', 'C2', 'merged'};
    disc = rv.excluded;   % cs_excluded_edfs.m
    O.mouse = struct('mouse', {}, 'group', {}, 'scope', {}, 'rec', {}, 'gaps', {}, 'dx', {}, 'T_eeg', {}, 'T_sz', {});
    edf_rows = {};
    flags = cell(0, 3);
    for m = 1:numel(mice)
        mouse = mice{m};
        grp = units.group{find(strcmp(units.mouse, mouse), 1)};
        dxm = rv.dx(strcmp(rv.dx.mouse, mouse), :);
        dx_all = [dxm.start_posix dxm.end_posix];
        dx_vid = [dxm.start_posix(dxm.video_reviewed) dxm.end_posix(dxm.video_reviewed)];
        for s = 1:numel(scopes)
            sc = scopes{s};
            camps = units.campaign(strcmp(units.mouse, mouse));
            if ~strcmp(sc, 'merged')
                camps = intersect(camps, {sc});
            end
            if isempty(camps)
                continue;
            end
            E = edf_tbl(strcmp(edf_tbl.mouse, mouse) & ismember(edf_tbl.campaign, camps) & edf_tbl.analyzed, :);
            is_disc = ismember(E.edf, disc.edf(strcmp(disc.mouse, mouse)));
            E = E(~is_disc, :);
            spans = [E.start_posix E.end_posix];
            if strcmp(sc, 'merged') && size(spans, 1) > 1
                ov = sum(spans(:, 2) - spans(:, 1)) - cs_iv('length', spans);
                if ov > 1
                    flags(end+1, :) = {'revisar', 'EDFs solapados', sprintf('%s: %.0f s de EDFs que se solapan en el tiempo (unidos, no sumados)', mouse, ov)}; %#ok<AGROW>
                end
            end
            gaps = zeros(0, 2);
            for i = 1:height(E)
                gaps = [gaps; E.start_posix(i) + E.gaps{i}]; %#ok<AGROW>
            end
            rec = cs_iv('union', spans);
            T_eeg = cs_iv('subtract', cs_iv('subtract', rec, gaps), dx_all);
            T_sz = cs_iv('union', T_eeg, cs_iv('intersect', dx_vid, rec));
            O.mouse(end+1) = struct('mouse', mouse, 'group', grp, 'scope', sc, 'rec', rec, 'gaps', gaps, ...
                'dx', cs_iv('intersect', dx_all, rec), 'T_eeg', T_eeg, 'T_sz', T_sz);
            if strcmp(sc, 'merged')
                for i = 1:height(E)
                    sp = [E.start_posix(i) E.end_posix(i)];
                    [~, d_rec] = cs_clock_split(sp, tz, [0 1440]);
                    [~, d_eeg] = cs_clock_split(cs_iv('intersect', T_eeg, sp), tz, [0 1440]);
                    [~, d_sz] = cs_clock_split(cs_iv('intersect', T_sz, sp), tz, [0 1440]);
                    [~, d_dx] = cs_clock_split(cs_iv('intersect', dx_all, sp), tz, [0 1440]);
                    for k = 1:height(d_rec)
                        dd = d_rec.date(k);
                        edf_rows(end+1, :) = {mouse, grp, E.campaign{i}, E.edf{i}, dd, d_rec.total_s(k) / 3600, ...
                            pick(d_dx, dd) / 3600, pick(d_eeg, dd) / 3600, pick(d_sz, dd) / 3600}; %#ok<AGROW>
                    end
                end
            end
        end
    end
    O.edf_day = cell2table(edf_rows, 'VariableNames', {'mouse', 'group', 'campaign', 'edf', 'date', 'rec_h', ...
        'dx_h', 'T_eeg_h', 'T_seizure_obs_h'});

    % mouse x day (union across EDFs) and clock exposures, merged scope
    day_rows = {};
    for k = 1:numel(O.mouse)
        o = O.mouse(k);
        [O.mouse(k).hour_eeg_s, ~] = cs_clock_split(o.T_eeg, tz, 0:60:1440);
        [O.mouse(k).hour_sz_s, dsz] = cs_clock_split(o.T_sz, tz, 0:60:1440);
        [O.mouse(k).phase_sz_s, ~] = cs_clock_split(o.T_sz, tz, cfg.phase_edges_min);
        [O.mouse(k).phase_eeg_s, ~] = cs_clock_split(o.T_eeg, tz, cfg.phase_edges_min);
        if strcmp(o.scope, 'merged')
            [~, deeg] = cs_clock_split(o.T_eeg, tz, [0 1440]);
            [~, drec] = cs_clock_split(o.rec, tz, [0 1440]);
            for i = 1:height(dsz)
                day_rows(end+1, :) = {o.mouse, o.group, dsz.date(i), pick(drec, dsz.date(i)) / 3600, ...
                    pick(deeg, dsz.date(i)) / 3600, dsz.total_s(i) / 3600}; %#ok<AGROW>
            end
        end
    end
    O.day = cell2table(day_rows, 'VariableNames', {'mouse', 'group', 'date', 'rec_h', 'T_eeg_h', 'T_seizure_obs_h'});
    O.day.valid_day = O.day.T_seizure_obs_h >= cfg.min_valid_h_day;
    O.flags = flags;
end

function v = pick(dtbl, d)
    j = find(dtbl.date == d, 1);
    if isempty(j)
        v = 0;
    else
        v = dtbl.total_s(j);
    end
end
