function K = cs_iid_channel_intervals(C, S)
% CS_IID_CHANNEL_INTERVALS  Per channel of the main IID run (C, variant
% 'principal'), in absolute time (POSIX s): the analyzed intervals (the
% channel span minus detect_iid's exclusion intervals) and the onsets of
% the re-detected complexes. Used by the day windows, the daily
% association and the peri-ictal analysis (v2, section 6).
%
%   K: campaign, mouse, group, edf, region, start_posix, duration_s, inc (Nx2), on_abs (Mx1)

    Cm = C(strcmp(C.variant, 'principal'), :);
    n = height(Cm);
    inc = cell(n, 1);
    on_abs = cell(n, 1);
    grp = cell(n, 1);
    for i = 1:n
        inc{i} = cs_iv('subtract', [0 Cm.duration_s(i)], Cm.excl{i}) + Cm.start_posix(i);
        on_abs{i} = Cm.complexes{i}(:, 1) + Cm.start_posix(i);
        grp{i} = S.units.group{find(strcmp(S.units.mouse, Cm.mouse{i}), 1)};
    end
    K = table(Cm.campaign, Cm.mouse, grp, Cm.edf, Cm.region, Cm.start_posix, Cm.duration_s, inc, on_abs, ...
        'VariableNames', {'campaign', 'mouse', 'group', 'edf', 'region', 'start_posix', 'duration_s', 'inc', 'on_abs'});
end
