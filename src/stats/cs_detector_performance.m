function P = cs_detector_performance(cfg, S, O, rv)
% CS_DETECTOR_PERFORMANCE  Section 6.2, per campaign x mouse (in-scope
% units) plus campaign x group totals. Event level (a bilateral event =
% one detection). TP = detected events labeled TP (crisette excluded and
% reported apart), FP = events labeled FP, FN = missed seizures with EEG
% (source missed_eeg; video seizures are NOT false negatives: no EEG).
% FP/24 h over T_eeg of that campaign; FP inside a disconnection (any
% overlap with a DX interval of the mouse) reported apart, and
% FP_outside_dx/24 h over the same T_eeg.

    cats = {'IIDs', 'grooming/scratching', 'other', 'disconnect'};
    units = S.units(S.units.included, :);
    ev = S.events;
    rows = {};
    for i = 1:height(units)
        c = units.campaign{i}; m = units.mouse{i};
        sel = ev.usable & strcmp(ev.campaign, c) & strcmp(ev.mouse, m);
        tp = sum(sel & ismember(ev.label, {'TP', 'CONFLICT'}));
        fpe = ev(sel & strcmp(ev.label, 'FP'), :);
        ncr = sum(strcmp(S.crisette.campaign, c) & strcmp(S.crisette.mouse, m));
        fn = sum(strcmp(S.seizures.mouse, m) & strcmp(S.seizures.source, 'missed_eeg') & contains(S.seizures.campaign, c));
        o = O.mouse(strcmp({O.mouse.mouse}, m) & strcmp({O.mouse.scope}, c));
        Teeg = cs_iv('length', o.T_eeg) / 3600;
        dxm = rv.dx(strcmp(rv.dx.mouse, m), :);
        in_dx = cs_iv('overlap', [fpe.start_posix fpe.end_posix], [dxm.start_posix dxm.end_posix]) > 0;
        byc = cellfun(@(k) sum(strcmp(fpe.comments, k)), cats);
        n_other_cat = height(fpe) - sum(byc);
        rows(end+1, :) = {c, m, units.group{i}, tp, height(fpe), fn, ncr, tp / (tp + fn), tp / (tp + height(fpe)), ...
            Teeg, height(fpe) / Teeg * 24, sum(in_dx), sum(~in_dx), sum(~in_dx) / Teeg * 24, ...
            byc(1), byc(2), byc(3), byc(4), n_other_cat}; %#ok<AGROW>
    end
    names = {'campaign', 'mouse', 'group', 'TP', 'FP', 'FN', 'n_crisette', 'sensitivity', 'PPV', 'T_eeg_h', ...
        'FP_per_24h', 'FP_in_dx', 'FP_outside_dx', 'FP_outside_dx_per_24h', 'FP_IIDs', 'FP_grooming_scratching', ...
        'FP_other', 'FP_disconnect', 'FP_mixed_or_unlabeled_AF'};
    T = cell2table(rows, 'VariableNames', names);
    tot = {};
    for c = {'C1', 'C2', 'C1+C2'}
        for g = {'epileptic', 'control', 'all'}
            sel = true(height(T), 1);
            if ~strcmp(c{1}, 'C1+C2'), sel = sel & strcmp(T.campaign, c{1}); end
            if ~strcmp(g{1}, 'all'), sel = sel & strcmp(T.group, g{1}); end
            t = T(sel, :);
            if isempty(t), continue; end
            tp = sum(t.TP); fp = sum(t.FP); fn = sum(t.FN); te = sum(t.T_eeg_h);
            tot(end+1, :) = {c{1}, ['TOTAL:' g{1}], g{1}, tp, fp, fn, sum(t.n_crisette), tp / (tp + fn), tp / (tp + fp), te, ...
                fp / te * 24, sum(t.FP_in_dx), sum(t.FP_outside_dx), sum(t.FP_outside_dx) / te * 24, sum(t.FP_IIDs), ...
                sum(t.FP_grooming_scratching), sum(t.FP_other), sum(t.FP_disconnect), sum(t.FP_mixed_or_unlabeled_AF)}; %#ok<AGROW>
        end
    end
    P = [T; cell2table(tot, 'VariableNames', names)];
end
