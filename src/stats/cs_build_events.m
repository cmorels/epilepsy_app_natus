function ev = cs_build_events(rows)
% CS_BUILD_EVENTS  One row per detected event = (campaign, mouse, EDF,
% seizure_id). seizure_id is the bilateral id shared by the HPCr/HPCl rows
% of bilateral_reconcile.m; an event may carry several rows per region
% (fragments, e.g. 110-s seizure 15). Span = [min start, max end] over all
% its rows; detected hemispheres = regions with an accepted/rescued row
% (imputed rows are copies from the other hemisphere).
%
% Event label: the non-empty AE values of its rows (an empty row takes the
% label of its sibling rows); TP and FP together -> 'CONFLICT'.
% Event comment set: union of the rows' normalized AF values.

    [g, ev] = findgroups(rows(:, {'campaign', 'mouse', 'edf', 'seizure_id'}));
    n = height(ev);
    ev.event_id = strcat(ev.campaign, '|', ev.mouse, '|', ev.edf, '|', arrayfun(@num2str, ev.seizure_id, 'UniformOutput', false));
    ev.start_posix = splitapply(@min, rows.start_posix, g);
    ev.end_posix = splitapply(@max, rows.end_posix, g);
    ev.start_s = splitapply(@min, rows.start_s, g);
    ev.end_s = splitapply(@max, rows.end_s, g);
    ev.duration_s = ev.end_posix - ev.start_posix;
    ev.n_rows = splitapply(@numel, rows.start_s, g);
    ev.max_ll_ratio = splitapply(@max, rows.ll_ratio, g);
    ev.max_peak_energy_ratio = splitapply(@max, rows.peak_energy_ratio, g);
    ev.label = repmat({''}, n, 1);
    ev.n_rows_unlabeled = zeros(n, 1);
    ev.comments = repmat({''}, n, 1);
    ev.laterality = repmat({''}, n, 1);
    ev.statuses = repmat({''}, n, 1);
    for j = 1:n
        r = rows(g == j, :);
        labs = unique(r.label(~cellfun(@isempty, r.label)));
        if isempty(labs)
            ev.label{j} = '';
        elseif numel(labs) == 1
            ev.label{j} = labs{1};
        else
            ev.label{j} = 'CONFLICT';
        end
        ev.n_rows_unlabeled(j) = sum(cellfun(@isempty, r.label));
        cm = unique(r.comment(~cellfun(@isempty, r.comment)));
        ev.comments{j} = strjoin(cm', ';');
        det = unique(r.region(ismember(r.detection_status, {'accepted', 'rescued'})));
        if all(ismember({'HPCr', 'HPCl'}, det))
            ev.laterality{j} = 'bilateral';
        elseif isequal(det, {'HPCr'})
            ev.laterality{j} = 'HPCr_only';
        elseif isequal(det, {'HPCl'})
            ev.laterality{j} = 'HPCl_only';
        else
            ev.laterality{j} = 'unknown';
        end
        ev.statuses{j} = strjoin(strcat(r.region, ':', r.detection_status)', ';');
    end
end
