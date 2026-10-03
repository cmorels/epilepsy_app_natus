function [events, notes] = bilateral_events(seizures_events, cfg)
% BILATERAL_EVENTS  Group robust-branch detections of the same animal and
% session across channels into EVENTS, give each a shared event_id and a
% category, without touching any detection.
%
%   [events, notes] = bilateral_events(seizures_events, cfg)
%
% Input: seizures_events rows (the table run_pipeline_edf.m builds per
% recording with build_seizure_event_rows), with an ll_status column
% ('accepted' | 'in_band'; legacy rows are always 'accepted'). Output: the
% same rows, same order, with event_columns('seizures_events') appended.
% notes: cellstr of log messages (single-channel animals, unknown regions).
% Called after all channels of a recording are detected and before any
% summary is written or any event figure is drawn.
%
% GROUPING (robust rows only -- the review band and the categories are a
% robust-branch concept; legacy rows get the column defaults and no event):
%  - key = (subject_id, session_start). NOT the EDF file: in
%    cfg.edf.channels.mode='log' one EDF holds several animals.
%  - rows of DIFFERENT regions belong to one event if [start_s, end_s]
%    overlap or are less than cfg.seizure_robust.bilateral_tol_s apart;
%    transitive (union-find). Rows of the SAME region are never linked
%    directly; if they end up in one event they stay separate rows with the
%    same event_id and fragmented = true.
%  - event_id: 1..K per (subject_id, session_start), chronological by the
%    group's earliest start, across all categories. seizure_id (the
%    channel-local index) is left untouched.
%  - ref_start_s / ref_end_s: union of the group's rows.
%
% CATEGORY, counting CHANNELS (distinct regions) of the event:
%    n_accepted >= 2                  -> 'Crisis'
%    n_accepted == 1                  -> 'Candidates'
%    n_accepted == 0, n_in_band >= 1  -> 'Candidates_in_band'
% 'Candidates_in_band' is material the binary pipeline discarded outright:
% never mix it into Crisis counts or rates. A single-channel animal can at
% most reach 'Candidates' (noted in notes).
%
% EXTENSION POINT (cross-hemisphere rescue): rows added later for a
% hemisphere that did not accept an event (e.g. bilateral_reconcile.m's
% detection_status 'rescued' / 'imputed') can go through this same
% function: counts_toward_category() below is the single place deciding
% which rows count for n_accepted / n_in_band, and every other column
% (event_id, ref window, fragmented, hemisphere) applies to them unchanged.

    notes = {};
    events = seizures_events;
    n = height(events);
    events = add_event_columns(events, 'seizures_events');
    if n == 0
        return;
    end
    events.hemisphere = hemisphere_of_region(events.region, cfg);
    unknown = unique(events.region(strcmp(events.hemisphere, 'unknown')));
    for i = 1:numel(unknown)
        notes{end+1} = sprintf('region "%s" is in neither cfg.output.hemisphere list: hemisphere=unknown (extra joint-figure column)', unknown{i}); %#ok<AGROW>
    end

    robust = strcmp(events.seizure_mode, 'robust');
    if ~any(robust)
        return;
    end
    tol = cfg.seizure_robust.bilateral_tol_s;
    keys = strcat(events.subject_id, '|', cellstr(string(events.session_start, 'yyyy-MM-dd''T''HH:mm:ss.SSS')));

    for key = unique(keys(robust), 'stable')'
        rows = find(robust & strcmp(keys, key{1}));
        m = numel(rows);
        parent = 1:m;
        for a = 1:m
            for b = a+1:m
                ra = rows(a); rb = rows(b);
                if ~strcmp(events.region{ra}, events.region{rb}) && ...
                        max(events.start_s(ra), events.start_s(rb)) - min(events.end_s(ra), events.end_s(rb)) < tol
                    parent = uf_union(parent, a, b);
                end
            end
        end
        roots = arrayfun(@(i) uf_find(parent, i), 1:m);
        [~, ~, g] = unique(roots);
        g = g(:);
        nG = max(g);
        g_start = accumarray(g, events.start_s(rows), [nG 1], @min);
        g_end = accumarray(g, events.end_s(rows), [nG 1], @max);
        [~, order] = sortrows([g_start, g_end]);
        eid = zeros(nG, 1);
        eid(order) = 1:nG;

        regions_in_session = unique(events.region(rows), 'stable');
        if numel(regions_in_session) < 2
            notes{end+1} = sprintf('%s: only one robust channel with detections (%s) -- events can be at most ''Candidates''', ...
                key{1}, strjoin(regions_in_session, ',')); %#ok<AGROW>
        end

        for k = 1:nG
            idx = rows(g == k);
            regs = events.region(idx);
            counts = counts_toward_category(events(idx, :));
            acc_regs = unique(regs(counts & strcmp(events.ll_status(idx), 'accepted')), 'stable');
            band_regs = unique(regs(counts & strcmp(events.ll_status(idx), 'in_band')), 'stable');
            n_acc = numel(acc_regs);
            n_band = numel(setdiff(band_regs, acc_regs));
            if n_acc >= 2
                cat = 'Crisis';
            elseif n_acc == 1
                cat = 'Candidates';
            else
                cat = 'Candidates_in_band';
            end
            for i = idx(:)'
                events.event_id(i) = eid(k);
                events.category{i} = cat;
                events.n_accepted(i) = n_acc;
                events.n_in_band(i) = n_band;
                events.n_channels_in_group(i) = numel(unique(regs));
                events.fragmented(i) = nnz(strcmp(regs, events.region{i})) > 1;
                events.accepted_in_regions{i} = strjoin(acc_regs, ';');
                events.ref_start_s(i) = g_start(k);
                events.ref_end_s(i) = g_end(k);
            end
        end
    end
end

function tf = counts_toward_category(rows)
% Which rows of an event count for n_accepted / n_in_band. Today: every row
% (the detector's own). If rescued/imputed rows are added (see EXTENSION
% POINT above), only rows the channel produced on its own merit count.
    if ismember('detection_status', rows.Properties.VariableNames)
        tf = ismember(rows.detection_status, {'accepted', 'not_reconciled'});
    else
        tf = true(height(rows), 1);
    end
end

function r = uf_find(parent, i)
    while parent(i) ~= i
        i = parent(i);
    end
    r = i;
end

function parent = uf_union(parent, a, b)
    ra = uf_find(parent, a);
    rb = uf_find(parent, b);
    if ra ~= rb
        parent(max(ra, rb)) = min(ra, rb);
    end
end
