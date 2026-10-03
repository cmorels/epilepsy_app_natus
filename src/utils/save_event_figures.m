function [events, notes] = save_event_figures(events, channels, seizures_dir, cfg)
% SAVE_EVENT_FIGURES  Per-event figures for the rows bilateral_events.m
% grouped, organized by category:
%
%   03_seizures/<Crisis|Candidates|Candidates_in_band>/individual/
%       {subject}_{session}_event{ID:02d}_{region}.png/.fig   (one per channel)
%   03_seizures/<category>/joint/
%       {subject}_{session}_event{ID:02d}_joint.png/.fig     (4 rows x N channels)
%
%   [events, notes] = save_event_figures(events, channels, seizures_dir, cfg)
%
% events   : bilateral_events.m output for ONE recording; its
%            figure_individual_path / figure_joint_path are filled in.
% channels : struct array, one per ROBUST channel of the recording, with
%            .region .subject_id .session_start .label .t_rel .signal
%            .trace (energy_full, bp_full, ll_full, ll_median_global, threshold)
%            in manifest order. Every channel gets a column / an individual
%            figure for every event, detected there or not: seeing the
%            hemisphere that did not fire is the point of the joint figure.
%
% Joint figure: columns = left-hemisphere regions, then right, then any
% unknown region (manifest order) -- resolved ONLY from the region name via
% hemisphere_of_region.m. Rows = the four traces drawn by plot_trace_voltage
% / _bandpassed / _energy / _linelength (the SAME functions draw the
% individual figures, so both views of one event can never diverge). Same
% time window on every panel (reference window +/- cfg.seizure.zoom_margin_s,
% x-axes linked) and the SAME y-limits across columns within each row; the
% energy row switches to a log scale on EVERY column when the columns'
% maxima differ by more than 100x (never per-column scales). Each axes is
% tagged 'r<row>_c<col>' with UserData.region/.hemisphere for programmatic
% checks.
%
% Housekeeping, scoped to this subject + session only (prefix
% {subject}_{session}_event): any existing file with that prefix in any of
% the three category folders that is not one of the files this run
% produces (event gone, category changed, channel gone) is deleted, so an
% event never lives in two folders. With cfg.general.overwrite = false an
% expected file that already exists is not redrawn.

    notes = {};
    cats = {'Crisis', 'Candidates', 'Candidates_in_band'};
    robust = strcmp(events.seizure_mode, 'robust') & ~isnan(events.event_id);
    if isempty(channels) || ~any(robust)
        return;
    end
    dirs = make_event_dirs(seizures_dir, cats);
    fmts = cellstr(cfg.output.joint_figure_format);
    [ll_accept, ll_reject] = resolve_ll_band(cfg);

    subj = channels(1).subject_id;
    ss = channels(1).session_start;
    stamp = char(string(ss, 'yyyyMMdd_HHmmss'));
    prefix = sprintf('%s_%s_event', subj, stamp);

    % column order: left, right, unknown (manifest order within each)
    hemi = arrayfun(@(c) hemisphere_of_region(c.region, cfg), channels, 'UniformOutput', false);
    order = [find(strcmp(hemi, 'left')), find(strcmp(hemi, 'right')), find(strcmp(hemi, 'unknown'))];
    channels = channels(order);
    hemi = hemi(order);
    if any(strcmp(hemi, 'unknown'))
        notes{end+1} = sprintf('%s: region(s) %s not in cfg.output.hemisphere -> extra column(s) on the right', ...
            prefix, strjoin(arrayfun(@(c) c.region, channels(strcmp(hemi, 'unknown')), 'UniformOutput', false), ',')); %#ok<AGROW>
    end

    ev_ids = unique(events.event_id(robust))';
    expected = {};
    plan = struct('eid', {}, 'cat', {}, 'rows', {}, 'ind', {}, 'joint', {});
    for eid = ev_ids
        r = find(robust & events.event_id == eid);
        cat = events.category{r(1)};
        ind = cell(1, numel(channels));
        for c = 1:numel(channels)
            name = sprintf('%s%02d_%s', prefix, eid, channels(c).region);
            ind{c} = fullfile(dirs.(cat).individual, name);
            expected = [expected, strcat(cat, '|individual|', name, '.', fmts)]; %#ok<AGROW>
        end
        joint = '';
        if cfg.output.joint_figures
            name = sprintf('%s%02d_joint', prefix, eid);
            joint = fullfile(dirs.(cat).joint, name);
            expected = [expected, strcat(cat, '|joint|', name, '.', fmts)]; %#ok<AGROW>
        end
        plan(end+1) = struct('eid', eid, 'cat', cat, 'rows', r, 'ind', {ind}, 'joint', joint); %#ok<AGROW>
    end

    n_deleted = cleanup_stale(dirs, cats, prefix, expected);
    if n_deleted > 0
        notes{end+1} = sprintf('%s: removed %d stale event figure file(s) (old event ids / categories)', prefix, n_deleted);
    end

    for p = plan
        rows = events(p.rows, :);
        ref = [rows.ref_start_s(1), rows.ref_end_s(1)];
        win = [ref(1) - cfg.seizure.zoom_margin_s, ref(2) + cfg.seizure.zoom_margin_s];
        info = column_info(rows, channels, hemi, ref, ll_accept, ll_reject);
        ttl = sprintf('%s | session %s | event %02d | %s | start %s | reference window %.1f s', ...
            subj, stamp, p.eid, p.cat, char(string(ss + seconds(ref(1)), 'yyyy-MM-dd HH:mm:ss')), ref(2) - ref(1));

        for c = 1:numel(channels)
            if ~all_exist(p.ind{c}, fmts) || cfg.general.overwrite
                draw_grid(channels(c), info(c), win, ref, cfg, ttl, p.ind{c}, fmts);
            end
            events.figure_individual_path(p.rows(strcmp(rows.region, channels(c).region))) = {[p.ind{c} '.png']};
        end
        if ~isempty(p.joint)
            if ~all_exist(p.joint, fmts) || cfg.general.overwrite
                draw_grid(channels, info, win, ref, cfg, ttl, p.joint, fmts);
            end
            events.figure_joint_path(p.rows) = {[p.joint '.png']};
        end
    end
end

%% ======================================================================
function dirs = make_event_dirs(seizures_dir, cats)
% The category tree under 03_seizures/ (sibling of run_pipeline_edf.m's
% make_output_dirs: created only when an event figure is actually written).
    for i = 1:numel(cats)
        for sub = {'individual', 'joint'}
            d = fullfile(seizures_dir, cats{i}, sub{1});
            if ~isfolder(d)
                mkdir(d);
            end
            dirs.(cats{i}).(sub{1}) = d;
        end
    end
end

function n = cleanup_stale(dirs, cats, prefix, expected)
% expected: 'category|sub|filename' keys (never full paths: the same folder
% can be spelled several ways -- / vs \, 8.3 short names -- on Windows).
    n = 0;
    for i = 1:numel(cats)
        for sub = {'individual', 'joint'}
            d = dirs.(cats{i}).(sub{1});
            f = dir(fullfile(d, [prefix '*']));
            for k = 1:numel(f)
                if ~any(strcmpi(sprintf('%s|%s|%s', cats{i}, sub{1}, f(k).name), expected))
                    delete(fullfile(d, f(k).name));
                    n = n + 1;
                end
            end
        end
    end
end

function tf = all_exist(base, fmts)
    tf = all(cellfun(@(e) exist([base '.' e], 'file') == 2, fmts));
end

function info = column_info(rows, channels, hemi, ref, ll_accept, ll_reject)
    hemi_es = containers.Map({'left', 'right', 'unknown'}, {'izquierdo', 'derecho', 'desconocido'});
    for c = numel(channels):-1:1
        ch = channels(c);
        mine = rows(strcmp(rows.region, ch.region), :);
        if any(strcmp(mine.ll_status, 'accepted'))
            status = 'accepted';
        elseif any(strcmp(mine.ll_status, 'in_band'))
            status = 'in_band';
        else
            status = 'no_detection';
        end
        if ~isempty(mine)
            ll = max(mine.ll_ratio);
            spans = [mine.start_s, mine.end_s];
        else
            gs = round(ref(1) * ch.fs) + 1; ge = round(ref(2) * ch.fs) + 1;
            ll = median(ch.trace.ll_full(gs:ge), 'omitnan') / ch.trace.ll_median_global;
            spans = zeros(0, 2);
        end
        [d, which] = min(abs(ll - [ll_accept, ll_reject]));
        names = {'accept', 'reject'};
        label = regexprep(ch.label, '^EEG\s*', '');
        if strcmp(status, 'no_detection')
            head = sprintf('%s (%s, %s) - sin detección - ll_ratio en la ventana %.2f (a %.2f de %s)', ...
                ch.region, hemi_es(hemi{c}), label, ll, d, names{which});
        else
            head = sprintf('%s (%s, %s) - %s - ll_ratio %.2f (a %.2f de %s)', ...
                ch.region, hemi_es(hemi{c}), label, status, ll, d, names{which});
        end
        info(c) = struct('status', status, 'head', head, 'spans', spans, 'hemi', hemi{c});
    end
end

function draw_grid(channels, info, win, ref, cfg, ttl, base, fmts)
% One column per channel (1 column = individual figure, N = joint figure),
% rows = the four shared trace functions. Common y-limits per row.
    N = numel(channels);
    fig = figure('Visible', 'off', 'Position', [50, 50, 650 * N + 150, 1150]);
    tl = tiledlayout(fig, 4, N, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, ttl, 'Interpreter', 'none', 'FontWeight', 'bold');
    ax = gobjects(4, N);
    rng = nan(4, N, 2);
    emax = nan(1, N); emin_pos = nan(1, N);
    for c = 1:N
        ch = channels(c);
        t = ch.t_rel(:);
        ll_norm = ch.trace.ll_full / ch.trace.ll_median_global;
        for r = 1:4
            ax(r, c) = nexttile(tl, (r - 1) * N + c);
            ax(r, c).Tag = sprintf('r%d_c%d', r, c);
            ax(r, c).UserData = struct('region', ch.region, 'hemisphere', info(c).hemi);
            shade(ax(r, c), ref, info(c));
            switch r
                case 1, yr = plot_trace_voltage(ax(r, c), t, ch.signal(:), win);
                case 2, yr = plot_trace_bandpassed(ax(r, c), t, ch.trace.bp_full, win, cfg);
                case 3
                    yr = plot_trace_energy(ax(r, c), t, ch.trace.energy_full, win, ch.trace.threshold);
                    iz = t >= win(1) & t <= win(2);
                    e = ch.trace.energy_full(iz);
                    emax(c) = max(e, [], 'omitnan');
                    ep = e(e > 0);
                    if ~isempty(ep), emin_pos(c) = min(ep); end
                case 4, yr = plot_trace_linelength(ax(r, c), t, ll_norm, win, cfg);
            end
            rng(r, c, 1) = yr(1);
            rng(r, c, 2) = yr(2);
        end
        title(ax(1, c), info(c).head, 'Interpreter', 'none', 'FontSize', 9);
    end
    % same x on every panel; same y within each row
    linkaxes(ax(:), 'x');
    xlim(ax(1, 1), win);
    use_log = N > 1 && max(emax) / max(min(emax(emax > 0)), eps) > 100;
    for r = 1:4
        lo = min(rng(r, :, 1), [], 'omitnan'); hi = max(rng(r, :, 2), [], 'omitnan');
        if r == 3 && use_log
            lo = min(emin_pos, [], 'omitnan');
            set(ax(r, :), 'YScale', 'log');
            yl = [lo, hi * 1.2];
        else
            pad = 0.05 * max(hi - lo, eps);
            yl = [lo - pad, hi + pad];
        end
        if all(isfinite(yl)) && yl(2) > yl(1)
            set(ax(r, :), 'YLim', yl);
        end
    end
    status_legend(ax(1, 1));
    formats_save(fig, base, fmts);
    close(fig);
end

function shade(ax, ref, info)
    st = ll_status_style(info.status);
    xregion(ax, ref(1), ref(2), 'FaceColor', st.face, 'FaceAlpha', 0.35, 'EdgeColor', st.edge, ...
        'LineWidth', 1.5, 'LineStyle', st.line);
    for k = 1:size(info.spans, 1)
        xline(ax, info.spans(k, 1), '--', 'Color', st.edge, 'LineWidth', 1.2);
        xline(ax, info.spans(k, 2), '--', 'Color', st.edge, 'LineWidth', 1.2);
    end
end

function status_legend(ax)
    names = {'accepted', 'in_band', 'no_detection'};
    h = gobjects(1, 3);
    labels = cell(1, 3);
    for i = 1:3
        st = ll_status_style(names{i});
        labels{i} = st.label;
        h(i) = patch(ax, NaN, NaN, st.face, 'EdgeColor', st.edge, 'LineStyle', st.line, 'LineWidth', 1.5, 'FaceAlpha', 0.35);
    end
    legend(ax, h, labels, 'Location', 'northeast', 'AutoUpdate', 'off', 'FontSize', 7);
end

function formats_save(fig, base, fmts)
    for i = 1:numel(fmts)
        switch lower(fmts{i})
            case 'png', exportgraphics(fig, [base '.png'], 'Resolution', 150);
            case 'fig', savefig(fig, [base '.fig']);
            otherwise, error('save_event_figures:BadFormat', 'Unsupported figure format ''%s''.', fmts{i});
        end
    end
end
