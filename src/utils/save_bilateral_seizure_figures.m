function files = save_bilateral_seizure_figures(out_dir, base, cfg, data, tr, seizures, seizure_mode, region)
% SAVE_BILATERAL_SEIZURE_FIGURES  Seizure figures for a reconciled channel
% (see bilateral_reconcile.m). Replaces the detectors' own figures when
% cfg.bilateral.rescue_mode ~= 'off' -- the detectors' figures are indexed
% by their per-channel id and know nothing about other channels, and
% neither detector file may be modified.
%
%   files = save_bilateral_seizure_figures(out_dir, base, cfg, data, tr, seizures, seizure_mode, region)
%
% Writes, into out_dir:
%   {base}_seizures.fig/.png            panorama; rescued/imputed events
%                                        outlined with a DASHED edge
%   {base}_seizure{seizure_id:02d}.fig/.png  one zoom per shared seizure_id
%                                        (all of this channel's rows for that
%                                        id, so a fragmented event is one
%                                        figure), titled with its status
% Colours: accepted = red, rescued = orange, imputed = blue (legend on every
% figure). The reference window (union of accepted rows across channels)
% is drawn as black dotted lines on the zoom figures.
%
% Any stale {base}_seizure* file from a previous run is deleted first, so
% the folder never shows more events than seizures_events.csv.

    if ~isfolder(out_dir)
        mkdir(out_dir);
    end
    delete(fullfile(out_dir, [base '_seizure*']));
    files = struct('panorama_fig', '', 'panorama_png', '', 'zoom', {{}});
    if height(seizures) == 0
        return;
    end

    t_rel = data.t_rel(:);
    signal = data.signal(:);
    ll_norm = tr.ll_full / tr.ll_median_global;
    gap_t = gap_segments_seconds(data.valid_mask(:), t_rel);
    sty = status_styles();

    %% panorama
    target_points = 20000;
    fig = figure('Position', [50, 50, 1400, 1100], 'Visible', 'off');
    traces = {signal, tr.bp_full, tr.energy_full, ll_norm};
    colors = {[0.3 0.3 0.3], [0.2 0.4 0.7], [0.2 0.6 0.3], [0.6 0.3 0.6]};
    n_acc = nnz(strcmp(seizures.detection_status, 'accepted'));
    n_ev = numel(unique(seizures.seizure_id));
    titles = {sprintf('%s - Cleaned LFP (%d events reported, %d accepted here; %s branch, bilateral %s)', ...
                  region, n_ev, n_acc, seizure_mode, cfg.bilateral.rescue_mode), ...
              sprintf('Band-passed [%g-%g Hz]', cfg.seizure.bandpass_band(1), cfg.seizure.bandpass_band(2)), ...
              'Energy Metric', 'Line-Length Ratio'};
    ylabels = {'uV', 'Normalized', 'Energy (a.u.)', 'll / ll\_median\_global'};
    for p = 1:4
        subplot(4, 1, p);
        shade_gaps(gap_t); hold on;
        [tp, yp] = decimate_minmax(t_rel, traces{p}, target_points);
        plot(tp, yp, 'Color', colors{p}, 'LineWidth', 0.5);
        add_threshold_line(p, tr, cfg, seizure_mode);
        for k = 1:height(seizures)
            shade_event(seizures.start_s(k), seizures.end_s(k), sty.(seizures.detection_status{k}));
        end
        if p == 1
            status_legend(sty);
        end
        title(titles{p}, 'Interpreter', 'none');
        xlabel('Time (s)'); ylabel(ylabels{p}); grid on; hold off;
    end
    files.panorama_fig = fullfile(out_dir, [base '_seizures.fig']);
    files.panorama_png = fullfile(out_dir, [base '_seizures.png']);
    savefig(fig, files.panorama_fig);
    saveas(fig, files.panorama_png, 'png');
    close(fig);

    %% one zoom per shared seizure_id
    for sid = unique(seizures.seizure_id)'
        rows = seizures(seizures.seizure_id == sid, :);
        status = rows.detection_status{1};
        z0 = max(min([rows.start_s; rows.ref_start_s(1)]) - cfg.seizure.zoom_margin_s, t_rel(1));
        z1 = min(max([rows.end_s; rows.ref_end_s(1)]) + cfg.seizure.zoom_margin_s, t_rel(end));
        iz = t_rel >= z0 & t_rel <= z1;

        fig = figure('Position', [50, 50, 1400, 1100], 'Visible', 'off');
        for p = 1:4
            subplot(4, 1, p);
            shade_gaps(gap_t); hold on;
            plot(t_rel(iz), traces{p}(iz), 'k', 'LineWidth', 0.8);
            add_threshold_line(p, tr, cfg, seizure_mode);
            for k = 1:height(rows)
                shade_event(rows.start_s(k), rows.end_s(k), sty.(rows.detection_status{k}));
                xline(rows.start_s(k), '--', 'Color', sty.(status).edge, 'LineWidth', 1.5);
                xline(rows.end_s(k), '--', 'Color', sty.(status).edge, 'LineWidth', 1.5);
            end
            xline(rows.ref_start_s(1), ':k', 'LineWidth', 1.2);
            xline(rows.ref_end_s(1), ':k', 'LineWidth', 1.2);
            if p == 1
                title({zoom_title(sid, region, rows, cfg), ...
                    sprintf('accepted in %d/%d channel(s): %s | dotted = reference window [%.1f, %.1f] s', ...
                    rows.accepted_in_n_channels(1), rows.n_channels_total(1), empty_as_dash(rows.accepted_in_regions{1}), ...
                    rows.ref_start_s(1), rows.ref_end_s(1))}, 'Interpreter', 'none');
                status_legend(sty);
            elseif p == 4
                title(sprintf('Line-Length Ratio: ll=%.2f | peak\\_energy=%.1f | hf\\_db=%.1f | env\\_cv=%.2f', ...
                    rows.ll_ratio(1), rows.peak_energy_ratio(1), rows.hf_ratio_db(1), rows.envelope_cv(1)), 'Interpreter', 'tex');
            else
                title(titles{p}, 'Interpreter', 'none');
            end
            xlabel('Time (s)'); ylabel(ylabels{p}); grid on; xlim([z0 z1]); hold off;
        end
        zf = fullfile(out_dir, sprintf('%s_seizure%02d.fig', base, sid));
        zp = fullfile(out_dir, sprintf('%s_seizure%02d.png', base, sid));
        savefig(fig, zf);
        saveas(fig, zp, 'png');
        close(fig);
        files.zoom{end+1} = zp;
    end
end

%% ======================================================================
function s = zoom_title(sid, region, rows, cfg)
    status = rows.detection_status{1};
    switch status
        case 'accepted'
            detail = 'accepted';
            if rows.fragmented(1)
                detail = sprintf('accepted, fragmented in %d rows', height(rows));
            end
        case 'rescued'
            switch rows.rejected_by{1}
                case 'll_ratio'
                    detail = sprintf('rescued, ll_ratio=%.2f < %.2f', rows.ll_ratio(1), cfg.seizure_robust.ll_threshold);
                case 'min_duration'
                    detail = sprintf('rescued, duration=%.1f s < %g s', rows.duration_s(1), cfg.seizure.min_seizure_duration);
                otherwise
                    detail = sprintf('rescued (%s)', rows.rejected_by{1});
            end
        case 'imputed'
            if strcmp(rows.rejected_by{1}, 'no_valid_signal')
                detail = 'imputed: no valid signal in this channel';
            else
                detail = 'imputed: no candidate in this channel';
            end
        otherwise
            detail = status;
    end
    s = sprintf('Crisis %d (%s, %s)', sid, region, detail);
end

function sty = status_styles()
    sty.accepted = struct('face', [1 0.7 0.7], 'edge', [0.85 0 0], 'line', '-', 'label', 'accepted (this channel)');
    sty.rescued = struct('face', [1 0.85 0.55], 'edge', [0.9 0.45 0], 'line', '--', 'label', 'rescued (rejected candidate)');
    sty.imputed = struct('face', [0.75 0.85 1], 'edge', [0.1 0.35 0.9], 'line', '--', 'label', 'imputed (no candidate)');
end

function shade_event(s, e, st)
    xregion(s, e, 'FaceColor', st.face, 'FaceAlpha', 0.45, 'EdgeColor', st.edge, 'LineWidth', 2, 'LineStyle', st.line);
end

function status_legend(sty)
    names = fieldnames(sty);
    h = gobjects(numel(names), 1);
    for i = 1:numel(names)
        st = sty.(names{i});
        h(i) = patch(NaN, NaN, st.face, 'EdgeColor', st.edge, 'LineStyle', st.line, 'LineWidth', 2, 'FaceAlpha', 0.45);
    end
    legend(h, cellfun(@(n) sty.(n).label, names, 'UniformOutput', false), 'Location', 'northeast', 'AutoUpdate', 'off');
end

function add_threshold_line(p, tr, cfg, seizure_mode)
    if p == 3
        yline(tr.threshold, 'r--', sprintf('Threshold (%.2e)', tr.threshold), 'LineWidth', 1.5);
    elseif p == 4
        lbl = sprintf('ll\\_threshold (%.2f)', cfg.seizure_robust.ll_threshold);
        if ~strcmp(seizure_mode, 'robust')
            lbl = [lbl ' - informative, legacy branch'];
        end
        yline(cfg.seizure_robust.ll_threshold, 'r--', lbl, 'LineWidth', 1.5);
    end
end

function gap_t = gap_segments_seconds(valid_mask, t_rel)
    gap_idx = mask_to_segments(~valid_mask);
    if isempty(gap_idx)
        gap_t = zeros(0, 2);
    else
        gap_t = [t_rel(gap_idx(:, 1)), t_rel(gap_idx(:, 2))];
    end
end

function shade_gaps(gap_t)
    for i = 1:size(gap_t, 1)
        xregion(gap_t(i, 1), gap_t(i, 2), 'FaceColor', [0.85 0.85 0.85], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    end
end

function s = empty_as_dash(s)
    if isempty(s)
        s = '-';
    end
end
