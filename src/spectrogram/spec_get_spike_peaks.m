function pk = spec_get_spike_peaks(ev, r, seg, a, b, seizures, cfg)
% SPEC_GET_SPIKE_PEAKS  Spike times of one IID interval [a, b] (t_rel, s)
% of one channel, from the source chosen in cfg.spec.iid.spike_source.
%
%   pk = spec_get_spike_peaks(ev, r, seg, a, b, seizures, cfg)
%
%   'pipeline'     : <...>_IID_results.mat of module 2 (found by
%                    spec_resolve_files): iid_results.individual_peaks_table
%                    .Time_s inside [a, b]; exclusion zones =
%                    iid_results.exclusion_table (seizures + buffer, gaps,
%                    manual). Those peaks come from |band-pass 15-70 Hz|:
%                    they carry no polarity and are not the exact extremum
%                    of the analysed signal (spec_iid_interval re-aligns).
%                    detect_iid writes no .mat for a channel with zero
%                    spike complexes: then this is an error for the event.
%   'detect_local' : spec_detect_spikes_local.m on the interval.
%   'manual'       : ev.peaks (t_rel numbers, or times in ev.time_ref).
%   For 'detect_local' / 'manual' the exclusion zones are the seizures of
%   the same channel in cfg.spec.seizures_file (if given) +- exclusion_buffer_s.
%
% pk: t (column, t_rel), source (text), zones (table start_s, end_s,
%     reason), detect (threshold info or []), note

    p = cfg.spec.iid;
    pk = struct('t', zeros(0, 1), 'source', p.spike_source, 'zones', empty_zones(), 'detect', [], 'note', '');
    switch p.spike_source
        case 'pipeline'
            if isempty(r.iid_mat)
                error('spec_get_spike_peaks:NoIidMat', ['spike_source = ''pipeline'' but there is no %s_%s_%s_IID_results.mat ' ...
                    'under the subject folder (module 2 not run, or 0 spike complexes in that channel). ' ...
                    'Use spike_source = ''detect_local'' or ''manual'' for this event.'], r.subject, r.session_key, r.region);
            end
            S = load(r.iid_mat, 'iid_results');
            ir = S.iid_results;
            t = ir.individual_peaks_table.Time_s;
            pk.t = t(t >= a & t <= b);
            ex = ir.exclusion_table;
            if height(ex) > 0
                pk.zones = table(ex.start_s, ex.end_s, cellstr(string(ex.reason)), ...
                    'VariableNames', {'start_s', 'end_s', 'reason'});
            end
            pk.source = ['pipeline: ' r.iid_mat];
        case 'detect_local'
            pk.zones = seizure_zones(seizures, r, p.exclusion_buffer_s);
            include = seg.valid & seg.t_rel >= a & seg.t_rel <= b & ~in_zones(seg.t_rel, pk.zones);
            [pk.t, pk.detect] = spec_detect_spikes_local(seg.x, seg.fs, seg.t_rel, include, p.detect);
            pk.source = sprintf('detect_local (%s, threshold %.1f uV)', p.detect.threshold_mode, pk.detect.threshold_uV);
        case 'manual'
            v = ev.peaks;
            if isempty(v)
                error('spec_get_spike_peaks:NoManualPeaks', 'spike_source = ''manual'' but event %s has no peaks.', ev.id);
            end
            if isnumeric(v)
                t = v(:);
            else
                v = cellstr(v);
                t = zeros(numel(v), 1);
                for k = 1:numel(v)
                    t(k) = spec_resolve_times('to_t_rel', v{k}, ev.time_ref, r);
                end
            end
            out = t < a | t > b;
            if any(out)
                pk.note = sprintf('%d manual peak(s) outside the interval ignored', sum(out));
            end
            pk.t = t(~out);
            pk.zones = seizure_zones(seizures, r, p.exclusion_buffer_s);
            pk.source = sprintf('manual (%d)', numel(pk.t));
        otherwise
            error('spec_get_spike_peaks:BadSource', 'spike_source must be pipeline, detect_local or manual (got ''%s'').', p.spike_source);
    end
    pk.t = sort(pk.t(:));
end

function Z = empty_zones()
    Z = table(zeros(0, 1), zeros(0, 1), cell(0, 1), 'VariableNames', {'start_s', 'end_s', 'reason'});
end

function Z = seizure_zones(seizures, r, buffer_s)
    Z = empty_zones();
    if isempty(seizures)
        return;
    end
    [~, stem] = fileparts(r.source_file);
    stems = cellfun(@stem_of, seizures.source_file, 'UniformOutput', false);
    m = strcmp(seizures.subject_id, r.subject) & strcmpi(stems, stem) & strcmp(seizures.region, r.region);
    if any(m)
        Z = table(seizures.start_s(m) - buffer_s, seizures.end_s(m) + buffer_s, repmat({'seizure'}, sum(m), 1), ...
            'VariableNames', {'start_s', 'end_s', 'reason'});
    end
end

function s = stem_of(f)
    [~, s] = fileparts(f);
end

function m = in_zones(t, Z)
    m = false(size(t));
    for k = 1:height(Z)
        m = m | (t >= Z.start_s(k) & t <= Z.end_s(k));
    end
end
