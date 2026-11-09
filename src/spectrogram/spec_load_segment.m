function seg = spec_load_segment(r, t0, t1, cfg, cache_dir)
% SPEC_LOAD_SEGMENT  Samples of one channel between t_rel t0 and t1 (s),
% from the binary cache (spec_signal_cache.m), never the whole recording.
%
%   seg = spec_load_segment(r, t0, t1, cfg, cache_dir)
%
%   r   : struct from spec_resolve_files.m
%   cfg : cfg from spectrogram_config.m (cfg.spec.*)
%
% Samples outside the recording come back NaN. The signal is in uV:
%   'clean' : the *_clean.txt as written by clean_lfp.m (gain, outlier
%             interpolation and notch already applied -- not re-applied).
%   'raw'   : the 01_txt * gain_applied (clean header); then the same notch
%             (apply_notch_blocks.m with the clean header's notch_freqs /
%             halfwidth / order) only if the clean header says
%             notch_applied = true. No outlier interpolation.
%
% Outliers: when cfg.spec.count_outliers and the 01_txt exists, counts the
% samples of [t0, t1] whose raw * gain falls outside the clean header's
% [lower_threshold_uV, upper_threshold_uV], i.e. those clean_lfp.m
% interpolated. Without the 01_txt the count is NaN (unknown).
%
% seg: x, t_rel, i0, i1, fs, valid, source, gain_applied, notch_applied,
%      case_applied, n_outliers_segment, outlier_note, txt, load_info

    s = cfg.spec;
    fs = r.fs;
    i0 = floor(t0 * fs) + 1;
    i1 = ceil(t1 * fs) + 1;

    main = spec_signal_cache('open', r.txt, cache_dir, s.cache_chunk);
    x = spec_signal_cache('read', main, i0, i1);

    ch = r.clean_hdr;
    gain = 1;
    notch_on = false;
    case_applied = '';
    lower_uV = NaN; upper_uV = NaN;
    if ~isempty(ch)
        gain = num_field(ch.meta, 'gain_applied', 1);
        notch_on = strcmpi(text_field(ch.meta, 'notch_applied'), 'true');
        case_applied = text_field(ch.meta, 'case_applied');
        lower_uV = num_field(ch.meta, 'lower_threshold_uV', NaN);
        upper_uV = num_field(ch.meta, 'upper_threshold_uV', NaN);
    end

    seg = struct();
    if strcmp(s.signal_source, 'raw')
        x = x * gain;
        if notch_on
            ncfg = struct();
            ncfg.precondition.notch_freqs = sscanf(text_field(ch.meta, 'notch_freqs'), '%f')';
            ncfg.precondition.notch_harmonics = [];
            ncfg.precondition.notch_halfwidth_hz = num_field(ch.meta, 'notch_halfwidth_hz', 1);
            ncfg.precondition.notch_order = num_field(ch.meta, 'notch_order', 4);
            ncfg.seizure.bandpass_band = [0 fs / 2];   % let every configured notch frequency through
            ncfg.iid.bandpass_band = [0 fs / 2];
            x = apply_notch_blocks(x, ~isnan(x), fs, ncfg);
        end
    end

    n_out = NaN;
    out_note = '';
    if s.count_outliers
        if ~isempty(r.raw) && ~isnan(lower_uV) && ~isnan(upper_uV)
            if strcmp(s.signal_source, 'raw')
                xr = spec_signal_cache('read', main, i0, i1) * gain;   % before notch, as clean_lfp compares
            else
                rawinfo = spec_signal_cache('open', r.raw, cache_dir, s.cache_chunk);
                xr = spec_signal_cache('read', rawinfo, i0, i1) * gain;
            end
            n_out = sum(xr < lower_uV | xr > upper_uV);
            if n_out > 0
                out_note = sprintf(['%d sample(s) of this window are outside [%g, %g] uV and were linearly ' ...
                    'interpolated by clean_lfp'], n_out, lower_uV, upper_uV);
                if strcmp(s.signal_source, 'clean')
                    out_note = [out_note ' -- consider signal_source = ''raw'''];
                end
            end
        elseif isempty(r.raw)
            n_hdr = NaN;
            if ~isempty(ch)
                n_hdr = num_field(ch.meta, 'n_outliers', NaN);
            end
            out_note = sprintf(['outliers in this window unknown (no 01_txt); the clean header reports %g ' ...
                'interpolated sample(s) in the whole session'], n_hdr);
        end
    end

    seg.x = x;
    seg.i0 = i0;
    seg.i1 = i1;
    seg.t_rel = ((i0:i1)' - 1) / fs;
    seg.fs = fs;
    seg.valid = ~isnan(x);
    seg.source = s.signal_source;
    seg.gain_applied = gain;
    seg.notch_applied = notch_on;
    seg.case_applied = case_applied;
    seg.n_outliers_segment = n_out;
    seg.outlier_note = out_note;
    seg.txt = r.txt;
    seg.load_info = rmfield(main, {'header_pairs'});
end

function v = num_field(meta, name, default)
    v = default;
    if isfield(meta, name)
        t = str2double(meta.(name));
        if ~isnan(t)
            v = t;
        end
    end
end

function s = text_field(meta, name)
    s = '';
    if isfield(meta, name)
        s = strtrim(char(meta.(name)));
    end
end
