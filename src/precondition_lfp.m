function cond = precondition_lfp(data, q, case_spec, cfg)
% PRECONDITION_LFP  Applies GAIN only (notch happens later, inside
% clean_lfp.m -- see README.md for the order-of-operations rationale).
%
%   cond = precondition_lfp(data, q, case_spec, cfg)
%
%   data      : struct from load_lfp_txt.m (raw signal, in real uV)
%   q         : struct from signal_quality.m, for THIS data
%   case_spec : struct from resolve_case.m (gain_mode, gain, notch_mode --
%               notch_mode is carried through untouched for clean_lfp.m)
%   cfg       : struct from pipeline_config.m (uses cfg.precondition.*,
%               cfg.quality.reference_sigma_uV)
%
% OUTPUT (struct cond): a copy of data with .signal replaced (or NOT even
% reassigned when gain==1 exactly, to guarantee bit-identical output in
% that case -- see gain_mode='off' and the dead band below), plus
% gain_applied, gain_estimate_raw, gain_source ('off'|'explicit'|'auto'),
% reference_used, reference_source, sensitivity_equivalent_uV_per_mm
% (100/gain_applied -- e.g. gain 10 -> 10 uV/mm, gain 20 -> 5 uV/mm, the
% strongest external validation available: it's checked against the
% sensitivity that was actually needed in Natus), and .quality = q.

    cond = data;

    switch case_spec.gain_mode
        case 'off'
            gain_applied = 1;
            gain_estimate_raw = NaN;
            gain_source = 'off';
            reference_used = NaN;
            reference_source = 'off';

        case 'explicit'
            gain_applied = case_spec.gain;
            gain_estimate_raw = case_spec.gain;
            gain_source = 'explicit';
            reference_used = NaN;
            reference_source = 'off';
            [gain_applied, clamped] = clamp_gain(gain_applied, cfg.precondition.gain_max);
            if clamped
                warning('precondition_lfp:GainClamped', ...
                    '%s: explicit gain %.3g clamped to cfg.precondition.gain_max = %.3g.', ...
                    data.file, gain_estimate_raw, cfg.precondition.gain_max);
            end

        case 'auto'
            region = field_or(data.meta, 'region', '');
            [reference_used, reference_source] = resolve_reference_value(region, cfg.quality.reference_sigma_uV);
            if isnan(reference_used)
                error('precondition_lfp:NoReference', ...
                    ['%s (region "%s") is declared as a gain case but cfg.quality.reference_sigma_uV has no ' ...
                     'value for this region. Run estimate_reference_sigma.m and set cfg.quality.reference_sigma_uV first.'], ...
                    data.file, region);
            end

            gain_estimate_raw = reference_used / q.sigma_band_uV;
            gain_source = 'auto';

            deadband = cfg.precondition.gain_deadband;
            if gain_estimate_raw >= 1/deadband && gain_estimate_raw <= deadband
                gain_applied = 1;
            else
                [gain_applied, clamped] = clamp_gain(gain_estimate_raw, cfg.precondition.gain_max);
                if clamped
                    warning('precondition_lfp:GainClamped', ...
                        '%s: auto gain estimate %.3g clamped to cfg.precondition.gain_max = %.3g.', ...
                        data.file, gain_estimate_raw, cfg.precondition.gain_max);
                end
            end

        otherwise
            error('precondition_lfp:BadGainMode', 'case_spec.gain_mode must be ''off'', ''explicit'', or ''auto'' (got ''%s'').', case_spec.gain_mode);
    end

    if gain_applied == 1
        cond.signal = data.signal;  % no multiplication at all -- bit-identical passthrough
    else
        cond.signal = data.signal * gain_applied;
    end

    cond.gain_applied = gain_applied;
    cond.gain_estimate_raw = gain_estimate_raw;
    cond.gain_source = gain_source;
    cond.reference_used = reference_used;
    cond.reference_source = reference_source;
    cond.sensitivity_equivalent_uV_per_mm = 100 / gain_applied;
    cond.quality = q;
end

function [gain_out, clamped] = clamp_gain(gain_in, gain_max)
    clamped = gain_in > gain_max;
    if clamped
        gain_out = gain_max;
    else
        gain_out = gain_in;
    end
end

function v = field_or(s, name, default)
    if isfield(s, name)
        v = s.(name);
    else
        v = default;
    end
end

function [ref_val, ref_source] = resolve_reference_value(region, cfg_ref)
    if isa(cfg_ref, 'containers.Map')
        if isKey(cfg_ref, region)
            ref_val = cfg_ref(region);
            ref_source = 'config_map';
            return;
        end
        ref_val = NaN;
        ref_source = 'none';
        return;
    end
    if isnumeric(cfg_ref) && isscalar(cfg_ref) && ~isnan(cfg_ref)
        ref_val = cfg_ref;
        ref_source = 'config_scalar';
        return;
    end
    ref_val = NaN;
    ref_source = 'none';
end
