function [ll_accept, ll_reject, band_on] = resolve_ll_band(cfg)
% RESOLVE_LL_BAND  The robust branch's ll_ratio review band.
%
%   [ll_accept, ll_reject, band_on] = resolve_ll_band(cfg)
%
% band_on = false (and ll_accept = ll_reject = ll_threshold) when either
% cfg.seizure_robust.ll_accept or .ll_reject is missing, [] or NaN: the
% band machinery is then OFF and every output is byte-identical to the
% binary ll_threshold pipeline. band_on = true otherwise -- including
% ll_accept == ll_reject (zero-width band), which reproduces the binary
% classification while still producing ll_status / event categories.
% Errors if ll_reject > ll_accept.

    sr = cfg.seizure_robust;
    a = field_or_empty(sr, 'll_accept');
    r = field_or_empty(sr, 'll_reject');
    if isempty(a) || isempty(r) || isnan(a) || isnan(r)
        ll_accept = sr.ll_threshold;
        ll_reject = sr.ll_threshold;
        band_on = false;
        return;
    end
    if ~(isnumeric(a) && isscalar(a) && isnumeric(r) && isscalar(r))
        error('resolve_ll_band:BadValue', 'cfg.seizure_robust.ll_accept / ll_reject must be numeric scalars.');
    end
    if r > a
        error('resolve_ll_band:BadBand', ...
            'cfg.seizure_robust.ll_reject (%.4g) must be <= cfg.seizure_robust.ll_accept (%.4g).', r, a);
    end
    ll_accept = a;
    ll_reject = r;
    band_on = true;
end

function v = field_or_empty(s, name)
    if isfield(s, name)
        v = s.(name);
    else
        v = [];
    end
end
