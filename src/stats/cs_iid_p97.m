function q = cs_iid_p97(data, exclusion_table, cap_uV)
% CS_IID_P97  The real 97th percentile of |signal| over the samples
% detect_iid included -- the value its quantized baseline search estimates
% on the cfg.iid.baseline grid (100:5:130 uV, so the baseline saturates at
% 125 uV whenever P97 >= 130 uV).
%
% The inclusion mask is rebuilt as detect_iid's build_exclusion_zones
% does: valid samples minus every zone of its exclusion_table, both bounds
% inclusive on t_rel (boundary samples located with the same >= / <=
% comparisons on data.t_rel, so the mask is exact).
%
%   q.p97_uV             97th percentile of |signal| (uV), included samples
%   q.frac_below_cap     fraction of included samples with |signal| < cap_uV
%                        (detect_iid's grid top; >= 0.97 means no saturation)
%   q.analyzed_min       sum(mask) / fs / 60 (must equal detect_iid's
%                        analyzed_duration_min)

    t = data.t_rel(:);
    n = numel(t);
    fs = data.fs;
    mask = data.valid_mask(:);
    for i = 1:height(exclusion_table)
        i1 = first_ge(t, exclusion_table.start_s(i), fs, n);
        i2 = last_le(t, exclusion_table.end_s(i), fs, n);
        if i1 <= i2
            mask(i1:i2) = false;
        end
    end
    x = abs(data.signal(mask));
    q.p97_uV = prctile(x, 97);
    q.frac_below_cap = mean(x < cap_uV);
    q.analyzed_min = sum(mask) / fs / 60;
end

function k = first_ge(t, s, fs, n)
% first index with t(k) >= s (n+1 if none)
    k = min(max(floor(s * fs) - 1, 1), n);
    while k > 1 && t(k - 1) >= s
        k = k - 1;
    end
    while k <= n && t(k) < s
        k = k + 1;
    end
end

function k = last_le(t, e, fs, n)
% last index with t(k) <= e (0 if none)
    k = min(max(floor(e * fs) + 3, 1), n);
    while k < n && t(k + 1) <= e
        k = k + 1;
    end
    while k >= 1 && t(k) > e
        k = k - 1;
    end
end
