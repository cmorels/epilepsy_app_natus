function [Z, label, used] = spec_normalize(P, mode, Pbl)
% SPEC_NORMALIZE  Normalize a power map (f x t, linear units).
%
%   [Z, label, used] = spec_normalize(P, mode, Pbl)
%
%   mode : 'db'           10 log10(P)
%          'db_baseline'  10 log10(P ./ mean(Pbl, 2))           (dB re baseline)
%          'zscore'       (10 log10 P - mean_bl) ./ std_bl, both of 10 log10(Pbl)
%          'percent'      100 (P - mean(Pbl)) ./ mean(Pbl)
%   Pbl  : baseline power columns (f x n_baseline), same rows as P; NaN
%          columns ignored. If it is empty / all NaN, a baseline mode falls
%          back to 'db' (used says which mode was applied).
%   label: colorbar label

    used = mode;
    if ~strcmp(mode, 'db')
        if nargin < 3 || isempty(Pbl) || all(isnan(Pbl(:)))
            used = 'db';
        end
    end
    switch used
        case 'db'
            Z = 10 * log10(P);
            label = 'power (dB re 1 uV^2)';
        case 'db_baseline'
            m = mean(Pbl, 2, 'omitnan');
            Z = 10 * log10(P ./ m);
            label = 'dB re baseline';
        case 'zscore'
            L = 10 * log10(Pbl);
            Z = (10 * log10(P) - mean(L, 2, 'omitnan')) ./ std(L, 0, 2, 'omitnan');
            label = 'z (dB re baseline)';
        case 'percent'
            m = mean(Pbl, 2, 'omitnan');
            Z = 100 * (P - m) ./ m;
            label = '% change re baseline';
        otherwise
            error('spec_normalize:BadMode', 'norm must be ''db'', ''db_baseline'', ''zscore'' or ''percent'' (got ''%s'').', mode);
    end
end
