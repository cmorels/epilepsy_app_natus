function [fex, txt] = spec_notch_band(r, notch_applied, extra_hz)
% SPEC_NOTCH_BAND  Frequencies removed by the pipeline's mains notch, to
% keep them out of peak / dominant-frequency searches: the notch dip is
% deeper in the baseline than in a seizure, so a ratio can "peak" there.
%
%   [fex, txt] = spec_notch_band(r, notch_applied, extra_hz)
%
%   r             : spec_resolve_files struct (clean header: notch_freqs,
%                   notch_halfwidth_hz)
%   notch_applied : whether the analysed signal carries the notch
%   extra_hz      : added on each side (e.g. half the multitaper resolution)
%   fex           : @(f) logical, true where f is within
%                   notch_freq +- (halfwidth + 1 + extra_hz)
%   txt           : e.g. '46.5-53.5 Hz (50 Hz notch)', '' when none

    fex = @(f) false(size(f));
    txt = '';
    if ~notch_applied || isempty(r.clean_hdr) || ~isfield(r.clean_hdr.meta, 'notch_freqs')
        return;
    end
    nf = sscanf(char(r.clean_hdr.meta.notch_freqs), '%f');
    hw = 1;
    if isfield(r.clean_hdr.meta, 'notch_halfwidth_hz')
        v = str2double(r.clean_hdr.meta.notch_halfwidth_hz);
        if ~isnan(v)
            hw = v;
        end
    end
    if isempty(nf)
        return;
    end
    half = hw + 1 + extra_hz;     % +1 Hz: the band-stop's skirts
    lo = nf(:)' - half;
    hi = nf(:)' + half;
    fex = @(f) any(f(:) >= lo & f(:) <= hi, 2);
    parts = arrayfun(@(a, b, c) sprintf('%.1f-%.1f Hz (%g Hz notch)', a, b, c), lo, hi, nf(:)', 'UniformOutput', false);
    txt = strjoin(parts, '; ');
end
