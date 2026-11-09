function y = rp_edf_phys(x, s)
% RP_EDF_PHYS  Digital EDF samples -> physical values, as edfread does:
% pmin + (x - dmin) * (pmax - pmin) / (dmax - dmin). s = one element of
% rp_edf_read's E.sig.
    y = s.pmin + (double(x) - s.dmin) * (s.pmax - s.pmin) / (s.dmax - s.dmin);
end
