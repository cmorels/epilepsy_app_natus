function tf = cs_in_iv(t, iv)
% CS_IN_IV  Which points t fall inside the union of the intervals iv
% (Nx2 [start end], any order or overlap). One discretize over the
% interval edges, so it is fast for many points.

    tf = false(size(t));
    iv = cs_iv('union', iv);
    if isempty(iv) || isempty(t)
        return;
    end
    b = discretize(t, reshape(iv', 1, []));
    tf = mod(b, 2) == 1;   % odd bins lie inside an interval; NaN (outside every edge) -> false
end
