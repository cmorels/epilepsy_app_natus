function out = cs_iv(op, a, b)
% CS_IV  Interval arithmetic on Nx2 [start end] matrices (any numeric time
% base, here POSIX seconds). Every result is sorted and non-overlapping.
%
%   U = cs_iv('union', A)          merge overlapping/touching intervals
%   U = cs_iv('union', A, B)       union of A and B
%   D = cs_iv('subtract', A, B)    A minus B
%   I = cs_iv('intersect', A, B)   A intersect B
%   L = cs_iv('length', A)         total length after union (scalar)
%   F = cs_iv('overlap', A, B)     per-row overlap length of A(i,:) with union(B)

    if nargin < 3
        b = zeros(0, 2);
    end
    a = clean(a);
    b = clean(b);
    switch op
        case 'union'
            out = merge([a; b]);
        case 'subtract'
            out = subtract(merge(a), merge(b));
        case 'intersect'
            out = intersect_iv(merge(a), merge(b));
        case 'length'
            u = merge(a);
            out = sum(u(:, 2) - u(:, 1));
        case 'overlap'
            ub = merge(b);
            out = zeros(size(a, 1), 1);
            for i = 1:size(a, 1)
                x = intersect_iv(a(i, :), ub);
                out(i) = sum(x(:, 2) - x(:, 1));
            end
        otherwise
            error('cs_iv:BadOp', 'Unknown op "%s".', op);
    end
end

function a = clean(a)
    if isempty(a)
        a = zeros(0, 2);
        return;
    end
    a = double(a);
    a = a(all(isfinite(a), 2) & a(:, 2) > a(:, 1), :);
end

function u = merge(a)
    if isempty(a)
        u = zeros(0, 2);
        return;
    end
    a = sortrows(a, 1);
    u = a(1, :);
    for i = 2:size(a, 1)
        if a(i, 1) <= u(end, 2)
            u(end, 2) = max(u(end, 2), a(i, 2));
        else
            u(end+1, :) = a(i, :); %#ok<AGROW>
        end
    end
end

function d = subtract(a, b)
    d = zeros(0, 2);
    for i = 1:size(a, 1)
        segs = a(i, :);
        for j = 1:size(b, 1)
            next = zeros(0, 2);
            for k = 1:size(segs, 1)
                s = segs(k, :);
                if b(j, 2) <= s(1) || b(j, 1) >= s(2)
                    next(end+1, :) = s; %#ok<AGROW>
                    continue;
                end
                if b(j, 1) > s(1)
                    next(end+1, :) = [s(1) b(j, 1)]; %#ok<AGROW>
                end
                if b(j, 2) < s(2)
                    next(end+1, :) = [b(j, 2) s(2)]; %#ok<AGROW>
                end
            end
            segs = next;
            if isempty(segs)
                break;
            end
        end
        d = [d; segs]; %#ok<AGROW>
    end
end

function x = intersect_iv(a, b)
    x = zeros(0, 2);
    i = 1; j = 1;
    while i <= size(a, 1) && j <= size(b, 1)
        lo = max(a(i, 1), b(j, 1));
        hi = min(a(i, 2), b(j, 2));
        if hi > lo
            x(end+1, :) = [lo hi]; %#ok<AGROW>
        end
        if a(i, 2) < b(j, 2)
            i = i + 1;
        else
            j = j + 1;
        end
    end
end
