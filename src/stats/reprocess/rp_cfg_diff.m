function d = rp_cfg_diff(a, b, ignore)
% RP_CFG_DIFF  Dotted paths where two pipeline cfg structs differ
% (containers.Map compared by keys and values). ignore: paths not
% compared (default {'paths.output_root'}).

    if nargin < 3, ignore = {'paths.output_root'}; end
    fa = flat(a, ''); fb = flat(b, '');
    k = union(fieldnames(fa), fieldnames(fb));
    d = {};
    for i = 1:numel(k)
        p = strrep(k{i}, '__', '.');
        if ismember(p, ignore), continue; end
        if ~isfield(fa, k{i}) || ~isfield(fb, k{i}) || ~isequaln(fa.(k{i}), fb.(k{i}))
            d{end+1} = p; %#ok<AGROW>
        end
    end
end

function f = flat(v, pre)
    f = struct();
    if isstruct(v) && isscalar(v)
        for n = fieldnames(v)'
            g = flat(v.(n{1}), [pre '__' n{1}]);
            for m = fieldnames(g)', f.(m{1}) = g.(m{1}); end
        end
    else
        if isa(v, 'containers.Map')
            v = [keys(v); values(v)];
        end
        f.(regexprep(pre, '^__', '')) = v;
    end
end
