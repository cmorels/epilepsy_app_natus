function [b, kind] = spec_baseline_interval(spec, is_t_rel, time_ref, r)
% SPEC_BASELINE_INTERVAL  Baseline setting -> 'relative' offsets or an
% 'absolute' [t1 t2] in t_rel.
%
%   [b, kind] = spec_baseline_interval(spec, is_t_rel, time_ref, r)
%
%   spec : [a b] numbers -> relative to the event anchor (onset / each
%          peak), unless is_t_rel (then absolute t_rel);
%          {'dd/MM/yyyy HH:mm:ss.SSS', '...'} -> absolute, in time_ref
%          ('natus' | 'pipeline')
%   r    : spec_resolve_files.m struct of the session

    if isnumeric(spec)
        if numel(spec) ~= 2 || spec(2) <= spec(1)
            error('spec_baseline_interval:Bad', 'baseline must be [a b] with b > a.');
        end
        b = spec(:)';
        if is_t_rel
            kind = 'absolute';
        else
            kind = 'relative';
        end
        return;
    end
    spec = cellstr(spec);
    if numel(spec) ~= 2
        error('spec_baseline_interval:Bad', 'An absolute baseline needs two times (start, end).');
    end
    b = [spec_resolve_times('to_t_rel', spec{1}, time_ref, r), spec_resolve_times('to_t_rel', spec{2}, time_ref, r)];
    if b(2) <= b(1)
        error('spec_baseline_interval:Bad', 'Baseline end must be after its start.');
    end
    kind = 'absolute';
end
