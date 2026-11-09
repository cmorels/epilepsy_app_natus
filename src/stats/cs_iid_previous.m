function [rr, src, why] = cs_iid_previous(dirs, ch, sz, man, prov)
% CS_IID_PREVIOUS  An earlier IID re-detection of this channel that can be
% reused: same zone INTERVALS (seizures-argument zones and manual zones,
% compared as sorted [start end] rows with ==; zone ids are labels and may
% be renumbered between runs) and same code (SHA-1 of detect_iid.m,
% cs_iid_run_channel.m, cs_iid_p97.m).
%
%   dirs : output roots searched in order (this run's first, then earlier
%          runs, e.g. v1), each laid out <root>/<mouse>/<EDF>/<base>_IID_rerun.mat
%   ch   : one channel row (mouse, edf, file)
%   sz, man : the zones cs_iid_channel_zones.m gives for this run
%   prov : this run's provenance (cs_provenance.m)
%   rr   : the stored result ([] if none reusable); src = its file
%   why  : 'igual', or why nothing was reused (first root that had a file)

    rr = [];
    src = '';
    why = 'sin resultado previo';
    [~, nm] = fileparts(ch.file{1});
    base = [strrep(nm, '_clean', '') '_IID_rerun.mat'];
    a_sz = sortrows(reshape(sz{:, {'start_s', 'end_s'}}, [], 2));
    a_man = sortrows(reshape(man, [], 2));
    first_why = '';
    for d = 1:numel(dirs)
        f = fullfile(dirs{d}, ch.mouse{1}, ch.edf{1}, base);
        if ~isfile(f)
            continue;
        end
        q = load(f, 'rr');
        b_sz = sortrows(reshape(q.rr.seizure_zones{:, {'start_s', 'end_s'}}, [], 2));
        b_man = sortrows(reshape(q.rr.manual_zones, [], 2));
        if ~isfield(q.rr, 'p97_uV') || ~strcmp(code_key(q.rr.provenance), code_key(prov))
            w = 'codigo distinto';
        elseif isequal(a_sz, b_sz) && isequal(a_man, b_man)
            rr = q.rr;
            src = f;
            why = 'igual';
            return;
        else
            w = sprintf('zonas distintas (crisis +%d/-%d, manuales +%d/-%d)', ...
                sum(~ismember(a_sz, b_sz, 'rows')), sum(~ismember(b_sz, a_sz, 'rows')), ...
                sum(~ismember(a_man, b_man, 'rows')), sum(~ismember(b_man, a_man, 'rows')));
        end
        if isempty(first_why)
            first_why = w;
        end
    end
    if ~isempty(first_why)
        why = first_why;
    end
end

function k = code_key(prov)
    k = [prov.iid_code{strcmp(prov.iid_code(:, 1), 'detect_iid.m'), 2}, ...
        prov.stats_code{strcmp(prov.stats_code(:, 1), 'cs_iid_run_channel.m'), 2}, ...
        prov.stats_code{strcmp(prov.stats_code(:, 1), 'cs_iid_p97.m'), 2}];
end
