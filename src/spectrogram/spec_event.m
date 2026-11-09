function ev = spec_event(id, tipo, subject, sesion, region, inicio, fin, time_ref, overrides, origin)
% SPEC_EVENT  One requested event in the module's canonical form (what
% spec_read_event_list.m and spec_select_seizures.m produce and spec_run.m
% consumes).
%
%   id       : text, unique within subject + session (output file name)
%   tipo     : 'seizure' | 'iid'
%   subject  : exact ('005-s' ~= '005')
%   sesion   : 'yyyyMMdd_HHmmss' or the EDF file name
%   region   : 'HPCr' | 'HPCl' | 'both'
%   inicio, fin : t_rel number, datetime, or text 'dd/MM/yyyy HH:mm:ss.SSS'
%   time_ref : 'natus' | 'pipeline' | 't_rel'
%   overrides: struct of per-event cfg.spec overrides (see spec_run.m), or []
%   origin   : where the event came from (for the preview / log)
%
% Extra fields filled by the selectors: peaks (manual spike times), sz
% (seizures_events rows behind a detected seizure, one per region).

    if nargin < 9 || isempty(overrides)
        overrides = struct();
    end
    if nargin < 10
        origin = 'manual';
    end
    tipo = lower(strtrim(char(tipo)));
    if ~any(strcmp(tipo, {'seizure', 'iid'}))
        error('spec_event:BadTipo', 'Event %s: tipo must be ''seizure'' or ''iid'' (got ''%s'').', char(id), tipo);
    end
    ev = struct('id', strtrim(char(id)), 'tipo', tipo, 'subject', strtrim(char(subject)), ...
        'sesion', strtrim(char(sesion)), 'region', strtrim(char(region)), 'inicio', {inicio}, 'fin', {fin}, ...
        'time_ref', lower(strtrim(char(time_ref))), 'overrides', overrides, 'origin', char(origin), ...
        'peaks', [], 'sz', struct([]));
    if isempty(ev.id)
        error('spec_event:NoId', 'Every event needs an id.');
    end
end
