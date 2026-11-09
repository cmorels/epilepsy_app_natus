function [A, stop, flags] = cs_phaseA_checks(cfg, edf_tbl, ev, rv, S, LC, RD)
% CS_PHASEA_CHECKS  Section 3 (phase A) of the v2 task: snapshot vs v1,
% 084-s completeness, 048-s 20260622, 128-s 20260729 (exclusion, DX-21,
% channel A2 in recording_log, procedure for a future single-EDF run) and
% the channels to re-detect. Returns one row per check with whether it
% triggers a STOP condition of the task.
%
%   A    : table {punto, resultado, detalle, stop}
%   stop : any STOP condition met
%   flags: {severity, topic, detail}

    tz = cfg.tz;
    loc = @(p) char(datetime(round(p), 'ConvertFrom', 'posixtime', 'TimeZone', tz, 'Format', 'yyyy-MM-dd HH:mm:ss'));   % display only
    rows = cell(0, 4);
    flags = cell(0, 3);

    % ---- A.1 snapshot ----------------------------------------------------
    used = LC(LC.read_by_analysis, :);
    ch = used(used.changed_vs_v1, :);
    d = strjoin(arrayfun(@(i) sprintf('%s %s (%d filas; EDF: %s; modificado %s)', ch.campaign{i}, ch.mouse{i}, ch.n_rows_changed(i), ...
        ch.edfs_changed{i}, ch.modified{i}), (1:height(ch))', 'UniformOutput', false), ' || ');
    if isempty(ch)
        recent = find(datetime(used.modified, 'InputFormat', 'yyyy-MM-dd HH:mm:ss') > datetime(2026, 10, 28));
        d = sprintf(['ninguno: ningun archivo cambio desde que los leyo la ultima corrida de v1 (%d modificados despues de ' ...
            'su inicio). Modificados desde el 28/10: %s'], sum(used.modified_after_v1), strjoin(arrayfun(@(i) sprintf('%s %s %s', ...
            used.campaign{i}, used.mouse{i}, used.modified{i}), recent', 'UniformOutput', false), '; '));
    end
    rows(end+1, :) = {'A.1 snapshot de etiquetas', sprintf('%d archivos copiados (%s); %d con AE/AF distinto de v1', height(LC), ...
        LC.snapshot_taken{1}, height(ch)), d, false};

    % ---- A.2 084-s C2 -------------------------------------------------------
    x = rv.excluded(strcmp(rv.excluded.mouse, '084-s') & strcmp(rv.excluded.kind, 'sin AE completo'), :);
    E = edf_tbl(strcmp(edf_tbl.campaign, 'C2') & strcmp(edf_tbl.mouse, '084-s'), :);
    det = strjoin(arrayfun(@(i) sprintf('%s: %d eventos, %d sin AE', E.edf{i}, ...
        sum(strcmp(ev.campaign, 'C2') & strcmp(ev.mouse, '084-s') & strcmp(ev.edf, E.edf{i})), ...
        sum(strcmp(ev.campaign, 'C2') & strcmp(ev.mouse, '084-s') & strcmp(ev.edf, E.edf{i}) & cellfun(@isempty, ev.label))), ...
        (1:height(E))', 'UniformOutput', false), ' | ');
    u = S.units(strcmp(S.units.campaign, 'C2') & strcmp(S.units.mouse, '084-s'), :);
    rows(end+1, :) = {'A.2 084-s C2: todos los eventos con AE', sprintf('%d de %d EDF con eventos sin AE -> excluidos; unidad %s (incluida=%d)', ...
        height(x), height(E), u.status{1}, u.included(1)), [det ' (la crisette O-13 no requiere AE)'], false};

    % ---- A.3 048-s 20260622 ---------------------------------------------
    e0 = '20260622~ TRAP_fd13b728-0df8-45ad-8f8c-2b79a894d3ef';
    j0 = find(strcmp(edf_tbl.mouse, '048-s') & strcmp(edf_tbl.edf, e0));
    t1 = posixtime(datetime(2026, 6, 22, 0, 0, 0, 'TimeZone', tz));
    t2 = posixtime(datetime(2026, 6, 25, 0, 0, 0, 'TimeZone', tz));
    W = sortrows(edf_tbl(strcmp(edf_tbl.mouse, '048-s') & edf_tbl.start_posix < t2 & edf_tbl.end_posix > t1, :), 'start_posix');
    cover = strjoin(arrayfun(@(i) sprintf('%s %s -> %s', W.edf{i}, loc(W.start_posix(i)), loc(W.end_posix(i))), ...
        (1:height(W))', 'UniformOutput', false), ' | ');
    q1 = load(fullfile(cfg.v1_dir, 'phase1_results.mat'), 'O', 'S');
    q2 = load(fullfile(cfg.v1_dir, 'phase2_results.mat'), 'Z');
    o = q1.O.mouse(strcmp({q1.O.mouse.mouse}, '048-s') & strcmp({q1.O.mouse.scope}, 'merged'));
    span = [edf_tbl.start_posix(j0) edf_tbl.end_posix(j0)];
    h_rec = cs_iv('length', cs_iv('intersect', o.rec, span)) / 3600;
    h_sz = cs_iv('length', cs_iv('intersect', o.T_sz, span)) / 3600;
    n_sz = sum(strcmp(q1.S.seizures.mouse, '048-s') & contains(q1.S.seizures.edf, e0));
    pz = q2.Z.channels(strcmp(q2.Z.channels.mouse, '048-s') & strcmp(q2.Z.channels.edf, e0), :);
    ev0 = q1.S.events(strcmp(q1.S.events.mouse, '048-s') & strcmp(q1.S.events.edf, e0), :);
    green = W(W.start_posix >= span(2), :);
    stop3 = ~isempty(green);
    rows(end+1, :) = {'A.3 048-s: EDF que cubre el tramo 22/06-24/06', ...
        sprintf('%s cubre solo %s -> %s (%.1f h); el tramo posterior lo cubren %d EDF distintos', e0, loc(span(1)), loc(span(2)), ...
        diff(span) / 3600, height(green)), ['EDF de 048-s entre 22/06 y 24/06: ' cover], stop3};
    rows(end+1, :) = {'A.3 048-s: exclusion de 20260622 en v1', sprintf(['v1 lo excluyo bien: %.2f h en el registro, %.2f h en T_seizure_obs, ' ...
        '%d crisis, %d eventos (%d TP) marcados descartados, IIDs: %d canales corridos'], h_rec, h_sz, n_sz, height(ev0), ...
        sum(strcmp(ev0.label, 'TP')), sum(pz.run)), ['El raster v1 dibujaba solo T_seizure_obs (verde) y desconexiones dentro del ' ...
        'registro (rojo); el EDF excluido no quedaba en el registro, asi que no se dibujaba nada en su lugar. Lo verde de 22/06 ' ...
        '15:01 -> 24/06 18:14 son 20260622a~ TRA_51144d20, 20260622a~ TRA_0ca6ac7f y 20260623b~ TRA_e6bb2f2c. Impacto en v1: ninguno.'], false};

    % ---- A.4 128-s 20260729 ---------------------------------------------
    e7 = '20260729~ TRAP_d5154696-610f-41eb-91af-d50e5658bcd3';
    x7 = rv.excluded(strcmp(rv.excluded.mouse, '128-s') & strcmp(rv.excluded.edf, e7), :);
    dx21 = rv.dx(rv.dx.id == "DX-21", :);
    rp = isfield(cfg, 'include_128s_0729_reprocessed') && cfg.include_128s_0729_reprocessed;
    if (isempty(x7) && ~rp) || isempty(dx21)
        error('cs_phaseA_checks:128s', '128-s %s is not excluded or DX-21 is missing', e7);
    end
    if rp   % reprocessed A1 / A2 EDFs replace the original (run_campaign_stats.m, rp128_integrate.m)
        n = cellfun(@(s) sum(strcmp(ev.mouse, '128-s') & strcmp(ev.edf, s)), cfg.rp128.edfs);
        left = ~isempty(x7) || any(strcmp(ev.mouse, '128-s') & strcmp(ev.edf, e7));   % 083-s / 124-s share the EDF
        rows(end+1, :) = {'A.4a 128-s C2 20260729 reprocesado (A1 / A2)', sprintf('ED-17 reemplazado por %s', ...
            strjoin(cfg.rp128.edfs, ' + ')), sprintf('%s evento(s) detectado(s); el original de 128-s sigue presente: %d', ...
            mat2str(n), left), left};
    else
        ev7 = ev(strcmp(ev.mouse, '128-s') & strcmp(ev.edf, e7), :);
        rows(end+1, :) = {'A.4a 128-s C2 20260729 excluido de todo', sprintf('%s (%s): %s -> %s, %.1f h fuera', strjoin(x7.id', '/'), ...
            strjoin(x7.kind', '/'), loc(x7.start_posix(1)), loc(x7.end_posix(1)), (x7.end_posix(1) - x7.start_posix(1)) / 3600), ...
            sprintf('%d evento(s) detectado(s) dentro, AE=%s; sin tiempo, crisis, IIDs ni desempeno del detector', height(ev7), ...
            strjoin(unique(ev7.label)', ',')), false};
    end
    rows(end+1, :) = {'A.4b DX-21 registrado', sprintf('%s -> %s (%.0f s)', loc(dx21.start_posix), loc(dx21.end_posix), dx21.duration_s), ...
        'Queda en 2_Desconexiones / zonas R8; sin efecto mientras el EDF este excluido; se usara con el EDF nuevo', false};
    L = load_recording_log(cfg.recording_log);
    l7 = L(strcmp(regexprep(L.filename, '\.EDF$', '', 'ignorecase'), e7), :);
    a2 = l7(strcmp(l7.port, 'A2'), :);
    other = a2(~strcmp(a2.animal_id, '128-s'), :);
    rows(end+1, :) = {'A.4c recording_log: canal A2 en 20260729', ifelse(isempty(a2), 'A2 no asignado a ningun raton', ...
        sprintf('A2 asignado a %s', strjoin(a2.animal_id', ','))), ['Asignaciones del EDF: ' strjoin(strcat(l7.animal_id, '=', l7.port)', ', ')], ...
        ~isempty(other)};
    rows(end+1, :) = {'A.4d procedimiento futuro (no ejecutado)', 'ver detalle', procedure_text(cfg), false};

    % ---- A.5 channels to re-detect ----------------------------------------
    rd = RD(strcmp(RD.accion, 're-detectar'), :);
    ne = sum(~rd.explicado);
    [g, k] = findgroups(rd(:, {'mouse', 'motivo'}));
    summ = strjoin(arrayfun(@(i) sprintf('%s: %d (%s)', k.mouse{i}, sum(g == i), k.motivo{i}), (1:height(k))', 'UniformOutput', false), ' | ');
    gone = RD(strcmp(RD.accion, 'ya no se corre'), :);
    ru = erase(RD.accion(startsWith(RD.accion, 'reutilizar ')), 'reutilizar ');
    [u, ~, gu] = unique(ru);
    reuse = strjoin(arrayfun(@(i) sprintf('%s: %d', u{i}, sum(gu == i)), (1:numel(u))', 'UniformOutput', false), ', ');
    rows(end+1, :) = {'A.5 canales a re-detectar', sprintf(['%d corridas a re-detectar (%d sin explicacion), %d reutilizan una corrida ' ...
        'anterior (%s), %d de v1 ya no se corren'], height(rd), ne, numel(ru), reuse, height(gone)), [summ ' || ya no se corren: ' ...
        strjoin(strcat(gone.mouse, {' '}, gone.edf, {' '}, gone.region, {' ('}, gone.motivo, {')'})', '; ')], ne > 20};

    A = cell2table(rows, 'VariableNames', {'punto', 'resultado', 'detalle', 'stop'});
    stop = any(A.stop);
    for i = find(A.stop)'
        flags(end+1, :) = {'STOP', A.punto{i}, [A.resultado{i} '. ' A.detalle{i}]}; %#ok<AGROW>
    end
end

function s = ifelse(c, a, b)
    if c, s = a; else, s = b; end
end

function t = procedure_text(cfg)
    R = readtable(fullfile(cfg.campaigns(2).root, 'reference_per_region.csv'), 'TextType', 'char');
    ref = strjoin(arrayfun(@(i) sprintf('%s %.4f uV', R.region{i}, R.reference_uV(i)), (1:height(R))', 'UniformOutput', false), ', ');
    t = sprintf([ ...
        '1) Dejar el EDF combinado (A1 hasta 17:18:39 + A2 desde 17:23:55, decimado a 512 Hz como los demas) en EDF export\\decimated\\128-s\\. ' ...
        'Si conserva el nombre 20260729~ TRAP_d5154696..., apartar el original. ' ...
        '2) EEG_recording_log.xlsx: fila del EDF nuevo para 128-s con el Port y los canales HPCr/HPCl (C2/C4) que tenga el archivo combinado; ' ...
        'Attenuation vacio -> caso excel_clean (notch 50 Hz, sin ganancia). ' ...
        '3) Punto de entrada: un script como pipeline_output_campaign_20261014\\rerun_077s_campaign_20261014.m (pasos 3-4 de run_all_subjects.m): ' ...
        'cfg = campaign_config(); cfg.edf.channels.log_file = recording_log; la referencia de ganancia GUARDADA de la C2 ' ...
        '(pipeline_output_campaign_20261014\\reference_per_region.csv: %s) en cfg.quality.reference_sigma_uV (containers.Map por region) ' ...
        'y reference_sigma_fallback_uV. NO correr compute_run_reference.m: recalibraria con otro conjunto de canales y cambiaria la ' ...
        'ganancia respecto de la cohorte. cfg.edf.subject_id = ''128-s''; cfg.paths.output_root = una carpeta NUEVA (p. ej. ' ...
        'pipeline_output_campaign_20261014_128s_reproc\\128-s): no escribir en la C2 (con cfg.edf.reuse_import = true y el mismo nombre ' ...
        'de EDF se reutilizaria la importacion vieja de 01_txt). Llamar run_pipeline_edf({ruta del EDF nuevo}, cfg). ' ...
        '4) No correr merge_pipeline_runs sobre la C2 (reescribe su merged\\). ' ...
        '5) Etiquetar AE/AF del seizures_events nuevo y sumar esa carpeta como fuente en run_campaign_stats.m (reemplaza a ED-17); ' ...
        'DX-21 cubre el hueco. Los demas ratones no cambian: la referencia es la guardada y la deteccion de crisis es por canal y EDF.'], ref);
end
