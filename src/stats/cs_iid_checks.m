function [chk, flags, CRI, fpw] = cs_iid_checks(cfg, Z, C, I, O)
% CS_IID_CHECKS  Section-8 checks of the IID re-detection, plus the
% crisette_focal count (IIDs that would fall inside those zones, from the
% runs without them).
%
%   chk  : {control, pasa, detalle}
%   flags: {severity, topic, detail}
%   CRI  : per crisette x channel: complexes inside the (buffered) zone in
%          the run without crisette zones, and in the main run (must be 0)
%   fpw  : per FP-IIDs window: re-detected complexes inside it (both channels)

    buf = cfg.iid_excl_buffer_s;
    chk = cell(0, 3);
    flags = cell(0, 3);
    Cm = C(strcmp(C.variant, 'principal'), :);
    Zz = Z.zones;

    % 1. no re-detected IID inside a zone (seizure, crisette, non-IIDs FP, disconnection)
    bad = 0; badch = {};
    for i = 1:height(Cm)
        z = Zz(strcmp(Zz.mouse, Cm.mouse{i}), :);
        iv = [z.start_posix - buf, z.end_posix + buf] - Cm.start_posix(i);
        cx = Cm.complexes{i};
        k = 0;
        for r = 1:size(iv, 1)
            k = k + sum(cx(:, 2) >= iv(r, 1) & cx(:, 1) <= iv(r, 2));
        end
        if k > 0
            bad = bad + k;
            badch{end+1} = sprintf('%s %s %s %s (%d)', Cm.campaign{i}, Cm.mouse{i}, Cm.edf{i}(1:10), Cm.region{i}, k); %#ok<AGROW>
        end
    end
    chk(end+1, :) = {'ningun IID re-detectado dentro de crisis final, crisette, FP no-IIDs o desconexion (+/- 5 s)', bad == 0, ...
        sprintf('%d complejos en zona; %s', bad, strjoin(badch, ', '))};

    % 2. FP-IIDs windows have IIDs in the re-detection
    W = Z.fp_iids;
    nW = zeros(height(W), 1);
    in_other = false(height(W), 1);
    for w = 1:height(W)
        cs = find(strcmp(Cm.mouse, W.mouse{w}));
        for c = cs'
            t = Cm.complexes{c}(:, 1) + Cm.start_posix(c);
            nW(w) = nW(w) + sum(t >= W.start_posix(w) & t <= W.end_posix(w));
        end
        z = Zz(strcmp(Zz.mouse, W.mouse{w}), :);
        in_other(w) = any(z.start_posix - buf < W.end_posix(w) & z.end_posix + buf > W.start_posix(w));
    end
    fpw = W(:, {'event_id', 'mouse', 'campaign', 'edf', 'seizure_id', 'duration_s'});
    fpw.start_local = arrayfun(@(p) char(datetime(p, 'ConvertFrom', 'posixtime', 'TimeZone', cfg.tz, ...
        'Format', 'yyyy-MM-dd HH:mm:ss')), W.start_posix, 'UniformOutput', false);
    fpw.n_complejos_redetectados = nW;
    fpw.solapa_otra_zona = in_other;
    z0 = nW == 0;
    chk(end+1, :) = {'las ventanas FP-IIDs tienen IIDs en la re-deteccion', ~any(z0), sprintf(['%d/%d ventanas con >= 1 complejo ' ...
        '(mediana %g); %d con 0 (%d de ellas solapan otra zona de exclusion)'], sum(~z0), numel(nW), median(nW), sum(z0), ...
        sum(z0 & in_other))};
    if any(z0)
        x = fpw(z0, :);
        flags(end+1, :) = {'revisar', 'ventanas FP-IIDs sin IIDs re-detectados', sprintf('%d: %s', height(x), ...
            strjoin(strcat(x.event_id, ' (', x.start_local, ')')', '; '))};
    end

    % 3. channels with disconnections: baseline must not rise
    ba = I.before_after;
    d = ba(ba.disconnected_h > 0, :);
    up = d(d.baseline_uV > d.baseline_uV_orig, :);
    chk(end+1, :) = {'canales con desconexiones: el baseline IID re-detectado no sube', isempty(up) || height(up) == 0, ...
        sprintf('%d canales con desconexiones; baseline sube en %d, baja en %d, igual en %d', height(d), height(up), ...
        sum(d.baseline_uV < d.baseline_uV_orig), sum(d.baseline_uV == d.baseline_uV_orig))};

    % 4. T_iid <= T_eeg <= recording, per mouse x EDF; 5. T_iid = T_eeg - IID zones (1 s)
    E = O.edf_day;
    [g, ek] = findgroups(E(:, {'mouse', 'campaign', 'edf'}));
    ek.rec_h = splitapply(@sum, E.rec_h, g);
    ek.T_eeg_h = splitapply(@sum, E.T_eeg_h, g);
    [tf, j] = ismember(strcat(Cm.mouse, '|', Cm.campaign, '|', Cm.edf), strcat(ek.mouse, '|', ek.campaign, '|', ek.edf));
    T_iid = Cm.analyzed_min / 60;
    T_eeg = NaN(height(Cm), 1); rec = NaN(height(Cm), 1);
    T_eeg(tf) = ek.T_eeg_h(j(tf));
    rec(tf) = ek.rec_h(j(tf));
    tol = 1 / 3600;
    ok4 = tf & T_iid <= T_eeg + tol & T_eeg <= rec + tol;
    chk(end+1, :) = {'T_iid <= T_eeg <= registro, por raton x EDF (cada canal)', all(ok4), ...
        sprintf('%d/%d canales; %d sin fila de T_eeg', sum(ok4), height(Cm), sum(~tf))};
    P = Z.channels;
    [tf2, j2] = ismember(strcat(Cm.campaign, '|', Cm.mouse, '|', Cm.edf, '|', Cm.region), ...
        strcat(P.campaign, '|', P.mouse, '|', P.edf, '|', P.region));
    dif_s = abs(T_iid(tf2) - P.expected_T_iid_h(j2(tf2))) * 3600;
    chk(end+1, :) = {'T_iid (analyzed_duration) = T_eeg - zonas IID, tolerancia 1 s por EDF', all(tf2) && all(dif_s <= 1), ...
        sprintf('max |diferencia| = %.3f s en %d canales', max(dif_s), numel(dif_s))};

    % crisette counts
    Cs = C(strcmp(C.variant, 'sin_crisette'), :);
    cr = Zz(strcmp(Zz.rule, 'crisette'), :);
    rows = {};
    for i = 1:height(Cs)
        mi = find(strcmp(Cm.campaign, Cs.campaign{i}) & strcmp(Cm.mouse, Cs.mouse{i}) & strcmp(Cm.edf, Cs.edf{i}) & ...
            strcmp(Cm.region, Cs.region{i}));
        z = cr(strcmp(cr.mouse, Cs.mouse{i}) & cr.end_posix + buf > Cs.start_posix(i) & ...
            cr.start_posix - buf < Cs.start_posix(i) + Cs.duration_s(i), :);
        for r = 1:height(z)
            iv = [z.start_posix(r) - buf, z.end_posix(r) + buf] - Cs.start_posix(i);
            n_sin = sum(Cs.complexes{i}(:, 2) >= iv(1) & Cs.complexes{i}(:, 1) <= iv(2));
            n_con = sum(Cm.complexes{mi}(:, 2) >= iv(1) & Cm.complexes{mi}(:, 1) <= iv(2));
            rows(end+1, :) = {z.origin{r}, Cs.mouse{i}, Cs.campaign{i}, Cs.edf{i}, Cs.region{i}, z.start_local{r}, z.duration_s(r), ...
                n_sin, n_con, Cs.baseline_uV(i), Cm.baseline_uV(mi), Cs.n_complexes(i) - Cm.n_complexes(mi)}; %#ok<AGROW>
        end
    end
    CRI = cell2table(rows, 'VariableNames', {'obs_id', 'mouse', 'campaign', 'edf', 'region', 'inicio', 'duracion_s', ...
        'n_complejos_en_zona_sin_excluir', 'n_complejos_en_zona_corrida_principal', 'baseline_sin_zona_uV', ...
        'baseline_con_zona_uV', 'diferencia_total_complejos'});
end
