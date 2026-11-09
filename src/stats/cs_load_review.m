function rv = cs_load_review(cfg, edf_tbl)
% CS_LOAD_REVIEW  Sheets 1_Crisis_no_detectadas (ND), 2_Desconexiones (DX)
% and 3_EDF_estado (ED) of the review workbook (header row 4). Columns are
% taken by position and their header text is asserted. Workbook datetimes
% are local wall-clock (typed from Natus) and get TimeZone cfg.tz.
%
% DX open limits ('inicio del registro', 'fin del registro', 'EDF completo')
% are filled from that mouse's EDF span (*_clean.txt header: session start,
% start + n_samples/fs). The named EDF is used ONLY for that; intervals are
% later matched to recordings by mouse + absolute time (DX-08 is named
% 20260807 but lies inside 20260806).

    f = cfg.review_xlsx;
    nd = readtable(f, 'Sheet', '1_Crisis_no_detectadas', 'Range', 'A4', 'VariableNamingRule', 'preserve');
    dx = readtable(f, 'Sheet', '2_Desconexiones', 'Range', 'A4', 'VariableNamingRule', 'preserve');
    ed = readtable(f, 'Sheet', '3_EDF_estado', 'Range', 'A4', 'VariableNamingRule', 'preserve');
    assert_hdr(nd, [1 2 4 6 7 9 10 11], {'ID', 'Rat', 'EDF', 'Inicio', 'Fin', 'Tipo', 'Analizable', 'Incluir'});
    assert_hdr(dx, [1 2 4 5 6 7 8 15], {'ID', 'Rat', 'EDF', 'mite inicio', 'Inicio', 'Fin', 'mite fin', 'Video revisado'});

    nd = nd(startsWith(string(nd{:, 1}), 'ND-'), :);
    rv.nd = table(string(nd{:, 1}), string(nd{:, 2}), string(nd{:, 4}), to_posix(nd{:, 6}, cfg.tz), ...
        to_posix(nd{:, 7}, cfg.tz), string(nd{:, 9}), string(nd{:, 10}), string(nd{:, 11}), ...
        'VariableNames', {'id', 'mouse', 'edf', 'start_posix', 'end_posix', 'tipo', 'analizable', 'incluir'});

    dx = dx(startsWith(string(dx{:, 1}), 'DX-'), :);
    n = height(dx);
    st = to_posix(dx{:, 6}, cfg.tz);
    en = to_posix(dx{:, 7}, cfg.tz);
    lim_s = lower(strtrim(string(dx{:, 5})));
    lim_e = lower(strtrim(string(dx{:, 8})));
    filled = strings(n, 1);
    for i = 1:n
        mouse = string(dx{i, 2});
        edf = string(dx{i, 4});
        j = find(strcmp(edf_tbl.mouse, mouse) & strcmp(edf_tbl.edf, edf), 1);
        need = lim_s(i) ~= "valor" || lim_e(i) ~= "valor";
        if need && isempty(j)
            error('cs_load_review:NoEdf', '%s: EDF %s of %s not found among analyzed EDFs', dx{i, 1}{1}, edf, mouse);
        end
        if lim_s(i) == "inicio del registro" || lim_s(i) == "edf completo"
            st(i) = edf_tbl.start_posix(j);
            filled(i) = filled(i) + "start=EDF start; ";
        end
        if lim_e(i) == "fin del registro" || lim_e(i) == "edf completo"
            en(i) = edf_tbl.end_posix(j);
            filled(i) = filled(i) + "end=EDF end; ";
        end
    end
    video = startsWith(lower(strtrim(string(dx{:, 15}))), "s");
    rv.dx = table(string(dx{:, 1}), string(dx{:, 2}), string(dx{:, 4}), st, en, en - st, lim_s, lim_e, ...
        video, strtrim(filled), lim_s == "edf completo" & lim_e == "edf completo", ...
        'VariableNames', {'id', 'mouse', 'edf', 'start_posix', 'end_posix', 'duration_s', 'limit_start', ...
        'limit_end', 'video_reviewed', 'filled', 'whole_edf'});
    if any(~isfinite(rv.dx.start_posix) | ~isfinite(rv.dx.end_posix) | rv.dx.duration_s <= 0)
        error('cs_load_review:BadDx', 'A DX interval is incomplete after filling open limits.');
    end

    ed = ed(startsWith(string(ed{:, 1}), 'ED-'), :);
    rv.ed = table(string(ed{:, 1}), string(ed{:, 2}), string(ed{:, 4}), string(ed{:, 6}), ...
        'VariableNames', {'id', 'mouse', 'edf', 'estado'});
end

function p = to_posix(x, tz)
    if isdatetime(x)
        x.TimeZone = tz;
        p = posixtime(x);
    else
        p = NaN(numel(x), 1);
    end
    p = p(:);
end

function assert_hdr(T, idx, frag)
    names = T.Properties.VariableNames;
    for k = 1:numel(idx)
        if ~contains(names{idx(k)}, frag{k})
            error('cs_load_review:Header', 'Column %d is "%s", expected to contain "%s"', idx(k), names{idx(k)}, frag{k});
        end
    end
end
