function report = verify_notch(raw_txt, clean_txt, t_start_s, t_end_s)
% VERIFY_NOTCH  Diagnostic only: checks whether the 50 Hz notch was really
% applied, by comparing the raw txt (01_txt/) against the clean txt
% (*_clean.txt) over one time window (e.g. a seizure +/- 15 s).
%
%   report = verify_notch(raw_txt, clean_txt, t_start_s, t_end_s)
%
% Self-contained (MATLAB base + Signal Processing Toolbox): reads only the
% header and the requested window of each txt, never the pipeline modules.
% Figure saved to tools/verification/<clean name>_<t0>-<t1>s.fig/.png.

    report = struct();
    fprintf('\n==== verify_notch ====\nRAW   : %s\nCLEAN : %s\nWindow: %.2f - %.2f s\n', ...
        raw_txt, clean_txt, t_start_s, t_end_s);
    for p = {raw_txt, clean_txt}
        if exist(p{1}, 'file') ~= 2
            error('verify_notch:FileNotFound', 'File not found: %s', p{1});
        end
    end

    %% 0. File dates (suspect 1: stale clean file)
    dr = dir(raw_txt); dc = dir(clean_txt);
    fprintf('\nFile dates:\n  raw   %s\n  clean %s\n', datestr(dr.datenum), datestr(dc.datenum));
    if dc.datenum < dr.datenum
        fprintf('  WARNING: clean is OLDER than raw -> it may come from a previous run.\n');
    end

    %% 1. Headers
    [hr, ~] = read_header(raw_txt);
    [hc, nhc] = read_header(clean_txt);
    fs = str2double(hr.fs);
    if str2double(hc.fs) ~= fs
        error('verify_notch:FsMismatch', 'fs differs: raw %s vs clean %s', hr.fs, hc.fs);
    end

    fields = {'case_applied', 'case_source', 'notch_mode', 'notch_applied', 'notch_freqs', ...
        'notch_halfwidth_hz', 'notch_order', 'notch_blocks_skipped', 'line_ratio_db', ...
        'line_ratio_p95_db', 'pct_time_line_high', 'quality_class', 'gain_applied'};
    fprintf('\nClean header (%s):\n', clean_txt);
    missing = {};
    for i = 1:numel(fields)
        if isfield(hc, fields{i})
            fprintf('  %-22s = %s\n', fields{i}, hc.(fields{i}));
        else
            fprintf('  %-22s = <MISSING>\n', fields{i});
            missing{end+1} = fields{i}; %#ok<AGROW>
        end
    end
    if ~isempty(missing)
        fprintf('  FINDING: clean header lacks %s -- clean_lfp traceability incomplete.\n', strjoin(missing, ', '));
    end
    gain = str2double(field_or(hc, 'gain_applied', 'NaN'));
    header_says_notch = strcmpi(field_or(hc, 'notch_applied', ''), 'true');

    %% Load window
    i0 = max(1, floor(t_start_s * fs) + 1);
    i1 = floor(t_end_s * fs) + 1;
    n = i1 - i0 + 1;
    [~, nhr] = read_header(raw_txt);
    xr = read_window(raw_txt, nhr, i0, n);
    xc = read_window(clean_txt, nhc, i0, n);
    t = (i0 - 1 + (0:n-1)') / fs;
    if numel(xr) ~= n || numel(xc) ~= n
        error('verify_notch:ShortRead', 'Window exceeds file length (read %d/%d raw, %d/%d clean).', ...
            numel(xr), n, numel(xc), n);
    end
    ok = ~isnan(xr) & ~isnan(xc);
    fprintf('\nWindow: %d samples at %g Hz (%d NaN). Raw samples beyond +/-2500 uV (outlier-cleaned): %d\n', ...
        n, fs, nnz(~ok), nnz(abs(xr) > 2500));
    xr0 = fillmissing(xr, 'linear'); xc0 = fillmissing(xc, 'linear');

    %% 2. Spectra
    win = round(4 * fs);
    [Pr, f] = pwelch(xr0 - mean(xr0), hann(win), round(win/2), [], fs);
    Pc = pwelch(xc0 - mean(xc0), hann(win), round(win/2), [], fs);

    prom = @(P, f0) 10*log10(mean(P(f >= f0-1 & f <= f0+1)) / ...
        median(P((f >= f0-5 & f <= f0-2) | (f >= f0+2 & f <= f0+5))));
    pr50 = prom(Pr, 50); pc50 = prom(Pc, 50); att50 = pr50 - pc50;
    if att50 < 10
        v50 = 'NO APLICADO';
    elseif att50 > 25
        v50 = 'APLICADO';
    else
        v50 = 'INSUFICIENTE';
    end

    harm = [100 150 200];
    harm = harm(harm < fs/2 - 5);
    harm_att = arrayfun(@(h) prom(Pr, h) - prom(Pc, h), harm);
    harm_filtered = harm_att > 10;

    %% 3. Analysis band and phase
    ratio_db = 10*log10(Pc ./ Pr);
    in_band = (f >= 5 & f <= 45) | (f >= 55 & f <= 75);
    dev_lo = max(abs(ratio_db(f >= 5 & f <= 45)));
    dev_hi = max(abs(ratio_db(f >= 55 & f <= 75)));
    dev = max(abs(ratio_db(in_band)));

    br = bandpass(xr0, [5 40], fs); bc = bandpass(xc0, [5 40], fs);
    maxlag = round(0.5 * fs);
    [xcf, lags] = xcorr(br - mean(br), bc - mean(bc), maxlag, 'coeff');
    [pk, k] = max(xcf);
    lag = lags(k);

    %% 4. Figure
    out_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'tools', 'verification');
    if ~isfolder(out_dir), mkdir(out_dir); end
    [~, cname] = fileparts(clean_txt);
    tag = sprintf('%s_%.0f-%.0fs', cname, t_start_s, t_end_s);

    fig = figure('Visible', 'off', 'Position', [50 50 1400 1500]);
    tl = tiledlayout(fig, 5, 1, 'TileSpacing', 'compact');
    title(tl, sprintf('verify\\_notch: %s', strrep(tag, '_', '\_')));
    ylim_all = [min([xr0; xc0]) max([xr0; xc0])];

    ax1 = nexttile(tl); plot(ax1, t, xr0, 'Color', [0.8 0.2 0.2]); ylim(ax1, ylim_all);
    title(ax1, sprintf('CRUDA  (std %.1f \\muV, rango [%.0f, %.0f] \\muV)', std(xr0), min(xr0), max(xr0)));
    ylabel(ax1, '\muV');
    ax2 = nexttile(tl); plot(ax2, t, xc0, 'Color', [0.1 0.3 0.8]); ylim(ax2, ylim_all);
    title(ax2, sprintf('LIMPIA  (std %.1f \\muV, rango [%.0f, %.0f] \\muV)', std(xc0), min(xc0), max(xc0)));
    ylabel(ax2, '\muV');
    ax3 = nexttile(tl); hold(ax3, 'on');
    plot(ax3, t, xr0, 'Color', [0.8 0.2 0.2]); plot(ax3, t, xc0, 'Color', [0.1 0.3 0.8]);
    legend(ax3, {'cruda', 'limpia'}); title(ax3, 'Superpuestas'); ylabel(ax3, '\muV'); xlabel(ax3, 't (s)');
    linkaxes([ax1 ax2 ax3], 'x'); xlim(ax1, [t(1) t(end)]);

    ax4 = nexttile(tl);
    loglog(ax4, f, Pr, 'Color', [0.8 0.2 0.2]); hold(ax4, 'on');
    loglog(ax4, f, Pc, 'Color', [0.1 0.3 0.8]);
    xlim(ax4, [1 fs/2]);
    for fx = [50 harm], xline(ax4, fx, ':k'); end
    legend(ax4, {sprintf('cruda (50 Hz: %+.1f dB)', pr50), sprintf('limpia (50 Hz: %+.1f dB)', pc50)}, 'Location', 'southwest');
    title(ax4, sprintf('PSD (pwelch 4 s, 50%% solape) -- atenuacion 50 Hz = %.1f dB', att50));
    xlabel(ax4, 'Hz'); ylabel(ax4, '\muV^2/Hz');

    ax5 = nexttile(tl);
    tz0 = (t_start_s + t_end_s) / 2 - 0.25;
    z = t >= tz0 & t <= tz0 + 0.5;
    plot(ax5, t(z), xr0(z), '.-', 'Color', [0.8 0.2 0.2]); hold(ax5, 'on');
    plot(ax5, t(z), xc0(z), '.-', 'Color', [0.1 0.3 0.8]);
    legend(ax5, {'cruda', 'limpia'}); xlabel(ax5, 't (s)'); ylabel(ax5, '\muV');
    title(ax5, sprintf('Zoom 0.5 s (%.2f-%.2f s): 50 Hz = 25 ciclos en esta ventana', tz0, tz0 + 0.5));

    savefig(fig, fullfile(out_dir, [tag '.fig']));
    exportgraphics(fig, fullfile(out_dir, [tag '.png']), 'Resolution', 110);
    close(fig);

    %% 5. Verdict
    coherent = header_says_notch == ~strcmp(v50, 'NO APLICADO');
    fprintf('\n---- VEREDICTO ----\n');
    fprintf('NOTCH 50 Hz : de %.1f dB a %.1f dB (atenuacion %.1f dB)  -> %s\n', pr50, pc50, att50, v50);
    fprintf('ARMONICOS   : ');
    for i = 1:numel(harm)
        fprintf('%d Hz %s (%.1f dB)  ', harm(i), yesno(harm_filtered(i)), harm_att(i));
    end
    fprintf('(esperado: NO)\n');
    fprintf('BANDA UTIL  : desviacion maxima fuera de 45-55 Hz = %.2f dB (5-45: %.2f, 55-75: %.2f)  -> %s\n', ...
        dev, dev_lo, dev_hi, okbad(dev < 1, 'OK', 'DEFORMADA'));
    fprintf('FASE        : retardo del pico de correlacion = %d muestras (r=%.4f)  -> %s\n', ...
        lag, pk, okbad(lag == 0, 'OK', 'DESPLAZADA'));
    fprintf('GANANCIA    : gain_applied = %g (esperado 1)%s\n', gain, okbad(gain == 1, '', '  <-- NO ESPERADO'));
    fprintf('CABECERA vs SENAL : %s (cabecera notch_applied=%s)\n', ...
        okbad(coherent, 'COHERENTE', 'INCONSISTENTE'), field_or(hc, 'notch_applied', '<missing>'));
    fprintf('Figura: %s.png/.fig\n', fullfile(out_dir, tag));

    report = struct('prom_raw_db', pr50, 'prom_clean_db', pc50, 'att_db', att50, 'verdict', v50, ...
        'harm_hz', harm, 'harm_att_db', harm_att, 'band_dev_db', dev, 'lag_samples', lag, ...
        'gain_applied', gain, 'header_coherent', coherent, 'missing_fields', {missing}, ...
        'figure', fullfile(out_dir, [tag '.png']));
end

function [h, nlines] = read_header(path)
    fid = fopen(path, 'rt');
    c = onCleanup(@() fclose(fid));
    h = struct(); nlines = 0;
    while true
        pos = ftell(fid);
        ln = fgetl(fid);
        if ~ischar(ln), break; end
        if ~startsWith(ln, '#')
            fseek(fid, pos, 'bof');
            break;
        end
        nlines = nlines + 1;
        tok = regexp(ln, '^#\s*([^=]+?)\s*=\s?(.*)$', 'tokens', 'once');
        if ~isempty(tok)
            h.(matlab.lang.makeValidName(tok{1})) = strtrim(tok{2});
        end
    end
end

function x = read_window(path, n_header, i0, n)
    fid = fopen(path, 'rt');
    c = onCleanup(@() fclose(fid));
    C = textscan(fid, '%f', n, 'HeaderLines', n_header + i0 - 1, 'TreatAsEmpty', {'NaN', 'nan'});
    x = C{1};
end

function v = field_or(s, name, default)
    if isfield(s, name), v = s.(name); else, v = default; end
end

function s = yesno(tf)
    if tf, s = 'SI'; else, s = 'NO'; end
end

function s = okbad(tf, a, b)
    if tf, s = a; else, s = b; end
end
