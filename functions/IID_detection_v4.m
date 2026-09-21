%% --- 3.X IID DETECTION: Walsh pre-detector + Envelope confirmation ---

% Assumes clean_signal_uV (µV), t (s), fs (Hz) already exist

%% (A) Band-pass filter
bp_band = [15 70];              % Hz típico para IIDs
lfp_bp  = bandpass(clean_signal_uV, bp_band, fs, 'Steepness', 0.5);

trimN   = round(fs);            % recorta 1 s para bordes de filtro
lfp_bp  = lfp_bp(trimN+1:end-trimN);
t_bp    = t(trimN+1:end-trimN);
t_bp    = t_bp(:);

%% (B) Walsh pre-detector
walsh_raw = Walsh(lfp_bp(:));   % misma longitud que lfp_bp

% Suavizado corto para estabilidad (5–15 ms)
winW   = max(1, round(0.010 * fs));    % 10 ms
walsh  = movmean(walsh_raw, winW);

% Umbral Walsh: mediana + kW * MAD  (por grabación)
kW     = 1;                               % más bajo = más sensible
medW   = median(walsh);
madW   = mad(walsh, 1);                   % MAD (L1)
thrW   = medW + kW * madW;

is_walsh_high = walsh > thrW;

%% (C) Dilatar candidatos Walsh para no cortar inicio/fin del IID
dilate_sec  = 0.03;                       % 30 ms
dilate_samp = max(1, round(dilate_sec*fs));

is_walsh_cand = is_walsh_high;
if any(is_walsh_high)
    % dilatación binaria 1D simple
    kern = true(1, 2*dilate_samp+1);
    is_walsh_cand = conv(double(is_walsh_high), double(kern), 'same') > 0;
end

%% (D) Envolvente (Hilbert) + global threshold basado en controles

% Hilbert envelope
xa    = hilbert(lfp_bp);
env   = abs(xa);

% Smoothing window (debe coincidir con build_env_threshold_from_controls)
winE_sec = 0.020;                      
winE     = max(1, round(winE_sec * fs));
env_s    = movmean(env, winE);

% Cargar stats globales de controles
load('envelope_ctrl_threshold.mat', ...
     'medE_ctrl_global', 'madE_ctrl_global', 'bp_band', 'winE_sec');

kE   = 13;     % factor de umbral basado en controles
thrE = medE_ctrl_global + kE * madE_ctrl_global;

is_env_high = env_s > thrE;

%% (E) Candidatos finales = intersección (precisión) o unión (sensibilidad)
use_intersection = true;
if use_intersection
    iid_mask = is_walsh_cand & is_env_high;
else
    iid_mask = is_walsh_cand | is_env_high;
end

%% (F) Post-procesado: duración mínima + fusión de huecos
min_dur_ms   = 20;                         % duración mínima IID
merge_gap_ms = 80;                         % fusionar huecos cortos

min_samp     = max(1, round(min_dur_ms/1000 * fs));
merge_gap    = max(1, round(merge_gap_ms/1000 * fs));

d  = diff([0; iid_mask(:); 0]);
st = find(d==1);
en = find(d==-1)-1;

% eliminar eventos muy cortos
dur = en - st + 1;
keep = dur >= min_samp;
st   = st(keep);
en   = en(keep);

% fusionar intervalos separados por ≤ merge_gap
if ~isempty(st)
    st2 = st(1); en2 = en(1);
    ST = []; EN = [];
    for k = 2:numel(st)
        if st(k) - en2 - 1 <= merge_gap
            en2 = en(k);                 % extend
        else
            ST(end+1) = st2; EN(end+1) = en2; %#ok<AGROW>
            st2 = st(k); en2 = en(k);
        end
    end
    ST(end+1) = st2; EN(end+1) = en2;     % último
else
    ST = []; EN = [];
end

n_iid = numel(ST);

% --- Contar cuántos mini-segmentos (st/en) hay dentro de cada IID fusionado ---
if ~isempty(st)
    st2 = st(1); en2 = en(1);
    ST = []; EN = [];
    for k = 2:numel(st)
        if st(k) - en2 - 1 <= merge_gap
            en2 = en(k);                 % extend
        else
            ST(end+1) = st2; EN(end+1) = en2; %#ok<AGROW>
            st2 = st(k); en2 = en(k);
        end
    end
    ST(end+1) = st2; EN(end+1) = en2;     % último
else
    ST = []; EN = [];
end

n_iid = numel(ST);

%% === SINGLE vs POLYSPIKE (conteo de picos por IID, sin findpeaks) ===

num_peaks = zeros(n_iid,1);

for m = 1:n_iid

    % segmento de señal banda-pasada correspondiente a este IID
    seg = lfp_bp(ST(m):EN(m));
    seg_abs = abs(seg(:));            % columna, valor absoluto

    if isempty(seg_abs) || max(seg_abs) <= 0
        num_peaks(m) = 0;
        continue;
    end

    % --- parámetros de "qué es un pico" (ajustables) ---
    min_pk_dist = round(0.010 * fs);    % 10 ms entre picos
    min_pk_height   = 0.4 * max(seg_abs);   % ≥ 40% del máximo del IID

    x = seg_abs;

    % 1) máximos locales simples: x(i) > x(i-1) y x(i) >= x(i+1)
    cand = find( x(2:end-1) >  x(1:end-2) & ...
                 x(2:end-1) >= x(3:end) ) + 1;

    % 2) filtrar por altura mínima
    cand = cand( x(cand) >= min_pk_height );

    % 3) imponer distancia mínima entre picos
    if isempty(cand)
        num_peaks(m) = 0;
    else
        keep = true(size(cand));
        last = cand(1);
        for j = 2:numel(cand)
            if cand(j) - last < min_pk_dist
                keep(j) = false;          % demasiado cerca del pico anterior
            else
                last = cand(j);
            end
        end
        cand = cand(keep);
        num_peaks(m) = numel(cand);
    end
end

Type   = repmat("single", n_iid, 1);
Type(num_peaks >= 2) = "polyspike";
IsPoly = num_peaks >= 2;


%% (G) Máscaras lógicas para stats
mask_env_only   = is_env_high & ~is_walsh_cand;
mask_walsh_only = is_walsh_cand & ~is_env_high;
mask_both       = is_walsh_cand &  is_env_high;

%% (H) Tabla de resultados + impresión por pantalla
if n_iid > 0
    iid_times = table( (1:n_iid)', ...
        t_bp(ST(:)), t_bp(EN(:)), ...
        (EN(:)-ST(:)+1)/fs, ...
        Type, num_peaks, ...
        'VariableNames', {'Id','Start_s','End_s','Duration_s', ...
                          'Type','num_peaks'});
else
    iid_times = table([],[],[],[],[],[],[], ...
        'VariableNames', {'Id','Start_s','End_s','Duration_s', ...
                          'Type','num_peaks'});
end

disp(iid_times);

%% (I) Guardar resumen en Excel (1 fila por ejecución/animal)
perc_env_only   = 100 * mean(mask_env_only);
perc_walsh_only = 100 * mean(mask_walsh_only);
perc_both       = 100 * mean(mask_both);

% Si tienes una etiqueta de animal/grupo, ponla aquí:
mouse_label = "unknown";   % e.g. "control", "epileptic", "mouse01"

summary = table( ...
    mouse_label, ...
    n_iid, ...
    bp_band(1), bp_band(2), ...
    kE, kW, ...
    perc_env_only, perc_walsh_only, perc_both, ...
    'VariableNames', { ...
        'MouseLabel', ...
        'Total_IIDs', ...
        'BP_low_Hz', 'BP_high_Hz', ...
        'kE', 'kW', ...
        'Perc_env_only', ...
        'Perc_walsh_only', ...
        'Perc_both' ...
    } ...
);

out_file = 'IID_detection_summary.xlsx';

if isfile(out_file)
    writetable(summary, out_file, 'WriteMode', 'append');
else
    writetable(summary, out_file);
end

%% (J) Plot compacto: LFP band-passed + marcadores IID
figure; tiledlayout(1,1);

nexttile;
plot(t_bp, lfp_bp);
ylabel('\muV');
title('Band-passed LFP');
hold on;

% Asteriscos de inicio/fin
for k = 1:n_iid
    idx_start = ST(k);
    idx_end   = EN(k);
    plot(t_bp(idx_start), lfp_bp(idx_start), '*g', 'MarkerSize', 8); % inicio
    plot(t_bp(idx_end),   lfp_bp(idx_end),   '*r', 'MarkerSize', 8); % fin
end

% Texto con parámetros
mask_env_only   = is_env_high & ~is_walsh_cand;
mask_walsh_only = is_walsh_cand & ~is_env_high;
mask_both       = is_walsh_cand &  is_env_high;

fprintf('\n===== IID DETECTION SUMMARY =====\n');
fprintf('Total IIDs detected : %d\n', n_iid);
fprintf('Band-pass (Hz)      : [%.1f – %.1f]\n', bp_band(1), bp_band(2));
fprintf('kE (envelope thr)   : %.2f\n', kE);
fprintf('kW (Walsh thr)      : %.2f\n', kW);
fprintf('env only   : %.3f %% of samples\n', 100*mean(mask_env_only));
fprintf('walsh only : %.3f %% of samples\n', 100*mean(mask_walsh_only));
fprintf('both       : %.3f %% of samples\n', 100*mean(mask_both));
param_str = sprintf('BP [%.0f–%.0f] Hz | kE = %.1f | kW = %.1f | IIDs = %d', ...
                    bp_band(1), bp_band(2), kE, kW, n_iid);

text(0.01, 0.97, param_str, ...
     'Units','normalized', ...
     'HorizontalAlignment','left', ...
     'VerticalAlignment','top', ...
     'FontSize',12, ...
     'BackgroundColor','w', ...
     'Margin',2);

hold off;

%% WALSH FUNCTION (local function at end of script)

function walsh_metric = Walsh(x)
% WALSH  Walsh-based burstiness metric.
%   - Conv 'same' para mantener misma longitud que x
%   - Usa 3 escalas (4, 8, 16 puntos) de Hadamard/Walsh

    x = x(:);  % columna

    % Hadamard matrices
    H4  = hadamard(4);
    H8  = hadamard(8);
    H16 = hadamard(16);

    % Filas seleccionadas (patrones ±1)
    w4_1  = H4(3,:);
    w8_1  = H8(5,:);
    w16_1 = H16(9,:);

    % Convoluciones centradas (misma longitud)
    W4  = conv(x,  w4_1 , 'same');
    W8  = conv(x,  w8_1 , 'same');
    W16 = conv(x, w16_1 , 'same');

    % Métrica de burstiness (multi-escala)
    walsh_metric = abs(W4 + W8 + W16);
end
