% remove_outliers_epilepsy.m
% Detecta outliers técnicos (artefactos) sin eliminar eventos fisiológicos reales
%
% Estrategia:
% - Umbrales más permisivos (k_factor alto)
% - Solo corrige artefactos extremos (fuera de rango fisiológico)


%% 1. USER SETTINGS 
file_path = '';             % '' -> usa el primer .txt de la carpeta

% === PARÁMETROS ===
k_factor  = 10;              % Para crisis epilépticas: usar 15-20
                            % (valores típicos sin crisis: 5-10)
                           

use_fixed_thresholds = true; % true: usar umbrales fijos absolutos
                              % false: usar umbrales automáticos (MAD)
                              
                              % YO EVITARIA UTILIZAR UMBRALES AUTOMATICOS
                              % EN SENALES INTERICTALES (MUY AGRESIVO)

% Umbrales fijos basados en rango fisiológico de LFP hipocampal
fixed_lower_uV = -2500;     
fixed_upper_uV = +2500;    

% NOTA: Crisis epilépticas típicas en hipocampo:
% - Baseline normal: ±200 µV
% - Crisis moderada: ±500-1000 µV
% - Crisis severa: ±1000-2500 µV
% - Artefacto técnico: > ±2500 µV (mov. cable, etc.)

%% 2. LOAD DATA
data = load_LFP_intan_txt(file_path);
fs   = data.fs;             
signal_uV = data.signal;

%% 3. INFORMACIÓN DEL ARCHIVO
fprintf('\n=== Información del archivo ===\n');
fprintf('Archivo:                %s\n', data.file);
fprintf('Frecuencia de muestreo: %.2f Hz\n', fs);
fprintf('Duración:               %.3f segundos (%.1f min)\n', ...
    max(data.time_seconds), max(data.time_seconds)/60);
fprintf('Número de muestras:     %d\n', numel(signal_uV));

%% 4. ESTADÍSTICAS DE LA SEÑAL
baseline_uV = median(signal_uV, 'omitnan');
noise_uV    = mad(signal_uV, 1);

if noise_uV == 0
    warning('MAD = 0; usando un valor mínimo ficticio de ruido (1 µV).');
    noise_uV = 1;
end

% Umbrales automáticos (basados en MAD)
min_auto_uV = baseline_uV - k_factor * noise_uV;
max_auto_uV = baseline_uV + k_factor * noise_uV;

% Decidir qué umbrales usar
if use_fixed_thresholds
    min_allowed_uV = fixed_lower_uV;
    max_allowed_uV = fixed_upper_uV;
    threshold_type = 'FIXED';
else
    min_allowed_uV = min_auto_uV;
    max_allowed_uV = max_auto_uV;
    threshold_type = sprintf('AUTOMATIC (k=%.1f)', k_factor);
end

fprintf('\n=== Estadísticas de la señal ===\n');
fprintf('Baseline (mediana):     %.1f µV\n', baseline_uV);
fprintf('Ruido (MAD):            %.1f µV\n', noise_uV);
fprintf('Rango observado:        [%.1f, %.1f] µV\n', min(signal_uV), max(signal_uV));
fprintf('Rango dinámico:         %.1f µV\n', range(signal_uV));

% Estimar si hay posibles crisis
signal_std = std(signal_uV);
if signal_std > 200  % típico de señales con crisis
    fprintf('\n SEÑAL CON ALTA VARIABILIDAD (std = %.1f µV)\n', signal_std);
    fprintf('    Posibles crisis epilépticas presentes.\n');
    fprintf('    Usando umbrales permisivos para preservar eventos fisiológicos.\n');
end

fprintf('Umbrales de detección\n');
fprintf('Tipo:                   %s\n', threshold_type);
fprintf('Umbral INFERIOR:        %.1f µV\n', min_allowed_uV);
fprintf('Umbral SUPERIOR:        %.1f µV\n', max_allowed_uV);

if use_fixed_thresholds
    fprintf('\n Usando umbrales FIJOS.\n');
    fprintf('   Corrigiendo artefactos técnicos extremos.\n');
else
    fprintf('\nUmbral AUTO inferior:  %.1f µV (k=%.1f)\n', min_auto_uV, k_factor);
    fprintf('Umbral AUTO superior:  %.1f µV (k=%.1f)\n', max_auto_uV, k_factor);
end

%% 5. REMOVER OUTLIERS
[clean_signal_uV, outlier_info] = remove_outliers_amplitude( ...
    signal_uV, min_allowed_uV, max_allowed_uV);

num_outliers = size(outlier_info, 1);
fprintf('\n=== Resultados ===\n');
fprintf('Outliers detectados:    %d (%.4f%% de la señal)\n', ...
    num_outliers, 100 * num_outliers / numel(signal_uV));

if num_outliers > 0
    fprintf('Valores outliers:       [%.1f, %.1f] µV\n', ...
        min(outlier_info(:,2)), max(outlier_info(:,2)));
    
    % Advertencia si hay muchos outliers
    outlier_percent = 100 * num_outliers / numel(signal_uV);
    if outlier_percent > 1.0
        fprintf('\n  ADVERTENCIA: %.2f%% de outliers es alto.\n', outlier_percent);
        fprintf('    Considera:\n');
        fprintf('    - Aumentar k_factor (actualmente %.1f)\n', k_factor);
        fprintf('    - Verificar calidad de grabación\n');
        fprintf('    - Revisar umbrales fijos (min: %.0f, max: %.0f)\n', ...
            fixed_lower_uV, fixed_upper_uV);
     % Advertencia moderada: posibles interictal spikes
    
    elseif outlier_percent > 0.3
        fprintf('\n ADVERTENCIA: %.2f%% de outliers detectados.\n', outlier_percent);
        fprintf('    POSIBLE REMOCIÓN DE INTERICTAL SPIKES (IIS)\n');
        fprintf('    \n');
        fprintf('    En ratones con epilepsia crónica SIN crisis activas,\n');
        fprintf('    los "outliers" pueden ser IIS legítimos, NO artefactos.\n');
        fprintf('    \n');
        fprintf('      RECOMENDACIÓN para actividad interictal:\n');
        fprintf('       → Usar umbrales FIJOS:\n');
        fprintf('          Umbral bajo = -2500 uV\n');
        fprintf('          Umbral alto = +2500 uV\n');
        fprintf('       → O aumentar k_factor (valor actual=%.1f');
        fprintf('    \n');
        fprintf('    **  IIS típicos: 200-800 µV, duración <100 ms\n');
        fprintf('        Artefactos: >2500 µV, sin patrón fisiológico\n');
        
    elseif outlier_percent > 0.1
        % Info: pocos outliers (probablemente OK)
        fprintf('\n  %.2f%% de outliers (nivel aceptable).\n', outlier_percent);
        fprintf('    Probablemente solo artefactos técnicos aislados.\n');
    end
end


% Vector lógico para marcar outliers
is_outlier = false(size(signal_uV));
if num_outliers > 0
    idx = outlier_info(:,1);
    idx = idx(idx >= 1 & idx <= numel(signal_uV));
    is_outlier(idx) = true;
end

%% 6. GRAFICAR 
t = data.time_seconds;

fig = figure('Position', [50, 50, 1400, 900]);
fig.Name = sprintf('Outlier Removal: %s', data.file);

% --- Panel 1: Señal original con outliers y umbrales ---
subplot(3,1,1)
plot(t, signal_uV, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5); hold on

if num_outliers > 0
    plot(t(is_outlier), signal_uV(is_outlier), '.r', 'MarkerSize', 8);
end

yline(min_allowed_uV, '--r', sprintf('Min: %.0f µV', min_allowed_uV), ...
    'LabelHorizontalAlignment', 'left', 'LineWidth', 2);
yline(max_allowed_uV, '--r', sprintf('Max: %.0f µV', max_allowed_uV), ...
    'LabelHorizontalAlignment', 'left', 'LineWidth', 2);
yline(baseline_uV, '--k', sprintf('Baseline: %.0f µV', baseline_uV), ...
    'LabelHorizontalAlignment', 'right', 'LineWidth', 1);

% Marcar zona de crisis típica 
yline(1000, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
yline(-1000, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
text(max(t)*0.98, 1000, 'Typical seizure range', 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom', 'Color', [0.5 0.5 0.5]);
x_area = [min(t), max(t), max(t), min(t)];
y_area = [-1000, -1000, 1000, 1000];
fill(x_area, y_area, [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.3);

title_str = sprintf('Señal Original – %d outliers (%.4f%%)', ...
    num_outliers, 100*num_outliers/numel(signal_uV));
title(title_str, 'Interpreter', 'none')
xlabel('Tiempo (s)')
ylabel('Voltaje (µV)')
grid on
if num_outliers > 0
    legend({'Señal', 'Outliers'}, 'Location', 'best')
end
hold off

% --- Panel 2: Señal limpia con zoom en eventos grandes ---
subplot(3,1,2)
plot(t, clean_signal_uV, 'b', 'LineWidth', 0.5)
yline(baseline_uV, '--k', sprintf('Baseline: %.0f µV', baseline_uV), ...
    'LabelHorizontalAlignment', 'right', 'LineWidth', 1);
title('Señal Después de Remover Outliers')
xlabel('Tiempo (s)')
ylabel('Voltaje (µV)')
grid on

% --- Panel 3: Diferencia (solo artefactos corregidos) ---
subplot(3,1,3)
diff_signal = signal_uV - clean_signal_uV;
plot(t, diff_signal, 'r', 'LineWidth', 0.5)
yline(0, '--k', 'LineWidth', 1);
title('Diferencia (Original - Limpia)')
xlabel('Tiempo (s)')
ylabel('Diferencia (µV)')
grid on

% Si no hay diferencias, anotar
if num_outliers == 0
    text(max(t)/2, 0, 'Sin correcciones', ...
        'HorizontalAlignment', 'center', 'FontSize', 14, ...
        'Color', [0 0.5 0], 'FontWeight', 'bold');
end

% Guardar figura
try
    [~, fname, ~] = fileparts(data.file);
    fig_name = sprintf('%s_outliers.fig', fname);
    saveas(gcf, fig_name);
    fprintf('\nFigura guardada: %s\n', fig_name);
catch ME
    warning('No se pudo guardar la figura: %s', ME.message);
end

%% 7. ANÁLISIS ADICIONAL: DISTRIBUCIÓN DE AMPLITUDES
figure('Position', [100, 100, 1200, 400]);
figure_name = sprintf('Histograma: %s', data.file);
set(gcf, 'Name', figure_name);

subplot(1,2,1)
histogram(signal_uV, 100, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
hold on
xline(min_allowed_uV, '--r', 'Min', 'LineWidth', 2);
xline(max_allowed_uV, '--r', 'Max', 'LineWidth', 2);
xline(baseline_uV, '--k', 'Baseline', 'LineWidth', 1);
title('Distribución de Amplitudes Original')
xlabel('Amplitud (µV)')
ylabel('Frecuencia')
grid on
hold off

subplot(1,2,2)
histogram(clean_signal_uV, 100, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
hold on
xline(baseline_uV, '--k', 'Baseline', 'LineWidth', 1);
title('Distribución de Amplitudes Limpia')
xlabel('Amplitud (µV)')
ylabel('Frecuencia')
grid on
hold off

try
    [~, fname, ~] = fileparts(data.file);
    hist_name = sprintf('%s_histogram.fig', fname);
    saveas(gcf, hist_name);
    fprintf('Histograma guardado: %s\n', hist_name);
catch
    warning('No se pudo guardar el histograma.');
end

%% 8. GUARDAR SEÑAL LIMPIA AUTOMÁTICAMENTE
fprintf('\n=== Guardando señal limpia ===\n');

[~, fname, ~] = fileparts(data.file);
out_name = sprintf('%s_clean.txt', fname);

fid = fopen(out_name, 'w');
if fid == -1
    error('No se pudo crear el archivo: %s', out_name);
end

% Header
fprintf(fid, '# Original file: %s\n', data.file);
fprintf(fid, '# fs = %.10f\n', fs);
fprintf(fid, '# time_unit = seconds\n');
fprintf(fid, '# offset_applied = none (original signal)\n');
fprintf(fid, '# outliers_removed = %d samples (%.4f%%)\n', ...
    num_outliers, 100*num_outliers/numel(signal_uV));
fprintf(fid, '# removal_method = Amplitude thresholds\n');
fprintf(fid, '# k_factor = %.1f\n', k_factor);
fprintf(fid, '# threshold_type = %s\n', threshold_type);
fprintf(fid, '# lower_threshold = %.1f microvolts\n', min_allowed_uV);
fprintf(fid, '# upper_threshold = %.1f microvolts\n', max_allowed_uV);
fprintf(fid, '# signal_range_original = [%.2f, %.2f] microvolts\n', ...
    min(signal_uV), max(signal_uV));
fprintf(fid, '# signal_range_clean = [%.2f, %.2f] microvolts\n', ...
    min(clean_signal_uV), max(clean_signal_uV));
fprintf(fid, '# columns = amplitude_microvolts\n');

% Datos
fprintf(fid, '%.6f\n', clean_signal_uV);
fclose(fid);

fprintf(' Señal limpia guardada: %s\n', out_name);
fprintf('   (El archivo se sobrescribirá si ya existe)\n');
fprintf('\n Procesamiento completado.\n');

%% LOCAL FUNCTION
function [clean_signal, outlier_info] = remove_outliers_amplitude(signal, min_allowed, max_allowed)
    signal = signal(:);
    
    is_out = signal < min_allowed | signal > max_allowed | isnan(signal);
    
    out_idx   = find(is_out);
    out_vals  = signal(is_out);
    
    outlier_info = [out_idx, out_vals];
    clean_signal = signal;
    
    if isempty(out_idx)
        return;
    end
    
    good_idx  = find(~is_out);
    
    if numel(good_idx) < 2
        warning('Muy pocos puntos válidos para interpolar. Se deja la señal sin cambios.');
        return;
    end
    
    clean_signal(is_out) = interp1( ...
        good_idx, signal(good_idx), out_idx, 'linear', 'extrap');
end