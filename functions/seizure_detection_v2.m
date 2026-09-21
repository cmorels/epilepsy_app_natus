% seizure_detection.m
% Detección de crisis epilépticas con hysteresis + gap-fill + merge
% VERSIÓN MEJORADA: Carga automáticamente archivo *_clean_epilepsy.txt

clear; clc;

%% 1. CARGAR ARCHIVO LIMPIO
file_path = '';  % '' -> busca automáticamente el primer *_clean.txt

if isempty(file_path)
    files = dir('*_clean.txt');
    if isempty(files)
        error('No se encontró archivo limpio (*_clean.txt). Ejecuta remove_outliers primero.');
    end
    file_path = files(1).name;
    fprintf('Archivo seleccionado automáticamente: %s\n', file_path);
else
    fprintf('Cargando archivo especificado: %s\n', file_path);
end

% Cargar datos
data = load_LFP_intan_txt(file_path);
lfp_clean = data.signal;     % µV
fs = data.fs;
t = data.time_seconds;

fprintf('Duración: %.2f s (%.2f min)\n', max(t), max(t)/60);
fprintf('Fs: %.2f Hz\n\n', fs);

%% 2. PARÁMETROS DE DETECCIÓN
bandpass_band      = [5 75];     % Hz
window_sec_energy  = 2;          % moving-average window (s)
min_seizure_sec    = 15;         % minimum duration to keep (s)

% Robustness knobs
enter_factor = 8;      % median-energy multiplier to ENTER a seizure
exit_ratio   = 0.50;   % fraction of enter threshold to EXIT a seizure (hysteresis)
gap_sec      = 2;      % fill holes inside seizures up to this length (s)
merge_sec    = 2;      % merge adjacent segments separated by <= this gap (s)

%% 3. NORMALIZAR Y FILTRAR
% Normalize 0-1 (like the app) and band-pass 5-75 Hz
min_val   = min(lfp_clean);
max_val   = max(lfp_clean);
lfp_norm  = (lfp_clean - min_val) / (max_val - min_val);

lfp_bp = bandpass(lfp_norm, bandpass_band, fs, 'Steepness', 0.5);
trim_samples = round(fs);                       % avoid filter edge effects
lfp_bp = lfp_bp(trim_samples+1:end-trim_samples);
t_bp   = t(trim_samples+1:end-trim_samples);

%% 4. CALCULAR ENERGÍA
% Rectify and emphasize bursts (energy proxy)
lfp_pos    = max(lfp_bp, 0);
lfp_power  = lfp_pos .^ 4;

% Smooth energy
win_samples   = max(1, round(window_sec_energy * fs));
energy_metric = movmean(lfp_power, win_samples);
t_energy      = t_bp;

%% 5. UMBRALES AUTOMÁTICOS + HYSTERESIS
median_energy = median(energy_metric);
thr_enter     = enter_factor * median_energy;
thr_exit      = exit_ratio   * thr_enter;

% Hysteresis pass: once "in seizure", stay in until energy < thr_exit
is_enter = energy_metric > thr_enter;
is_exit  = energy_metric > thr_exit;

in_sz = false(size(energy_metric));
state = false;
for i = 1:numel(energy_metric)
    if ~state && is_enter(i)
        state = true;
    elseif state && ~is_exit(i)
        state = false;
    end
    in_sz(i) = state;
end

%% 6. FILL SHORT HOLES (≤ gap_sec)
gap_samples   = round(gap_sec * fs);
d             = diff([0; in_sz; 0]);
starts        = find(d==1);
ends          = find(d==-1)-1;

if ~isempty(starts)
    in_sz_filled = in_sz;
    for k = 1:numel(starts)-1
        hole = starts(k+1) - ends(k) - 1;
        if hole > 0 && hole <= gap_samples
            in_sz_filled(ends(k)+1 : starts(k+1)-1) = true;
        end
    end
else
    in_sz_filled = in_sz;
end

%% 7. MERGE SEGMENTS (≤ merge_sec)
merge_samples = round(merge_sec * fs);
d2      = diff([0; in_sz_filled; 0]);
st      = find(d2==1);
en      = find(d2==-1)-1;

keep_st = [];
keep_en = [];
if ~isempty(st)
    cur_st = st(1);
    cur_en = en(1);
    for k = 2:numel(st)
        if st(k) - cur_en - 1 <= merge_samples
            % extend current segment
            cur_en = en(k);
        else
            % save and start new
            keep_st(end+1) = cur_st; %#ok<AGROW>
            keep_en(end+1) = cur_en; %#ok<AGROW>
            cur_st = st(k);
            cur_en = en(k);
        end
    end
    % save last
    keep_st(end+1) = cur_st;
    keep_en(end+1) = cur_en;
end

%% 8. APLICAR DURACIÓN MÍNIMA
dur_samples = keep_en - keep_st + 1;
dur_seconds = dur_samples / fs;
mask        = dur_seconds >= min_seizure_sec;
keep_st     = keep_st(mask);
keep_en     = keep_en(mask);
dur_seconds = dur_seconds(mask);

%% 9. BUILD RESULTS TABLE
n_seizures = numel(keep_st);
if n_seizures > 0
    seiz_start_time = t_energy(keep_st);
    seiz_end_time   = t_energy(keep_en);
    seizures = table( (1:n_seizures)', ...
                      seiz_start_time(:), ...
                      seiz_end_time(:), ...
                      dur_seconds(:), ...
                      'VariableNames', {'Id','start_time_s','end_time_s','Duration_s'});
else
    seizures = table([],[],[],[], ...
                     'VariableNames', {'Id','start_time_s','end_time_s','Duration_s'});
end

%% 10. REPORT
fprintf('\n=== SEIZURE DETECTION ===\n');
fprintf('Archivo:         %s\n', file_path);
fprintf('Median energy:   %.3e\n', median_energy);
fprintf('Enter thr:       %.3e  (x%.1f)\n', thr_enter, enter_factor);
fprintf('Exit thr:        %.3e  (%.0f%% of enter)\n', thr_exit, 100*exit_ratio);
fprintf('Hole fill <=:    %.1f s\n', gap_sec);
fprintf('Merge gaps <=:   %.1f s\n', merge_sec);
fprintf('Min duration >=: %.1f s\n', min_seizure_sec);
display_seizures(seizures);

%% 11. PLOTS
fig = figure('Position', [50, 50, 1400, 900]);
fig.Name = sprintf('Seizure Detection: %s', file_path);

tiledlayout(3,1)

% (1) Cleaned LFP
nexttile
plot(t, lfp_clean, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5); 
title('Cleaned LFP (\muV)'); 
xlabel('Time (s)'); 
ylabel('\muV'); 
hold on
for k = 1:n_seizures
    xline(seizures.start_time_s(k),'r','LineWidth',1.5);
    xline(seizures.end_time_s(k)  ,'r','LineWidth',1.5);
end
grid on
hold off

% (2) Band-passed
nexttile
plot(t_bp, lfp_bp, 'Color', [0.2 0.4 0.7], 'LineWidth', 0.5); 
title('Band-passed (5-75 Hz)'); 
xlabel('Time (s)'); 
ylabel('Norm.'); 
hold on
for k = 1:n_seizures
    xline(seizures.start_time_s(k),'r','LineWidth',1.5);
    xline(seizures.end_time_s(k)  ,'r','LineWidth',1.5);
end
grid on
hold off

% (3) Energy + thresholds
nexttile
plot(t_energy, energy_metric, 'Color', [0.2 0.6 0.3], 'LineWidth', 0.5); 
hold on
yline(thr_enter,'r--','Enter','LineWidth',2);
yline(thr_exit,'k--','Exit','LineWidth',1.5);
title('Energy metric'); 
xlabel('Time (s)'); 
ylabel('a.u.');
for k = 1:n_seizures
    xline(seizures.start_time_s(k),'r','LineWidth',1.5);
    xline(seizures.end_time_s(k)  ,'r','LineWidth',1.5);
end
grid on
hold off

%% 12. GUARDAR FIGURA
try
    [~, fname, ~] = fileparts(file_path);
    fig_name = sprintf('%s_seizure_detection.fig', fname);
    saveas(fig, fig_name);
    fprintf('\nFigura guardada: %s\n', fig_name);
catch ME
    warning('No se pudo guardar la figura: %s', ME.message);
end

%% 13. GUARDAR RESULTADOS
if n_seizures > 0
    [~, fname, ~] = fileparts(file_path);
    results_file = sprintf('%s_seizures.mat', fname);
    save(results_file, 'seizures', 'file_path', 'fs', ...
         'enter_factor', 'exit_ratio', 'gap_sec', 'merge_sec', 'min_seizure_sec');
    fprintf('Resultados guardados: %s\n', results_file);
end

fprintf('\n*** Detección completada ***\n\n');

%% ========================================================================
%% FUNCIÓN LOCAL
%% ========================================================================
function display_seizures(tbl)
    n = height(tbl);
    
    if n == 0
        fprintf('\n--- No seizures detected ---\n\n');
        return;
    end
    
    % Definir anchos de columnas
    w_id = 5;
    w_start = 12;
    w_end = 12;
    w_dur = 12;
    
    % Header
    fprintf('\nDetected Seizures: %d\n', n);
    fprintf('%-*s  %-*s  %-*s  %-*s\n', ...
        w_id, 'Id', w_start, 'Start(s)', w_end, 'End(s)', w_dur, 'Duration(s)');
    
    % Línea separadora automática
    fprintf('%s  %s  %s  %s\n', ...
        repmat('-', 1, w_id), ...
        repmat('-', 1, w_start), ...
        repmat('-', 1, w_end), ...
        repmat('-', 1, w_dur));
    
    % Datos con 3 decimales fijos
    for i = 1:n
        fprintf('%-*d  %*.3f  %*.3f  %*.3f\n', ...
            w_id, tbl.Id(i), ...
            w_start, tbl.start_time_s(i), ...
            w_end, tbl.end_time_s(i), ...
            w_dur, tbl.Duration_s(i));
    end
    fprintf('\n');
end