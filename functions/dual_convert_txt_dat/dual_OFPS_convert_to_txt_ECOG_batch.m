% dual_OFPS_convert_to_txt_ECOG_v3.m
% Author: Esther Bliard + Camila Morel + ajustes
% "ONE FILE PER SIGNAL" con encabezado fs y recorte por ventana [t0 t1].
% VERSIÓN ACTUALIZADA: Nueva estructura de carpetas 
% MODIFICACIÓN: Guarda archivos en carpeta actual Y en batch_analyses

clear; clc;
read_Intan_RHD2000_file;

%% --- USUARIO: ventana a exportar [t0 t1] en segundos --------------------
% Ejemplos:
% export_window_sec = [120 300];   % de 120 s a 300 s
% export_window_sec = [];          % exporta TODO
export_window_sec = [];
% -------------------------------------------------------------------------

%% --- Configurar carpeta de destino batch_analyses ---
batch_folder = 'C:\Users\camila.morel\Documents\Analyses\ECoG\batch_analyses';

% Verificar que existe la carpeta batch_analyses
if ~exist(batch_folder, 'dir')
    warning('La carpeta batch_analyses no existe. Se creará automáticamente.');
    mkdir(batch_folder);
end

%% --- Parámetros de sesión / nombres robustos
% Nueva estructura de carpetas:
% C:\...\ECoG\{ID_A}-{ID_B}\{yyyymmdd}\{ID_A}_{ID_B}_{hhmmss}

current_path = pwd;
path_parts = split(current_path, filesep);

% Valores por defecto
mouse_id_A   = 'A';
mouse_id_B   = 'B';
session_date = 'unknownDate';
session_time = 'unknownTime';
date_time = datestr(now, 'dd-mmm-yyyy HH:MM:SS');

% Intentar extraer de la estructura de carpetas
try
    % Última carpeta: {ID_A}_{ID_B}_{hhmmss}
    session_folder = path_parts{end};
    session_parts = split(session_folder, '_');
    
    if numel(session_parts) >= 3
        mouse_id_A   = session_parts{1};  % Primer ID (puerto A)
        mouse_id_B   = session_parts{2};  % Segundo ID (puerto B)
        session_time = session_parts{3};  % hhmmss
    end
    
    % Penúltima carpeta: {yyyymmdd}
    if numel(path_parts) >= 2
        date_folder = path_parts{end-1};
        if numel(date_folder) == 8 && all(isstrprop(date_folder, 'digit'))
            session_date = date_folder;  % yyyymmdd
            
            % Construir datetime completo: yyyymmdd + hhmmss
            datetime_str = strcat(session_date, session_time);
            session_datetime = datetime(datetime_str, 'InputFormat', 'yyyyMMddHHmmss');
            date_time = datestr(session_datetime, 'dd-mmm-yyyy HH:MM:SS');
        end
    end
    
    fprintf('\n=== Información de Sesión ===\n');
    fprintf('Ratón en puerto A:  %s\n', mouse_id_A);
    fprintf('Ratón en puerto B:  %s\n', mouse_id_B);

    if ~strcmp(session_date, 'unknownDate')
        formatted_date = datestr(session_datetime, 'dd-mm-yyyy');
    else
        formatted_date = 'unknownDate';
    end
    fprintf('Fecha:           %s\n', formatted_date);
    
    % Formatear hora como HH:MM:SS
    if ~strcmp(session_time, 'unknownTime') && numel(session_time) == 6
        formatted_time = sprintf('%s:%s:%s', session_time(1:2), session_time(3:4), session_time(5:6));
    else
        formatted_time = 'unknownTime';
    end
    fprintf('Hora de inicio:  %s\n', formatted_time);
    fprintf('\n');
    
catch ME
    warning('No se pudo extraer información de la ruta. Usando valores por defecto.');
    warning('Error: %s', ME.message);
end

% No necesitamos lfp_filename_root, usamos directamente session_date y session_time

%% --- fs desde Intan
fs = frequency_parameters.amplifier_sample_rate;

%% --- Cargar multicanal desde amplifier.dat
num_channels = numel(amplifier_channels);
fileinfo = dir('amplifier.dat');
if isempty(fileinfo)
    error('No se encontró amplifier.dat en la carpeta actual.');
end
num_samples = fileinfo.bytes / (num_channels * 2); % int16 = 2 bytes
fid = fopen('amplifier.dat', 'r');
v = fread(fid, [num_channels, num_samples], 'int16');
fclose(fid);
v = v * 0.195; % -> microvoltios

%% --- Cargar vector de tiempo (en segundos)
tinfo = dir('time.dat');
if isempty(tinfo)
    error('No se encontró time.dat en la carpeta actual.');
end
fid = fopen('time.dat', 'r');
time = fread(fid, tinfo.bytes/4, 'int32');
fclose(fid);
time = double(time) / fs; % segundos

% Calcular duración del registro
duration_seconds = time(end) - time(1);
duration_hours = floor(duration_seconds / 3600);
duration_minutes = floor(mod(duration_seconds, 3600) / 60);
duration_secs = floor(mod(duration_seconds, 60));

fprintf('Duración:        %02d:%02d:%02d\n\n', ...
    duration_hours, duration_minutes, duration_secs);

%% --- Mapeo región por sufijo de canal (ajústalo a tu headstage)
region_mapping = containers.Map( ...
    {'016','023'}, ...
    {'HPCleft','HPCright'} );

% Armar listas por puerto (A/B) a partir de la info en amplifier_channels
channels_A = struct('idx', {}, 'region', {}, 'native_name', {});
channels_B = struct('idx', {}, 'region', {}, 'native_name', {});

for i = 1:num_channels
    native = amplifier_channels(i).native_channel_name;  % ej. 'A-016'
    if startsWith(native, 'A-') || startsWith(native, 'B-')
        suffix = extractAfter(native, '-');
        if isKey(region_mapping, suffix)
            region = region_mapping(suffix);
            ch_info = struct('idx', i, 'region', region, 'native_name', native);
            if startsWith(native, 'A-')
                channels_A(end+1) = ch_info; %#ok<SAGROW>
            else
                channels_B(end+1) = ch_info; %#ok<SAGROW>
            end
        end
    end
end

%% --- Exportar por puerto (un archivo por señal en la ventana [t0 t1])
procesar_puerto_window(channels_A, 'A', mouse_id_A, v, time, fs, session_date, session_time, date_time, export_window_sec, batch_folder);
procesar_puerto_window(channels_B, 'B', mouse_id_B, v, time, fs, session_date, session_time, date_time, export_window_sec, batch_folder);

%% FUNCION LOCAL
function procesar_puerto_window(channels, port_label, mouse_id, v, time, fs, session_date, session_time, date_time, export_window_sec, batch_folder)
% Genera un .txt por señal (región) con header (incluye fs) y recorte exacto [t0 t1].
% MODIFICACIÓN: Exporta señal SIN offset (+10000)
% También grafica la señal recortada para verificación rápida.
% NUEVA: Guarda en carpeta actual Y en batch_analyses

    if isempty(channels), return; end

    % --- calcular índices de recorte globales ---
    ns_total = min(size(v,2), numel(time));
    t = time(1:ns_total);

    if isempty(export_window_sec)
        idx_start = 1;
        idx_end   = ns_total;
        tag_str   = '';
    else
        t0 = export_window_sec(1);
        t1 = export_window_sec(2);
        if t1 <= t0
            error('export_window_sec debe ser [t0 t1] con t1 > t0.');
        end
        idx_start = max(1, floor(t0*fs) + 1);
        idx_end   = min(ns_total, floor(t1*fs) + 1);
        if idx_start >= idx_end
            error('La ventana solicitada no contiene muestras válidas.');
        end
        tag_str = sprintf('_%ds-%ds', round(t(idx_start)), round(t(idx_end)));
    end

    fig = figure('Name', sprintf('Puerto %s (ventana)', port_label), 'NumberTitle', 'off');
    hold on; legend_entries = {};

    for i = 1:numel(channels)
        ch = channels(i);
        fprintf("Puerto %s: %s (%s)\n", port_label, ch.region, ch.native_name);

        lfp_full = v(ch.idx, 1:ns_total);   % µV
        lfp_seg  = lfp_full(idx_start:idx_end);

        % *** CAMBIO PRINCIPAL: NO aplicar offset ***
        % Versión original: lfp_seg_positive = lfp_seg + 10000;
        lfp_seg_export = lfp_seg;  % Señal original sin modificar

        % tiempos del segmento
        t_seg = t(idx_start:idx_end);

        % nombre de archivo de salida (un archivo por señal)
        % Formato: {ID_raton}_{yyyymmdd}_{hhmmss}_{region}.txt
        % Ejemplo: 39921_20250429_143247_HPCleft.txt
        
        if isempty(tag_str)
            out_name = sprintf('%s_%s_%s_%s.txt', ...
                mouse_id, session_date, session_time, ch.region);
        else
            out_name = sprintf('%s_%s_%s_%s%s.txt', ...
                mouse_id, session_date, session_time, ch.region, tag_str);
        end

        % ---- Preparar header común ----
        header_lines = {
            sprintf('# mouse_id = %s\n', mouse_id)
            sprintf('# fs = %.5f\n', fs)
            sprintf('# time_unit = seconds\n')
            sprintf('# port = %s\n', port_label)
            sprintf('# region = %s\n', ch.region)
            sprintf('# native_name = %s\n', ch.native_name)
            sprintf('# session_time = %s\n', date_time)
            sprintf('# exported_window = [%0.6f %0.6f] seconds\n', t_seg(1), t_seg(end))
             sprintf('# columns = amplitude_microvolts\n')
            sprintf('# signal_range = [%.2f, %.2f] microvolts\n', min(lfp_seg_export), max(lfp_seg_export))
        };

        % ---- Guardar en carpeta ACTUAL ----
        fid = fopen(out_name, 'w');
        if fid == -1
            warning('No se pudo abrir para escritura (carpeta actual): %s', out_name);
        else
            % Escribir header
            for h = 1:numel(header_lines)
                fprintf(fid, '%s', header_lines{h});
            end
            % Escribir datos
            fprintf(fid, '%.6f\n', lfp_seg_export);
            fclose(fid);
            fprintf(" Guardado (actual):      %s (%d muestras)\n", ...
                    out_name, numel(lfp_seg_export));
        end

        % ---- Guardar en carpeta BATCH_ANALYSES ----
        batch_filepath = fullfile(batch_folder, out_name);
        fid_batch = fopen(batch_filepath, 'w');
        if fid_batch == -1
            warning('No se pudo abrir para escritura (batch_analyses): %s', batch_filepath);
        else
            % Escribir header
            for h = 1:numel(header_lines)
                fprintf(fid_batch, '%s', header_lines{h});
            end
            % Escribir datos
            fprintf(fid_batch, '%.6f\n', lfp_seg_export);
            fclose(fid_batch);
            fprintf(" Guardado (batch):       %s\n", batch_filepath);
        end

        % plot del segmento (SIN OFFSET)
        color = 'b'; if strcmp(ch.region,'HPCright'), color = 'r'; end
        plot(t_seg, lfp_seg_export, color);
        legend_entries{end+1} = sprintf('%s', ch.region); %#ok<AGROW>
    end

    title(sprintf('Recording of mouse %s in port %s (%s)', mouse_id, port_label, date_time), 'Interpreter','none');
    xlabel('Time (s)'); ylabel('(\muV)');
    legend(legend_entries, 'Location','best'); 
    yline(0,'--k','LineWidth',1.5);
    grid on;

    % guardar figura (con tag si aplica)
    try
        png_name = sprintf('Graph_%s.png', mouse_id);
        fig_name = sprintf('Graph_%s.fig', mouse_id);
        
        saveas(fig, png_name);
        savefig(fig, fig_name);
    catch ME
        warning("No se pudo guardar figura:\n%s", ME.message);
    end
     close(fig);
end