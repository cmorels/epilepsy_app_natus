% dual_OFPC_convert_to_txt_ECOG_v4_batch.m
% Author: Esther Bliard + Camila Morel + ajustes
% "ONE FILE PER CHANNEL" con encabezado completo y recorte por ventana [t0 t1].
% VERSIÓN ACTUALIZADA: Nueva estructura de carpetas, sin offset
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
    
catch ME
    warning('No se pudo extraer información de la ruta. Usando valores por defecto.');
    warning('Error: %s', ME.message);
end

%% --- fs desde Intan
fs = frequency_parameters.amplifier_sample_rate;

%% --- Cargar vector de tiempo (en segundos)
tinfo = dir('time.dat');
if isempty(tinfo)
    error('No se encontró time.dat en la carpeta actual.');
end
fid = fopen('time.dat', 'r');
time = fread(fid, tinfo.bytes/4, 'int32');
fclose(fid);
time = double(time) / fs;   % segundos
nTime = numel(time);

% Calcular duración del registro
duration_seconds = time(end) - time(1);
duration_hours = floor(duration_seconds / 3600);
duration_minutes = floor(mod(duration_seconds, 3600) / 60);
duration_secs = floor(mod(duration_seconds, 60));

fprintf('Duración:        %02d:%02d:%02d\n\n', ...
    duration_hours, duration_minutes, duration_secs);

%% --- Definir mapeo de sufijos a regiones (ajústalo a tu headstage)
region_mapping = containers.Map( ...
    {'016','023'}, ...
    {'HPCleft','HPCright'});  % Mapea el sufijo del canal a la región

%% --- Precalcular índices de ventana (globales) sobre el eje de tiempo
if isempty(export_window_sec)
    idx_start = 1;
    idx_end   = nTime;
    tag_str   = '';
else
    t0 = export_window_sec(1);
    t1 = export_window_sec(2);
    if t1 <= t0
        error('export_window_sec debe ser [t0 t1] con t1 > t0.');
    end
    idx_start = max(1, floor(t0*fs) + 1);
    idx_end   = min(nTime, floor(t1*fs) + 1);
    if idx_start >= idx_end
        error('La ventana solicitada no contiene muestras válidas.');
    end
    tag_str = sprintf('_%ds-%ds', round(time(idx_start)), round(time(idx_end)));
end

%% --- Buscar y procesar archivos amp-*.dat
channel_files = dir('amp-*.dat');
if isempty(channel_files)
    error('No se encontraron archivos amp-*.dat en la carpeta actual.');
end

% Figuras separadas por puerto
fig_A = figure('Name', 'Puerto A (segment)', 'NumberTitle', 'off');
hold on; legend_entries_A = {};

fig_B = figure('Name', 'Puerto B (segment)', 'NumberTitle', 'off');
hold on; legend_entries_B = {};

for i = 1:length(channel_files)
    fname = channel_files(i).name;
    parts = regexp(fname, 'amp-(A|B)-(\d+)\.dat', 'tokens');
    if isempty(parts), continue; end

    port      = parts{1}{1};   % 'A' o 'B'
    ch_suffix = parts{1}{2};   % '016', '023', etc.
    native_name = sprintf('%s-%s', port, ch_suffix);

    if ~isKey(region_mapping, ch_suffix)
        % si no está mapeado, lo saltamos
        continue;
    end
    region = region_mapping(ch_suffix);
    
    % Determinar mouse_id según puerto
    if strcmp(port, 'A')
        mouse_id = mouse_id_A;
    else
        mouse_id = mouse_id_B;
    end

    % Cargar señal LFP de ese archivo
    finfo       = dir(fname);
    num_samples = finfo.bytes / 2;      % int16
    fid         = fopen(fname, 'r');
    if fid == -1
        warning('No se pudo abrir %s. Se omite.', fname);
        continue;
    end
    lfp = fread(fid, [1, num_samples], 'int16');
    fclose(fid);

    % convertir a microvoltios
    lfp = double(lfp) * 0.195;

    % Alinear longitudes por seguridad
    ns = min(numel(lfp), nTime);
    lfp = lfp(1:ns);
    t   = time(1:ns);

    % Aplicar ventana [idx_start idx_end]
    idx_end_local = min(idx_end, ns);
    if idx_start >= idx_end_local
        warning('Canal %s (%s): se quedó sin muestras en la ventana solicitada. Se omite.', native_name, region);
        continue;
    end

    lfp_seg = lfp(idx_start:idx_end_local);
    t_seg   = t(idx_start:idx_end_local);

    % *** SIN OFFSET - señal original ***
    lfp_seg_export = lfp_seg;

    % Nombre de archivo de salida .txt
    % Formato: {ID_raton}_{yyyymmdd}_{hhmmss}_{region}.txt
    if isempty(tag_str)
        new_file = sprintf('%s_%s_%s_%s.txt', ...
            mouse_id, session_date, session_time, region);
    else
        new_file = sprintf('%s_%s_%s_%s%s.txt', ...
            mouse_id, session_date, session_time, region, tag_str);
    end

    % ---- Preparar header común ----
    header_lines = {
        sprintf('# mouse_id = %s\n', mouse_id)
        sprintf('# fs = %.5f\n', fs)
        '# time_unit = seconds\n'
        sprintf('# port = %s\n', port)
        sprintf('# region = %s\n', region)
        sprintf('# native_name = %s\n', native_name)
        sprintf('# session_time = %s\n', date_time)
        sprintf('# exported_window = [%0.6f %0.6f] seconds\n', t_seg(1), t_seg(end))
        '# columns = amplitude_microvolts\n'
        sprintf('# signal_range = [%.2f, %.2f] microvolts\n', min(lfp_seg_export), max(lfp_seg_export))
    };

    % ---- Guardar en carpeta ACTUAL ----
    fid = fopen(new_file, 'w');
    if fid == -1
        warning('No se pudo abrir el archivo para escritura (actual): %s', new_file);
    else
        % Escribir header
        for h = 1:numel(header_lines)
            fprintf(fid, '%s', header_lines{h});
        end
        % Escribir datos
        fprintf(fid, '%.6f\n', lfp_seg_export);
        fclose(fid);

        fprintf("✅ Guardado (actual):      %s (%d muestras, %.1f–%.1f s)\n", ...
                new_file, numel(lfp_seg_export), t_seg(1), t_seg(end));
    end

    % ---- Guardar en carpeta BATCH_ANALYSES ----
    batch_filepath = fullfile(batch_folder, new_file);
    fid_batch = fopen(batch_filepath, 'w');
    if fid_batch == -1
        warning('No se pudo abrir el archivo para escritura (batch): %s', batch_filepath);
    else
        % Escribir header
        for h = 1:numel(header_lines)
            fprintf(fid_batch, '%s', header_lines{h});
        end
        % Escribir datos
        fprintf(fid_batch, '%.6f\n', lfp_seg_export);
        fclose(fid_batch);
        
        fprintf("📦 Guardado (batch):       %s\n", batch_filepath);
    end

    % ---- Graficar segmento para verificación rápida ----
    color = 'b';
    if strcmp(region, 'HPCright'), color = 'r'; end
    
    if strcmp(port, 'A')
        figure(fig_A);
        plot(t_seg, lfp_seg_export, color);
        legend_entries_A{end+1} = sprintf('%s', region); %#ok<AGROW>
    else
        figure(fig_B);
        plot(t_seg, lfp_seg_export, color);
        legend_entries_B{end+1} = sprintf('%s', region); %#ok<AGROW>
    end
end

% Finalizar figura Puerto A
figure(fig_A);
title(sprintf('Recording of mouse %s in port A (%s)', mouse_id_A, date_time), 'Interpreter','none');
xlabel('Time (s)');
ylabel('(\muV)');
if ~isempty(legend_entries_A)
    legend(legend_entries_A, 'Location', 'best');
end
yline(0, '--k', 'LineWidth', 1.5);
grid on;

try
    png_name_A = sprintf('Graph_%s.png', mouse_id_A);
    fig_name_A = sprintf('Graph_%s.fig', mouse_id_A);
    saveas(fig_A, png_name_A);
    savefig(fig_A, fig_name_A);
    fprintf("📊 Figura guardada: %s\n", fig_name_A);
catch ME
    warning("No se pudo guardar figura puerto A:\n%s", ME.message);
end

% Finalizar figura Puerto B
figure(fig_B);
title(sprintf('Recording of mouse %s in port B (%s)', mouse_id_B, date_time), 'Interpreter','none');
xlabel('Time (s)');
ylabel('(\muV)');
if ~isempty(legend_entries_B)
    legend(legend_entries_B, 'Location', 'best');
end
yline(0, '--k', 'LineWidth', 1.5);
grid on;

try
    png_name_B = sprintf('Graph_%s.png', mouse_id_B);
    fig_name_B = sprintf('Graph_%s.fig', mouse_id_B);
    saveas(fig_B, png_name_B);
    savefig(fig_B, fig_name_B);
    fprintf("📊 Figura guardada: %s\n", fig_name_B);
catch ME
    warning("No se pudo guardar figura puerto B:\n%s", ME.message);
end

fprintf('\n*** CONVERSIÓN COMPLETA ***\n\n');