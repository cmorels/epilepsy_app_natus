function data = load_LFP_intan_txt(file_path, fs)
% 1-channel LFP from Intan-exported .txt
% 
% - Si no se entrega file_path, carga el primer .txt de la carpeta actual.
% - Si no se entrega fs, intenta leerlo del header (# fs = ...) del archivo.
%
% OUTPUT (struct data):
%   .signal       : vector columna con la señal
%   .time_seconds : vector tiempo en segundos (1/fs, 2/fs, ...)
%   .fs           : frecuencia de muestreo
%   .file         : nombre del archivo
%   .folder       : carpeta

    %% elegir archivo si no se entrega file_path
    if nargin < 1 || isempty(file_path)
        files = dir('*.txt');
        if isempty(files)
            error('load_LFP_intan_txt:NoFilesFound', ...
                  'No .txt files found in the current directory.');
        end
        file_path = fullfile(files(1).folder, files(1).name);
    end

    %% normalizar file_path y revisar extensión
    if ~(ischar(file_path) || isstring(file_path))
        error('load_LFP_intan_txt:BadInput', ...
              'file_path must be a char/string with full path to the .txt file');
    end
    file_path = char(file_path);  % ensure char

    [folder, file, ext] = fileparts(file_path);
    if ~strcmpi(ext, '.txt')
        error('load_LFP_intan_txt:InputError', ...
              'Expected a .txt file, got "%s".', ext);
    end

    %% abrir archivo
    fid = fopen(file_path, 'rt');
    if fid == -1
        error('load_LFP_intan_txt:FileOpenError', ...
              'Could not open file: %s', file_path);
    end

    %% leer header (líneas que empiezan con #) y buscar fs
    header_fs = [];
    header_mouse_id = ''; 
    header_region = '';
    header_session_time = '';
    pos_after_header = ftell(fid);

    line = fgetl(fid);
    while ischar(line) && startsWith(strtrim(line), '#')
        tokens = regexp(line, 'fs\s*=\s*([0-9.eE+\-]+)', 'tokens', 'once');
            if ~isempty(tokens)
                header_fs = str2double(tokens{1});
            end
    
        tokens = regexp(line, 'mouse_id\s*=\s*(.+)', 'tokens', 'once');
            if ~isempty(tokens)
                header_mouse_id = strtrim(tokens{1});
            end

        tokens = regexp(line, 'region\s*=\s*(.+)', 'tokens', 'once');
            if ~isempty(tokens)
                header_region = strtrim(tokens{1});
            end
        tokens = regexp(line, 'session_time\s*=\s*(.+)', 'tokens', 'once');
            if ~isempty(tokens)
                header_session_time = strtrim(tokens{1});
            end

        pos_after_header = ftell(fid);
        line = fgetl(fid);
    end

    % Volver justo al final del header para leer solo números
    fseek(fid, pos_after_header, 'bof');

    %% decidir fs (argumento > header)
    if nargin < 2 || isempty(fs)
        if ~isempty(header_fs)
            fs = header_fs;
        else
            fclose(fid);
            error('load_LFP_intan_txt:NoFsFound', ...
                  'No fs provided and no "fs = ..." found in header of %s', file_path);
        end
    end

    %% leer datos numéricos (señal)
    C = textscan(fid, '%f');  % solo números
    fclose(fid);

    signal = C{1};
    n      = numel(signal);

    %% construir vector de tiempo
    time_seconds = (1:n) / fs;   % 1/fs, 2/fs, ...

    %% empaquetar en struct
    data = struct();
    data.session_time = header_session_time;
    data.mouse_id     = header_mouse_id;   
    data.region       = header_region;     
    data.signal       = signal;
    data.time_seconds = time_seconds;
    data.fs           = fs;
    data.file         = [file ext];
    data.folder       = folder;
end

