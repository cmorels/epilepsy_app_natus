function out = spec_signal_cache(action, varargin)
% SPEC_SIGNAL_CACHE  Read only the needed window of a pipeline txt
% (01_txt / 02_clean) without holding the whole recording in memory.
%
%   info = spec_signal_cache('open', txt_path, cache_dir, chunk)
%   x    = spec_signal_cache('read', info, i0, i1)   % samples i0..i1 (1-based), NaN outside the file
%   spec_signal_cache('reset')                       % forget opened files (and delete temp caches)
%
% The first 'open' of a txt streams its data section in chunks (textscan
% '%f', exactly the values load_lfp_txt.m would return) into
% <cache_dir>\<txt name>.bin (float64), and its header into a .mat next to
% it. Later opens (same or later runs) reuse it while the txt keeps the
% same size and modification date; otherwise it is rebuilt. A
% 24-48 h txt holds ~45-85 M samples: load_lfp_txt.m would need ~2 GB of
% RAM for it, the cache needs one chunk.
%
% Header, fs, session_start and meta come from load_lfp_txt.m itself, run
% on a temporary copy of the header plus one sample, so they are parsed
% exactly as the pipeline parses them.
%
% info: path, bin, n, fs, session_start, meta, header_pairs, file, folder,
%       build_s (seconds spent building, 0 if reused), reused (logical).

    persistent opened temp_bins
    if isempty(opened)
        opened = containers.Map('KeyType', 'char', 'ValueType', 'any');
        temp_bins = {};
    end

    switch action
        case 'open'
            txt_path = char(varargin{1});
            cache_dir = varargin{2};
            chunk = 5e6;
            if numel(varargin) >= 3 && ~isempty(varargin{3})
                chunk = varargin{3};
            end
            temp = isempty(cache_dir);
            if temp
                cache_dir = fullfile(tempdir, 'spec_cache');
            end
            if isKey(opened, txt_path)
                out = opened(txt_path);
                out.build_s = 0;          % already opened in this run
                out.reused = true;
                return;
            end
            out = open_one(txt_path, cache_dir, chunk);
            opened(txt_path) = out;
            if temp && ~out.reused
                temp_bins{end+1} = out.bin;
            end

        case 'read'
            info = varargin{1};
            i0 = varargin{2};
            i1 = varargin{3};
            out = read_range(info, i0, i1);

        case 'reset'
            for k = 1:numel(temp_bins)
                if exist(temp_bins{k}, 'file') == 2
                    delete(temp_bins{k});
                end
                m = strrep(temp_bins{k}, '.bin', '_meta.mat');
                if exist(m, 'file') == 2
                    delete(m);
                end
            end
            opened = containers.Map('KeyType', 'char', 'ValueType', 'any');
            temp_bins = {};
            out = [];

        otherwise
            error('spec_signal_cache:BadAction', 'Unknown action ''%s''.', action);
    end
end

%% ======================================================================
function info = open_one(txt_path, cache_dir, chunk)
    if exist(txt_path, 'file') ~= 2
        error('spec_signal_cache:FileNotFound', 'File not found: %s', txt_path);
    end
    d = dir(txt_path);
    [~, name] = fileparts(txt_path);
    bin = fullfile(cache_dir, [name '.bin']);
    meta_file = fullfile(cache_dir, [name '_meta.mat']);

    if exist(meta_file, 'file') == 2 && exist(bin, 'file') == 2
        m = load(meta_file);
        b = dir(bin);
        if isfield(m, 'info') && strcmp(m.info.path, txt_path) && m.txt_bytes == d.bytes ...
                && abs(m.txt_datenum - d.datenum) < 1e-9 && b.bytes == 8 * m.info.n
            info = m.info;
            info.bin = bin;
            info.build_s = 0;
            info.reused = true;
            return;
        end
    end

    if ~isfolder(cache_dir)
        mkdir(cache_dir);
    end
    t0 = tic;
    fid = fopen(txt_path, 'rt');
    if fid < 0
        error('spec_signal_cache:FileOpenError', 'Could not open %s', txt_path);
    end
    closer = onCleanup(@() fclose(fid));
    parse_header(fid);   % leaves fid at the first data line
    h = spec_peek_header(txt_path);
    header_info = rmfield(h, 'n_lines');

    tmp_bin = [bin '.part'];
    fb = fopen(tmp_bin, 'w');
    if fb < 0
        error('spec_signal_cache:CacheWrite', 'Could not write the cache %s', tmp_bin);
    end
    n = 0;
    while true
        C = textscan(fid, '%f', chunk);
        v = C{1};
        if isempty(v)
            break;
        end
        fwrite(fb, v, 'double');
        n = n + numel(v);
        if numel(v) < chunk
            break;
        end
    end
    fclose(fb);
    clear closer
    movefile(tmp_bin, bin, 'f');

    info = header_info;
    info.path = txt_path;
    info.bin = bin;
    info.n = n;
    info.build_s = toc(t0);
    info.reused = false;
    if isfield(info.meta, 'n_samples')
        n_hdr = str2double(info.meta.n_samples);
        if ~isnan(n_hdr) && n_hdr ~= n
            warning('spec_signal_cache:SampleCount', '%s: header n_samples = %d but %d samples read.', ...
                txt_path, n_hdr, n);
        end
    end
    txt_bytes = d.bytes;
    txt_datenum = d.datenum;
    save(meta_file, 'info', 'txt_bytes', 'txt_datenum');
end

function x = read_range(info, i0, i1)
    i0 = round(i0);
    i1 = round(i1);
    x = NaN(i1 - i0 + 1, 1);
    a = max(i0, 1);
    b = min(i1, info.n);
    if b < a
        return;
    end
    fb = fopen(info.bin, 'r');
    if fb < 0
        error('spec_signal_cache:CacheRead', 'Could not read the cache %s', info.bin);
    end
    fseek(fb, 8 * (a - 1), 'bof');
    v = fread(fb, b - a + 1, 'double');
    fclose(fb);
    x((a - i0 + 1):(b - i0 + 1)) = v;
end
