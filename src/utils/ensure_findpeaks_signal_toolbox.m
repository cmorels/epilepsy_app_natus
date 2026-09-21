function ensure_findpeaks_signal_toolbox()
% ENSURE_FINDPEAKS_SIGNAL_TOOLBOX  Forces MATLAB to resolve `findpeaks` to
% the Signal Processing Toolbox implementation instead of a Chronux copy
% that may be earlier on the path (Chronux's findpeaks has an incompatible
% signature). Verbatim port of the path-fixing block from
% IID_detection_FINAL.m (lines 10-29) -- unchanged logic, only relocated
% into a function so it is called once instead of duplicated per file.
% Cheap to call repeatedly: it only touches the path when a Chronux
% findpeaks is actually found ahead of the toolbox one.

    sigtool = fullfile(matlabroot, 'toolbox', 'signal', 'signal');
    if exist(sigtool, 'dir')
        addpath(sigtool, '-begin');
    end

    fps = which('findpeaks', '-all');
    for i = 1:numel(fps)
        if contains(fps{i}, 'chronux', 'IgnoreCase', true)
            chronux_dir = fileparts(fps{i});
            prev = '';
            while ~isempty(chronux_dir) && ~strcmpi(get_last_dir(chronux_dir), 'chronux_2_12') && ~strcmp(prev, chronux_dir)
                prev = chronux_dir;
                chronux_dir = fileparts(chronux_dir);
            end
            if contains(chronux_dir, 'chronux_2_12', 'IgnoreCase', true)
                rmpath(genpath(chronux_dir));
            end
            break;
        end
    end
    rehash toolboxcache;
end

function d = get_last_dir(p)
    [~, d] = fileparts(p);
end
