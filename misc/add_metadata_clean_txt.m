% Add mouse_id, session_time, and region headers to _clean.txt files
% Modifies files in place without creating duplicates

% Get all files ending in _clean.txt
files = dir('*_clean.txt');

fprintf('Found %d files ending in _clean.txt\n', length(files));

for i = 1:length(files)
    filename = files(i).name;
    fprintf('\nProcessing: %s\n', filename);
    
    try
        % Parse filename to extract metadata
        [mouse_id, session_time, region] = parse_filename(filename);
        
        % Read the entire file
        fid = fopen(filename, 'r');
        if fid == -1
            error('Could not open file: %s', filename);
        end
        
        % Read all lines
        file_content = cell(0, 1);  % Initialize as column vector
        line_count = 0;
        while ~feof(fid)
            line = fgetl(fid);
            if ischar(line)
                line_count = line_count + 1;
                file_content{line_count, 1} = line;
            end
        end
        fclose(fid);
        
        % Check if headers already exist
        has_mouse_id = false;
        has_session_time = false;
        has_region = false;
        
        for j = 1:min(20, length(file_content))  % Check first 20 lines
            if contains(file_content{j}, '# mouse_id')
                has_mouse_id = true;
            end
            if contains(file_content{j}, '# session_time')
                has_session_time = true;
            end
            if contains(file_content{j}, '# region')
                has_region = true;
            end
        end
        
        % If all headers already exist, skip this file
        if has_mouse_id && has_session_time && has_region
            fprintf('  Headers already exist. Skipping.\n');
            continue;
        end
        
        % Find where to insert the new headers (after existing headers)
        insert_position = 0;
        for j = 1:length(file_content)
            if ~isempty(file_content{j}) && file_content{j}(1) == '#'
                insert_position = j;
            else
                break;
            end
        end
        
        % Create new header lines as column cell array
        new_headers = {
            sprintf('# mouse_id = %s', mouse_id);
            sprintf('# session_time = %s', session_time);
            sprintf('# region = %s', region)
        };
        
        % Insert new headers after existing headers
        if insert_position > 0
            file_content = [file_content(1:insert_position); new_headers; file_content(insert_position+1:end)];
        else
            % If no headers found, add at the beginning
            file_content = [new_headers; file_content];
        end
        
        % Write back to the same file
        fid = fopen(filename, 'w');
        if fid == -1
            error('Could not open file for writing: %s', filename);
        end
        
        for j = 1:length(file_content)
            fprintf(fid, '%s\n', file_content{j});
        end
        fclose(fid);
        
        fprintf('  Added headers successfully.\n');
        fprintf('    mouse_id: %s\n', mouse_id);
        fprintf('    session_time: %s\n', session_time);
        fprintf('    region: %s\n', region);
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
    end
end

fprintf('\nProcessing complete!\n');

%% Helper function to parse filename
function [mouse_id, session_time, region] = parse_filename(filename)
    % Remove .txt extension
    base = strrep(filename, '.txt', '');
    
    % Remove _clean suffix
    base = strrep(base, '_clean', '');
    
    % Split by underscore
    parts = strsplit(base, '_');
    
    % Extract mouse_id (first part)
    mouse_id = parts{1};
    
    % Extract date and time (second and third parts)
    date_str = parts{2};  % yyyymmdd
    time_str = parts{3};  % hhmmss
    
    % Parse date
    year = str2double(date_str(1:4));
    month = str2double(date_str(5:6));
    day = str2double(date_str(7:8));
    
    % Parse time
    hour = str2double(time_str(1:2));
    minute = str2double(time_str(3:4));
    second = str2double(time_str(5:6));
    
    % Format session_time
    session_time = sprintf('%02d %s %04d %02d:%02d:%02d', ...
        day, get_month_name(month), year, hour, minute, second);
    
    % Extract region (everything after time until 'part' or end)
    region_parts = {};
    for i = 4:length(parts)
        if strcmp(parts{i}, 'part')
            break;
        end
        region_parts{end+1} = parts{i};
    end
    region = strjoin(region_parts, '_');
end

%% Helper function to get month name
function month_name = get_month_name(month_num)
    months = {'january', 'february', 'march', 'april', 'may', 'june', ...
              'july', 'august', 'september', 'october', 'november', 'december'};
    month_name = months{month_num};
end