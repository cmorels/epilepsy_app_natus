%% SPLIT LONG RECORDINGS INTO 3-HOUR SEGMENTS
% Automatically splits any recording >3 hours into manageable chunks
% Output files maintain all metadata and are ready for your pipeline

clear; clc;

%% CONFIGURATION
folder_path = '';  % '' = current folder
segment_duration_hours = 3;

%% FIND FILES
if isempty(folder_path)
    folder_path = pwd;
end

% Find all .txt files (exclude already processed ones)
all_files = dir(fullfile(folder_path, '*.txt'));
files = all_files(~contains({all_files.name}, {'_clean.txt', '_segment', '_part'}));

if isempty(files)
    error('No .txt files found in: %s', folder_path);
end

fprintf('LONG RECORDING SPLITTER\n');
fprintf('Segment duration: %.1f hours\n', segment_duration_hours);
fprintf('Files found: %d\n\n', length(files));

%% PROCESS EACH FILE
for file_idx = 1:length(files)
    
    file_path = fullfile(files(file_idx).folder, files(file_idx).name);
    [~, fname, ~] = fileparts(files(file_idx).name);
    
    fprintf('Processing: %s\n', files(file_idx).name);
    
    try
        % Load data
        data = load_LFP_intan_txt(file_path);
        fs = data.fs;
        signal_uV = data.signal;
        t = data.time_seconds;
        
        total_duration_hours = max(t) / 3600;
        n_samples = length(signal_uV);
        
        fprintf('  Duration: %.2f hours (%.1f min)\n', total_duration_hours, total_duration_hours*60);
        fprintf('  Samples: %d\n', n_samples);
        fprintf('  Sampling rate: %.2f Hz\n', fs);
        
        % Check if splitting is needed
        if total_duration_hours <= segment_duration_hours
            fprintf('  → No splitting needed (< %.1f hours)\n\n', segment_duration_hours);
            continue;
        end
        
        % Calculate segments
        samples_per_segment = round(segment_duration_hours * 3600 * fs);
        num_segments = ceil(n_samples / samples_per_segment);
        
        fprintf('  → Splitting into %d segments\n', num_segments);
        
        % Create segments subfolder
        segments_folder = fullfile(folder_path, 'segments');
        if ~exist(segments_folder, 'dir')
            mkdir(segments_folder);
        end
        
        % Calculate original signal range
        original_signal_range = [min(signal_uV), max(signal_uV)];
        
        % Split and save
        for seg = 1:num_segments
            
            % Calculate sample range
            start_sample = (seg-1) * samples_per_segment + 1;
            end_sample = min(seg * samples_per_segment, n_samples);
            
            % Extract segment
            segment_signal = signal_uV(start_sample:end_sample);
            segment_time = t(start_sample:end_sample);
            
            % Calculate segment duration
            seg_duration_hours = length(segment_signal) / fs / 3600;
            seg_duration_min = seg_duration_hours * 60;
            
            % Time offset for this segment
            time_offset_hours = (seg-1) * segment_duration_hours;
            
            % Calculate segment signal range
            segment_signal_range = [min(segment_signal), max(segment_signal)];
            
            % Create filename
            segment_filename = sprintf('%s_part%d_of_%d.txt', fname, seg, num_segments);
            segment_path = fullfile(segments_folder, segment_filename);
            
            % Write file with metadata
            fid = fopen(segment_path, 'w');
            
            % ===== ORIGINAL METADATA (preserved from parent file) =====
            fprintf(fid, '# mouse_id = %s\n', data.mouse_id);
            fprintf(fid, '# fs = %.10f\n', fs);
            fprintf(fid, '# time_unit = seconds\n');
            
            if isfield(data, 'port')
                fprintf(fid, '# port = %s\n', data.port);
            end
            
            fprintf(fid, '# region = %s\n', data.region);
            
            if isfield(data, 'native_name')
                fprintf(fid, '# native_name = %s\n', data.native_name);
            end
            
            fprintf(fid, '# session_time = %s\n', data.session_time);
            fprintf(fid, '# signal_range = [%.2f, %.2f] microvolts\n', segment_signal_range(1), segment_signal_range(2));
            
            % ===== SEGMENTATION METADATA =====
            fprintf(fid, '# Original_file = %s\n', files(file_idx).name);
            fprintf(fid, '# Segment = %d of %d\n', seg, num_segments);
            fprintf(fid, '# Time_range = %.2f - %.2f hours from start\n', ...
                time_offset_hours, time_offset_hours + seg_duration_hours);
            fprintf(fid, '# Segment_duration = %.2f hours (%.1f min)\n', ...
                seg_duration_hours, seg_duration_min);
            fprintf(fid, '# Original_total_duration = %.2f hours (%.1f min)\n', ...
                total_duration_hours, total_duration_hours * 60);
            fprintf(fid, '# samples = %d\n', length(segment_signal));
            fprintf(fid, '# Original_signal_range = [%.2f, %.2f] microvolts\n', ...
                original_signal_range(1), original_signal_range(2));
            
            % Data column description
            fprintf(fid, '# columns = amplitude_microvolts\n');
            
            % Write signal data
            fprintf(fid, '%.6f\n', segment_signal);
            
            fclose(fid);
            
            fprintf('    Part %d: %.1f min [%.1f - %.1f hours] → %s\n', ...
                seg, seg_duration_min, time_offset_hours, ...
                time_offset_hours + seg_duration_hours, segment_filename);
        end
        
        fprintf('  ✓ Saved %d segments to: segments/\n\n', num_segments);
        
    catch ME
        fprintf('  ✗ ERROR: %s\n\n', ME.message);
        continue;
    end
end

fprintf('SPLITTING COMPLETE\n');
fprintf('Segmented files saved in: %s\n', fullfile(folder_path, 'segments'));
fprintf('\nYou can now run your main pipeline on the files in segments/\n');
