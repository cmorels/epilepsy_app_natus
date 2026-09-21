L_channels = [1,2];

path_name = 'C:\Users\camila.morel\Desktop\Codigos_Cami';
files = dir(fullfile(path_name, '*.edf'));

L_file_names = {files.name};

for k = 1:length(L_file_names)

    file_name = L_file_names{k};
    full_file_name = fullfile(path_name, file_name);
    
    % read only the info file for date, number of records, samples, annotations and labels
    info = edfinfo(full_file_name);
    disp('----------------info file imported----------------')
    display(file_name)

    date_str = char(info.StartDate);
    time_str = char(info.StartTime);

    fs = info.NumSamples(1) / seconds(info.DataRecordDuration);
    Record_duration = seconds(info.DataRecordDuration(1));
    nrec = info.NumDataRecords;

    total_time = duration(seconds(info.NumSamples(1) * info.NumDataRecords / fs));
    fprintf("total time of the edf file is %g hours\n", hours(total_time))

    % read EDF using relative times
    [hdr, ~] = edfread(full_file_name, 'TimeOutputType', 'duration');
    disp('signal imported')

    record_times = hdr.("Record Time");   % duration, relative to file start

    % differences between consecutive record start times
    dt = seconds(diff(record_times));
    
    % tolerance to handle small floating-point inconsistencies
    tol = 1e-3;

    % a new segment starts when the gap is larger than one record duration
    break_idx = find(abs(dt - Record_duration) > tol);
    diff_dt = dt - Record_duration;
    disp('---diff dt in discontinuous idx---')
    display(diff_dt(break_idx))
    display(dt(break_idx))
    

    % start/end record indices of each continuous segment
    segment_starts = [1; break_idx + 1];
    segment_ends   = [break_idx; nrec];
    display(segment_starts)
    display(segment_ends)

    fprintf('Found %d continuous segments\n', length(segment_starts))

    % absolute start time of EDF file
    base_datetime = datetime([date_str ' ' time_str], ...
        'InputFormat', 'dd.MM.yy HH.mm.ss', ...
        'TimeZone', 'Europe/Paris');
    
    destination_folder = fullfile(path_name, 'txt');
    if ~exist(destination_folder, 'dir')
        mkdir(destination_folder);
    end

    for seg = 1:length(segment_starts)

        start_rec = segment_starts(seg);
        end_rec   = segment_ends(seg);

        n_segment_records = end_rec - start_rec + 1;

        % segment start time relative to EDF file start
        segment_start_time = record_times(start_rec);
        %display(segment_start_time)

        % absolute start time, for synchronizing with video (not useful...)
        segment_start_datetime = base_datetime + segment_start_time;
        %display(segment_start_datetime)
        segment_unix_time = round(posixtime(segment_start_datetime));

        fprintf('Segment %d: rec %d -> %d, start at %s\n', ...
            seg, start_rec, end_rec, char(segment_start_datetime))

        % create one output matrix per segment
        % assumes same samples/record for all channels
        samples_per_record_ref = info.NumSamples(L_channels(1));
        % matrix (number of channels x total number of samples in the segment)
        L_y = zeros(length(L_channels), n_segment_records * samples_per_record_ref);

        for i = 1:length(L_channels)
            channel = L_channels(i);
            % Select variable name for the given channel index from header
            channel_name = hdr.Properties.VariableNames{channel};
            % Number of samples in this channel's record 
            samples_per_record = info.NumSamples(channel);

            % if channels have different NumSamples, reinitialize row length safely
            if samples_per_record ~= samples_per_record_ref
                error('Selected channels do not all have the same NumSamples per record.')
            end

            for record = start_rec:end_rec
            
                % Convert global record index (in the full EDF file) 
                % to a local index within the current segment (starting from 1)
                local_rec = record - start_rec + 1;
            
                % Determine where this record starts in the continuous signal vector
                start_idx = 1 + (local_rec - 1) * samples_per_record;
                
                % Determine where this record ends in the continuous signal vector
                end_idx   = local_rec * samples_per_record;
                
                % Get the raw signal samples (in µV) for this record and channel
                samples = hdr.(channel_name){record};
                
                % Place the samples into the output matrix at the correct position
                % Convert from µV to mV and ensure row orientation
                L_y(i, start_idx:end_idx) = transpose(samples) / 1000; % uV -> mV
                
            end
        end

        
        % save one txt per channel for this segment
        for i = 1:length(L_channels)
            channel = L_channels(i);
            channel_name = hdr.Properties.VariableNames{channel};

            % sprintf('%d_%d_segment%d_%s.txt', segment_unix_time, start_rec ,seg, channel_name);
            new_file_name = fullfile(destination_folder, ...
                sprintf('%d_segment%d_%s.txt', segment_unix_time ,seg, channel_name));

            N = size(L_y, 2);
            t_segment = seconds(segment_start_time) + (0:N-1)/fs;
            data_to_save = [t_segment(:), L_y(i,:)'];
            writematrix(data_to_save, new_file_name, 'delimiter', '\t');
        end
    end
end