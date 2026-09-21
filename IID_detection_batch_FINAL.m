% IID_detection_batch.m
% Batch processing: Improved interictal spike detection with polyspike complex identification
% NOW COMPATIBLE WITH PIPELINE BATCH OUTPUT
% - Processes ALL *_clean.txt files in the folder
% - No offset removal needed
% - Consistent data loading with pipeline batch

clear; clc;

%% Force 'findpeaks' from Signal Toolbox
sigtool = fullfile(matlabroot,'toolbox','signal','signal');
if exist(sigtool,'dir'); addpath(sigtool,'-begin'); end
fps = which('findpeaks','-all');
for i = 1:numel(fps)
    if contains(fps{i}, 'chronux', 'IgnoreCase', true)
        chronux_dir = fileparts(fps{i});
        prev = '';
        while ~isempty(chronux_dir) && ~strcmpi(getLastDir(chronux_dir),'chronux_2_12') && ~strcmp(prev,chronux_dir)
            prev = chronux_dir;
            chronux_dir = fileparts(chronux_dir);
        end
        if contains(chronux_dir,'chronux_2_12','IgnoreCase',true)
            rmpath(genpath(chronux_dir));
        end
        break
    end
end
rehash toolboxcache
disp('findpeaks in path (should list toolbox/signal first):'); which findpeaks -all

%% ----------------- 1) Find all clean files in folder -----------------
folder_path = '';  % '' = current folder

if isempty(folder_path)
    folder_path = pwd;
end

% Find all *_clean.txt files
files = dir(fullfile(folder_path, '*_clean.txt'));

if isempty(files)
    error('No *_clean.txt files found in: %s\nPlease run the pipeline batch first.', folder_path);
end

num_files = length(files);

fprintf('\n=== IID DETECTION BATCH PROCESSING ===\n');
fprintf('Folder: %s\n', folder_path);
fprintf('Files found: %d\n\n', num_files);

for i = 1:num_files
    fprintf('%2d. %s\n', i, files(i).name);
end

fprintf('\n');

% Create output subfolder
output_folder = fullfile(folder_path, 'IID_results');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
    fprintf('Created output folder: %s\n\n', output_folder);
else
    fprintf('Using existing output folder: %s\n\n', output_folder);
end

% Summary storage (for OUTPUT 1: overview CSV)
batch_summary = cell(num_files, 12);

%% ----------------- 2) Process each file -----------------
for file_idx = 1:num_files
    
    fprintf('\n');
    fprintf('==========================================================\n');
    fprintf('PROCESSING FILE %d/%d: %s\n', file_idx, num_files, files(file_idx).name);
    fprintf('==========================================================\n');
    
    file_path = fullfile(files(file_idx).folder, files(file_idx).name);
    [~, fname, ~] = fileparts(files(file_idx).name);
    fname_base = strrep(fname, '_clean', '');
    
    try
        
        %% LOAD DATA
        fprintf('\n[LOADING DATA]\n');
        data        = load_LFP_intan_txt(file_path);
        session     = data.session_time;
        mouse_id    = data.mouse_id;
        region      = data.region;
        fs          = data.fs;
        signal_uV   = data.signal;
        t           = data.time_seconds;
        
        total_recording_min = max(t) / 60;
        n_samples = numel(signal_uV);
        
        fprintf('  File:     %s\n', data.file);
        fprintf('  Session:  %s\n', session);
        fprintf('  Mouse ID: %s\n', mouse_id);
        fprintf('  Region:   %s\n', region);
        fprintf('  Duration: %.1f min\n', total_recording_min);
        fprintf('  Samples:  %d\n', n_samples);
        fprintf('  Fs:       %.2f Hz\n', fs);
        
        sig_uV = signal_uV;
        
        %% BASELINE calculation
        Voltage_mV      = (sig_uV / 1000);
        ABS_Voltage_mV  = abs(Voltage_mV);
        
        Start_Value_uV  = 100;
        Step_uV         = 5;
        TargetFrac      = 0.97;
        Max_uV          = 130;
        
        Percent = 0;
        sv = Start_Value_uV;
        
        while Percent < TargetFrac && sv <= Max_uV
            thr_mV = sv / 1000;
            Percent = mean(ABS_Voltage_mV < thr_mV);
            if Percent < TargetFrac
                sv = sv + Step_uV;
            end
        end
        
        if sv > Max_uV
            Baseline_uV = Max_uV - Step_uV;
        else
            Baseline_uV = sv - Step_uV;
        end
        
        Lower_Threshold_uV = 2.5 * Baseline_uV;
        Upper_Threshold_uV = 2000;
        Lower_Threshold_mV = Lower_Threshold_uV / 1000;
        Upper_Threshold_mV = Upper_Threshold_uV / 1000;
        
        fprintf('\n[BASELINE DETECTION]\n');
        fprintf('  Baseline:         %.0f µV\n', Baseline_uV);
        fprintf('  Lower threshold:  %.0f µV\n', Lower_Threshold_uV);
        fprintf('  Upper threshold:  %.0f µV\n', Upper_Threshold_uV);
        
        %% Band-pass filter 15–70 Hz
        x_mV  = (sig_uV - median(sig_uV,'omitnan'))/1000;
        x_bp  = bandpass(x_mV, [15 70], fs, 'Steepness', 0.5);
        abs_mV = abs(x_bp);
        
        %% Initial spike detection
        MinPeakDistance = max(1, round(0.030 * fs));
        MaxPeakWidth    = max(1, round(0.060 * fs));
        MinPeakProm     = 0.2;
        
        [spikes, pos, width, prom] = findpeaks( ...
            abs_mV, ...
            'MinPeakDistance',   MinPeakDistance, ...
            'MinPeakProminence', MinPeakProm, ...
            'MinPeakHeight',     Lower_Threshold_mV, ...
            'MaxPeakWidth',      MaxPeakWidth);
        
        keep  = spikes <= Upper_Threshold_mV;
        spikes = spikes(keep);
        pos    = pos(keep);
        width  = width(keep);
        prom   = prom(keep);
        
        initial_spike_count = numel(spikes);
        pos_s = (pos-1)/fs;
        
        fprintf('\n[INITIAL SPIKE DETECTION]\n');
        fprintf('  Raw peak count: %d\n', initial_spike_count);
        fprintf('  Detection rate: %.2f spikes/min\n', initial_spike_count / total_recording_min);
        
        %% GROUP SPIKES INTO COMPLEXES
        max_interspike_ms = 150;
        max_interspike_samples = round((max_interspike_ms/1000) * fs);
        
        if isempty(pos)
            Spike_Complex_Record = zeros(0,8);
            polyspike_indices = [];
        else
            ISI_samples = diff(pos);
            complex_id = 1;
            spike_complex_ids = ones(size(pos));
            
            for i = 1:numel(ISI_samples)
                if ISI_samples(i) > max_interspike_samples
                    complex_id = complex_id + 1;
                end
                spike_complex_ids(i+1) = complex_id;
            end
            
            num_complexes = max(spike_complex_ids);
            Spike_Complex_Record = zeros(num_complexes, 8);
            
            for c = 1:num_complexes
                complex_spikes = spike_complex_ids == c;
                complex_indices = find(complex_spikes);
                
                n_spikes_in_complex = sum(complex_spikes);
                first_spike_idx = complex_indices(1);
                last_spike_idx = complex_indices(end);
                
                start_time_s = pos_s(first_spike_idx);
                end_time_s = pos_s(last_spike_idx);
                duration_ms = (end_time_s - start_time_s) * 1000;
                
                max_amplitude = max(spikes(complex_spikes));
                mean_amplitude = mean(spikes(complex_spikes));
                is_polyspike = n_spikes_in_complex >= 2;
                
                Spike_Complex_Record(c,:) = [...
                    c, start_time_s, end_time_s, duration_ms, ...
                    n_spikes_in_complex, max_amplitude, mean_amplitude, is_polyspike];
            end
            
            polyspike_indices = find(Spike_Complex_Record(:,8) == 1);
        end
        
        num_spike_complexes = size(Spike_Complex_Record, 1);
        num_polyspikes = numel(polyspike_indices);
        num_single_spikes = num_spike_complexes - num_polyspikes;
        
        fprintf('\n[SPIKE COMPLEX ANALYSIS]\n');
        fprintf('  Total spike complexes:   %d\n', num_spike_complexes);
        fprintf('    - Single spikes:       %d (%.1f%%)\n', num_single_spikes, ...
            100*num_single_spikes/max(1,num_spike_complexes));
        fprintf('    - Polyspike complexes: %d (%.1f%%)\n', num_polyspikes, ...
            100*num_polyspikes/max(1,num_spike_complexes));
        fprintf('  Complexes/min:           %.2f\n', num_spike_complexes / total_recording_min);
        fprintf('  Polyspikes/min:          %.2f\n', num_polyspikes / total_recording_min);
        
        %% Create summary tables
        Individual_Peaks_Table = table(pos_s(:), spikes(:), width(:), prom(:), ...
            'VariableNames', {'Time_s','Amplitude_mV','Width_samples','Prominence_mV'});
        
        if ~isempty(Spike_Complex_Record)
            Spike_Complex_Table = table( ...
                Spike_Complex_Record(:,1), Spike_Complex_Record(:,2), ...
                Spike_Complex_Record(:,3), Spike_Complex_Record(:,4), ...
                Spike_Complex_Record(:,5), Spike_Complex_Record(:,6), ...
                Spike_Complex_Record(:,7), logical(Spike_Complex_Record(:,8)), ...
                'VariableNames', {'Complex_ID','Start_s','End_s','Duration_ms', ...
                                 'N_spikes','Max_amplitude_mV','Mean_amplitude_mV','Is_polyspike'});
            
            Single_Spike_Table = Spike_Complex_Table(~Spike_Complex_Table.Is_polyspike, :);
            Polyspike_Table = Spike_Complex_Table(Spike_Complex_Table.Is_polyspike, :);
        else
            Spike_Complex_Table = table([], [], [], [], [], [], [], [], ...
                'VariableNames', {'Complex_ID','Start_s','End_s','Duration_ms', ...
                                 'N_spikes','Max_amplitude_mV','Mean_amplitude_mV','Is_polyspike'});
            Single_Spike_Table = Spike_Complex_Table;
            Polyspike_Table = Spike_Complex_Table;
        end
        
        %% Group into bursts
        maxGap_s = 5;
        
        if isempty(Spike_Complex_Record)
            Burst_Record = zeros(0,5);
        else
            complex_times_s = Spike_Complex_Record(:,2);
            Interval_dt = diff(complex_times_s);
            Status = Interval_dt <= maxGap_s;
            
            D = [];
            complex_tally = 1;
            dur = 0;
            t_ini = complex_times_s(1);
            
            for k = 1:numel(Interval_dt)
                if Status(k)
                    complex_tally = complex_tally + 1;
                    dur = dur + Interval_dt(k);
                else
                    D = [D; complex_tally, dur, t_ini, t_ini+dur, 0]; %#ok<AGROW>
                    complex_tally = 1;
                    dur = 0;
                    t_ini = complex_times_s(k+1);
                end
            end
            D = [D; complex_tally, dur, t_ini, t_ini+dur, 0];
            Burst_Record = D;
        end
        
        if ~isempty(Burst_Record)
            mask = Burst_Record(:,1) > 3 & ...
                   Burst_Record(:,2) > 4 & ...
                   Burst_Record(:,2) < 40;
            Burst_Record = Burst_Record(mask,:);
        end
        
        burst_count = size(Burst_Record,1);
        fprintf('\n[BURST DETECTION]\n');
        fprintf('  Bursts detected: %d\n', burst_count);
        fprintf('  Bursts/hour:     %.2f\n', burst_count / (total_recording_min/60));
        
        if ~isempty(Burst_Record)
            Burst_Table = table( ...
                Burst_Record(:,3), Burst_Record(:,4), ...
                Burst_Record(:,1), Burst_Record(:,2), ...
                'VariableNames', {'Start_s','End_s','N_complexes','Duration_s'});
        else
            Burst_Table = table([],[],[],[], ...
                'VariableNames', {'Start_s','End_s','N_complexes','Duration_s'});
        end
        
        %% PLOTS
        Time_total = (numel(x_bp)-1)/fs;
        Time = (0:1/fs:Time_total).';
        
        % Figure 1: Full view
        fig1 = figure('Name',sprintf('%s - IID Detection', fname_base), ...
            'Position', [100 100 1400 600], 'Visible', 'off');
        plot(Time, x_bp, 'k', 'DisplayName', 'Signal'); hold on
        
        single_spike_times = [];
        single_spike_values = [];
        if ~isempty(Single_Spike_Table)
            for i = 1:height(Single_Spike_Table)
                idx = find(abs(pos_s - Single_Spike_Table.Start_s(i)) < 0.001, 1);
                if ~isempty(idx)
                    single_spike_times(end+1) = pos_s(idx); %#ok<AGROW>
                    single_spike_values(end+1) = x_bp(pos(idx)); %#ok<AGROW>
                end
            end
        end
        
        poly_spike_times = [];
        poly_spike_values = [];
        first_box = true;
        if ~isempty(Polyspike_Table)
            for i = 1:height(Polyspike_Table)
                t_start = Polyspike_Table.Start_s(i);
                t_end = Polyspike_Table.End_s(i);
                
                complex_peak_mask = (pos_s >= t_start) & (pos_s <= t_end);
                complex_peak_indices = find(complex_peak_mask);
                
                for j = complex_peak_indices'
                    poly_spike_times(end+1) = pos_s(j); %#ok<AGROW>
                    poly_spike_values(end+1) = x_bp(pos(j)); %#ok<AGROW>
                end
                
                y_range = [min(x_bp(pos(complex_peak_indices))) max(x_bp(pos(complex_peak_indices)))];
                y_margin = 0.1 * diff(y_range);
                rectangle('Position', [t_start, y_range(1)-y_margin, t_end-t_start, diff(y_range)+2*y_margin], ...
                          'EdgeColor', 'r', 'LineWidth', 1.5, 'LineStyle', '--');
                
                if first_box
                    plot(NaN, NaN, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Polyspike boundaries');
                    first_box = false;
                end
            end
        end
        
        if ~isempty(single_spike_times)
            plot(single_spike_times, single_spike_values, 'bo', ...
                'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'Single spikes');
        end
        
        if ~isempty(poly_spike_times)
            plot(poly_spike_times, poly_spike_values, 'ro', ...
                'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', 'Polyspikes');
        end
        
        if burst_count > 0
            plot([Burst_Record(1,3) Burst_Record(1,4)], [1 1]*1.2, ...
                'm', 'LineWidth', 3, 'DisplayName', 'Bursts');
            for i=2:burst_count
                plot([Burst_Record(i,3) Burst_Record(i,4)], [1 1]*1.2, ...
                    'm', 'LineWidth', 3, 'HandleVisibility', 'off');
            end
        end
        
        yline(Lower_Threshold_mV, 'r-', 'LineWidth', 1, 'HandleVisibility', 'off');
        yline(-Lower_Threshold_mV, 'r-', 'LineWidth', 1, 'HandleVisibility', 'off');
        xlabel('Time (s)');
        ylabel('Voltage (mV)');
        title(sprintf('%s - %s - %s | Blue = single spikes, Red = polyspikes', ...
            mouse_id, region, session), 'Interpreter', 'none');
        legend('Location', 'best');
        grid on;
        hold off;
        
        % Save figure
        savefig(fig1, fullfile(output_folder, sprintf('%s_IID_overview.fig', fname_base)));
        close(fig1);
        
        % Figure 2: Polyspike examples
        if ~isempty(Polyspike_Table) && height(Polyspike_Table) > 0
            n_examples = min(4, height(Polyspike_Table));
            fig2 = figure('Name',sprintf('%s - Polyspike Examples', fname_base), ...
                'Position', [200 200 1400 800], 'Visible', 'off');
            
            for i = 1:n_examples
                subplot(2, 2, i);
                
                t_center = (Polyspike_Table.Start_s(i) + Polyspike_Table.End_s(i)) / 2;
                window_s = 0.5;
                
                t_start_plot = max(0, t_center - window_s);
                t_end_plot = min(Time(end), t_center + window_s);
                
                idx_plot = (Time >= t_start_plot) & (Time <= t_end_plot);
                plot(Time(idx_plot), x_bp(idx_plot), 'k'); hold on;
                
                complex_peak_mask = (pos_s >= Polyspike_Table.Start_s(i)) & ...
                                   (pos_s <= Polyspike_Table.End_s(i));
                complex_peak_indices = find(complex_peak_mask);
                
                for j = complex_peak_indices'
                    plot(pos_s(j), x_bp(pos(j)), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
                end
                
                xline(Polyspike_Table.Start_s(i), '--r', 'LineWidth', 1.5);
                xline(Polyspike_Table.End_s(i), '--r', 'LineWidth', 1.5);
                
                yline(Lower_Threshold_mV, ':r');
                yline(-Lower_Threshold_mV, ':r');
                
                title(sprintf('Polyspike #%d: %d spikes, %.1f ms', ...
                    i, Polyspike_Table.N_spikes(i), Polyspike_Table.Duration_ms(i)));
                xlabel('Time (s)');
                ylabel('Voltage (mV)');
                grid on;
                hold off;
            end
            
            savefig(fig2, fullfile(output_folder, sprintf('%s_IID_polyspike_examples.fig', fname_base)));
            close(fig2);
        end
        
        %% Save results
        
        % ═══════════════════════════════════════════════════════════════
        % OUTPUT 2: Individual file CSV - ALL spike complexes (single + polyspikes)
        % Format: (original_filename)_IID.csv
        % ═══════════════════════════════════════════════════════════════
        
        % Create combined table with all spike complexes
        if ~isempty(Spike_Complex_Table)
            % Create the output table with required columns
            Individual_IID_Table = table( ...
                Spike_Complex_Table.Complex_ID, ...
                Spike_Complex_Table.Start_s, ...
                Spike_Complex_Table.End_s, ...
                Spike_Complex_Table.Duration_ms / 1000, ... % Convert to seconds
                Spike_Complex_Table.Max_amplitude_mV, ...
                Spike_Complex_Table.Mean_amplitude_mV, ...
                Spike_Complex_Table.Is_polyspike, ...
                Spike_Complex_Table.N_spikes, ...
                'VariableNames', {'complex_ID', 'start_s', 'end_s', 'duration_s', ...
                                 'max_amplitude_mV', 'mean_amplitude_mV', ...
                                 'is_polyspike', 'n_spikes'});
        else
            Individual_IID_Table = table([], [], [], [], [], [], [], [], ...
                'VariableNames', {'complex_ID', 'start_s', 'end_s', 'duration_s', ...
                                 'max_amplitude_mV', 'mean_amplitude_mV', ...
                                 'is_polyspike', 'n_spikes'});
        end
        
        % Save individual file CSV with original filename + _IID
        individual_csv_name = sprintf('%s_IID.csv', fname_base);
        writetable(Individual_IID_Table, fullfile(output_folder, individual_csv_name));
        
        % ═══════════════════════════════════════════════════════════════
        % Save .mat file for additional analysis
        % ═══════════════════════════════════════════════════════════════
        IID_results = struct();
        IID_results.mouse_id = mouse_id;
        IID_results.session = session;
        IID_results.region = region;
        IID_results.fs = fs;
        IID_results.total_recording_min = total_recording_min;
        IID_results.baseline_uV = Baseline_uV;
        IID_results.lower_threshold_uV = Lower_Threshold_uV;
        IID_results.upper_threshold_uV = Upper_Threshold_uV;
        IID_results.spike_complex_table = Spike_Complex_Table;
        IID_results.single_spike_table = Single_Spike_Table;
        IID_results.polyspike_table = Polyspike_Table;
        IID_results.burst_table = Burst_Table;
        IID_results.individual_peaks_table = Individual_Peaks_Table;
        
        save(fullfile(output_folder, sprintf('%s_IID_results.mat', fname_base)), 'IID_results');
        
        fprintf('\n[FILES SAVED]\n');
        fprintf('  - %s\n', individual_csv_name);
        fprintf('  - %s_IID_results.mat\n', fname_base);
        fprintf('  - %s_IID_overview.fig\n', fname_base);
        if ~isempty(Polyspike_Table)
            fprintf('  - %s_IID_polyspike_examples.fig\n', fname_base);
        end
        
        % Store summary (for OUTPUT 1: overview CSV)
        batch_summary{file_idx, 1} = files(file_idx).name;
        batch_summary{file_idx, 2} = session;
        batch_summary{file_idx, 3} = total_recording_min;
        batch_summary{file_idx, 4} = mouse_id;
        batch_summary{file_idx, 5} = region;
        batch_summary{file_idx, 6} = num_single_spikes;
        batch_summary{file_idx, 7} = num_spike_complexes;
        batch_summary{file_idx, 8} = num_spike_complexes / total_recording_min;  % spike_complexes_per_min
        batch_summary{file_idx, 9} = num_polyspikes;
        batch_summary{file_idx, 10} = num_polyspikes / total_recording_min;      % polyspikes_per_min
        batch_summary{file_idx, 11} = burst_count;
        batch_summary{file_idx, 12} = burst_count / (total_recording_min/60);    % bursts_per_hour
        
        fprintf('\n[COMPLETED]\n');
        
    catch ME
        fprintf('\n[ERROR]\n');
        fprintf('  Message: %s\n', ME.message);
        fprintf('  File: %s\n', ME.stack(1).file);
        fprintf('  Line: %d\n', ME.stack(1).line);
        
        % Store error
        batch_summary{file_idx, 1} = files(file_idx).name;
        batch_summary{file_idx, 2} = 'ERROR';
        batch_summary{file_idx, 3} = NaN;
        batch_summary{file_idx, 4} = 'ERROR';
        batch_summary{file_idx, 5} = 'ERROR';
        batch_summary{file_idx, 6} = NaN;
        batch_summary{file_idx, 7} = NaN;
        batch_summary{file_idx, 8} = NaN;
        batch_summary{file_idx, 9} = NaN;
        batch_summary{file_idx, 10} = NaN;
        batch_summary{file_idx, 11} = NaN;
        batch_summary{file_idx, 12} = NaN;
    end
end

%% FINAL SUMMARY
fprintf('\n');
fprintf('==========================================================\n');
fprintf('BATCH PROCESSING COMPLETE\n');
fprintf('==========================================================\n\n');

% ═══════════════════════════════════════════════════════════════
% OUTPUT 1: Summary CSV with all files
% Format: (mouse_ID)_IID_summary.csv
% ═══════════════════════════════════════════════════════════════

summary_table = cell2table(batch_summary, ...
    'VariableNames', {'filename', 'session', 'total_recording_min', 'mouse_id', 'region', ...
                      'total_single_spikes', 'total_spike_complexes', 'spike_complexes_per_min', ...
                      'total_polyspikes', 'polyspikes_per_min', 'total_bursts', 'bursts_per_hour'});

disp(summary_table);

% Determine output filename based on first successful file's mouse_ID
summary_filename = 'IID_summary.csv';  % Default
for i = 1:num_files
    if ~strcmp(batch_summary{i,2}, 'ERROR')
        first_mouse_id = batch_summary{i,4};
        summary_filename = sprintf('%s_IID_summary.csv', first_mouse_id);
        break;
    end
end

% Save OUTPUT 1
writetable(summary_table, fullfile(output_folder, summary_filename));
fprintf('\n═══════════════════════════════════════════════════════════\n');
fprintf('OUTPUT 1 saved: %s\n', summary_filename);
fprintf('═══════════════════════════════════════════════════════════\n');

% Statistics
successful = sum(~strcmp(batch_summary(:,2), 'ERROR'));
failed = num_files - successful;

fprintf('\nSTATISTICS:\n');
fprintf('  Files processed: %d\n', num_files);
fprintf('  Successful:      %d\n', successful);
fprintf('  Failed:          %d\n', failed);

if successful > 0
    valid_data = ~strcmp(batch_summary(:,2), 'ERROR');
    total_single = sum([batch_summary{valid_data,6}]);
    total_complexes = sum([batch_summary{valid_data,7}]);
    total_polyspikes = sum([batch_summary{valid_data,9}]);
    total_bursts = sum([batch_summary{valid_data,11}]);
    
    fprintf('\n  Total single spikes:       %d\n', total_single);
    fprintf('  Total spike complexes:     %d\n', total_complexes);
    fprintf('  Total polyspikes:          %d\n', total_polyspikes);
    fprintf('  Total bursts:              %d\n', total_bursts);
end

fprintf('\n*** IID BATCH PROCESSING COMPLETE ***\n\n');

%% Local functions
function d = getLastDir(p)
    [~,d] = fileparts(p);
end