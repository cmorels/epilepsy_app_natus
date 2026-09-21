% IID_detection_single_with_exclusions.m
% Single-file interictal spike detection with polyspike complex identification
% NOW WITH SEIZURE EXCLUSION ZONES
% - Manual file selection via dialog
% - Excludes user-defined time periods (e.g., seizure epochs)
% - Works with *_clean.txt files (post-outlier removal)

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

%% =====================================================================
%  USER CONFIGURATION: EXCLUSION ZONES (SEIZURE PERIODS)
%  =====================================================================
%  Define time periods to EXCLUDE from IID detection (e.g., seizures)
%  Format: [start_time_s, end_time_s] for each exclusion zone
%  
%  Examples:
%  exclusion_zones = [];                              % No exclusions
%  exclusion_zones = [100, 150];                      % Exclude 100-150s
%  exclusion_zones = [100, 150; 300, 350; 500, 520]; % Multiple zones
%  
%  You can also load from seizure detection results:
%  load('yourfile_seizures.mat', 'seizures');
%  exclusion_zones = [seizures.start_time_s, seizures.end_time_s];

exclusion_zones = [465.8, 510.3]; % <-- EDIT THIS

% Optional: Add buffer around exclusion zones (seconds)
exclusion_buffer_s = 5;  % Add 5s before/after each exclusion zone

%% =====================================================================

%% ----------------- 1) Select file manually -----------------
[file, path] = uigetfile('*.txt', 'Select a cleaned LFP file (*_clean.txt)');

if isequal(file, 0)
    fprintf('User canceled file selection.\n');
    return;
end

file_path = fullfile(path, file);
fprintf('\nSelected file: %s\n', file);

%% ----------------- 2) Load data -----------------
fprintf('\n=== LOADING DATA ===\n');
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

%% ----------------- 3) Process exclusion zones -----------------
if ~isempty(exclusion_zones)
    % Add buffer to exclusion zones
    if exclusion_buffer_s > 0
        exclusion_zones(:,1) = exclusion_zones(:,1) - exclusion_buffer_s;
        exclusion_zones(:,2) = exclusion_zones(:,2) + exclusion_buffer_s;
        % Clip to recording bounds
        exclusion_zones(:,1) = max(exclusion_zones(:,1), t(1));
        exclusion_zones(:,2) = min(exclusion_zones(:,2), t(end));
    end
    
    % Create inclusion mask (true = analyze, false = exclude)
    inclusion_mask = true(size(t));
    
    for i = 1:size(exclusion_zones, 1)
        zone_start = exclusion_zones(i, 1);
        zone_end = exclusion_zones(i, 2);
        
        % Mark samples in this zone as excluded
        in_zone = (t >= zone_start) & (t <= zone_end);
        inclusion_mask(in_zone) = false;
    end
    
    n_excluded_samples = sum(~inclusion_mask);
    excluded_duration_min = n_excluded_samples / fs / 60;
    
    fprintf('\n=== EXCLUSION ZONES ===\n');
    fprintf('  Number of zones:      %d\n', size(exclusion_zones, 1));
    fprintf('  Buffer applied:       %.1f s\n', exclusion_buffer_s);
    fprintf('  Excluded duration:    %.2f min (%.1f%% of recording)\n', ...
        excluded_duration_min, 100*excluded_duration_min/total_recording_min);
    fprintf('  Analyzed duration:    %.2f min\n', total_recording_min - excluded_duration_min);
    
    fprintf('\n  Exclusion zones:\n');
    for i = 1:size(exclusion_zones, 1)
        fprintf('    Zone %d: %.2f - %.2f s (%.1f s duration)\n', ...
            i, exclusion_zones(i,1), exclusion_zones(i,2), ...
            exclusion_zones(i,2) - exclusion_zones(i,1));
    end
else
    inclusion_mask = true(size(t));
    fprintf('\n=== NO EXCLUSION ZONES ===\n');
    fprintf('  Analyzing entire recording\n');
end

%% ----------------- 4) BASELINE calculation -----------------
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

fprintf('\n=== BASELINE DETECTION ===\n');
fprintf('  Baseline:         %.0f µV\n', Baseline_uV);
fprintf('  Lower threshold:  %.0f µV\n', Lower_Threshold_uV);
fprintf('  Upper threshold:  %.0f µV\n', Upper_Threshold_uV);

%% ----------------- 5) Band-pass filter 15–70 Hz -----------------
x_mV  = (sig_uV - median(sig_uV,'omitnan'))/1000;
x_bp  = bandpass(x_mV, [15 70], fs, 'Steepness', 0.5);
abs_mV = abs(x_bp);

%% ----------------- 6) Initial spike detection (individual peaks) -----------------
MinPeakDistance = max(1, round(0.030 * fs));
MaxPeakWidth    = max(1, round(0.060 * fs));
MinPeakProm     = 0.2;

[spikes, pos, width, prom] = findpeaks( ...
    abs_mV, ...
    'MinPeakDistance',   MinPeakDistance, ...
    'MinPeakProminence', MinPeakProm, ...
    'MinPeakHeight',     Lower_Threshold_mV, ...
    'MaxPeakWidth',      MaxPeakWidth);

% Filter by upper threshold
keep  = spikes <= Upper_Threshold_mV;
spikes = spikes(keep);
pos    = pos(keep);
width  = width(keep);
prom   = prom(keep);

% *** CRITICAL: APPLY EXCLUSION MASK TO DETECTED SPIKES ***
spikes_in_included_regions = inclusion_mask(pos);
spikes = spikes(spikes_in_included_regions);
pos = pos(spikes_in_included_regions);
width = width(spikes_in_included_regions);
prom = prom(spikes_in_included_regions);

initial_spike_count = numel(spikes);
pos_s = (pos-1)/fs;

% Calculate analyzed time for rate calculations
analyzed_duration_min = sum(inclusion_mask) / fs / 60;

fprintf('\n=== INITIAL SPIKE DETECTION ===\n');
fprintf('  Raw peak count:   %d\n', initial_spike_count);
fprintf('  Detection rate:   %.2f spikes/min (analyzed periods only)\n', ...
    initial_spike_count / analyzed_duration_min);

%% ----------------- 7) GROUP SPIKES INTO COMPLEXES -----------------
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

fprintf('\n=== SPIKE COMPLEX ANALYSIS ===\n');
fprintf('  Total spike complexes:   %d\n', num_spike_complexes);
fprintf('    - Single spikes:       %d (%.1f%%)\n', num_single_spikes, ...
    100*num_single_spikes/max(1,num_spike_complexes));
fprintf('    - Polyspike complexes: %d (%.1f%%)\n', num_polyspikes, ...
    100*num_polyspikes/max(1,num_spike_complexes));
fprintf('  Complexes/min:           %.2f (analyzed periods only)\n', ...
    num_spike_complexes / analyzed_duration_min);
fprintf('  Polyspikes/min:          %.2f (analyzed periods only)\n', ...
    num_polyspikes / analyzed_duration_min);

%% ----------------- 8) Create summary tables -----------------
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

fprintf('\n=== SPIKE COMPLEX TABLE (first 10) ===\n');
disp(Spike_Complex_Table(1:min(10,height(Spike_Complex_Table)),:));

if ~isempty(Polyspike_Table)
    fprintf('\n=== POLYSPIKE COMPLEXES ===\n');
    disp(Polyspike_Table);
    
    if height(Polyspike_Table) > 0
        fprintf('\n=== POLYSPIKE STATISTICS ===\n');
        fprintf('  Mean spikes per complex: %.1f\n', mean(Polyspike_Table.N_spikes));
        fprintf('  Max spikes in complex:   %d\n', max(Polyspike_Table.N_spikes));
        fprintf('  Mean duration:           %.1f ms\n', mean(Polyspike_Table.Duration_ms));
        fprintf('  Mean amplitude:          %.3f mV\n', mean(Polyspike_Table.Mean_amplitude_mV));
    end
end

%% ----------------- 9) Group into bursts -----------------
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
fprintf('\n=== BURST DETECTION ===\n');
fprintf('  Bursts detected: %d\n', burst_count);
fprintf('  Bursts/hour:     %.2f (analyzed periods only)\n', ...
    burst_count / (analyzed_duration_min/60));

if ~isempty(Burst_Record)
    Burst_Table = table( ...
        Burst_Record(:,3), Burst_Record(:,4), ...
        Burst_Record(:,1), Burst_Record(:,2), ...
        'VariableNames', {'Start_s','End_s','N_complexes','Duration_s'});
    disp(Burst_Table);
else
    Burst_Table = table([],[],[],[], ...
        'VariableNames', {'Start_s','End_s','N_complexes','Duration_s'});
end

%% ----------------- 10) PLOTS -----------------
Time_total = (numel(x_bp)-1)/fs;
Time = (0:1/fs:Time_total).';

% Figure 1: Full view with exclusion zones highlighted
figure('Name','IID Detection with Exclusion Zones', 'Position', [100 100 1400 600]);
plot(Time, x_bp, 'k', 'DisplayName', 'Signal'); hold on

% Highlight exclusion zones
if ~isempty(exclusion_zones)
    for i = 1:size(exclusion_zones, 1)
        xregion(exclusion_zones(i,1), exclusion_zones(i,2), ...
            'FaceColor', [0.9 0.9 0.9], 'FaceAlpha', 0.5, ...
            'EdgeColor', 'none');
    end
    % Add to legend (only once)
    plot(NaN, NaN, 's', 'Color', [0.9 0.9 0.9], 'MarkerFaceColor', [0.9 0.9 0.9], ...
        'MarkerSize', 10, 'DisplayName', 'Excluded zones');
end

% Plot single spike positions
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

% Plot polyspike complexes with boxes
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
title(sprintf('%s - %s - %s | Blue = single, Red = polyspikes, Gray = excluded', ...
    mouse_id, region, session), 'Interpreter', 'none');
legend('Location', 'best');
grid on;
hold off;

% Figure 2: Polyspike examples (zoomed view)
if ~isempty(Polyspike_Table) && height(Polyspike_Table) > 0
    n_examples = min(4, height(Polyspike_Table));
    figure('Name','Polyspike Complex Examples', 'Position', [200 200 1400 800]);
    
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
end

%% ----------------- 11) Save results -----------------
try
    [~, fname, ~] = fileparts(data.file);
    fname_base = strrep(fname, '_clean', '');
    
    % ═══════════════════════════════════════════════════════════════
    % OUTPUT 2: Individual file CSV - ALL spike complexes (single + polyspikes)
    % Format: (original_filename)_IID.csv
    % ═══════════════════════════════════════════════════════════════
    
    % Create combined table with all spike complexes
    if ~isempty(Spike_Complex_Table)
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
    
    % Save OUTPUT 2
    individual_csv_name = sprintf('%s_IID.csv', fname_base);
    writetable(Individual_IID_Table, individual_csv_name);
    
    % ═══════════════════════════════════════════════════════════════
    % OUTPUT 1: Summary CSV for this file
    % Format: (mouse_ID)_IID_summary.csv
    % ═══════════════════════════════════════════════════════════════
    
    % Calculate rates
    spike_complexes_per_min = num_spike_complexes / analyzed_duration_min;
    polyspikes_per_min = num_polyspikes / analyzed_duration_min;
    bursts_per_hour = burst_count / (analyzed_duration_min / 60);
    
    % Create summary table for this file
    summary_data = {
        data.file, ...                      % filename
        session, ...                        % session
        analyzed_duration_min, ...          % total_recording_min (analyzed only)
        mouse_id, ...                       % mouse_id
        region, ...                         % region
        num_single_spikes, ...              % total_single_spikes
        num_spike_complexes, ...            % total_spike_complexes
        spike_complexes_per_min, ...        % spike_complexes_per_min
        num_polyspikes, ...                 % total_polyspikes
        polyspikes_per_min, ...             % polyspikes_per_min
        burst_count, ...                    % total_bursts
        bursts_per_hour                     % bursts_per_hour
    };
    
    summary_table = cell2table(summary_data, ...
        'VariableNames', {'filename', 'session', 'total_recording_min', 'mouse_id', 'region', ...
                          'total_single_spikes', 'total_spike_complexes', 'spike_complexes_per_min', ...
                          'total_polyspikes', 'polyspikes_per_min', 'total_bursts', 'bursts_per_hour'});
    
    % Save OUTPUT 1
    summary_csv_name = sprintf('%s_IID_summary.csv', mouse_id);
    writetable(summary_table, summary_csv_name);
    
    % ═══════════════════════════════════════════════════════════════
    % Save exclusion zones info (if any)
    % ═══════════════════════════════════════════════════════════════
    if ~isempty(exclusion_zones)
        Exclusion_Table = table( ...
            (1:size(exclusion_zones,1))', ...
            exclusion_zones(:,1), ...
            exclusion_zones(:,2), ...
            exclusion_zones(:,2) - exclusion_zones(:,1), ...
            'VariableNames', {'Zone_ID','Start_s','End_s','Duration_s'});
        writetable(Exclusion_Table, sprintf('%s_IID_exclusion_zones.csv', fname_base));
    end
    
    % ═══════════════════════════════════════════════════════════════
    % Save .mat file for additional analysis
    % ═══════════════════════════════════════════════════════════════
    IID_results = struct();
    IID_results.mouse_id = mouse_id;
    IID_results.session = session;
    IID_results.region = region;
    IID_results.fs = fs;
    IID_results.total_recording_min = total_recording_min;
    IID_results.analyzed_duration_min = analyzed_duration_min;
    IID_results.excluded_duration_min = total_recording_min - analyzed_duration_min;
    IID_results.exclusion_zones = exclusion_zones;
    IID_results.exclusion_buffer_s = exclusion_buffer_s;
    IID_results.baseline_uV = Baseline_uV;
    IID_results.lower_threshold_uV = Lower_Threshold_uV;
    IID_results.upper_threshold_uV = Upper_Threshold_uV;
    IID_results.spike_complex_table = Spike_Complex_Table;
    IID_results.single_spike_table = Single_Spike_Table;
    IID_results.polyspike_table = Polyspike_Table;
    IID_results.burst_table = Burst_Table;
    IID_results.individual_peaks_table = Individual_Peaks_Table;
    
    save(sprintf('%s_IID_results.mat', fname_base), 'IID_results');
    
    fprintf('\n═══════════════════════════════════════════════════════════\n');
    fprintf('FILES SAVED\n');
    fprintf('═══════════════════════════════════════════════════════════\n');
    fprintf('  OUTPUT 1 (summary): %s\n', summary_csv_name);
    fprintf('  OUTPUT 2 (details): %s\n', individual_csv_name);
    if ~isempty(exclusion_zones)
        fprintf('  Exclusion zones:    %s_IID_exclusion_zones.csv\n', fname_base);
    end
    fprintf('  MAT file:           %s_IID_results.mat\n', fname_base);
    fprintf('═══════════════════════════════════════════════════════════\n');
    
catch ME
    warning('Could not save files: %s', ME.message);
end

fprintf('\n*** IID DETECTION COMPLETE ***\n\n');

%% ----------------- Local functions -----------------
function d = getLastDir(p)
    [~,d] = fileparts(p);
end