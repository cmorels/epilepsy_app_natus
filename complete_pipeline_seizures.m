%% PIPELINE BATCH: Load → Remove Outliers → Seizure Detection
% Procesa automáticamente todos los archivos .txt en la carpeta
% FIGURAS GUARDADAS EN FORMATO .FIG (MATLAB nativo)
% 
% Pipeline:
%   1. Load LFP from .txt
%   2. Remove outliers → *_clean.txt (CON ANÁLISIS DETALLADO)
%   3. Seizure detection → *_seizures.mat (SIMPLIFIED METHOD)

clear; clc;

%% CONFIGURATION

folder_path = '';  % '' -> carpeta actual

% OUTLIER REMOVAL SETTINGS
outlier_config = struct();
outlier_config.use_fixed_thresholds = true;
outlier_config.fixed_lower_uV = -2500;
outlier_config.fixed_upper_uV = +2500;
outlier_config.k_factor = 10;

% SEIZURE DETECTION SETTINGS 
seizure_config = struct();
seizure_config.bandpass_band = [5 75];
seizure_config.power_exponent = 4;            
seizure_config.window_sec = 2;                % smoothing window
seizure_config.median_factor = 10;           % threshold = median(energy) × factor
seizure_config.min_seizure_duration = 15;     % minimum sustained duration 

%% FIND FILES

if isempty(folder_path)
    folder_path = pwd;
end

% find/create batch_results subfolder
batch_results_folder = fullfile(folder_path, 'batch_results');
if ~exist(batch_results_folder, 'dir')
    mkdir(batch_results_folder);
end

% buscar archivos .txt 
all_files = dir(fullfile(folder_path, '*.txt'));
files = all_files(~contains({all_files.name}, '_clean.txt'));

if isempty(files)
    error('No .txt files found in: %s', folder_path);
end

num_files = length(files);

fprintf('BATCH PROCESSING PIPELINE------------------------------------\n');
fprintf('\nFolder: %s\n', folder_path);
fprintf('Files found: %d\n', num_files);
fprintf('\nPipeline:\n');
fprintf('  1. Load LFP data\n');
if outlier_config.use_fixed_thresholds
    fprintf('  2. Remove outliers: FIXED thresholds [%.0f, %.0f] µV\n', ...
        outlier_config.fixed_lower_uV, outlier_config.fixed_upper_uV);
else
    fprintf('  2. Remove outliers: AUTOMATIC thresholds (k=%.1f)\n', ...
        outlier_config.k_factor);
end
fprintf('  3. Detect seizures: median × %.1f, ≥%.0fs duration \n', ...
    seizure_config.median_factor, seizure_config.min_seizure_duration);
fprintf('\n');

for i = 1:num_files
    fprintf('%2d. %s\n', i, files(i).name);
end

% fprintf('\nPress any key to start processing...\n');
% pause;


%% BATCH PROCESSING
outlier_summary = cell(num_files, 12);  
seizure_summary = cell(0, 15);           

for file_idx = 1:num_files
    
    fprintf('\n');
    fprintf('---------------------------------------------------------\n');
    fprintf('FILE %d/%d: %s\n', file_idx, num_files, files(file_idx).name);
    
    file_path = fullfile(files(file_idx).folder, files(file_idx).name);
    [~, fname, ~] = fileparts(files(file_idx).name);
    
    try        
        %% STEP 1: LOAD DATA      
        fprintf('\n[1/3] LOADING DATA...\n');         
        data        = load_LFP_intan_txt(file_path);
        session     = data.session_time;
        mouse_id    = data.mouse_id;     
        region      = data.region;
        fs          = data.fs;
        signal_uV   = data.signal;
        t           = data.time_seconds;
        
        total_recording_min = max(t) / 60;
        n_samples = numel(signal_uV);
        
        fprintf('      Duration: %.1f min\n', total_recording_min);
        fprintf('      Samples: %d\n', n_samples);
        fprintf('      Fs: %.2f Hz\n', fs);
        
      
        %% STEP 2: REMOVE OUTLIERS
        
        fprintf('\n[2/3] REMOVING OUTLIERS...\n');
        
        % Calcular estadísticas de señal
        baseline_uV = median(signal_uV, 'omitnan');
        noise_uV = mad(signal_uV, 1);
        if noise_uV == 0
            warning('MAD = 0; usando valor mínimo de ruido (1 µV).');
            noise_uV = 1;
        end
        
        signal_std = std(signal_uV);
        signal_range = range(signal_uV);
        
        % Calcular umbrales automáticos (siempre, para comparación)
        min_auto_uV = baseline_uV - outlier_config.k_factor * noise_uV;
        max_auto_uV = baseline_uV + outlier_config.k_factor * noise_uV;
            
        % NUEVA LÓGICA: Decidir umbrales con auto-switch
        trigger_fixed = outlier_config.use_fixed_thresholds;
        threshold_type = '';
        switched_reason = '';
        
        if ~outlier_config.use_fixed_thresholds
            % Modo automático inicial: hacer detección preliminar
            [~, outlier_info_test] = remove_outliers_amplitude( ...
                signal_uV, min_auto_uV, max_auto_uV);
            
            outlier_pct_test = 100 * size(outlier_info_test, 1) / numel(signal_uV);
            
            % CRITERIO AUTO-SWITCH: si >0.1% outliers, trigger fixed
            if outlier_pct_test > 0.1
                trigger_fixed = true;
                switched_reason = sprintf('Auto-switched (%.2f%% > 0.1%%)', outlier_pct_test);
                
                fprintf('     AUTO-SWITCH TRIGGERED\n');
                fprintf('     Preliminary detection: %.2f%% outliers\n', outlier_pct_test);
                fprintf('     Switching to FIXED thresholds for this file.\n');
                fprintf('     Likely interictal spikes, not artifacts.\n');
            end
        end
        
        % Aplicar umbrales finales
        if trigger_fixed
            min_allowed = outlier_config.fixed_lower_uV;
            max_allowed = outlier_config.fixed_upper_uV;
            
            if isempty(switched_reason)
                threshold_type = 'FIXED (user config)';
            else
                threshold_type = sprintf('FIXED');
            end
        else
            min_allowed = min_auto_uV;
            max_allowed = max_auto_uV;
            threshold_type = sprintf('AUTOMATIC (k=%.1f)', outlier_config.k_factor);
        end
        
        % Mostrar estadísticas
        fprintf('\n      SIGNAL STATISTICS\n');
        fprintf('        Baseline (median):  %.1f µV\n', baseline_uV);
        fprintf('        Noise (MAD):        %.1f µV\n', noise_uV);
        fprintf('        Std deviation:      %.1f µV\n', signal_std);
        fprintf('        Range observed:     [%.1f, %.1f] µV\n', min(signal_uV), max(signal_uV));
        fprintf('        Dynamic range:      %.1f µV\n', signal_range);
        
        fprintf('\n      OUTLIER THRESHOLDS\n');
        fprintf('        Type: %s\n', threshold_type);
        fprintf('        Lower threshold: %.1f µV\n', min_allowed);
        fprintf('        Upper threshold: %.1f µV\n', max_allowed);
        
        if trigger_fixed
            fprintf('\n     Triggered FIXED thresholds.\n');
            fprintf('       Correcting technical artifacts only.\n');
            
            if ~outlier_config.use_fixed_thresholds
                fprintf('\n        (Auto thresholds values: [%.1f, %.1f] µV)\n', ...
                    min_auto_uV, max_auto_uV);
            end
        else
            fprintf('\n   Using AUTOMATIC thresholds (MAD-based).\n');       
        end
        
        [clean_signal_uV, outlier_info] = remove_outliers_amplitude( ...
            signal_uV, min_allowed, max_allowed);
        
        num_outliers = size(outlier_info, 1);
        outlier_pct = 100 * num_outliers / n_samples;
               
        fprintf('\n       OUTLIERS DETECTED:  %d (%.4f%%)\n', num_outliers, outlier_pct);   
       
        % Save clean signal
        clean_file = fullfile(batch_results_folder, sprintf('%s_clean.txt', fname));
        fid = fopen(clean_file, 'w');
        fprintf(fid, '# Original file: %s\n', files(file_idx).name);
        fprintf(fid, '# fs = %.10f\n', fs);
        fprintf(fid, '# time_unit = seconds\n');
        fprintf(fid, '# outliers_removed = %d samples (%.4f%%)\n', num_outliers, outlier_pct);
        fprintf(fid, '# removal_method = Amplitude thresholds\n');
        fprintf(fid, '# k_factor = %.1f\n', outlier_config.k_factor);
        fprintf(fid, '# threshold_type = %s\n', threshold_type);
        fprintf(fid, '# lower_threshold = %.1f microvolts\n', min_allowed);
        fprintf(fid, '# upper_threshold = %.1f microvolts\n', max_allowed);
        
        % Nueva línea: indicar si hubo auto-switch
        if ~isempty(switched_reason)
            fprintf(fid, '# auto_switch = true (reason: %s)\n', switched_reason);
        else
            fprintf(fid, '# auto_switch = false\n');
        end
        
        fprintf(fid, '# signal_range_original = [%.2f, %.2f] microvolts\n', ...
            min(signal_uV), max(signal_uV));
        fprintf(fid, '# signal_range_clean = [%.2f, %.2f] microvolts\n', ...
            min(clean_signal_uV), max(clean_signal_uV));
        fprintf(fid, '# columns = amplitude_microvolts\n');
        fprintf(fid, '%.6f\n', clean_signal_uV);
        fclose(fid);
        
        fprintf('\n      SAVED: %s_clean.txt\n', fname);
        
        
        %% STEP 3: SEIZURE DETECTION 
        
        fprintf('\n[3/3] SEIZURE DETECTION...\n');
        fprintf('       Method: Hilbert env. + energy metric + duration\n');  
        fprintf('       Config.: Power=%d, Window=%.1fs\n', ...
                       seizure_config.power_exponent, seizure_config.window_sec);
        fprintf('       Threshold: median × %.1f\n', seizure_config.median_factor);
        fprintf('       Minimum duration: ≥ %.0f s\n', seizure_config.min_seizure_duration);
 
        % Normalize and band-pass
        min_val = min(clean_signal_uV);
        max_val = max(clean_signal_uV);
        lfp_norm = (clean_signal_uV - min_val) / (max_val - min_val);
        
        lfp_bp = bandpass(lfp_norm, seizure_config.bandpass_band, fs, 'Steepness', 0.5);
        trim_samples = round(fs);
        lfp_bp = lfp_bp(trim_samples+1:end-trim_samples);
        t_bp = t(trim_samples+1:end-trim_samples);
        
        % Compute energy metric
        window_samples = round(seizure_config.window_sec * fs);
        envelope = abs(hilbert(lfp_bp));
       
            % Power → Smooth
            power = envelope .^ seizure_config.power_exponent;
            energy = movmean(power, window_samples);
   
        % Energy threshold
        threshold = median(energy) * seizure_config.median_factor;
        above_threshold = energy > threshold;
        
        n_samples_above = sum(above_threshold);
        pct_above = 100 * n_samples_above / numel(energy);
        
        fprintf('\n       Median energy:      %.2e\n', median(energy));
        fprintf('       Threshold value:    %.2e\n', threshold);
        fprintf('       Samples above thr:  %d (%.2f%%)\n', n_samples_above, pct_above);
        
        % Find continuous segments
        d = diff([0; above_threshold; 0]);
        starts = find(d == 1);
        ends = find(d == -1) - 1;
        
        durations_samples = ends - starts + 1;
        durations_sec = durations_samples / fs;
        
        fprintf('       Segments found:     %d\n', numel(starts));
        if ~isempty(durations_sec)
            fprintf('       Duration range:     %.2f - %.2f seconds\n', ...
                min(durations_sec), max(durations_sec));
        end
        
        % Duration filter
        if ~isempty(durations_sec)
            valid_mask = durations_sec >= seizure_config.min_seizure_duration;
            
            seizure_starts = starts(valid_mask);
            seizure_ends = ends(valid_mask);
            seizure_durations = durations_sec(valid_mask);
            
            rejected = sum(~valid_mask);
            fprintf('       Rejected (short):   %d segments\n', rejected);
            fprintf('       Valid seizures:     %d\n', sum(valid_mask));
        else
            seizure_starts = [];
            seizure_ends = [];
            seizure_durations = [];
        end
        
        % Build results table
        n_seizures = numel(seizure_starts);
        if n_seizures > 0
            seiz_start_time = t_bp(seizure_starts);
            seiz_end_time = t_bp(seizure_ends);
            seizures = table((1:n_seizures)', seiz_start_time(:), seiz_end_time(:), seizure_durations(:), ...
                'VariableNames', {'id','start_time_s','end_time_s','duration_s'});
            
            total_seizure_time = sum(seizure_durations);
        else
            seizures = table([],[],[],[], ...
                'VariableNames', {'id','start_time_s','end_time_s','duration_s'});
            total_seizure_time = 0;
        end
        
        fprintf('\n     SEIZURE DETECTION RESULTS    \n');
        fprintf('      Seizures detected:  %d\n', n_seizures);
        if n_seizures > 0
            fprintf('      Total seizure time: %.1f s (%.1f%% of recording)\n', ...
                total_seizure_time, 100*total_seizure_time/max(t_bp));
            fprintf('      Mean duration:      %.1f s\n', mean(seizure_durations));
            fprintf('      Seizure times:\n');
            for k = 1:n_seizures
                fprintf('       #%d: %.1f - %.1f s (%.1f s)\n', ...
                    k, seiz_start_time(k), seiz_end_time(k), seizure_durations(k));
            end
        end
        
      % save results 
        if n_seizures > 0
            results_file = fullfile(batch_results_folder, sprintf('%s_seizures.mat', fname));
            save(results_file, 'seizures', 'fs', 'seizure_config');
            fprintf('\n        SAVED: %s_seizures.mat\n', fname);
            
          
            fig_seizures = figure('Position', [50, 50, 1400, 900], 'Visible', 'off');
            
            % panel 1: cleaned LFP
            subplot(3,1,1)
            plot(t, clean_signal_uV, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5); hold on
            for k = 1:n_seizures
                start_t = seiz_start_time(k);
                end_t = seiz_end_time(k);
                xregion(start_t, end_t, 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
                    'EdgeColor', 'r', 'LineWidth', 2);
            end
            title(sprintf('%s - Cleaned LFP (%d seizures)', fname, n_seizures), 'Interpreter', 'none');
            xlabel('Time (s)'); ylabel('µV'); grid on; hold off;
            
            % panel 2: band-passed
            subplot(3,1,2)
            plot(t_bp, lfp_bp, 'Color', [0.2 0.4 0.7], 'LineWidth', 0.5); hold on
            for k = 1:n_seizures
                start_t = seiz_start_time(k);
                end_t = seiz_end_time(k);
                xregion(start_t, end_t, 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
                    'EdgeColor', 'r', 'LineWidth', 2);
            end
            title('Band-passed [5-75 Hz]');
            xlabel('Time (s)'); ylabel('Normalized'); grid on; hold off;
            
            % panel 3: energy metric with threshold
            subplot(3,1,3)
            plot(t_bp, energy, 'Color', [0.2 0.6 0.3], 'LineWidth', 0.5); hold on
            yline(threshold, 'r--', sprintf('Threshold (median × %.1f)', seizure_config.median_factor), ...
                'LineWidth', 2);
            plot(t_bp(above_threshold), energy(above_threshold), 'r.', 'MarkerSize', 3);
            
            % highlight detected seizures
            for k = 1:n_seizures
                start_t = seiz_start_time(k);
                end_t = seiz_end_time(k);
                xregion(start_t, end_t, 'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
                    'EdgeColor', 'r', 'LineWidth', 2);
            end
            
            title('Energy Metric');
            xlabel('Time (s)'); ylabel('Energy (a.u.)'); grid on; hold off;
            
            % save
            fig_seizures_file = fullfile(batch_results_folder, sprintf('%s_seizures.fig', fname));
            savefig(fig_seizures, fig_seizures_file);
            close(fig_seizures);

                % set(fig_seizures, 'Visible', 'on');
                % drawnow;
                % pause(0.1);
                % 
                % fig_seizures_file = fullfile(folder_path, sprintf('%s_seizures.fig', fname));
                % savefig(fig_seizures, fig_seizures_file);
                % close(fig_seizures);
            
            fprintf('      SAVED: %s_seizures.fig\n', fname);

        % zoom seizures
            for k = 1:n_seizures
                zoom_margin = 5; % s
                zoom_start = max(seiz_start_time(k) - zoom_margin, t_bp(1));
                zoom_end = min(seiz_end_time(k) + zoom_margin, t_bp(end));              
                idx_zoom = (t_bp >= zoom_start) & (t_bp <= zoom_end);
                idx_zoom_full = (t >= zoom_start) & (t <= zoom_end);

                fig_zoom = figure('Position', [50, 50, 1400, 900], 'Visible', 'off');
                
                % panel 1: full LFP 
                subplot(3,1,1)
                plot(t(idx_zoom_full), clean_signal_uV(idx_zoom_full), 'k', 'LineWidth', 0.8);
                hold on;  
                xregion(seiz_start_time(k), seiz_end_time(k), ...  % highlight seizure period
                    'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
                    'EdgeColor', 'r', 'LineWidth', 2);
                xline(seiz_start_time(k), 'r--', 'Start', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
                xline(seiz_end_time(k), 'r--', 'End', 'LineWidth', 2, 'LabelHorizontalAlignment', 'right');
                
                title(sprintf('Seizure #%d (full LFP)', k), 'Interpreter', 'none');
                xlabel('Time (s)');
                ylabel('Voltage (µV)');
                grid on;
                xlim([zoom_start, zoom_end]);
                hold off;
                
                % panel 2: band-passed signal 
                subplot(3,1,2)
                plot(t_bp(idx_zoom), lfp_bp(idx_zoom), 'k', 'LineWidth', 0.8);
                hold on;
                
                xregion(seiz_start_time(k), seiz_end_time(k), ...
                    'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
                    'EdgeColor', 'r', 'LineWidth', 2);
                
                xline(seiz_start_time(k), 'r--', 'LineWidth', 2);
                xline(seiz_end_time(k), 'r--', 'LineWidth', 2);
                
                title(sprintf('Seizure %d (band-passed [5-75 Hz])', k), 'Interpreter', 'none');
                xlabel('Time (s)');
                ylabel('Normalized Amplitude');
                grid on;
                xlim([zoom_start, zoom_end]);
                hold off;
                
                % panel 3: energy metric (zoomed)
                subplot(3,1,3)
                plot(t_bp(idx_zoom), energy(idx_zoom), 'Color', [0.2 0.6 0.3], 'LineWidth', 1.0);
                hold on;
                
                % threshold line
                yline(threshold, 'r--', sprintf('Threshold (%.2e)', threshold), 'LineWidth', 2);
                
                % highlight above-threshold region
                above_idx = idx_zoom & above_threshold';
                if any(above_idx)
                    plot(t_bp(above_idx), energy(above_idx), 'r.', 'MarkerSize', 6);
                end
                
                xregion(seiz_start_time(k), seiz_end_time(k), ...
                    'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
                    'EdgeColor', 'r', 'LineWidth', 2);
                
                xline(seiz_start_time(k), 'r--', 'LineWidth', 2);
                xline(seiz_end_time(k), 'r--', 'LineWidth', 2);
                
                title(sprintf('Hilbert envelope(^%d) over a %ds smoothed window', ...
                              seizure_config.power_exponent, seizure_config.window_sec));
                xlabel('Time (s)');
                ylabel('Energy (a.u.)');
                grid on;
                xlim([zoom_start, zoom_end]);
                hold off;
                
                % Save zoomed figure
                
                zoom_fig_file = fullfile(batch_results_folder, sprintf('%s_seizure%d.fig', fname, k));
                savefig(fig_zoom, zoom_fig_file);
                close(fig_zoom);

                    % set(fig_zoom, 'Visible', 'on');
                    % drawnow;
                    % pause(0.1);
                    % 
                    % zoom_fig_file = fullfile(folder_path, sprintf('%s_seizure%d.fig', fname, k));
                    % savefig(fig_zoom, zoom_fig_file);
                    % close(fig_zoom);
                
                fprintf('          Seizure #%d zoom saved: %s_seizure%d.fig\n', k, fname, k);
            end
        end
        
        %% SAVE TO SUMMARY
        
        % outliers 
        outlier_summary{file_idx, 1} = files(file_idx).name;
        outlier_summary{file_idx, 2} = session;
        outlier_summary{file_idx, 3} = mouse_id;
        outlier_summary{file_idx, 4} = region;
        outlier_summary{file_idx, 5} = total_recording_min;
        outlier_summary{file_idx, 6} = n_samples;
        outlier_summary{file_idx, 7} = baseline_uV;
        outlier_summary{file_idx, 8} = threshold_type;  
        outlier_summary{file_idx, 9} = min_allowed;     
        outlier_summary{file_idx, 10} = max_allowed;     
        outlier_summary{file_idx, 11} = num_outliers;
        outlier_summary{file_idx, 12} = outlier_pct;
        
        % seizures
        if n_seizures > 0
            for k = 1:n_seizures
                row_idx = size(seizure_summary, 1) + 1;  
                
                seizure_summary{row_idx, 1} = files(file_idx).name;
                seizure_summary{row_idx, 2} = session;
                seizure_summary{row_idx, 3} = mouse_id;
                seizure_summary{row_idx, 4} = region;
                seizure_summary{row_idx, 5} = total_recording_min;
                seizure_summary{row_idx, 6} = k;                           
                seizure_summary{row_idx, 7} = seizures.start_time_s(k);     
                seizure_summary{row_idx, 8} = seizures.end_time_s(k);       
                seizure_summary{row_idx, 9} = seizures.duration_s(k);

                % detection metrics (same for all seizures in this file)
                seizure_summary{row_idx, 10} = median(energy);
                seizure_summary{row_idx, 11} = seizure_config.median_factor; 
                seizure_summary{row_idx, 12} = threshold;               
                seizure_summary{row_idx, 13} = pct_above;                 
                seizure_summary{row_idx, 14} = numel(starts);             
                seizure_summary{row_idx, 15} = rejected;                  
            end
        end
        
        fprintf('\n       COMPLETED\n');
        
    catch ME
        fprintf('\n       ERROR: %s\n', ME.message);
        fprintf('         File: %s\n', ME.stack(1).file);
        fprintf('         Line: %d\n', ME.stack(1).line);
        
        % Store error in outlier summary
        outlier_summary{file_idx, 1} = files(file_idx).name;
        outlier_summary{file_idx, 2} = 'ERROR';
        outlier_summary{file_idx, 3} = 'ERROR';
        outlier_summary{file_idx, 4} = 'ERROR';
        outlier_summary{file_idx, 5} = NaN;
        outlier_summary{file_idx, 6} = NaN;
        outlier_summary{file_idx, 7} = NaN;
        outlier_summary{file_idx, 8} = ME.message;
        outlier_summary{file_idx, 9} = NaN;
        outlier_summary{file_idx, 10} = NaN;
        outlier_summary{file_idx, 11} = NaN;
        outlier_summary{file_idx, 12} = NaN;
    end
end




%% FINAL SUMMARY

fprintf('\nBATCH PROCESSING COMPLETE -------------------------------\n\n');

outlier_table = cell2table(outlier_summary, ...
    'VariableNames', {'filename', 'session', 'mouse_id', 'region', 'total_recording_min', ...
                      'n_samples', 'baseline_uV', 'threshold', 'lower_thr_uV', ...
                      'upper_thr_uV', 'num_outliers', 'outlier_pct'});

if ~isempty(seizure_summary)
    seizure_table = cell2table(seizure_summary, ...
        'VariableNames', {'filename', 'session', 'mouse_id', 'region', 'total_recording_min', ...
                          'seizure_number', 'start_time_s', 'end_time_s', 'duration_s', ...
                          'median_energy', 'median_factor', 'thr_value', 'pct_above_thr', ...
                          'total_segments', 'rejected_segments'});
end


fprintf('OUTLIER REMOVAL SUMMARY\n');

disp(outlier_table);
fprintf('\n');

if ~isempty(seizure_summary)
    fprintf('SEIZURE DETECTION SUMMARY\n');
    disp(seizure_table);
    fprintf('\n');
end

outlier_csv = fullfile(batch_results_folder, 'outlier_summary.csv');
writetable(outlier_table, outlier_csv);

if ~isempty(seizure_summary)
    seizure_csv = fullfile(batch_results_folder, 'seizure_summary.csv');
    writetable(seizure_table, seizure_csv);
end

fprintf('  CSV files saved:\n');
fprintf('    - outlier_summary.csv\n');
if ~isempty(seizure_summary)
    fprintf('    - seizure_summary.csv\n');
end

%% le files 
excel_file = fullfile(batch_results_folder, 'batch_summary.xlsx');

writetable(outlier_table, excel_file, 'Sheet', 'outliers');
if ~isempty(seizure_summary)
    writetable(seizure_table, excel_file, 'Sheet', 'seizures');
end

fprintf('  Excel file saved:\n');
fprintf('    - batch_summary.xlsx\n');
fprintf('      Sheets: Outliers');
if ~isempty(seizure_summary)
    fprintf(', Seizures');
end
fprintf('\n\n');

% statistics
successful = sum(~strcmp(outlier_summary(:,2), 'ERROR'));
failed = num_files - successful;
total_seizures = size(seizure_summary, 1);

fprintf('STATISTICS\n');
fprintf('  Files processed:    %d\n', num_files);
fprintf('  Successful:         %d\n', successful);
fprintf('  Failed:             %d\n', failed);
fprintf('  Total seizures:     %d\n', total_seizures);

if total_seizures > 0
    fprintf('  Total seizure time: %.1f s\n', sum(seizure_table.duration_s));
    fprintf('  Mean duration:      %.1f s\n', mean(seizure_table.duration_s));
    fprintf('  Min duration:       %.1f s\n', min(seizure_table.duration_s));
    fprintf('  Max duration:       %.1f s\n', max(seizure_table.duration_s));
end

fprintf('\n  Detection method: Hilbert envelope + energy metric + duration \n');
fprintf('  Parameters:\n');
fprintf('    - Band-pass:      [%d-%d] Hz\n', seizure_config.bandpass_band(1), seizure_config.bandpass_band(2));
fprintf('    - Power:          %d\n', seizure_config.power_exponent);
fprintf('    - Window:         %.1f s\n', seizure_config.window_sec);
fprintf('    - Threshold:      median × %.1f\n', seizure_config.median_factor);
fprintf('    - Min duration:   %.0f s\n', seizure_config.min_seizure_duration);

fprintf('\n  Output files generated:\n');
fprintf('    Per-file:\n');
fprintf('      - %d *_clean.txt files\n', successful);

% Count files with seizures
files_with_seizures = 0;
if ~isempty(seizure_summary)
    unique_files = unique(seizure_summary(:,1));
    unique_files = unique_files(~strcmp(unique_files, 'ERROR'));
    files_with_seizures = numel(unique_files);
end

if files_with_seizures > 0
    fprintf('      - %d *_seizures.mat files (only files with seizures)\n', files_with_seizures);
    fprintf('      - %d *_seizures.fig files (only files with seizures)\n', files_with_seizures);
else
    fprintf('      - 0 *_seizures.mat files (no seizures detected)\n');
    fprintf('      - 0 *_seizures.fig files (no seizures detected)\n');
end

fprintf('\n    Summaries:\n');
fprintf('      - outlier_summary.csv\n');
if total_seizures > 0
    fprintf('      - seizure_summary.csv\n');
end
fprintf('      - batch_summary.xlsx (Outliers');
if total_seizures > 0
    fprintf(', Seizures');
end
fprintf(')\n');

fprintf('\n*** PIPELINE COMPLETE ***\n\n');


%% LOCAL FUNCTIONS

function [clean_signal, outlier_info] = remove_outliers_amplitude(signal, min_allowed, max_allowed)
    signal = signal(:);
    is_out = signal < min_allowed | signal > max_allowed | isnan(signal);
    out_idx = find(is_out);
    out_vals = signal(is_out);
    outlier_info = [out_idx, out_vals];
    clean_signal = signal;
    
    if isempty(out_idx)
        return;
    end
    
    good_idx = find(~is_out);
    
    if numel(good_idx) < 2
        warning('Too few valid points to interpolate.');
        return;
    end
    
    clean_signal(is_out) = interp1(good_idx, signal(good_idx), out_idx, 'linear', 'extrap');
end