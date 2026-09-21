%% REGENERATE SEIZURE FIGURES FROM MODIFIED .mat FILE
% Quick script to recreate seizure visualizations after manual edits

clear; clc;

%% SELECT FILES
[txt_file, txt_path] = uigetfile('*.txt', 'Select cleaned LFP .txt file');
if txt_file == 0, error('No .txt file selected'); end

[mat_file, mat_path] = uigetfile('*.mat', 'Select modified seizures .mat file');
if mat_file == 0, error('No .mat file selected'); end

%% LOAD DATA
fprintf('Loading data...\n');

% Load LFP
data = load_LFP_intan_txt(fullfile(txt_path, txt_file));
fs = data.fs;
t = data.time_seconds;
clean_signal_uV = data.signal;

% Load seizures
load(fullfile(mat_path, mat_file), 'seizures', 'seizure_config');

% Extract base filename
[~, fname, ~] = fileparts(txt_file);

n_seizures = height(seizures);
fprintf('Found %d seizures in .mat file\n', n_seizures);

if n_seizures == 0
    fprintf('No seizures to plot.\n');
    return;
end

%% PREPROCESS (same as pipeline)
fprintf('Preprocessing signal...\n');

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
power = envelope .^ seizure_config.power_exponent;
energy = movmean(power, window_samples);

threshold = median(energy) * seizure_config.median_factor;
above_threshold = energy > threshold;

%% OVERVIEW FIGURE
fprintf('Creating overview figure...\n');

fig_overview = figure('Position', [50, 50, 1400, 900]);

% Panel 1: Cleaned LFP
subplot(3,1,1)
plot(t, clean_signal_uV, 'Color', [0.3 0.3 0.3], 'LineWidth', 0.5); hold on
for k = 1:n_seizures
    xregion(seizures.start_time_s(k), seizures.end_time_s(k), ...
        'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
        'EdgeColor', 'r', 'LineWidth', 2);
end
title(sprintf('%s - Cleaned LFP (%d seizures)', fname, n_seizures), 'Interpreter', 'none');
xlabel('Time (s)'); ylabel('µV'); grid on; hold off;

% Panel 2: Band-passed
subplot(3,1,2)
plot(t_bp, lfp_bp, 'Color', [0.2 0.4 0.7], 'LineWidth', 0.5); hold on
for k = 1:n_seizures
    xregion(seizures.start_time_s(k), seizures.end_time_s(k), ...
        'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
        'EdgeColor', 'r', 'LineWidth', 2);
end
title('Band-passed [5-75 Hz]');
xlabel('Time (s)'); ylabel('Normalized'); grid on; hold off;

% Panel 3: Energy metric
subplot(3,1,3)
plot(t_bp, energy, 'Color', [0.2 0.6 0.3], 'LineWidth', 0.5); hold on
yline(threshold, 'r--', sprintf('Threshold (median × %.1f)', seizure_config.median_factor), 'LineWidth', 2);
plot(t_bp(above_threshold), energy(above_threshold), 'r.', 'MarkerSize', 3);

for k = 1:n_seizures
    xregion(seizures.start_time_s(k), seizures.end_time_s(k), ...
        'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, ...
        'EdgeColor', 'r', 'LineWidth', 2);
end
title('Energy Metric');
xlabel('Time (s)'); ylabel('Energy (a.u.)'); grid on; hold off;

% Save overview
fig_file = fullfile(txt_path, sprintf('%s_seizures_REGENERATED.fig', fname));
savefig(fig_overview, fig_file);
fprintf('Saved: %s\n', fig_file);

%% ZOOM FIGURES FOR EACH SEIZURE
fprintf('Creating zoom figures...\n');

zoom_margin = 5; % seconds

for k = 1:n_seizures
    zoom_start = max(seizures.start_time_s(k) - zoom_margin, t_bp(1));
    zoom_end = min(seizures.end_time_s(k) + zoom_margin, t_bp(end));
    
    idx_zoom = (t_bp >= zoom_start) & (t_bp <= zoom_end);
    idx_zoom_full = (t >= zoom_start) & (t <= zoom_end);
    
    fig_zoom = figure('Position', [50, 50, 1400, 900]);
    
    % Panel 1: Full LFP
    subplot(3,1,1)
    plot(t(idx_zoom_full), clean_signal_uV(idx_zoom_full), 'k', 'LineWidth', 0.8); hold on;
    xregion(seizures.start_time_s(k), seizures.end_time_s(k), ...
        'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 2);
    xline(seizures.start_time_s(k), 'r--', 'Start', 'LineWidth', 2);
    xline(seizures.end_time_s(k), 'r--', 'End', 'LineWidth', 2);
    title(sprintf('Seizure #%d (full LFP)', k), 'Interpreter', 'none');
    xlabel('Time (s)'); ylabel('Voltage (µV)'); grid on; xlim([zoom_start, zoom_end]); hold off;
    
    % Panel 2: Band-passed
    subplot(3,1,2)
    plot(t_bp(idx_zoom), lfp_bp(idx_zoom), 'k', 'LineWidth', 0.8); hold on;
    xregion(seizures.start_time_s(k), seizures.end_time_s(k), ...
        'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 2);
    xline(seizures.start_time_s(k), 'r--', 'LineWidth', 2);
    xline(seizures.end_time_s(k), 'r--', 'LineWidth', 2);
    title(sprintf('Seizure %d (band-passed [5-75 Hz])', k), 'Interpreter', 'none');
    xlabel('Time (s)'); ylabel('Normalized Amplitude'); grid on; xlim([zoom_start, zoom_end]); hold off;
    
    % Panel 3: Energy metric
    subplot(3,1,3)
    plot(t_bp(idx_zoom), energy(idx_zoom), 'Color', [0.2 0.6 0.3], 'LineWidth', 1.0); hold on;
    yline(threshold, 'r--', sprintf('Threshold (%.2e)', threshold), 'LineWidth', 2);
    
    above_idx = idx_zoom & above_threshold';
    if any(above_idx)
        plot(t_bp(above_idx), energy(above_idx), 'r.', 'MarkerSize', 6);
    end
    
    xregion(seizures.start_time_s(k), seizures.end_time_s(k), ...
        'FaceColor', [1 0.7 0.7], 'FaceAlpha', 0.4, 'EdgeColor', 'r', 'LineWidth', 2);
    xline(seizures.start_time_s(k), 'r--', 'LineWidth', 2);
    xline(seizures.end_time_s(k), 'r--', 'LineWidth', 2);
    
    title(sprintf('Hilbert envelope(^%d) over a %ds smoothed window', ...
        seizure_config.power_exponent, seizure_config.window_sec));
    xlabel('Time (s)'); ylabel('Energy (a.u.)'); grid on; xlim([zoom_start, zoom_end]); hold off;
    
    % Save zoom figure
    zoom_fig_file = fullfile(txt_path, sprintf('%s_seizure%d_REGENERATED.fig', fname, k));
    savefig(fig_zoom, zoom_fig_file);
    close(fig_zoom);
    
    fprintf('  Seizure #%d saved: %s\n', k, zoom_fig_file);
end

fprintf('\nDone! All figures regenerated.\n');