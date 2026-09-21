% view_LFP_bandpassed.m

clear;
[filename, filepath] = uigetfile('*.txt', 'Select LFP .txt file');

if isequal(filename, 0)
    disp('cancelled file selection');
    return;
end

full_path = fullfile(filepath, filename);

data = load_LFP_intan_txt(full_path);

signal_uV = data.signal;
t = data.time_seconds;
fs = data.fs;

signal_bp = bandpass(signal_uV, [5 75], fs, 'Steepness', 0.5);


figure('Name', sprintf('%s - Bandpassed [5-75 Hz]', filename), ...
       'NumberTitle', 'off', 'Position', [100, 100, 1400, 600]);

plot(t, signal_bp, 'k', 'LineWidth', 0.5);

title(sprintf('%s - Bandpassed [5-75 Hz]', filename), 'Interpreter', 'none', 'FontSize', 12);
xlabel('Time (s)', 'FontSize', 11);
ylabel('Amplitude (µV)', 'FontSize', 11);
grid on;

fprintf('Signal loaded: %s\n', filename);
fprintf('Duration: %.2f minutes\n', max(t)/60);
fprintf('Bandpass filter: [5-75] Hz\n');