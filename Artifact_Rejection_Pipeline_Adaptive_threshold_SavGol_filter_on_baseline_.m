%% Adaptive Frequency-Dependent Artifact Rejection
% Description: Applies initial HP/BP pre-processing, then uses Savitzky-Golay 
% smoothing and local Median Absolute Deviation (MAD) on a Stim-OFF baseline 
% to dynamically detect and notch stimulation artifacts in a Stim-ON signal.

clear; clc; close all;

%% ========================================================================
%  1. USER PARAMETERS & CONFIGURATION
%  ========================================================================
% File Paths
cfg.dataPath     = fullfile(pwd, 'Data');
cfg.file_off     = fullfile(cfg.dataPath, 'SingleCh_StimOFF.mat');
cfg.file_on      = fullfile(cfg.dataPath, 'SingleCh_StimON.mat');

% Plotting Parameters
cfg.f_max_plot   = 600;       % Maximum frequency to display in plots (Hz)

% Pre-Processing Filters (HP & BP)
cfg.hp_cutoff    = 0.5;       % High-pass cutoff (Hz)
cfg.bp_cutoffs   = [20, 500]; % Band-pass edges (Hz)
cfg.prefilt_ord  = 3;         % General filter order for pre-processing

% Adaptive Detection Parameters
cfg.smooth_win   = 10;        % Window size for local MAD calculation (Hz)
cfg.k_mad        = 8;         % Multiplier for MAD thresholding (Stringency)
cfg.width_thr    = 4;         % Max peak width to be considered a stim artifact (Hz)
cfg.min_dist     = 5;         % Minimum distance between detected artifacts (Hz)

% Notch Filter Parameters
cfg.notch_width  = 1.0;       % Half-width of the notch filter (Hz)
cfg.notch_order  = 3;         % Butterworth notch filter order

%% ========================================================================
%  2. DATA INGESTION
%  ========================================================================
fprintf('Loading isolated EMG channel data...\n');
if ~exist(cfg.file_off, 'file') || ~exist(cfg.file_on, 'file')
    error('Data files missing. Ensure OpenSource_Export is in the current directory.');
end

d_off  = load(cfg.file_off);
d_on   = load(cfg.file_on);
Fs     = double(d_off.Fs);

raw_off = double(d_off.raw_off(:)');
raw_on  = double(d_on.raw_on(:)');

%% ========================================================================
%  3. PRE-PROCESSING: HIGH-PASS & BAND-PASS
%  ========================================================================
fprintf('Applying Initial High-Pass and Band-Pass Filters...\n');
[B_hp, A_hp] = butter(cfg.prefilt_ord, cfg.hp_cutoff / (Fs * 0.5), 'high');
[B_bp, A_bp] = butter(5, cfg.bp_cutoffs / (Fs * 0.5), 'bandpass'); % Order 5 for tight BP roll-off

% Filter OFF condition
temp_off = filtfilt(B_hp, A_hp, raw_off);
filt_off = filtfilt(B_bp, A_bp, temp_off);

% Filter ON condition
temp_on  = filtfilt(B_hp, A_hp, raw_on);
filt_on  = filtfilt(B_bp, A_bp, temp_on);

%% ========================================================================
%  4. BASELINE MODELING: STIM-OFF FFT & SMOOTHING
%  ========================================================================
fprintf('Computing baseline model and dynamic threshold...\n');
L_off = length(filt_off);

% STRICT INDEXING: Guarantees length match even if L_off is odd
f_off = Fs * (0:floor(L_off/2)) / L_off;

% Compute Single-Sided FFT on FILTERED data
Y_off = fft(filt_off);
P2_off = abs(Y_off/L_off);
P1_off = P2_off(1:floor(L_off/2)+1);
if length(P1_off) > 2
    P1_off(2:end-1) = 2 * P1_off(2:end-1);
end

% Mask for target frequency range
f_mask = f_off <= cfg.f_max_plot;
f_plot = f_off(f_mask);
P1_off_plot = P1_off(f_mask);

% Savitzky-Golay Smoothing (Log Domain for stable variance)
log_P1_off = log10(P1_off_plot + 1e-10);
sgolay_span = round(length(log_P1_off) / 8) * 2 + 1; % Force odd span
smooth_log = sgolayfilt(log_P1_off, 3, sgolay_span);
smooth_template = 10.^smooth_log;

% Compute Frequency-Dependent Threshold via Local MAD
mad_template = zeros(size(f_plot));
bin_hz = mean(diff(f_plot));
win_bins = round(cfg.smooth_win / bin_hz);

for i = 1:length(f_plot)
    idx1 = max(1, i - floor(win_bins/2));
    idx2 = min(length(f_plot), i + floor(win_bins/2));
    local_vals = P1_off_plot(idx1:idx2);
    mad_template(i) = median(abs(local_vals - median(local_vals)));
end
freq_dep_threshold = smooth_template + (cfg.k_mad * mad_template);

%% ========================================================================
%  5. ARTIFACT DETECTION: STIM-ON SIGNAL
%  ========================================================================
fprintf('Detecting artifacts in Stim-ON data...\n');
L_on = length(filt_on);
time_on = (0:L_on-1) / Fs;

% STRICT INDEXING: Guarantees length match even if L_on is odd
f_on = Fs * (0:floor(L_on/2)) / L_on;

% Compute Single-Sided FFT for FILTERED ON Signal
Y_on = fft(filt_on);
P2_on = abs(Y_on/L_on);
P1_on = P2_on(1:floor(L_on/2)+1);
if length(P1_on) > 2
    P1_on(2:end-1) = 2 * P1_on(2:end-1);
end

% Interpolate to match plotting/threshold grid safely
P1_on_plot = interp1(f_on, P1_on, f_plot, 'linear', 'extrap');

% Calculate difference spectrum (findpeaks requires scalar MinPeakHeight)
diff_spec = P1_on_plot - freq_dep_threshold; 

% Peak Detection via Adaptive Threshold
[~, pk_locs, pk_widths, ~] = findpeaks(diff_spec, f_plot, ...
    'MinPeakHeight', 0, ...
    'MinPeakDistance', cfg.min_dist, ...
    'WidthReference', 'halfprom');

% Isolate sharp peaks (artifacts) vs wide peaks (physiological)
sharp_mask = (pk_widths < cfg.width_thr) & (pk_locs > 5) & (pk_locs < (Fs/2 - 5));
artifact_freqs = pk_locs(sharp_mask);

fprintf('Detected %d narrow-band artifacts.\n', length(artifact_freqs));

%% ========================================================================
%  6. ARTIFACT REMOVAL: SEQUENTIAL NOTCH FILTERING
%  ========================================================================
fprintf('Applying targeted notch filters...\n');
filt_sig = filt_on; % Base notch filtering on the HP/BP filtered signal

for f_idx = 1:length(artifact_freqs)
    f0 = artifact_freqs(f_idx);
    Wn = [(f0 - cfg.notch_width)/(Fs/2), (f0 + cfg.notch_width)/(Fs/2)];
    
    % Edge-case guards
    Wn(1) = max(Wn(1), 1e-4);
    Wn(2) = min(Wn(2), 1 - 1e-4);
    
    [Bnot, Anot] = butter(cfg.notch_order, Wn, 'stop');
    filt_sig = filtfilt(Bnot, Anot, filt_sig);
end

% Compute Cleaned FFT
Y_filt = fft(filt_sig);
P1_filt = abs(Y_filt/L_on);
P1_filt = P1_filt(1:floor(L_on/2)+1);
P1_filt(2:end-1) = 2 * P1_filt(2:end-1);
P1_filt_plot = interp1(f_on, P1_filt, f_plot, 'linear', 'extrap');

%% ========================================================================
%  7. EXPLANATORY VISUALIZATIONS
%  ========================================================================
disp('Generating explanatory figures...');

% FIGURE 1: Baseline Modeling & Threshold Generation
figure('Name', 'Step 1: Baseline Thresholding', 'Position', [50, 500, 800, 350], 'Color', 'w');
plot(f_plot, P1_off_plot, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Filtered Stim-OFF'); hold on;
plot(f_plot, smooth_template, 'k', 'LineWidth', 1.5, 'DisplayName', 'Savitzky-Golay Baseline');
plot(f_plot, freq_dep_threshold, 'g--', 'LineWidth', 1.5, 'DisplayName', sprintf('Threshold (k=%d)', cfg.k_mad));
xlim([5 cfg.f_max_plot]);
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Step 1: Establishing the Adaptive Threshold from Baseline Noise');
legend('Location', 'northeast'); grid on;
make_pretty;
% FIGURE 2: Peak Detection
figure('Name', 'Step 2: Artifact Detection', 'Position', [50, 80, 800, 350], 'Color', 'w');
plot(f_plot, P1_on_plot, 'b', 'LineWidth', 1.0, 'DisplayName', 'Filtered Stim-ON'); hold on;
plot(f_plot, freq_dep_threshold, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Dynamic Threshold');
scatter(artifact_freqs, interp1(f_plot, P1_on_plot, artifact_freqs), 50, 'r', 'filled', 'v', 'DisplayName', 'Detected Artifacts');
xlim([5 cfg.f_max_plot]);
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Step 2: Dynamic Detection of Non-Physiological Peaks');
legend('Location', 'northeast');
make_pretty;
% FIGURE 3: Final Results (Time & Frequency)
figure('Name', 'Step 3: Before & After', 'Position', [900, 80, 800, 770], 'Color', 'w');

% Frequency Subplot
subplot(2,1,1);
plot(f_plot, P1_on_plot, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.0, 'DisplayName', 'Before Removal'); hold on;
plot(f_plot, P1_filt_plot, 'b', 'LineWidth', 1.2, 'DisplayName', 'After Removal');
xlim([5 cfg.f_max_plot]);
xlabel('Frequency (Hz)'); ylabel('Amplitude');
title('Frequency Domain: Artifact Suppression');
legend('Location', 'northeast');

% Time Subplot (Zoomed to first 2 seconds)
subplot(2,1,2);
plot(time_on, filt_on/1000, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Before Removal'); hold on;
plot(time_on, filt_sig/1000, 'b', 'DisplayName', 'After Removal');
xlabel('Time (s)'); ylabel('Amplitude (mV)');
title('Time Domain: Signal Integrity');
legend('Location', 'northeast');
fprintf('Execution complete. Ready for review.\n');
make_pretty;