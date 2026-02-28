%% EMG Artifact Removal Pipeline (Open-Source Version)
% Description: Demonstrates high-pass, band-pass, and dynamic notch 
% filtering to remove stimulation artifacts from a single EMG channel.

clear; clc; close all;

%% ========================================================================
%  1. CONFIGURATION & PARAMETERS
%  ========================================================================
% Target Data Path (Relative path for open-source users)
cfg.dataPath       = fullfile(pwd, 'Data');
cfg.file_off       = fullfile(cfg.dataPath, 'SingleCh_StimOFF.mat');
cfg.file_on        = fullfile(cfg.dataPath, 'SingleCh_StimON.mat');

% Hardware & Filtering Constants
cfg.hp_cutoff      = 0.5;           % High-pass cutoff (Hz)
cfg.bp_cutoffs     = [20, 500];     % Band-pass edges (Hz)
cfg.filt_order     = 3;             % General filter order 

% Artifact Removal Parameters (Notch)
cfg.notch_width    = 1.0;           % Hz to remove on either side of peak (+/-)
cfg.off_notch_freq = [60, 195, 300, 420]; % Known ambient line noise peaks
cfg.stim_pk_thresh = 1000;          % Minimum peak height threshold for stim detection
cfg.stim_min_delay = 20;            % Minimum peak distance/delay for stim detection

%% ========================================================================
%  2. DATA LOADING
%  ========================================================================
fprintf('Loading isolated EMG channel data...\n');
if ~exist(cfg.file_off, 'file') || ~exist(cfg.file_on, 'file')
    error('Data files not found. Ensure the Data folder is in the current directory.');
end

% Load data (Fs is loaded automatically from the .mat files)
data_off = load(cfg.file_off);
data_on  = load(cfg.file_on);

cfg.Fs = data_off.Fs; % Assuming same sampling rate for both
time_off = (0:length(data_off.raw_off)-1) / cfg.Fs;
time_on  = (0:length(data_on.raw_on)-1) / cfg.Fs;

%% ========================================================================
%  3. PREPROCESSING: HIGH-PASS & BAND-PASS FILTERING
%  ========================================================================
fprintf('Applying High-Pass and Band-Pass Filters...\n');

[B_hp, A_hp] = butter(cfg.filt_order, cfg.hp_cutoff / (cfg.Fs * 0.5), 'high');
[B_bp, A_bp] = butter(5, cfg.bp_cutoffs / (cfg.Fs * 0.5), 'bandpass'); 

% Apply to OFF condition
% temp_off = filtfilt(B_hp, A_hp, data_off.raw_off);
% filt_off = filtfilt(B_bp, A_bp, temp_off);
filt_off = filtfilt(B_hp, A_hp, data_off.raw_off);


% Apply to ON condition
filt_on = filtfilt(B_hp, A_hp, data_on.raw_on);
% temp_on = filtfilt(B_hp, A_hp, data_on.raw_on);
% filt_on = filtfilt(B_bp, A_bp, temp_on);

%% ========================================================================
%  4. ARTIFACT REJECTION & VISUALIZATION: STIM OFF
%  ========================================================================
fprintf('Processing STIM OFF Baseline...\n');
[f_vec_off, P1_off_raw] = compute_fft(filt_off, cfg.Fs);

% Apply Static Notch Filters for Ambient Noise
sig_off_clean = filt_off;
for i = 1:length(cfg.off_notch_freq)
    f0 = cfg.off_notch_freq(i);
    Wn = [(f0 - cfg.notch_width)/(cfg.Fs*0.5), (f0 + cfg.notch_width)/(cfg.Fs*0.5)];
    [Bnot, Anot] = butter(cfg.filt_order, Wn, 'stop');
    sig_off_clean = filtfilt(Bnot, Anot, sig_off_clean);
end
[~, P1_off_clean] = compute_fft(sig_off_clean, cfg.Fs);

% Plot STIM OFF Diagnostics
fig1 = figure('Name', 'Diagnostics: STIM OFF', 'Position', [100, 100, 1200, 400], 'Color', 'w');
sgtitle(sprintf('STIM OFF Processing | Notch Targets: %s Hz | Notch Width: ±%.1f Hz', ...
        num2str(cfg.off_notch_freq), cfg.notch_width), 'FontWeight', 'bold');

subplot(1,2,1);
plot(time_off, filt_off/1000, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Bandpass Only'); hold on;
plot(time_off, sig_off_clean/1000, 'b', 'DisplayName', 'Notch Filtered');
title('Time Domain'); xlabel('Time (s)'); ylabel('Amplitude (mV)');
legend;

% xlim([10, 30]);
% ylim([-1.2, 1.2]);


subplot(1,2,2);
plot(f_vec_off, P1_off_raw, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Before Notch'); hold on;
plot(f_vec_off, P1_off_clean, 'b', 'DisplayName', 'After Notch');
title('Frequency Domain'); xlabel('Frequency (Hz)'); ylabel('|P1(f)|');
make_pretty;
legend;

%% ========================================================================
%  5. ARTIFACT REJECTION & VISUALIZATION: STIM ON
%  ========================================================================
fprintf('Processing STIM ON Artifact Removal...\n');
[f_vec_on, P1_on_raw] = compute_fft(filt_on, cfg.Fs);

% Dynamic Peak Detection (Requires find_notch_peaks_for_stim_artifact.m in path)
notch_targets = find_notch_peaks_for_stim_artifact(P1_on_raw, f_vec_on, cfg.Fs, cfg.stim_pk_thresh, cfg.stim_min_delay);

% Apply Dynamic Notch Filters
sig_on_clean = filt_on;
for i = 1:length(notch_targets)
    f0 = notch_targets(i);
    Wn = [(f0 - cfg.notch_width)/(cfg.Fs*0.5), (f0 + cfg.notch_width)/(cfg.Fs*0.5)];
    [Bnot, Anot] = butter(cfg.filt_order, Wn, 'stop');
    sig_on_clean = filtfilt(Bnot, Anot, sig_on_clean);
end

sig_on_clean = filtfilt(B_bp, A_bp, sig_on_clean);
[~, P1_on_clean] = compute_fft(sig_on_clean, cfg.Fs);

% Plot STIM ON Diagnostics
fig2 = figure('Name', 'Diagnostics: STIM ON', 'Position', [100, 600, 1200, 400], 'Color', 'w');
sgtitle(sprintf('STIM ON Artifact Rejection | Thresh: %d | Min Delay: %d | Notch Width: ±%.1f Hz', ...
        cfg.stim_pk_thresh, cfg.stim_min_delay, cfg.notch_width), 'FontWeight', 'bold');

subplot(1,2,1);
plot(time_on, filt_on/1000, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Raw Stim Artifact'); hold on;
plot(time_on, sig_on_clean/1000, 'r', 'DisplayName', 'Artifact Cleaned');
title('Time Domain'); xlabel('Time (s)'); ylabel('Amplitude (mV)');
legend;
% xlim([61, 81]);
% ylim([-1.2, 1.2]);

subplot(1,2,2);
plot(f_vec_on, P1_on_raw, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Before Removal'); hold on;
plot(f_vec_on, P1_on_clean, 'r', 'DisplayName', 'After Removal');
title('Frequency Domain'); xlabel('Frequency (Hz)'); ylabel('|P1(f)|');
% xlim([0, 500]);
% ylim([0, 50]);
legend;

fprintf('Pipeline execution complete.\n');
make_pretty;

%% ========================================================================
%  HELPER FUNCTIONS
%  ========================================================================
function [f_vec, P1] = compute_fft(signal, Fs)
    % Computes the single-sided amplitude spectrum of a 1D signal
    L = length(signal);
    Y = fft(signal);
    P2 = abs(Y / L);
    P1 = P2(1:floor(L/2)+1);
    if length(P1) > 2
        P1(2:end-1) = 2 * P1(2:end-1);
    end
    f_vec = Fs * (0:floor(L/2)) / L;
end