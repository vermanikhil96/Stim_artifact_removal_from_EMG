function notch_freq = find_notch_peaks_for_stim_artifact(FFT,Freq,Fs,freq_cutoff,lower_freq_cutoff);


figure
plot(Freq,FFT) 
title("FFT (Select point for peak Threshold)")
xlabel("Freq")
xlim([0.1,freq_cutoff])
make_pretty;
pause;
[~,peak_thresh] = ginput(1);
close

figure
plot(FFT) 
title(" FFT (Select minimum peak distance)")
xlabel("samples")
ylim([0,peak_thresh*3])
make_pretty;
pause;
[x,y]= ginput(2);
min_peak_distance = abs(x(2)-x(1));
%% plot FFT with Notch centers

[~,locs] = findpeaks(FFT, "MinPeakHeight",peak_thresh,"MinPeakDistance", min_peak_distance);
figure
plot(Freq,FFT) 
hold on
scatter(Freq(locs(Freq(locs)<freq_cutoff)),FFT(locs(Freq(locs)<freq_cutoff)))
title(" FFT (Select minimum peak distance)")
xlabel("samples")
xlim([0.1,Fs/2])
make_pretty;

notch_freq = Freq(locs(Freq(locs)<freq_cutoff));
notch_freq = notch_freq(notch_freq>lower_freq_cutoff);

% figure
% plot(Freq,FFT) 
% hold on
% scatter(notch_freq,FFT(notch_freq))
% title(" FFT (Select minimum peak distance)")
% xlabel("samples")
% xlim([0.1,Fs/2])
% make_pretty;

end

