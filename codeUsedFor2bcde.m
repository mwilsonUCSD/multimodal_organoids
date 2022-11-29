% Author: Madison Wilson

clear;
clc;

% specify recording file
nickname = 'NOD20\0428_run03';
%  Manual load processed\video_runX.mat
%load 'C:\Users\Maddie\Documents\Martin Exps\Processed Data\NOD20 video\0309_run4.mat';
%  Manual load processed_EP\lpf_YYY_runX.mat
load(['D:\Maddie\Martin Exps\Processed Data\' nickname '.mat']);

%% check sampling frequency, impedance
fs = rhd.frequency_parameters.amplifier_sample_rate;
impedance = [rhd.amplifier_channels.electrode_impedance_magnitude];
t_amplifier = rhd.t_amplifier;

% mapping of array channels (Gold facing down onto cortex)
% [4, 3,  1,  16]
% [5, 6,  2,  15]
% [7, 10, 14, 13]
% [8, 9,  11, 12]
mapping = [4,3,1,16,5,6,2,15,7,10,14,13,8,9,11,12];

impedance_mapped = impedance(9:24); %cut out unused channels
good_channels = find(impedance_mapped < 35e6); % consider channels below an impedance threshold "good channels"
impedance_mapped = impedance_mapped(good_channels);

data_good = rhd.amplifier_data(9:24,:); % cut out unused channels
data_good = data_good(good_channels, :); % reorder channels to match array mapping

% find ADC channel marks (times of stimuli onset)
ind = find(diff(rhd.board_adc_data(1,:) > 0.3) == 1) + 1;
t_stim = t_amplifier(ind);

%% plot raw data
figure;
for i = 1:size(data_good,1)
    plot(t_amplifier, data_good(i,1:end)+500*i); hold on;
    %plot(t_amplifier(t_range(1:end)), data_good(i,1:end)+100*i); hold on;
end
title([nickname ' Raw Data'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');

%% Filter data below fs_down Hz (easier on ICA)
fs_down = 6000; %(6-20 kHz)

d_low_3000 = designfilt('lowpassiir','FilterOrder',8, ...
'PassbandFrequency',3000,'PassbandRipple',0.2, ...
'SampleRate',fs);

data_clean_filt = filtfilt(d_low_3000, data_good');

data_down = data_clean_filt(1:fs/fs_down:end,:); % fs_down Hz sampling rate
t_down = t_amplifier(1:fs/fs_down:end);

%% ICA to remove noise
[Aica,ica_source_activity] = ica_denoise(data_down, t_down, 1:length(good_channels));


%% Remove noisy ICA components
good_ind = setdiff(1:length(good_channels), [5]); % second input is channel to remove based on imaging artifact
%good_ind = setdiff(1:length(good_channels), [1]); % second input is channel to remove based on imaging artifact % MNW
ECoG_denoise = ica_source_activity(:,good_ind) * Aica(:,good_ind)';

figure;
for i = 1:size(data_good,1)
    plot(t_down, ECoG_denoise(:,i)+500*i); hold on;
    plot(t_down, data_down(:,i)+500*i); hold on;
    %plot(t_amplifier(t_range(1:end)), data_good(i,1:end)+100*i); hold on;
end
title([nickname ' data post-ICA'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');

%% Filter data into LFP
d_low_250 = designfilt('lowpassiir','FilterOrder',8, ...
'PassbandFrequency',250,'PassbandRipple',0.2, ...
'SampleRate',fs);

% LFP is filtered < 250 Hz
data_lfp = filtfilt(d_low_250, ECoG_denoise);


%% Plot LFP (figure 3b)
figure;
for i = 1:size(data_lfp,2)
%for i = 2:8
    plot(t_down, data_lfp(:,i)+500*i); hold on;
    %plot(0:1/fs_down:10, data_lfp(i,760540:760540+fs_down*10)+500*i, 'LineWidth', 2); hold on;
end
ylim([0 (i+1)*500]);
for i = 1:length(t_stim)
    plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
end
title('LFP');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');


%% calculate and plot LFP trial averages for each channel (figure 3c)
% User parameters, note these replace the ones in the section above
leadTime = -1; % -0.25
lenTrial = 6; %0.5;
badEpochs = []; % Manually remove epochs with noise 

[trialAvgLFP, t_trial_lfp] = f_plotSDTrialAverages2D(0.001, [-200 200], [0,0.1], data_lfp', t_down, fs_down, t_stim, leadTime, lenTrial, badEpochs, mapping, good_channels, 'Amplitude (\muV)');
sgtitle([nickname ' LFP'], 'Interpreter', 'none');


%% detect and plot peak amplitudes and time delays (figre 2d)
% Modify the color axis limits inside the function peakColorMap()
upperBoundNegPeakEst = 0.05;
upperBoundUpperPeakEst = 0.12;
%peakColorMap(trialAvgLFP, t_trial_lfp, fs_down, upperBoundNegPeakEst, upperBoundUpperPeakEst, leadTime, lenTrial, good_channels, mapping, nickname);%MNW
peakColorMap_220927(trialAvgLFP, t_trial_lfp, fs_down, upperBoundNegPeakEst, upperBoundUpperPeakEst, leadTime, lenTrial, good_channels, mapping, nickname);


%% Average of trial (Morlet) spectrograms for single channel (figure 2e and Supplementary figure 4)
% User parameters
channel = 12; % channel of interest
freq_range = [32,150]; % frequency range of interest%freq_range = [32,150]; % frequency range of interest%MNW
deltaF = 0.5; % frequency resolution [Hz]
morletParams = [20,80]; % Range of cycles (low number of cycles => better temporal resolution, high # cycles => better freq resolution)

% Change index parameter to match downsampled data
i = find(good_channels == channel);
ind_down = zeros(size(ind));
for j = 1:length(ind)
    ind_down(j) = find(t_down >= t_amplifier(ind(j)),1);
end

% Compute frequency range and make holder variable
num_frex = (freq_range(2)-freq_range(1)+1)/deltaF;
rangeF = linspace(freq_range(1), freq_range(2), num_frex);
holdingM = zeros(length(rangeF), (lenTrial-leadTime)*fs_down);

numTrials = length(t_stim);
lastUsedJ = 1;

for j = 1:numTrials
    % skip bad epochs (due to mouse motion or other noise contamination)
    if find(badEpochs == j)
        continue;
    end
    
    %abridge the epochs with no noise
    A= ECoG_denoise(ind_down(lastUsedJ)+leadTime*fs_down:ind_down(lastUsedJ)-1, i);
    B = ECoG_denoise(ind_down(j):ind_down(j)+(lenTrial*fs_down)-1, i);
    tempTrialData = [A;B];
    
    % Morlet-based spectrogram
    M = my_morlet(tempTrialData, fs_down, freq_range(1), freq_range(2), deltaF, morletParams);
    holdingM = holdingM + M;
    disp(lastUsedJ);
    lastUsedJ = j+1;
end

%average
div = numTrials - length(badEpochs);
holdingM = holdingM/div;

%subtract baseline
%baseline = median(holdingM(:,1:(leadTime*fs_down*-1)-1), 2);%MNW
baseline = median(holdingM(:,1:(leadTime*fs_down*-1)-0.2*fs_down), 2);
holdingM = 20*(log10(holdingM) - log10(baseline));
t_trial = leadTime:1/fs_down:lenTrial-1/fs_down;

figure;
%plot_matrix(holdingM(:,1:10:end)',t_trial(1:10:end),rangeF,'l');
%imagesc(t_trial(1:10:end),rangeF,holdingM(:,1:10:end)');axis xy; colorbar; title('Spectrogram');

h = pcolor(t_trial(1:10:end),rangeF, holdingM(:,1:10:end));
set(h,'EdgeColor','none');
%set(gca, 'YScale', 'log');
%set(gca, 'FontSize', 25);
c = colorbar;
c.Label.String = '10*log10(signal/baseline)';
colormap jet;
caxis([1,8]);%caxis([0,4]);%MNW

xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title([nickname 'Channel: ' num2str(channel) ' mParams: ' num2str(morletParams) ' deltaF: ' num2str(deltaF)], 'Interpreter', 'none');%title([nickname 'Channel: ' num2str(channel)], 'Interpreter', 'none');


%% Test spectrogram with other STFT methods (Supplementary figure 12)

% To match Morlet ([4,40] , 0.11 deltaf) high freq range (32-250 Hz) results, use:
    % Overlap = 95%
    % Time res = 0.2 s
    % Leakage = 0.85 (Kaiser window with Beta ~= 6, approximates Hanning window)
% For low freq matching to Morlet results, use:
    % Overlap = 95%
    % Time res = 0.55 s
    % Leakage = 0.65
    
tres = 0.2;
op = 95;
leak = 0.85; %0.85 approximates hanning window
fRange = [1 150];

lastUsedJ = 1;

for j = 1:numTrials
    % skip bad epochs (due to mouse motion or other noise contamination)
    if find(badEpochs == j)
        continue;
    end
    
    %abridge the epochs with no noise
    %A= ECoG_denoise_down(ind_down(lastUsedJ)+leadTime*fs_down:ind_down(lastUsedJ)-1, i);
    %B = ECoG_denoise_down(ind_down(j):ind_down(j)+(lenTrial*fs_down)-1, i);
%     A= ECoG_denoise(ind_down(lastUsedJ)+leadTime*fs_down:ind_down(lastUsedJ)-1, i);
%     B = ECoG_denoise(ind_down(j):ind_down(j)+(lenTrial*fs_down)-1, i);
%     tempTrialData = [A;B];
    
    %abridge the epochs with no noise
    A= ECoG_denoise(ind_down(lastUsedJ)+leadTime*fs_down:ind_down(lastUsedJ)-1, i);
    B = ECoG_denoise(ind_down(j):ind_down(j)+(lenTrial*fs_down)-1, i);
    tempTrialData = [A;B];
    
    % STFT-based spectrogram
    [p,fp,tp] = pspectrum(tempTrialData, t_trial,'spectrogram', 'FrequencyLimits',fRange,'OverlapPercent', op,'TimeResolution',tres,'Leakage', leak); %'FrequencyResolution', 0.5,   
    if j == 1
        holdingP = zeros(length(fp),length(tp));
    end
    holdingP = holdingP + p;
    disp(lastUsedJ);
    lastUsedJ = j+1;
end

holdingP = holdingP/numTrials;

%tp = tp + leadTime;
closest2zero = find(tp < -0.2);
closest2zero = closest2zero(end);
bline = median(holdingP(:,1:closest2zero), 2);

figure; %[p,fp,tp] = pspectrum(EEG.data,EEG.srate,'spectrogram','FrequencyLimits',[1 150],'TimeResolution',0.8);
h2 = pcolor(tp, fp, 10*log10(holdingP) - 10*log10(bline));
set(h2,'EdgeColor','none');
%imagesc(tp,fp,10*log10(p) - 10*log10(bline));
set(gca, 'FontSize', 15);
c = colorbar;
c.Label.String = '10*log10(signal/baseline)';
colormap jet;
caxis([0.5,8]);
set(gca, 'YDir', 'normal');%set(gca, 'YDir', 'normal', 'YScale', 'log');

xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title([nickname ' STFT Ch: ' num2str(channel) ' tres: ' num2str(tres) ' op: ' num2str(op) ' leak: ' num2str(leak)], 'Interpreter', 'none');


%% 2D plot of Morlet method spectrograms (figure 2e)
% User parameters
freq_range = [32,150]; % frequency range of interest
deltaF = 0.5; % frequency resolution [Hz]
morletParams = [4,80]; % Range of cycles (low number of cycles => better temporal resolution, high # cycles => better freq resolution)

% Calculate frequency range
num_frex = (freq_range(2)-freq_range(1)+1)/deltaF;
rangeF = linspace(freq_range(1), freq_range(2), num_frex);
numTrials = length(t_stim);
    
figure;
for i = 1:length(good_channels)
    holdingM = zeros(length(rangeF), (lenTrial-leadTime)*fs_down);

    lastUsedJ = 1;

    for j = 1:numTrials
        % skip bad epochs (due to mouse motion or other noise contamination)
        if find(badEpochs == j)
            continue;
        end

        %abridge the epochs with no noise
        A= ECoG_denoise(ind_down(lastUsedJ)+leadTime*fs_down:ind_down(lastUsedJ)-1, i);
        B = ECoG_denoise(ind_down(j):ind_down(j)+(lenTrial*fs_down)-1, i);
        tempTrialData = [A;B];

        % Morlet-based spectrogram
        M = my_morlet(tempTrialData, fs_down, freq_range(1), freq_range(2), deltaF, morletParams);
        holdingM = holdingM + M;
        lastUsedJ = j+1;
        %disp(lastUsedJ);
    end

    %average
    div = numTrials - length(badEpochs);
    holdingM = holdingM/div;

    %subtract baseline
    %baseline = median(holdingM(:,1:(leadTime*fs_down*-1)-1), 2);
    baseline = median(holdingM(:,1:(leadTime*fs_down*-1)-0.2), 2);
    holdingM = 20*(log10(holdingM) - log10(baseline));
    t_trial = leadTime:1/fs_down:lenTrial-1/fs_down;

    subplot(4,4,find(mapping == good_channels(i)));
    h = pcolor(t_trial(1:10:end),rangeF, holdingM(:,1:10:end));
    set(h,'EdgeColor','none');
    %set(gca, 'YScale', 'log');
    %set(gca, 'FontSize', 25);
    c = colorbar;
    c.Label.String = '10*log10(signal/baseline)';
    colormap jet;
    caxis([1,8]);%caxis([0,4]); %MNW
    set(c, 'visible', 'off');
    title(['Channel: ' num2str(good_channels(i))]);
    xlabel('Time (seconds)');
    ylabel('Frequency (Hz)');
    
    disp(i);
end
sgtitle(nickname, 'Interpreter', 'none');





