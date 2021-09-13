% Author: Madison Wilson

clear;
clc;

% specify recording file
nickname = 'NOD13\0302_run3';
%  Manual load processed\video_runX.mat
%load 'C:\Users\Maddie\Documents\Martin Exps\Processed Data\NOD20 video\0309_run4.mat';
%  Manual load processed_EP\lpf_YYY_runX.mat
load(['C:\Users\Maddie\Documents\Martin Exps\Processed Data\' nickname '.mat']);

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
for i = 1:length(t_stim)
    plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
end
title([nickname ' Raw Data'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');


%% Evaluate ICA
%[Aica,ica_source_activity] = ica_denoise(data_good', t_amplifier, 1:length(good_channels));
f_plot_ica(real(Aica),real(ica_source_activity), good_channels, mapping, t_amplifier);


%% Remove noisy ICA components
ica_remove = [1,2];
good_ind = setdiff(1:length(good_channels), ica_remove); % second input is channel to remove based on imaging artifact
ECoG_denoise = ica_source_activity(:,good_ind) * Aica(:,good_ind)';

%Clear helper variables
clear good_ind;

%% Plot new and old data
figure;
for i = 1:size(data_good,1)
    %plot(t_amplifier, data_good(i,1:end)+500*i); hold on;
    plot(t_down, ECoG_denoise(:,i)+500*i); hold on;
end
for i = 1:length(t_stim)
    plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
end
title([nickname ' data post-ICA'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');


%% Filter data into LFP and MUA
fs_down = 1000;

tic;
d_low_250 = designfilt('lowpassiir','FilterOrder',8, ...
'PassbandFrequency',250,'PassbandRipple',0.2, ...
'SampleRate',fs);
d_low = designfilt('bandpassiir','FilterOrder',6, ...
'HalfPowerFrequency1',1,'HalfPowerFrequency2',250, ...
'SampleRate',fs);
d_low_100 = designfilt('lowpassiir','FilterOrder',8, ...
'PassbandFrequency',100,'PassbandRipple',0.2, ...
'SampleRate',fs);
d_bpass = designfilt('bandpassiir','FilterOrder',6, ...
'HalfPowerFrequency1',500,'HalfPowerFrequency2',3000, ...
'SampleRate',fs);
toc; tic;

% LFP is just filtered < 250 Hz
data_lfp = filtfilt(d_low_250, ECoG_denoise);
data_lfp = data_lfp(1:fs/fs_down:end,:);
t_down = t_amplifier(1:fs/fs_down:end);

% MUA is filtered, full-wave rectified, and LPF for smoothing
data_band_filt = filtfilt(d_bpass, ECoG_denoise);
data_rect = data_band_filt.^2;
data_mua = filtfilt(d_low_100, data_rect);
%data_mua = data_band_filt;
toc;

%Clear helper variables
%clear {d_low_250, d_low_100, d_bpass};
%% Plot LFP and MUA
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

figure;
for i = 1:size(data_lfp,2)
    plot(t_amplifier, data_band_filt(:,i)+500*i); hold on;
end
ylim([0 (i+1)*500]);
for i = 1:length(t_stim)
    plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
end
title('MUA');

%data_mua = data_mua';
%data_lfp = data_lfp';
%ECoG_denoise = ECoG_denoise';

%% calculate and plot MUA trial averages for each channel (May take awhile to run)
leadTime = -1;
lenTrial = 8; % seconds
badEpochs = []; % Manually remove epochs with noise 

% f_plotSDTrialAverages2D() may take awhile to run. It's done when the figure appears.
[trialAverages, t_trial] = f_plotSDTrialAverages2D(0.05, [0 200], [0,4], data_mua', t_amplifier, fs, t_stim, leadTime, lenTrial, badEpochs, mapping, good_channels, 'Amplitude (\muV^2)');
sgtitle([nickname ' MUA ICA removed: ' num2str(ica_remove(:)')], 'Interpreter', 'none');


%% calculate and plot LFP trial averages for each channel (figure 2c)
% User parameters, note these replace the ones in the section above

[trialAvgLFP, t_trial_lfp] = f_plotSDTrialAverages2D(0.001, [-200 200], [0,4], data_lfp', t_down, fs_down, t_stim, leadTime, lenTrial, badEpochs, mapping, good_channels, 'Amplitude (\muV)');
sgtitle([nickname ' LFP ICA removed: ' num2str(ica_remove(:)')], 'Interpreter', 'none');


%% Save LFP and MUA
root = 'C:\Users\Maddie\Documents\Martin Exps\Processed Data\';
filename = [root, nickname];
save(filename, 'data_band_filt', 'data_lfp', 'ica_remove', '-append');

%% detect and plot peak amplitudes and time delays (figre 2d)
% Modify the color axis limits inside the function peakColorMap()
upperBoundNegPeakEst = 0.055;
upperBoundUpperPeakEst = 0.12;
peakColorMap(trialAvgLFP, t_trial_lfp, fs_down, upperBoundNegPeakEst, upperBoundUpperPeakEst, leadTime, lenTrial, good_channels, mapping, nickname);


%% Average of trial (Morlet) spectrograms for single channel
% User parameters
channel = 12; % channel of interest
freq_range = [0,10]; % frequency range of interest
deltaF = 0.11; % frequency resolution [Hz]
morletParams = [4,40]; % Range of cycles (low number of cycles => better temporal resolution, high # cycles => better freq resolution)

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
    A= ECoG_denoise_down(ind_down(lastUsedJ)+leadTime*fs_down:ind_down(lastUsedJ)-1, i);
    B = ECoG_denoise_down(ind_down(j):ind_down(j)+(lenTrial*fs_down)-1, i);
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
%baseline = median(holdingM(:,1:(leadTime*fs_down*-1)-1), 2);
baseline = median(holdingM(:,1:(leadTime*fs_down*-1)-0.2*fs_down), 2);
holdingM = 10*(log10(holdingM) - log10(baseline));
t_trial = leadTime:1/fs_down:lenTrial-1/fs_down;

figure;
%plot_matrix(holdingM(:,1:10:end)',t_trial(1:10:end),rangeF,'l');
%imagesc(t_trial(1:10:end),rangeF,holdingM(:,1:10:end)');axis xy; colorbar; title('Spectrogram');

h = pcolor(t_trial(1:10:end),rangeF, holdingM(:,1:10:end));
set(h,'EdgeColor','none');
%set(gca, 'YScale', 'log');
set(gca, 'FontSize', 25);
c = colorbar;
c.Label.String = '10*log10(signal/baseline)';
colormap jet;
caxis([0.5,4]);

xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title([nickname 'Channel: ' num2str(channel)], 'Interpreter', 'none');


%% 2D plot of Morlet method spectrograms (figure 2e)
% User parameters
freq_range = [0,10]; % frequency range of interest
deltaF = 0.1; % frequency resolution [Hz]
morletParams = [4,40]; % Range of cycles (low number of cycles => better temporal resolution, high # cycles => better freq resolution)

ECoG_denoise_down = ECoG_denoise(1:fs/fs_down:end,:);

ind_down = zeros(size(ind));
for j = 1:length(ind)
    ind_down(j) = find(t_down >= t_amplifier(ind(j)),1);
end

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
        A= ECoG_denoise_down(ind_down(lastUsedJ)+leadTime*fs_down:ind_down(lastUsedJ)-1, i);
        B = ECoG_denoise_down(ind_down(j):ind_down(j)+(lenTrial*fs_down)-1, i);
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
    holdingM = 10*(log10(holdingM) - log10(baseline));
    t_trial = leadTime:1/fs_down:lenTrial-1/fs_down;

    subplot(4,4,find(mapping == good_channels(i)));
    h = pcolor(t_trial(1:10:end),rangeF, holdingM(:,1:10:end));
    set(h,'EdgeColor','none');
    %set(gca, 'YScale', 'log');
    %set(gca, 'FontSize', 25);
    c = colorbar;
    c.Label.String = '10*log10(signal/baseline)';
    colormap jet;
    caxis([0.5,4]);
    %set(c, 'visible', 'on');
    title(['Channel: ' num2str(good_channels(i))]);
    xlabel('Time (seconds)');
    ylabel('Frequency (Hz)');
    
    disp(i);
end
sgtitle([nickname ' Params: [' num2str(morletParams(1)) ',' num2str(morletParams(2)) ']' ], 'Interpreter', 'none');


%% Multitaper Spectrogram (Chronux toolkit) (not used in figures, but another method to make spectrograms if we need to match the MUA and LFP analysis in future)
channel = 7;
i = find(good_channels == channel);

params.tapers = [10 19];
params.pad = 0;
params.Fs = fs_down;
params.fpass = [32 200];
params.err = 0;
params.trialave = 0;
winsize = 0.5;
winstep = 0.1;
movingwin = [winsize winstep];

% num_frex = (params.fpass(2)-params.fpass(1)+1)/deltaF;
% rangeF = linspace(params.fpass(1), params.fpass(2), num_frex);
% holdingM = zeros(length(rangeF), (lenTrial-leadTime)*fs_down);

numTrials = length(t_stim);
lastUsedJ = 1;

for j = 1:numTrials
    % skip bad epochs (due to mouse motion contamination)
    if find(badEpochs == j)%j ~= 6
        continue;
    end
    
    %add the epochs with no motion
    %tempTrialData = data_good(i, ind(lastUsedJ)+leadTime*fs:(ind(j)+(lenTrial*fs)-1));
    A= ECoG_denoise(ind_down(lastUsedJ)+leadTime*fs_down:ind_down(lastUsedJ)-1, i);
    B = ECoG_denoise(ind_down(j):ind_down(j)+(lenTrial*fs_down)-1, i);
    tempTrialData = [A;B];
    
    % lower [] better in low freq
    [M,t,f]=mtspecgramc(tempTrialData,movingwin,params);
    if j == 1
        holdingM = zeros(length(t),length(f));
    end
    %M = my_morlet(tempTrialData, fs_down, params.fpass(1), params.fpass(2), deltaF, morletParams);
    holdingM = holdingM + M;
    lastUsedJ = j+1;
    disp(lastUsedJ);
end

%average
div = numTrials - length(badEpochs);
holdingM = holdingM/div;

%subtract baseline
%baseline = median(holdingM(:,1:(leadTime*fs_down*-1)-1), 2);
t = t-2;
baseline = median(holdingM(1:find(t >= -1,1),:), 1);
holdingM = 10*(log10(holdingM) - log10(baseline));
t_trial = leadTime:1/fs_down:lenTrial-1/fs_down;

figure;
%plot_matrix(holdingM,t, f);
imagesc(t,f,holdingM');axis xy; colorbar; title('Spectrogram');

%h = pcolor(t,f, holdingM');
%set(h,'EdgeColor','none');
%set(gca, 'YScale', 'log');
set(gca, 'FontSize', 25);
c = colorbar;
c.Label.String = '10*log10(signal/baseline)';
colormap jet;
caxis([0,4]);

xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
title([nickname 'Channel: ' num2str(channel)], 'Interpreter', 'none');























%% %%%%%%%%%%%%%%%% Below is old code I tested but didn't use for figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot 2D map by hand

t_trial = t_down(1:fs_down*(lenTrial-leadTime)) + leadTime;
trialAverages = zeros(16, fs_down*(lenTrial-leadTime));
numTrials = length(ind_down);

% Sum all the trials into a single matrix. The order is the electrode
% mapping order.
figure;
for i = 1:16
    subplot(4,4,i);
    for j = 1:numTrials
        if find(badEpochs == j)
            continue;
        end
        trialAverages(i,:) = trialAverages(i,:) + data_mua(i, ind_down(j)+(leadTime*fs_down):(ind_down(j)+(lenTrial*fs_down)-1));
        plot(t_trial, data_mua(i, ind_down(j)+(leadTime*fs_down):(ind_down(j)+(lenTrial*fs_down)-1)), 'Color', [0.65 0.65 0.65]); hold on;
    end
    trialAverages(i,:) = trialAverages(i,:)/(numTrials-length(badEpochs));
    %d = smoothdata(trialAverages(i,:), 'gaussian', 60);
    plot(t_trial, trialAverages(i,:), 'Color', 'k', 'LineWidth', 2);
end

% Divide matrix by number of trials to get trial averages
% trialAverages = trialAverages/numTrials;
% 
% d = smoothdata(trialAverages(electrode,:), 'gaussian', 200);
% plot(t_trial, d, 'Color', 'k', 'LineWidth', 2);
% xlim([-2 10]);
% xlabel('Time (seconds)');
% ylabel('Amplitude (\muV)');
% set(gca, 'FontSize', 25);


%% Calculate individual flashes
individualFlashes2 = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5];

[m, n] = f_plotdetailedTrialAverages2D([0 200], trialAverages, t_trial, fs_down, individualFlashes2, -0.05, 0.2, []);


%% plot trial averages for each electrode (1D map/list)
%t_trial = t_amplifier(1:fs*(lenTrial-leadTime)) + leadTime; 

onsets = zeros(16,1);
onsetidx = zeros(16,1);
range = [find(t_trial == fs*0) find(t_trial >= 0.05, 1)];
figure;
for i = 1:16
    yOffset = 2400 - 150*(i-1);
    if impedance_mapped(i) > 4e6
        plot(t_trial, trialAverages(i,:) + yOffset, 'r'); hold on;
        %plot(t_trial, trialAverages(i,:), 'r'); hold on;
    else
        plot(t_trial, trialAverages(i,:) + yOffset, 'b'); hold on;
        %plot(t_trial, trialAverages(i,:), 'b'); hold on;
        
        %find minimum/onset response
        [onsets(i,1), onsetidx(i,1)] = min(trialAverages(i,range(1):range(2)));
    end
    text(t_trial(end) + 0.02, yOffset, ['Ch ' num2str(mapping(i))]);
end

% Label graph
title('Trial Averages', 'FontSize', 16);
ylabel('Voltage (\muV)', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16);
ylim([0 2500]);
%xlim([-0.05 0.25]);

% Plot onset delay time for each electrode
% onsetidx = t_trial(onsetidx+range(1));
% onsets2D = reshape(onsetidx, [4 4])';
% figure;
% im = imagesc(onsets2D, [0.04 0.05]);
[m, n] = f_plotdetailedTrialAverages2D([0 200], trialAverages, 15, t_trial, fs, individualFlashes2, -0.05, 0.2, [], impedance_mapped, mapping);


%% plot onsets of stimulus for one channel
% figure;
% channel = 14;
% i = find(mapping == channel);
% leadTime = -1;
% lenTrial = 10; % seconds
% trialAverages = zeros(16, fs*(lenTrial-leadTime));
% numTrials = length(ind);
% avgTrial = [];
% 
% for j = 1:numTrials
%     
%     plot(t_trial, data_good(i, ind(j)+(leadTime*fs):(ind(j)+(lenTrial*fs)-1))); hold on;
% end
% 
% plot(t_trial, mean(data_good(:, ind(j)+(leadTime*fs):(ind(j)+(lenTrial*fs)-1)), 2), 'k', 'LineWidth', 2);

%% Average of trial spectrograms for single channel
channel = 14;
i = find(mapping == channel);
freq_range = [1,100];
deltaF = 0.5;
morletParams = [4,50];

badEpochs = [];

%abridged_data = data_good(:,ind(badEpochs):ind(badEpochs)
num_frex = (freq_range(2)-freq_range(1)+1)/deltaF;
rangeF = linspace(freq_range(1), freq_range(2), num_frex);
holdingM = zeros(length(rangeF), (lenTrial-leadTime)*fs);
lastUsedJ = 1;

for j = 1:numTrials
    % skip bad epochs (due to mouse motion contamination)
    if find(badEpochs == j)%j ~= 6
        continue;
    end
    
    %add the epochs with no motion
    %tempTrialData = data_good(i, ind(lastUsedJ)+leadTime*fs:(ind(j)+(lenTrial*fs)-1));
    A= data_mua(i, ind(lastUsedJ)+leadTime*fs:ind(lastUsedJ)-1);
    B = data_mua(i, ind(j):ind(j)+(lenTrial*fs)-1);
    tempTrialData = [A,B];
    
    % lower [] better in low freq
    M = my_morlet(tempTrialData, fs, freq_range(1), freq_range(2), deltaF, morletParams);
    holdingM = holdingM + M;
    lastUsedJ = j+1;
    disp(lastUsedJ);
end

%average
div = numTrials - length(badEpochs);
holdingM = holdingM/div;

%subtract baseline
baseline = median(holdingM(:,1:(leadTime*fs*-1)-1), 2);
holdingM = 10*(log10(holdingM) - log10(baseline));

figure;
h = pcolor(t_trial(1:10:end),rangeF, holdingM(:,1:10:end));
set(h,'EdgeColor','none');
set(gca, 'FontSize', 25);
c = colorbar;
c.Label.String = '10*log10(signal/baseline)';
colormap jet;
caxis([0,3.5]);

xlabel('Time (seconds)');
ylabel('Frequency (Hz)');
%title({[nickname]; ['Channel ' num2str(mapping(i)) ' Impedance: ' num2str(impedance_mapped(i)/1.0e6) ' MOhms'];['Morlet parameters: [' num2str(morletParams(1)) ',' num2str(morletParams(2)) ']']}, 'Interpreter', 'none');


%[m, n] = f_plotdetailedTrialAverages2D([0 200], trialAverages, 1, t_trial, fs_down, individualFlashes5, -0.05, 0.2, [], impedance_mapped, mapping);
