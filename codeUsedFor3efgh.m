 % Author: Madison Wilson
clear;
clc;

nickname = 'NOD18\0316_run03';
%  Manual load processed\video_runX.mat
%load 'C:\Users\Maddie\Documents\Martin Exps\Processed Data\NOD20 video\0309_run4.mat';
%  Manual load processed_EP\lpf_YYY_runX.mat
%load(['C:\Users\Maddie\Documents\Martin Exps\Processed Data\' nickname '.mat']);
load(['D:\Maddie\Martin Exps\Processed Data\' nickname '.mat']);

%% Gather time & frequency info from RHD file, map channels to array, and filter raw data into LFP and MUA
fs = rhd.frequency_parameters.amplifier_sample_rate;
impedance = [rhd.amplifier_channels.electrode_impedance_magnitude];

% mapping of array channels (Gold facing down onto cortex)
% [4, 3,  1,  16]
% [5, 6,  2,  15]
% [7, 10, 14, 13]
% [8, 9,  11, 12]
mapping = [4,3,1,16,5,6,2,15,7,10,14,13,8,9,11,12];

impedance_mapped = impedance(9:24); %cut out unused channels
good_channels = find(impedance_mapped < 10e6); % consider channels below an impedance threshold "good channels"
impedance_mapped = impedance_mapped(good_channels);

data_good = rhd.amplifier_data(9:24,:); % cut out unused channels
data_good = data_good(good_channels, :); % reorder channels to match array mapping

% Filter data below fs_down Hz (easier on ICA)
% d_low_3000 = designfilt('lowpassiir','FilterOrder',8, ...
% 'PassbandFrequency',3000,'PassbandRipple',0.2, ...
% 'SampleRate',fs);

%data_clean_filt = filtfilt(d_low_3000, data_good');
%data_down = data_clean_filt(1:fs/fs_down:end,:); % fs_down Hz sampling rate

% Create time axes for data and stimuli
t_amplifier = rhd.t_amplifier;
ind = find(diff(rhd.board_adc_data(1,:) > 0.3) == 1) + 1;
t_stim = t_amplifier(ind);

%Clear helper variables
%clear impedance;

%% plot raw data
figure;
for i = 1:size(data_good,1)
    plot(t_amplifier, data_good(i,1:end)+500*i); hold on;
end
for i = 1:length(t_stim)
    plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
end
title([nickname ' Raw Data'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');


%% Evaluate ICA
%[Aica,ica_source_activity] = ica_denoise(data_good', t_amplifier, 1:length(good_channels));
f_plot_ica(Aica,ica_source_activity, good_channels, mapping, t_amplifier);
for i = 1:length(t_stim)
    plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1);
end

%Clear helper variables
%clear data_good;

%% Remove noisy ICA components
%ica_remove = [1,2,3,4]; %Hardcode which ICA components to remove based on ICA plots in previous section
good_ind = setdiff(1:length(good_channels), ica_remove); % second input is channel to remove based on imaging artifact
ECoG_denoise = ica_source_activity(:,good_ind) * Aica(:,good_ind)';


%% Create LPF and MUA data and downsample
fs_down = 20000; % Choose downsample frequency

d_low = designfilt('lowpassiir','FilterOrder',8, ...
'PassbandFrequency',250,'PassbandRipple',0.2, ...
'SampleRate',fs);
d_bpass = designfilt('bandpassiir','FilterOrder',6, ...
'HalfPowerFrequency1',500,'HalfPowerFrequency2',3000, ...
'SampleRate',fs);

data_lfp = filtfilt(d_low, ECoG_denoise);
data_bpf = filtfilt(d_bpass, ECoG_denoise);

data_lfp_down = data_lfp(1:fs/fs_down:end,:); % fs_down Hz sampling rate
data_bpf_down = data_bpf(1:fs/fs_down:end,:); % fs_down Hz sampling rate
%data_bpf_down = data_bpf_down - median(data_bpf_down,2);
t_down = t_amplifier(1:fs/fs_down:end);


%% Plot all LFP and MUA (good) channels
figure;
for i = 1:size(data_lfp_down,2)
    plot(t_down, data_lfp_down(:,i)+250*i); hold on;
end
for i = 1:length(t_stim)
    plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1);
end
title(['LFP for ' nickname], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');

figure;
for i = 1:size(data_bpf_down,2)
    plot(t_down, data_bpf_down(:,i)+100*i); hold on;
end
for i = 1:length(t_stim)
    plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1);
end
title(['MUA for ' nickname], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');


%% MUA Trial Averages (For Figure 3e)
leadTime = -0.2;
lenTrial = 4;
badEpochs = [2,5,7,10];
lightLength = 5;

d_low_100 = designfilt('lowpassiir','FilterOrder',8, ...
'PassbandFrequency',100,'PassbandRipple',0.2, ...
'SampleRate',fs_down);

data_mua = filtfilt(d_low_100, data_bpf_down.^2);
%data_mua = sqrt(data_mua);

%data_mua = data_mua(1:fs_down/500:end,:);
%t_down = t_down(1:fs_down/500:end);

ind_down = zeros(size(ind));
for i = 1:length(ind_down)
    ind_down(i) = find(t_down >= t_stim(i),1);
end

t_trial = leadTime:1/fs_down:lenTrial-1/fs_down;

s = zeros(length(good_channels),1);
trialAverages_mua = zeros(round((lenTrial - leadTime)*fs_down), length(good_channels));
figure;
for j = 1:length(good_channels)
    temp = zeros(round((lenTrial - leadTime)*fs_down), length(t_stim)-length(badEpochs));
    lastUsedI = 1;
    count = 1;
    for i = 1:length(t_stim)
        if find(badEpochs == i)
            continue;
        end
        %ind = find(t_down >= t_stim(i),1);
        A = data_mua(ind_down(i)+leadTime*fs_down:ind_down(i)-1, j);
        %A = data_mua(ind_down(lastUsedI)+leadTime*fs_down:ind_down(lastUsedI)-1, j);
        B = data_mua(ind_down(i):ind_down(i)+(lenTrial*fs_down)-1, j);
        tempTrialData = [A;B];

        %temp(:,count) = tempTrialData;
        temp(:,count) = smoothdata(tempTrialData, 'gaussian', fs_down*0.0001);
        %temp(:,lastUsedJ) = smoothdata(data_mua(ind+leadTime*fs_down:ind+lenTrial*fs_down-1/fs_down, j),'gaussian', fs_down*0.01);
        lastUsedI = i + 1;
        count = count + 1;
    end
    %trialAverages_mua(:,j) = mean(temp,2);
    trialAverages_mua(:,j) = smoothdata(mean(temp,2), 'gaussian', fs_down*0.015); %mean(temp,2);
    
%     signal = trialAverages_mua(find(t_trial >= 0,1):(find(t_trial >= 0,1)-leadTime*fs_down));
%     noise = trialAverages_mua((find(t_trial >= 0,1)+leadTime*fs_down):find(t_trial >= 0,1));
%     s = snr(signal, noise);
    data_chunk = mean(temp(find(t_trial >= 0,1):(find(t_trial >= 0,1)-leadTime*fs_down),:),2);
    signal = findpeaks(data_chunk, 'SortStr','descend','NPeaks',1);
    noise = mean(trialAverages_mua((find(t_trial >= 0,1)+leadTime*fs_down):find(t_trial >= 0,1), j));
    s(j) = 20*log(signal/noise);
    
    %trialAverages_mua(:,j) = sqrt(trialAverages_mua(:,j));
    sd = smoothdata(std(temp,0,2), 'gaussian', fs_down*0.015) ./ sqrt(size(temp,2));
    curve1 = trialAverages_mua(:,j) + sd;
    curve2 = trialAverages_mua(:,j) - sd;
    subplot(4,4,find(mapping == good_channels(j)));
    fill([t_trial, fliplr(t_trial)], [curve1; flipud(curve2)], [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on;
    plot(t_trial, trialAverages_mua(:,j), 'k', 'LineWidth', 2); hold on;
    ylim([20 150]);
    xlim([-0.2 0.4]);%([leadTime, lenTrial]);
    fill([0,lightLength,lightLength,0], [0,0,2,2], [1 0.5 0.25], 'EdgeColor', 'none');
    title(['Channel ' num2str(good_channels(j)) ', SNR = ' num2str(s(j))]);
    disp(j);
end
sgtitle([nickname ' Average MUA' ], 'Interpreter', 'none');

%Clear helper variables
clear A; clear B; clear tempTrialData; clear temp;


%% MUA single LED pulse trial averages
pulseTimes = 0:0.5:3.9;%0:0.2:4.8;
leadTimeSingle = -0.1;
lenTrialSingle = 0.2;

trialAverages_mua = zeros(round((lenTrialSingle - leadTimeSingle)*fs_down), length(good_channels));
figure;
for j = 1:length(good_channels)
    temp = zeros(round((lenTrialSingle - leadTimeSingle)*fs_down), (length(t_stim)-length(badEpochs))*length(pulseTimes));
    lastUsedI = 1;
    count = 1;
    for i = 1:length(t_stim)
        if find(badEpochs == i)
            continue;
        end
        %ind = find(t_down >= t_stim(i),1);
        for k = 1:length(pulseTimes)
            A = data_mua(pulseTimes(k)*fs_down+(ind_down(i)+leadTimeSingle*fs_down:ind_down(i)-1), j);
            %A = data_mua(pulseTimes(k)*fs_down+(ind_down(lastUsedI)+leadTimeSingle*fs_down:ind_down(lastUsedI)-1), j);
            B = data_mua(pulseTimes(k)*fs_down+(ind_down(i):ind_down(i)+(lenTrialSingle*fs_down)-1), j);
            tempTrialData = [A;B];
            %temp(:,count) = tempTrialData;
            temp(:,count) = smoothdata(tempTrialData, 'gaussian', fs_down*0.0001);
            count = count + 1;
        end
        lastUsedI = i + 1;
    end
    trialAverages_mua(:,j) = smoothdata(mean(temp,2), 'gaussian', fs_down*0.001);
    % trialAverages_mua(:,j) = mean(temp,2); %
    subplot(4,4,find(mapping == good_channels(j)));
    plot(leadTimeSingle:1/fs_down:lenTrialSingle-1/fs_down, trialAverages_mua(:,j), 'k'); hold on;
    ylim([50 150]);
    xlim([leadTimeSingle, lenTrialSingle]);
    %fill([0,lightLength,lightLength,0], [0,0,1,1], [1 0.5 0.25], 'EdgeColor', 'none');
    title(['Channel ' num2str(good_channels(j))]);
end
sgtitle([nickname ' Average MUA' ], 'Interpreter', 'none');

%Clear helper variables
clear A; clear B; clear tempTrialData; clear temp;


%% Plot select trial averages
figure;
ch = [12,14,6,4]; %Select which channels to plot
for i = 1:length(ch)
    plot(leadTime:1/fs_down:lenTrial-1/fs_down, trialAverages_mua(:,find(good_channels == ch(i)))+50*i, 'k'); hold on;
end
fill([0,0.1,0.1,0], [20,20,21,21], [1 0.5 0.25], 'EdgeColor', 'none');
xlim([-0.2 0.4]);
xlabel('Latency (s)');
ylabel ('Amplitude (\muV^2)');


%% Spike detection using findpeaks.m
% Set epochs that have noise to zero so they aren't picked up by the threshold
badEpochs = [2,5,7,10];%[2,4,5,6,7,10]; % Define here if not already

data_bpf_noise_nan = data_bpf_down;
for i = 1:length(badEpochs)
    idx = find(t_down >= t_stim(badEpochs(i)),1);
    data_bpf_noise_nan(idx:idx+fs_down*10,:) = nan;
end
data_bpf_noise_nan(130*fs_down:end,:) = nan;
data_bpf_noise_nan(4*fs_down:9.15*fs_down,:) = nan;
data_bpf_noise_nan(27.06*fs_down:29.03*fs_down,:) = nan;
data_bpf_noise_nan(53.5*fs_down:55.95*fs_down,:) = nan;
data_bpf_noise_nan(73.5*fs_down:75.19*fs_down,:) = nan;
data_bpf_noise_nan(109.8*fs_down:111.7*fs_down,:) = nan;

std_multiplier = 3.25;
thresh = std_multiplier*nanstd(data_bpf_noise_nan); % modify SD multiplier based on SNR (3, 3.5, 4, or 4.5 usually good)

spikes = zeros(size(data_bpf_noise_nan));
for i = 1:length(good_channels)
    [pks, loc] = findpeaks(-1*data_bpf_noise_nan(:,i), 'MinPeakHeight', thresh(i), 'MinPeakDistance', fs_down*0.001);
    spikes(loc,i) = 1;
end

% Bin spikes and sum all channels
spikes_sum = sum(spikes,2);

% Plot spike raster and spiking vector
figure;
for i = 1:length(good_channels)
    temp = (spikes(:,i) ~= 0);
    scatter(t_down(temp), spikes(temp,i)*(i+1), '.'); hold on;
end
for i = 1:length(t_stim)
    plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
end
plot(t_down, spikes_sum);
%xlim([0 60]);
xlabel('Time (seconds)');
ylabel('Spike Rate (Spikes/sec) for population - light');
title(nickname, 'Interpreter', 'none');

%Clear helper variables
clear spikes_sum; clear pks; clear loc; clear temp;


%% plot select electrode channels, showing their MUA, MUA spikes, and thresholds (Figure 3b)
ch = [12,14,6,4]; % channel #, not location in array
y_offset = [0,100,200,300]; % Change this based on MUA data to evenly space MUA traces while keeping same y scale
    
spikes_plot = spikes;
spikes_plot(spikes_plot == 0) = NaN;

figure;
for j = 1:length(ch)
    i = find(good_channels == ch(j));
    
    plot([0 t_down(end)], [-1*thresh(i) -1*thresh(i)] + y_offset(j), 'Color', 'r', 'LineWidth', 2); hold on;
    plot(t_down, data_bpf_noise_nan(:,i) + y_offset(j), 'LineWidth', 1, 'Color', 'k'); hold on;
    scatter(t_down, spikes_plot(:,i)*(-1*thresh(i)) + y_offset(j) - 10, 'MarkerFaceColor', 'k'); hold on;
end

ylim([-50 350]);
xlim([10 20]);%([10 10.5]);
ylabel('Amplitude (\muV)');
xlabel('Time (seconds)');
title([nickname ' Channel: ' num2str(ch) ' SD: ' num2str(std_multiplier) ' ICA: ' num2str(ica_remove)], 'Interpreter', 'none');
set(gca, 'FontSize', 25);


%% Define periods of light and dark (for categorizing MUA spikes into ones happening during "light" or "dark" periods)
lenStim = 4; %time that the light was flashing (e.g. for 5 Hz, 2 sec stimuli => lenStim = 2)
numTrials = length(t_stim);
%badEpochs = [8,11,12,16,19]; %This is determined by-hand: looked at Martin's video code and choosing epochs/trials that had lots of motion detected

goodEpochs = setdiff(1:numTrials, badEpochs);
t_dark = zeros(length(goodEpochs)+1,2);
for i = 1:length(goodEpochs)+1
    if i == 1
        t_dark(i,:) = [1,find(t_down >= t_stim(goodEpochs(1)),1)];
    elseif i == length(goodEpochs)+1
        t_dark(i,:) = [find(t_down >= t_stim(goodEpochs(end))+lenStim,1), length(t_down)];
    else
        t_dark(i,:) = [find(t_down >= t_stim(goodEpochs(i-1))+lenStim,1), find(t_down >= t_stim(goodEpochs(i)),1)];
    end   
end

%t_light = [find(t_down >= t_stim(1),1),find(t_down >= t_stim(numTrials)+lenStim,1)]; %Define dark as before and after stimuli
t_light = [find(t_down >= t_stim(1),1), length(t_down)]; %Define dark as period just before stimuli onset
%t_light = []; %Define dark including time between light pulses

% double check values
disp(t_down(t_dark(:,1)));
disp(t_down(t_dark(:,2)));

%% Determine light and dark spike times
tLightOnset = 0:0.5:3.5;
lightLength = 0.1;

spikes_light = spikes; %ones(size(spikes));
%spikes_light_select = spikes;
spikes_dark = spikes;

for i = 1:size(t_dark,1)
    spikes_light(t_dark(i,1):t_dark(i,2),:) = 0; %set spikes that occur during dark to 0
end
%Remove spikes between light flashes
for i = 1:length(t_stim)
    for j = 1:length(tLightOnset)-1
        spikes_light((t_stim(i)+(tLightOnset(j)+lightLength))*fs_down:(t_stim(i)+tLightOnset(j+1))*fs_down,:) = 0;
    end
end
spikes_dark(t_light(1):t_light(2),:) = 0; %set spikes that occur during light to 0
%spikes_dark(1:214000,:) = 0; %(optional) Make same length as light trial time

%% Plot light & dark spike raster and spiking vector
figure;
for i = 1:length(good_channels)
    temp = (spikes_light(:,i) ~= 0);
    scatter(t_down(temp), spikes_light(temp,i)*i, 10, 'b', 'filled'); hold on;
    temp2 = (spikes_dark(:,i) ~= 0);
    scatter(t_down(temp2), spikes_dark(temp2,i)*i, 10, 'r', 'filled'); hold on;
    %disp(sum(temp2));
end
for i = 1:length(t_stim)
    plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1);
end
%xlim([0 60]);
xlabel('Time (seconds)');
ylabel('Spikes');
legend('Light', 'Dark');
title(nickname, 'Interpreter', 'none');

disp(['light spike count = ' num2str(sum(spikes_light,1))]);

%% PLV spectrogram (only at spike time points to reduce computation time) using Xin/Chronux multitaper method
% User parameters
tapers = [5,9];
params = struct('tapers', tapers,'pad',0,'Fs',fs_down,'fpass',[0,150]);
win = 1*fs_down; % window size around each MUA spike: win =(seconds)*fs_down

%Compute frequency range/resolution
nfft=max(2^(nextpow2(win)+params.pad),win);
[f,~]=getfgrid(params.Fs,nfft,params.fpass);
nFreq = length(f);

% Combine spike times across all channels so only need to compute
% spectrogram once, then will sort spikes by channel later
ind_light_all = find(sum(spikes_light,2) >= 1);
ind_dark_all = find(sum(spikes_dark,2) >= 1);
ind_dark_all(ind_dark_all < win/2) = []; %get rid of spikes too close to beginning of recording
ind_dark_all(ind_dark_all > length(t_down) - win/2) = []; % get rid of spikes too close to end of recording

tic;
phase_specgram_light = f_plv_XL(data_lfp_down, ind_light_all, win, nFreq, params); % Takes ~1 minute
toc; tic;
phase_specgram_dark = f_plv_XL(data_lfp_down, ind_dark_all, win, nFreq, params); %Takes ~2 minutes
toc;

%% Determine phases at MUA spikes for each channel, calculate PLV, and plot PLV vs. frequency (Figure 3d)
% This section can take awhile to run depending on number of spikes for each channel.
% The time grows exponentially: 50 spikes => 5 seconds, 170 spikes => 30 seconds, 320 spikes => 400 seconds.

% User parameter
nBoots = 2000; % Bootstrapping method, "generate nBoot bootstrapped samples"

[phi_light, phi_dark, pval, is_significant] = f_plot_PLV_vs_freq_221126(phase_specgram_light, phase_specgram_dark, spikes_light, spikes_dark, ...
        ind_light_all, ind_dark_all, t_down, nBoots, nFreq, win, f, mapping, good_channels, params, nickname);


%% Calculate PLV of MUA spikes to the LFP of each channel, not just the same channel. For 2D color plotting. (Takes a few minutes)
% User parameters
freq = [4,6]; % frequency band of interest (e.g. delta => [1,4], theta => [5,9], beta => [12,32])
nBoots = 2000; % number of bootstrap samples to take

%Calculate 2D PLV values
[plv_light_final_2D, plv_dark_final_2D, pval_color] = f_2D_PLV(phase_specgram_light, phase_specgram_dark, spikes_light, spikes_dark, ...
    ind_light_all, ind_dark_all, win, t_down, freq, f, nBoots, params);


%% Plot color maps for all frequencies (Figure 3h)
% User parameter
clim = [0.15, 0.35]; % choose color scale that best shows difference between light and dark

% Color map for light
figure;
for i = 1:length(good_channels)
    temp = NaN(16,1);
    temp(good_channels) = mean(plv_light_final_2D(i,:,:),3); %Average PLV across frequencies
    itpcMap = temp;
    itpcMap = itpcMap(mapping);
    itpc2D = reshape(itpcMap, [4 4])';
    itpc2D = fillmissing(itpc2D, 'linear',1, 'EndValues','nearest');
    itpc2D = flipud(itpc2D);
    itpc2D = interp2(itpc2D,2);

    subplot(4,4,find(mapping == good_channels(i)));
    cmap1 = pcolor(itpc2D);
    cmap1.FaceColor = 'interp';
    cmap1.EdgeColor = 'none';
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    c = colorbar;
    c.Label.String = 'PLV';
    caxis(clim);
    %title(['PLV of ' nickname], 'Interpreter', 'none');
    %title(['n = ' num2str(length(hAngles))]);
end
sgtitle(['Light PLV of ' nickname ' for ' num2str(freq(1)) '-' num2str(freq(2)) ' Hz'], 'Interpreter', 'none');

% Color map for dark
figure;
for i = 1:length(good_channels)
    temp = NaN(16,1);
    temp(good_channels) = mean(plv_dark_final_2D(i,:,:),3); %Average PLV across frequencies
    itpcMap = temp;
    itpcMap = itpcMap(mapping);
    itpc2D = reshape(itpcMap, [4 4])';
    itpc2D = fillmissing(itpc2D, 'linear',1, 'EndValues','nearest');
    itpc2D = flipud(itpc2D);
    itpc2D = interp2(itpc2D,2);

    subplot(4,4,find(mapping == good_channels(i)));
    cmap1 = pcolor(itpc2D);
    cmap1.FaceColor = 'interp';
    cmap1.EdgeColor = 'none';
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    c = colorbar;
    c.Label.String = 'PLV';
    caxis(clim);
    %title(['PLV of ' nickname], 'Interpreter', 'none');
    %title(['n = ' num2str(length(hAngles))]);
end
sgtitle(['Dark PLV of ' nickname ' for ' num2str(freq(1)) '-' num2str(freq(2)) ' Hz'], 'Interpreter', 'none');


%% Polar histograms of the PLV (figure 3g)
% User parameter
freq_single = 4.5; % choose a single frequency to look at (best to use one in the middle of previously used frequency range)

f_plv_histograms_221125(phase_specgram_light, phase_specgram_dark, spikes_light, spikes_dark, ind_light_all, ind_dark_all, win, t_down, freq_single, f, mapping, good_channels);






