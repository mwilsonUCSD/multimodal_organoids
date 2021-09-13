% Author: Madison Wilson
clear;
clc;

nickname = 'NOD18\0302_run1';

load(['C:\Users\Maddie\Documents\Martin Exps\Processed Data\' nickname '.mat']);

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

if isfield(rhd, 'board_adc_data')
    disp("ADC found");
    ind = find(diff(rhd.board_adc_data(1,:) > 0.3) == 1) + 1;
    t_stim = t_amplifier(ind);
end


%% plot raw data
figure;
for i = 1:size(data_good,1)
    plot(t_amplifier, ECoG_denoise(1:end,i)+500*i); hold on;
end
if isfield(rhd, 'board_adc_data')
    for i = 1:length(t_stim)
        plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
    end
end

title([nickname ' Raw Data'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');


%% Evaluate ICA
%[Aica,ica_source_activity] = ica_denoise(data_good', t_amplifier, 1:length(good_channels));
f_plot_ica(Aica,ica_source_activity, good_channels, mapping, t_amplifier);
if isfield(rhd, 'board_adc_data')
    for i = 1:length(t_stim)
        plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
    end
end


%% Remove noisy ICA components
%ica_remove = [];
good_ind = setdiff(1:length(good_channels), ica_remove); % second input is channel to remove based on imaging artifact
ECoG_denoise = ica_source_activity(:,good_ind) * Aica(:,good_ind)';


%% Filter and plot certain range (look for theta spindles)
f_bpass = [5,8];

theta_pass = designfilt('bandpassiir','FilterOrder',6, ...
'HalfPowerFrequency1',f_bpass(1),'HalfPowerFrequency2',f_bpass(2), ...
'SampleRate',fs);

data_theta = filtfilt(theta_pass, ECoG_denoise);
data_theta_down = data_theta(1:fs/fs_down_lfp:end,:);

figure;
for i = 1:size(data_theta_down,2)
    plot(t_down_lfp, data_theta_down(:,i)+500*i); hold on;
end
if isfield(rhd, 'board_adc_data')
    for i = 1:length(t_stim)
        plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
    end
end
title(['[' num2str(f_bpass) '] Hz filtered signal for ' nickname ' ICA: [' num2str(ica_remove) ']'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');


%% Envelope analysis
f_bpass = [13,30]; % Hz
thresh = 3; % # of stds for crossing threshold
t = [1,20]; % seconds
step = 0.5; % seconds
doSave = 0;

f_plot_envelope_map(ECoG_denoise, f_bpass, thresh, t, step, fs, t_amplifier, nickname, good_channels, mapping, ica_remove, doSave);


%% Create LPF and MUA data and downsample
fs_down = 6000; % Choose downsample frequency
fs_down_lfp = 500;

d_low = designfilt('lowpassiir','FilterOrder',8, ...
'PassbandFrequency',250,'PassbandRipple',0.2, ...
'SampleRate',fs);
d_bpass = designfilt('bandpassiir','FilterOrder',6, ...
'HalfPowerFrequency1',500,'HalfPowerFrequency2',3000, ...
'SampleRate',fs);
% d_60 = designfilt('bandstopiir','FilterOrder',8, ...
% 'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
% 'SampleRate',fs);

data_lfp = filtfilt(d_low, ECoG_denoise);
% data_lfp = filter(d_60, data_lfp);
data_bpf = filtfilt(d_bpass, ECoG_denoise);

%data_lfp_down_long = data_lfp(1:fs/fs_down:end,:);
data_lfp_down = data_lfp(1:fs/fs_down_lfp:end,:); % fs_down Hz sampling rate
% data_bpf_down = data_bpf(1:fs/fs_down:end,:); % fs_down Hz sampling rate
%data_bpf_down = data_bpf_down - median(data_bpf_down,2);
% t_down = t_amplifier(1:fs/fs_down:end);
t_down_lfp = t_amplifier(1:fs/fs_down_lfp:end);



%% Envelope Analysis loop
%doSave = 1;

oC = [0.8500, 0.3250, 0.0980];
cC = [0, 0.4470, 0.7410];
colors = [cC; cC; oC; oC; oC; oC; oC; cC; cC; oC; cC; cC; cC; cC; cC; cC];

% Plot general LFP
%thresh = 3*std(data_lfp,[],1);
figure;
for i = 1:size(data_good,1)
    %plot(t_amplifier, ECoG_denoise(:,i)+500*i, 'Color', colors(good_channels(i),:)); hold on;
    plot(t_down_lfp, abs(data_lfp_down(:,i)) + 1000*i, 'Color', colors(good_channels(i),:)); hold on;
    %plot(t_amplifier, data_bpf(:,i)+100*i, 'Color', colors(good_channels(i),:)); hold on;
    %plot(get(gca,'xlim'), 100*i-[thresh(i), thresh(i)]);
end
if isfield(rhd, 'board_adc_data')
    for i = 1:length(t_stim)
        plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
    end
end

title([nickname ' LFP'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');
set(gcf, 'Position', [1,41,1920,963]);
set(gca, 'FontSize', 15);

ylim([0,9000]);
%xlim([1 60]);

%% Burst movie (supplementary figure 5)
maxval = 500;
data_lfp_abs_mapped = reshape(abs(data_lfp_down(:,mapping)), [], 4, 4);

v = VideoWriter('bursts.avi');
open(v);

%pc = zeros(5,5);
%alphadata = zeros(20,20);
%alphadata([3,8,13,18],[3,8,13,18]) = 1;
map = [1:-0.01:0; 1:-0.01:0; ones(1,101)]';

%f = figure('visible','off');
figure;
I = imread('NOD18_brightfield_RGB.jpg');
imshow(I); hold on;
c = colorbar('FontSize', 15); colormap(map); caxis([0 maxval]); c.Label.String = '|Amplitude| (\muV)';
for i = 6.47*fs_down_lfp:0.004*fs_down_lfp:6.55*fs_down_lfp%length(t_down_lfp)
    %pc(2:5,1:4) = squeeze(abs(data_lfp_mapped(i,:,:)))';
    %pc = flipud(pc);
    %P = pcolor(pc);
    %set(P, 'AlphaData', alphadata);
    %set(h, 'EdgeColor', 'none');
    %caxis([0 500]);
    %pause(0.01);
    %c = colorbar;
    pc = squeeze(data_lfp_abs_mapped(i,:,:))';
    rarray = f_drawrectangles(pc, maxval);
    %c.Label.String = '10*log10(signal/baseline)';
    %colormap jet;
    %title([num2str(i/fs_down_lfp) ' s']);
    %text(50,1550,[num2str(i/fs_down_lfp) ' s'], 'FontSize', 30);
    frame = getframe(gca);
    %writeVideo(v,frame);
    imwrite(frame.cdata, ['burst_' num2str(i) '.jpg']);
    delete(rarray);
end

close(v);

% saveas(gcf, ['_210614_burst_analysis\10s_snippets\' nickname '__LFP'], 'jpeg');
% 
% if doSave == 1
%     saveas(gcf, ['_210614_burst_analysis\' nickname '__LFP_'], 'jpeg');
% end
% 
% 
% % Heat map analysis
% freq_bands = [1,4;5,7;8,12;13,30;31,59;61,100];
% %f_bpass = [13,30]; % Hz
% thresh = 3; % # of stds for crossing threshold
% t = [1,20]; % seconds
% step = 0.5; % seconds
% t_omit = [];
% dt = [-1, 1];
% 
% for i = 1:6 % loop through all freq bands
%     f_plot_envelope_map(ECoG_denoise, freq_bands(i,:), thresh, dt, t, step, fs, t_amplifier, nickname, good_channels, mapping, ica_remove, doSave, t_omit);
%     close all;
% end
% 
% disp('Done!');

%% Create mask of bursts

d_bpass = designfilt('bandpassiir','FilterOrder',6, ...
'HalfPowerFrequency1',12,'HalfPowerFrequency2',32, ...
'SampleRate',fs);

data_beta = filtfilt(d_bpass, ECoG_denoise);
%data_beta = data_beta(1:fs/fs_down_lfp:end,:);

% figure;
% plot(t_amplifier, smoothdata(data_beta.^2));
% xlim([17 35]);

thresh = 300;
mask = smoothdata(data_beta.^2) > thresh;

mask = smoothdata(sum(mask,2));
mask(mask >= 1) = 1;
mask(mask < 1) = 0;
%figure; plot(t_down_lfp, mask); xlim([17 35]); ylim([0 2]);
figure; plot(t_amplifier, mask); xlim([17 35]); ylim([0 2]);


%% Frequency analysis
snippet = [1,1.7;2.3,3.1;3.7,4.5;5.6,6.4]*fs_down_lfp; %[1.7,2.3;3.1,3.7;4.5,5.6;6.4,7.2;8.8,9.7]*fs_down_lfp; %[1,10]*fs_down_lfp; %[3.3,4; 5,5.8; 7.7,8.7;10.2,11;12.4,13.1]*fs_down_lfp; %
win_length = 2 * fs_down_lfp;
win_overlap = 1 * fs_down_lfp;
nfft = [];

%[pxx2,f,pxxc] = pwelch(data_lfp_down(snippet(1):snippet(2),:), win_length, win_overlap, nfft, fs_down_lfp, 'ConfidenceLevel', 0.95);

%[pxx,f] = pwelch(data_lfp_down(snippet(1):snippet(2),1), [], [], [], fs_down_lfp);

figure; hold on;
for i = 1:size(data_good,1)
    %[pxx,f,pxxc] = pwelch(data_lfp_down(snippet(1):snippet(2),i), win_length, win_overlap, nfft, fs_down_lfp, 'psd');
    pxx_tot = zeros(129,1);
    for j = 1:4
        %[pxx,f,pxxc] = pwelch(data_lfp_down(snippet(j,1):snippet(j,2),i), win_length, win_overlap, nfft, fs_down_lfp, 'psd');
        [pxx,f,pxxc] = pwelch(data_lfp_down(snippet(j,1):snippet(j,2),i),[],[],[],fs_down_lfp);
        pxx_tot = pxx+pxx_tot;
    end
    pxx_tot = pxx_tot / 4;
    %[pxx,f,pxxc] = pwelch(data_lfp_down(snippet(1):snippet(2),i),[],[],[],fs_down_lfp);
    plot(f, 10*log10(pxx_tot), 'Color', colors(good_channels(i),:));
    %[pxx,f,pxxc] = pwelch(data_lfp_down(snippet(j,1):snippet(j,2),i));
    %loglog(f, pxx, 'Color', colors(good_channels(i),:));
    %plot(f, (10*log10(pxx)).*f, 'Color', colors(good_channels(i),:));
    %plot(f, 10*log10(pxx2(:,i)./pxx(:,i)), 'Color', colors(good_channels(i),:));
    %plot(f, 10*log10(pxxc)', '-.', 'Color', colors(good_channels(i),:));
end
%set(gca, 'YScale', 'log');
%set(gca, 'XScale', 'log');
%ylim([-20 35]);
%xlim([0 100]);
title([nickname ' LFP time: ' num2str(snippet/fs_down_lfp)], 'Interpreter', 'none');


%%
figure;
for i = 1:size(data_good,1)
    %plot(t_amplifier, ECoG_denoise(:,i)+500*i, 'Color', colors(good_channels(i),:)); hold on;
    plot(t_down_lfp, data_lfp_down(:,i)+500*i, 'Color', colors(good_channels(i),:)); hold on;
end

title([nickname ' LFP'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');
set(gcf, 'Position', [1,41,1920,963]);
set(gca, 'FontSize', 15);


%% Plot all LFP and MUA (good) channels
figure;
hold on;
for i = 1:size(data_lfp_down,2)
    %plot(t_down, data_lfp_down_long(:,i)+500*i); hold on;
    plot(t_down_lfp, data_lfp_down(:,i)+500*i);
end
if isfield(rhd, 'board_adc_data')
    for i = 1:length(t_stim)
        plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1);
    end
end
title(['LFP for ' nickname ' ICA: [' num2str(ica_remove) ']'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');

figure;
hold on;
for i = 1:size(data_bpf_down,2)
    plot(t_down, data_bpf_down(:,i)+100*i);
end
if isfield(rhd, 'board_adc_data')
    for i = 1:length(t_stim)
        plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1);
    end
end
title(['MUA for ' nickname ' ICA: [' num2str(ica_remove) ']'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');

%% Chronux/custom frequency analysis
params.tapers = [3 5];
params.pad = 2;
params.Fs = fs_down;
params.fpass = [0 150];
params.err = 0;
params.trialave = 0;
winsize = 0.1;
winstep = 0.001;
movingwin = [winsize winstep];

tic;
[~,phi,allphi,t,f]=mtspecgramc_mnw(data_lfp_down,movingwin,params); % This takes around 1 minute, modified from Chronux toolbox's mtspecgramc() function
toc;

figure; plot_matrix(phi(:,:,1), t, f, 'n'); hold on;


%% Morlet frequency analysis (figure 4d)
% User parameters
do_mask = 0; % 0 = no mask, 1 = bursts, 2 = silent periods
time_range = [26,28]; % seconds
freq_range = [4,10]; % frequency range of interest
deltaF = 0.05; % frequency resolution [Hz]
morletParams = [2,10]; % Range of cycles (low number of cycles => better temporal resolution, high # cycles => better freq resolution)
cbar = [2,10];

% Calculate frequency range
num_frex = (freq_range(2)-freq_range(1)+1)/deltaF;
rangeF = linspace(freq_range(1), freq_range(2), num_frex);

time_range = time_range*fs_down_lfp;
if time_range(1) == 0
    time_range(1) = 1;
end
time_range = time_range(1):time_range(2);

if do_mask
    mask_trunc = boolean(mask(time_range));
end
M_tot = zeros(length(good_channels),length(rangeF));

%spectrum = zeros(16,750);

figure;
for i = 1:length(good_channels)
    
    % Morlet-based spectrogram
    M = my_morlet(data_lfp_down(time_range,i)', fs_down_lfp, freq_range(1), freq_range(2), deltaF, morletParams);
    
    %subtract baseline
    %baseline = median(holdingM(:,1:(leadTime*fs_down*-1)-1), 2);
    baseline =  mean(M(:,1:250), 2); %mean(M, 2); %
    M = 10*(log10(M) - log10(baseline));
    %M = 10*(log10(M));

    subplot(4,4,find(mapping == good_channels(i)));
    h = pcolor(t_down_lfp(time_range(1:10:end)),rangeF, M(:,1:10:end));
    set(h,'EdgeColor','none');
    %set(gca, 'YScale', 'log');
    %set(gca, 'FontSize', 25);
    c = colorbar;
    c.Label.String = '10*log10(signal/baseline)';
    colormap jet;
    caxis(cbar);
    %set(c, 'visible', 'on');
    title(['Channel: ' num2str(good_channels(i))]);
    xlabel('Time (seconds)');
    ylabel('Frequency (Hz)');
    
%     if do_mask
%         if do_mask == 1
%             M(:,~mask_trunc) = 0;
%         else
%             M(:,mask_trunc) = 0;
%         end
%         plot(rangeF, sum(M,2)./sum(mask_trunc), 'Color', colors(good_channels(i),:)); hold on;
%         M_tot(i,:) = sum(M,2)./sum(mask_trunc);
%     else
%         plot(rangeF, mean(M,2), 'Color', colors(good_channels(i),:)); hold on;
%         M_tot(i,:) = mean(M,2);
%     end
    
    %Average morlet spectrum during spikes
    %spectrum(i,:) = mean(M(:,boolean(binned_cor(:,i))),2);
    
    disp(i);
end
sgtitle([nickname ' Params: [' num2str(morletParams(1)) ',' num2str(morletParams(2)) '] ICA: [' num2str(ica_remove) ']'], 'Interpreter', 'none');

% Plot average PSD traces
%Cortex
%plot(rangeF,mean(M_tot([8,9,11,12,13,14,15,16],:),1), 'b', 'LineWidth', 2); hold on;
% Organoid
%plot(rangeF,mean(M_tot([3,4,5,6,7,10],:),1), 'r', 'LineWidth',2);

%xlabel('Frequency (Hz)');
%ylabel('Power (uV^2)');

%% Remove labels from 2D plots
for i = 1:16
    subplot(4,4,i);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'xlabel', []);
    set(gca, 'ylabel', []);
    set(gca, 'title', []);
end

%% Using M_tot, compute avg frequency powers for each freq band and compare cortex:org
d_band_power = mean(M_tot(:, (rangeF <= 4)),2);
t_band_power = mean(M_tot(:, (4 < rangeF) & (rangeF <= 7)),2);
a_band_power = mean(M_tot(:, (7 < rangeF) & (rangeF <= 12)),2);
b_band_power = mean(M_tot(:, (12 < rangeF) & (rangeF <= 32)),2);
lg_band_power = mean(M_tot(:, (32 < rangeF) & (rangeF <= 59)),2);
hg_band_power = mean(M_tot(:, (61 < rangeF) & (rangeF <= 150)),2);

combined_band_power_awake = [d_band_power, t_band_power, a_band_power, b_band_power, lg_band_power, hg_band_power];
orgChannels = [3,4,5,6,7,10];
cortexChannels = [8,9,11,12,13,14,15,16];
% meanCortex = mean(combined_band_power_awake([1,2,8,9,11,12,13,14,15,16],:),1);
% stdCortex = std(combined_band_power_awake([1,2,8,9,11,12,13,14,15,16],:),[],1);
% meanOrg = mean(combined_band_power_awake([3,4,5,6,7,10],:),1);
% stdOrg = std(combined_band_power_awake([3,4,5,6,7,10],:),[],1);

%Pair and compute ratios
% pairings_awake = zeros(length(cortexChannels)*length(orgChannels),6);
% for i = 1:length(cortexChannels)
%     for j = 1:length(orgChannels)
%         %disp((i-1)*length(orgChannels)+j);
%         ratio = combined_band_power(cortexChannels(i),:) ./ combined_band_power(orgChannels(j),:);
%         pairings_awake((i-1)*length(orgChannels)+j,:) = 10*log10(ratio);
%     end
% end
% 
% pop_mean_a = mean(pairings_a,1);
% 
% % Plot mean computed both ways (point pairings and mean/mean)
% figure;
% scatter(reshape(repmat(1:6,length(cortexChannels)*length(orgChannels),1),1,[]),reshape(pairings,1,[])); hold on;
% scatter(1:6,pop_mean, 'filled'); hold on;
% scatter(1:6, 10*log10(meanCortex./meanOrg), 'filled');


% figure;
% bar(10*log10(meanCortex./meanOrg));
% ylabel('10*log10(Cor/Org)');
%bar([meanCortex; meanOrg].');


%% Compute Cortex vs. organoid statistics for awake and asleep (figure 4b)
%Awake
for i = 1:6
    [h,p] = ttest2(combined_band_power_awake(cortexChannels,i),  combined_band_power_awake(orgChannels,i));
    disp(['i = ' num2str(i) ' h = ' num2str(h) ' p = ' num2str(p)]);
end
figure;
A = mean(combined_band_power_awake(cortexChannels,:),1)';
B = mean(combined_band_power_awake(orgChannels,:),1)';
bar([A./A, B./A]); hold on;
er1 = errorbar((1:6)'-0.125, A./A, std(combined_band_power_awake(cortexChannels,:),[],1)./A', 'k');
er2 = errorbar((1:6)'+0.125, B./A, std(combined_band_power_awake(orgChannels,:),[],1)./A', 'k');
er1.LineStyle = 'none';
er2.LineStyle = 'none';
xlabel('Frequency Bands');
ylabel('Power (uV^2)');
legend('Cortex', 'Organoid');

%Anesthetized
for i = 1:6
    [h,p] = ttest2(combined_band_power(cortexChannels,i),  combined_band_power(orgChannels,i));
    disp(['i = ' num2str(i) ' h = ' num2str(h) ' p = ' num2str(p)]);
end
figure;
A = mean(combined_band_power(cortexChannels,:),1)';
B = mean(combined_band_power(orgChannels,:),1)';
bar([A./A, B./A]); hold on;
er1 = errorbar((1:6)'-0.125, A./A, std(combined_band_power(cortexChannels,:),[],1)./A', 'k');
er2 = errorbar((1:6)'+0.125, B./A, std(combined_band_power(orgChannels,:),[],1)./A', 'k');
er1.LineStyle = 'none';
er2.LineStyle = 'none';
xlabel('Frequency Bands');
ylabel('Power (uV^2)');
legend('Cortex', 'Organoid');

%% Compute Students two-sample t-test (figure 4b)
%Pair and compute ratios
pairings_awake = zeros(length(cortexChannels)*length(orgChannels),6);
for i = 1:length(cortexChannels)
    for j = 1:length(orgChannels)
        %disp((i-1)*length(orgChannels)+j);
        ratio = combined_band_power_awake(cortexChannels(i),:) ./ combined_band_power_awake(orgChannels(j),:);
        pairings_awake((i-1)*length(orgChannels)+j,:) = 10*log10(ratio);
    end
end

%Pair and compute ratios
pairings_anes = zeros(length(cortexChannels)*length(orgChannels),6);
for i = 1:length(cortexChannels)
    for j = 1:length(orgChannels)
        %disp((i-1)*length(orgChannels)+j);
        ratio = combined_band_power_anes(cortexChannels(i),:) ./ combined_band_power_anes(orgChannels(j),:);
        pairings_anes((i-1)*length(orgChannels)+j,:) = 10*log10(ratio);
    end
end

for i = 1:6
    [h,p] = ttest2(pairings_anes(:,i), pairings_awake(:,i));
    disp(['i = ' num2str(i) ' h = ' num2str(h) ' p = ' num2str(p)]);
end

stds = [std(pairings_anes,[],1)',std(pairings_awake,[],1)'];
pop_mean_anes = mean(pairings_anes,1);
pop_mean_awake = mean(pairings_awake,1);

figure;
bar([pop_mean_anes', pop_mean_awake'], ['g'; 'c']); hold on;
%er = errorbar(repmat((1:6)',1,2), [pop_mean', pop_mean_a'], stds);
er1 = errorbar((1:6)'-0.125, pop_mean_anes', stds(:,1), 'k');
er2 = errorbar((1:6)'+0.125, pop_mean_awake', stds(:,2), 'k');
er1.LineStyle = 'none';
er2.LineStyle = 'none';
xlabel('Frequency Bands');
ylabel('10*log10(Cor/Org)');
legend('Anesthesia', 'Awake');

%% Spike detection using findpeaks.m
time_of_interest = [1,60]; % seconds

toi_mua = time_of_interest * fs;
toi_lfp = time_of_interest * fs_down_lfp;
toi_mua = toi_mua(1):toi_mua(2);
toi_lfp = toi_lfp(1):toi_lfp(2);
if toi_mua(1) == 0
    toi_mua(1) = 1;
    toi_lfp(1) = 1;
end
data_mua_trunc = data_bpf(toi_mua,:);
t_spikes = t_amplifier(toi_mua);

thresh_mua = 4*std(data_mua_trunc,[],1); % modify SD multiplier based on SNR (3, 3.5, 4, or 4.5 usually good)

spikes = zeros(size(data_mua_trunc));
for i = 1:length(good_channels)
    [pks, loc] = findpeaks(-1*data_mua_trunc(:,i), 'MinPeakHeight', thresh_mua(i), 'MinPeakDistance', fs*0.001);
    spikes(loc,i) = 1;
end

%sum_spikes_cor = sum(spikes(:,cortexChannels),2);
%sum_spikes_org = sum(spikes(:,orgChannels),2);

% Plot lfp, mua, and spikes
figure;
for i = 1:length(good_channels)
    temp = (spikes(:,i) ~= 0);
    %plot(t_spikes, data_mua_trunc(:,i)+100*i, 'Color', colors(good_channels(i),:)); hold on;
    %plot(t_down_lfp(toi_lfp), data_lfp_down(toi_lfp,i)+100*i, 'k'); hold on;
    plot(get(gca,'xlim'), -1*[thresh_mua(i), thresh_mua(i)]+100*i); hold on;
    scatter(t_spikes(temp), spikes(temp,i)*(100*i)-20, 100, '.', 'MarkerEdgeColor', colors(good_channels(i),:)); hold on;
%     if i == 1 || i == 2
%         continue;
%     elseif find(i==cortexChannels)
%         scatter(t_spikes(temp), spikes(temp,i)+778, 100, '.', 'MarkerEdgeColor', colors(good_channels(i),:)); hold on;
%     else
%         scatter(t_spikes(temp), spikes(temp,i), 100, '.', 'MarkerEdgeColor', colors(good_channels(i),:)); hold on;
%     end
end

plot(t_spikes, mask_spikes*1700);
% smoothby = 0.01; %Smooth by __ seconds
% binned_cor = movsum(sum_spikes_cor, [0,smoothby*fs],1);
% binned_cor = binned_cor(1:smoothby*fs:end) / smoothby;
% plot(t_spikes(1:200:end), binned_cor); hold on;
% plot(t_spikes, smoothdata(sum_spikes_cor, 'gaussian', fs*0.1)*2000); hold on;
% plot(t_spikes, smoothdata(sum_spikes_org, 'gaussian', fs*0.1)*2000); hold on;

% xlabel('Time (seconds)');
% ylabel('Spike Rate (Spikes/sec) for population - light');
% title(nickname, 'Interpreter', 'none');

% binned_cor = movsum(spikes, [0,40],1);
% binned_cor = binned_cor(1:fs/fs_down_lfp:end,:);

mask_spikes = mask(toi_mua);
count_tot = sum(spikes,1);
count_burst = sum(spikes.*mask_spikes,1);
count_silent = sum(spikes.*(~mask_spikes),1);
percent = count_burst./count_tot * 100;
nullHypothesisMean = (sum(mask_spikes) / length(mask_spikes)) * 100;
figure;
scatter(ones(1,length(cortexChannels)), percent(cortexChannels)); hold on;
scatter(2*ones(1,length(orgChannels)), percent(orgChannels)); hold on;
errorbar(0.75, mean(percent(cortexChannels)), std(percent(cortexChannels)), '-bo', 'MarkerFaceColor','b');
errorbar(1.75, mean(percent(orgChannels)), std(percent(orgChannels)), '-ro', 'MarkerFaceColor','r');
xlim([0 3]); ylim([0 100]);
plot(get(gca, 'xlim'), [nullHypothesisMean, nullHypothesisMean], 'k--');
legend('Cortex', 'Organoid');
ylabel('% overlap with bursts');

[h,p] = ttest2(percent(cortexChannels), percent(orgChannels));



%% Rolling frequency
rollFreq = movsum(spikes, [0,fs], 1);
figure;
for i = 1:length(good_channels)
    plot(t_spikes, smoothdata(rollFreq(:,i), 'gaussian', fs*0.5), 'Color', colors(i,:)); hold on;
end

temp = mean(rollFreq(1:15*fs,:), 1);
temp = mean(rollFreq(45*fs:end,:), 1);

figure;
scatter(ones(1,length(cortexChannels)), temp(cortexChannels)); hold on;
scatter(2*ones(1,length(orgChannels)), temp(orgChannels)); hold on;
errorbar(0.75, mean(temp(cortexChannels)), std(temp(cortexChannels)), '-bo', 'MarkerFaceColor','b');
errorbar(1.75, mean(temp(orgChannels)), std(temp(orgChannels)), '-ro', 'MarkerFaceColor','r');
xlim([0.5 2.5]);

[h,p] = ttest2(temp(cortexChannels), temp(orgChannels));


%% Compute average MUA frequency during bursts
mask_spikes = mask(toi_mua);
kStart = strfind(mask_spikes', [0 1]);
kEnd = strfind(mask_spikes', [1 0]);
% Get rid of bursts cutoff by beginning or end of recording
if kStart(1) > kEnd(1)
    kEnd(1) = [];
end
if kStart(end) > kEnd(end)
    kStart(end) = [];
end

average_mua_frequencies = nan(2,length(good_channels)); % Top row is burst, bottom row is silent periods

for i = 1:length(good_channels)
    mua_frequencies = zeros(1,length(kStart));
    for j = 1:length(kStart)
        % Burst calculation
        len = (kEnd(j) - kStart(j)) / fs; %Time of burst in seconds
        nSpikes = sum(spikes(kStart(j):kEnd(j),i)); % Num spikes during burst
        if len < 0.1
            mua_frequencies(1,j) = nan;
        else
            mua_frequencies(1,j) = nSpikes/len; % Frequency in Hz
        end
        % Silent calculation
        if j > 1
            len = (kStart(j) - kEnd(j-1)) / fs; %Time of burst in seconds
            nSpikes = sum(spikes(kEnd(j-1):kStart(j),i)); % Num spikes during burst
            if len < 0.1
                mua_frequencies(2,j) = nan;
            else
                mua_frequencies(2,j) = nSpikes/len; % Frequency in Hz
            end
        end
    end
    average_mua_frequencies(:,i) = nanmean(mua_frequencies,2);
end

figure;
scatter(ones(1,length(cortexChannels)), average_mua_frequencies(1,cortexChannels), 'b'); hold on;
scatter(1.25*ones(1,length(cortexChannels)), average_mua_frequencies(2,cortexChannels), 'b'); hold on;
plot(repmat([1,1.25], length(cortexChannels),1)', [average_mua_frequencies(1,cortexChannels); average_mua_frequencies(2,cortexChannels)], 'b'); hold on;

scatter(2*ones(1,length(orgChannels)), average_mua_frequencies(1,orgChannels), 'r'); hold on;
scatter(2.25*ones(1,length(orgChannels)), average_mua_frequencies(2,orgChannels), 'r'); hold on;
plot(repmat([2,2.25], length(orgChannels),1)', [average_mua_frequencies(1,orgChannels); average_mua_frequencies(2,orgChannels)], 'r'); hold on;

errorbar(0.75, mean(average_mua_frequencies(1,cortexChannels)), std(average_mua_frequencies(1,cortexChannels)), '-bo', 'MarkerFaceColor','b');
errorbar(1.5, mean(average_mua_frequencies(2,cortexChannels)), std(average_mua_frequencies(2,cortexChannels)), '-bo', 'MarkerFaceColor','b');
errorbar(1.75, mean(average_mua_frequencies(1,orgChannels)), std(average_mua_frequencies(1,orgChannels)), '-ro', 'MarkerFaceColor','r');
errorbar(2.5, mean(average_mua_frequencies(2,orgChannels)), std(average_mua_frequencies(2,orgChannels)), '-ro', 'MarkerFaceColor','r');
xlim([0 3]); ylim([0 10]);
legend('Cortex', 'Organoid');
ylabel('Average MUA event rate');

[h,p] = ttest2(average_mua_frequencies(1,cortexChannels), average_mua_frequencies(1,orgChannels)); % Burst cor to org
[h,p] = ttest2(average_mua_frequencies(2,cortexChannels), average_mua_frequencies(2,orgChannels)); % Silent cor to org
[h,p] = ttest2(average_mua_frequencies(1,cortexChannels), average_mua_frequencies(2,cortexChannels)); % Burst cor to silent cor
[h,p] = ttest2(average_mua_frequencies(1,orgChannels), average_mua_frequencies(2,orgChannels)); % Burst org to silent org

%% ISI
figure;
for i = 1:length(good_channels)
    isi_array = find(spikes(:,i));% .* (~mask_spikes));
    isi_array = isi_array(2:end) - isi_array(1:end-1);
    isi_array = isi_array ./ fs;
    isi_array(isi_array > 1) = [];
    subplot(4,4,find(mapping == good_channels(i)));
    histogram(isi_array, 'BinWidth', 0.01, 'Normalization', 'probability');
    title(['nspikes = ' num2str(length(isi_array) + 1)]);
    xlim([-0.1 1]);
    ylim([0 0.4]);
end


%% Compare awake and asleep spikes
figure;
scatter(ones(1,length(cortexChannels)), average_mua_frequencies(1,cortexChannels), 'b'); hold on;
scatter(1.25*ones(1,length(cortexChannels)), temp(cortexChannels), 'b'); hold on;
plot(repmat([1,1.25], length(cortexChannels),1)', [average_mua_frequencies(1,cortexChannels); temp(cortexChannels)], 'b'); hold on;

scatter(2*ones(1,length(orgChannels)), average_mua_frequencies(2,orgChannels), 'r'); hold on;
scatter(2.25*ones(1,length(orgChannels)), temp(orgChannels), 'r'); hold on;
plot(repmat([2,2.25], length(orgChannels),1)', [average_mua_frequencies(2,orgChannels); temp(orgChannels)], 'r'); hold on;

errorbar(0.75, mean(average_mua_frequencies(2,cortexChannels)), std(average_mua_frequencies(2,cortexChannels)), '-bo', 'MarkerFaceColor','b');
errorbar(1.5, mean(temp(cortexChannels)), std(temp(cortexChannels)), '-bo', 'MarkerFaceColor','b');
errorbar(1.75, mean(average_mua_frequencies(2,orgChannels)), std(average_mua_frequencies(2,orgChannels)), '-ro', 'MarkerFaceColor','r');
errorbar(2.5, mean(temp(orgChannels)), std(temp(orgChannels)), '-ro', 'MarkerFaceColor','r');
xlim([0 3]); ylim([0 15]);
legend('Cortex', 'Organoid');
ylabel('Average MUA event rate');

[h,p] = ttest2(average_mua_frequencies(1,cortexChannels), temp(cortexChannels)); % Burst cor to awake cor
[h,p] = ttest2(average_mua_frequencies(1,orgChannels), temp(orgChannels)); % Burst org to awake org
[h,p] = ttest2(average_mua_frequencies(1,cortexChannels), average_mua_frequencies(2,cortexChannels)); % Burst cor to silent cor
[h,p] = ttest2(average_mua_frequencies(1,orgChannels), average_mua_frequencies(2,orgChannels)); % Burst org to silent org

%% Single Channel Morlet frequency analysis
% User parameters
channel = 13;

freq_range = [0,150]; % frequency range of interest
deltaF = 0.5; % frequency resolution [Hz]
morletParams = [20,40]; % Range of cycles (low number of cycles => better temporal resolution, high # cycles => better freq resolution)

ch = find(good_channels == channel);

% Calculate frequency range
num_frex = (freq_range(2)-freq_range(1)+1)/deltaF;
rangeF = linspace(freq_range(1), freq_range(2), num_frex);

figure;

% Morlet-based spectrogram
M = my_morlet(data_lfp_down(time_range,ch), fs_down_lfp, freq_range(1), freq_range(2), deltaF, morletParams);

%subtract baseline
%baseline = median(holdingM(:,1:(leadTime*fs_down*-1)-1), 2);
baseline = median(M, 2);
M = 10*(log10(M) - log10(baseline));

h = pcolor(t_down_lfp(time_range(1:10:end)),rangeF, M(:,1:10:end));
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
    
title([nickname ' Params: [' num2str(morletParams(1)) ',' num2str(morletParams(2)) ']' ], 'Interpreter', 'none');


%% Coherency
% User params
time_range = [0 60] *fs_down_lfp;

params.tapers = [10 19];
params.pad = 0;
params.Fs = fs_down_lfp;
params.fpass = [0 150];
params.err = [1 0.05];
params.trialave = 0;
colors = ['y','m','c','r','g','b','k'];

if time_range(1) == 0
    time_range(1) = 1;
end
time_range = time_range(1):time_range(2);

num = length(good_channels);

% Instantiate heat maps
delta_heat = zeros(num,num);
ab_heat = zeros(num,num);
lowgamma_heat = zeros(num,num);
highgamma_heat = zeros(num,num);

% Instantiate frequency bands
N = length(time_range);
nfft=max(2^(nextpow2(N)+params.pad),N);
[f,~]=getfgrid(params.Fs,nfft,params.fpass);

delta_band = find(f>=2 & f<=4);
ab_band = find(f>=7 & f<=30);
lowgamma_band = find(f>=31 & f<=59);
highgamma_band = find(f>=61 & f<=100);

figure;
for i = 1:num
    for j = 1:i
        [C,phi,S12,S1,S2,f,confC,phistd] = coherencyc(data_lfp_down(time_range,j), data_lfp_down(time_range,i), params);
        
        delta_heat(i,j) = mean(C(delta_band));
        ab_heat(i,j) = mean(C(ab_band));
        lowgamma_heat(i,j) = mean(C(lowgamma_band));
        highgamma_heat(i,j) = mean(C(highgamma_band));
        
        subplot(num, num, (i-1)*num + j);
        %plot_vector(C,f,'n',confC,colors(mod(i,length(colors))+1));
        plot_vector(C,f,'n',confC,'k');
        set(gca, 'YLabel', []);
        set(gca, 'XLabel', []);
        title(['i = ' num2str(i) ' j = ' num2str(j)]);
        %plot_vector(smoothdata(C,'g',5),f,'n',[],colors(mod(i,length(colors))+1));
    end
    disp(['i = ' num2str(i)]);
end

%% Heat map of coherency at different frequency bands
figure;
temp = zeros(16,16);
temp(good_channels,good_channels) = delta_heat;
subplot(1,4,1);
imagesc(temp);
title('Delta (2-4 Hz)');

temp = zeros(16,16);
temp(good_channels,good_channels) = ab_heat;
subplot(1,4,2);
imagesc(temp);
title('Alpha/Beta (7-30 Hz)');

temp = zeros(16,16);
temp(good_channels,good_channels) = lowgamma_heat;
subplot(1,4,3);
imagesc(temp);
title('Low Gamma (31-59 Hz)');

temp = zeros(16,16);
temp(good_channels,good_channels) = highgamma_heat;
subplot(1,4,4);
imagesc(temp);
title('High Gamma (61-100 Hz)');
c = colorbar;
c.Label.String = 'Coherency';

sgtitle(nickname, 'Interpreter', 'none');



%% Magnitude Squared Coherence
% Test
fs_test = 200;
freq = 2*pi*20;
t = 0:1/fs_test:10;
x = cos(freq*t) + rand(1,length(t));
y = cos(freq*t);

figure; hold on;
plot(x);
plot(y);
[cxy,f] = mscohere(x,y,[],[],[],fs_test);

figure; plot(f,cxy);

% 












