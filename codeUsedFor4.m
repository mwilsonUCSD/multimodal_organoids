% Author: Madison Wilson
clear;
clc;

nickname = 'NOD18\0302_run3';

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

if isfield(rhd, 'board_adc_data')
    disp("ADC found");
    ind = find(diff(rhd.board_adc_data(1,:) > 0.3) == 1) + 1;
    t_stim = t_amplifier(ind);
end


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

figure;
for i = 1:size(data_theta,2)
    plot(t_amplifier, data_theta(:,i)+500*i); hold on;
end
if isfield(rhd, 'board_adc_data')
    for i = 1:length(t_stim)
        plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
    end
end
title(['[' num2str(f_bpass) '] Hz filtered signal for ' nickname ' ICA: [' num2str(ica_remove) ']'], 'Interpreter', 'none');
ylabel('Voltage (\muV)');
xlabel('Time (s)');


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


%% Morlet frequency analysis (figure 4d)
% User parameters
do_mask = 0; % 0 = no mask, 1 = bursts, 2 = silent periods
time_range = [26,28]; % seconds
freq_range = [1,150];%[10,32]; % frequency range of interest
deltaF = 0.5;%0.05; % frequency resolution [Hz]
morletParams = [4,40];%[5,20]; % Range of cycles (low number of cycles => better temporal resolution, high # cycles => better freq resolution)
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
    M = 10*(log10(M) - log10(baseline)); % 10*log(S/B) for power and 20*log(S/B) for amplitude
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
    
    if do_mask
        if do_mask == 1
            M(:,~mask_trunc) = 0;
        else
            M(:,mask_trunc) = 0;
        end
        plot(rangeF, sum(M,2)./sum(mask_trunc), 'Color', colors(good_channels(i),:)); hold on;
        M_tot(i,:) = sum(M,2)./sum(mask_trunc);
    else
        plot(rangeF, mean(M,2), 'Color', colors(good_channels(i),:)); hold on;
        M_tot(i,:) = mean(M,2);
    end
    
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
lC = length(cortexChannels);
lO = length(orgChannels);
scatter(repelem(1:6,lC)-0.125, reshape(combined_band_power_awake(cortexChannels,:)./repmat(A',lC,1), [1, lC*6]), 'ok', 'filled');
scatter(repelem(1:6,lO)+0.125, reshape(combined_band_power_awake(orgChannels,:)./repmat(A',lO,1), [1, lO*6]), 'ok', 'filled');
xlabel('Frequency Bands');
ylabel('Power (uV^2)');
legend('Cortex', 'Organoid');
title('Awake');

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
scatter(repelem(1:6,lC)-0.125, reshape(combined_band_power(cortexChannels,:)./repmat(A',lC,1), [1, lC*6]), 'ok', 'filled');
scatter(repelem(1:6,lO)+0.125, reshape(combined_band_power(orgChannels,:)./repmat(A',lO,1), [1, lO*6]), 'ok', 'filled');
xlabel('Frequency Bands');
ylabel('Power (uV^2)');
legend('Cortex', 'Organoid');
title('Anesthetized');




