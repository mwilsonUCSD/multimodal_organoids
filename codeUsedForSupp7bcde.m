%% Analyze correlation in MUA events between channels
% Author: Madison Wilson
% Algorithm proposed by Xin Liu
clear;
clc;

nickname = 'NOD18\0302_run1';
load(['D:\Maddie\Martin Exps\Processed Data\' nickname '.mat']);


% Extract variables
fs = rhd.frequency_parameters.amplifier_sample_rate;
impedance = [rhd.amplifier_channels.electrode_impedance_magnitude];
t_amplifier = rhd.t_amplifier;

mapping = [4,3,1,16,5,6,2,15,7,10,14,13,8,9,11,12];

impedance_mapped = impedance(9:24); %cut out unused channels
good_channels = find(impedance_mapped < 10e6); % consider channels below an impedance threshold "good channels"
impedance_mapped = impedance_mapped(good_channels);

data_good = rhd.amplifier_data(9:24,:); % cut out unused channels
data_good = data_good(good_channels, :); % reorder channels to match array mapping

%ind = find(diff(rhd.board_adc_data(1,:) > 0.3) == 1) + 1;
%t_stim = t_amplifier(ind);


% plot raw data
% figure;
% for i = 1:size(data_good,1)
%     plot(t_amplifier, data_good(i,1:end)+500*i); hold on;
%     %plot(t_amplifier(t_range(1:end)), data_good(i,1:end)+100*i); hold on;
% end
% title([nickname ' Raw Data'], 'Interpreter', 'none');
% ylabel('Voltage (\muV)');
% xlabel('Time (s)');


% Filter data below fs_down Hz (easier on ICA)
fs_down = 10000;

d_low_3000 = designfilt('lowpassiir','FilterOrder',8, ...
'PassbandFrequency',3000,'PassbandRipple',0.2, ...
'SampleRate',fs);

data_clean_filt = filtfilt(d_low_3000, data_good');

data_down = data_clean_filt(1:fs/fs_down:end,:); % fs_down Hz sampling rate
t_down = t_amplifier(1:fs/fs_down:end);


% ICA to remove noise
%[Aica,ica_source_activity] = ica_denoise(data_good', t_amplifier, 1:length(good_channels));
%f_plot_ica(Aica,ica_source_activity, good_channels, mapping, t_amplifier);
% for i = 1:length(t_stim)
%     plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1);
% end


% Remove noisy ICA components
good_ind = setdiff(1:length(good_channels), ica_remove); % second input is channel to remove based on imaging artifact
ECoG_denoise = ica_source_activity(:,good_ind) * Aica(:,good_ind)';

% Plot data
% figure;
% for i = 1:size(data_good,1)
%     plot(t_amplifier, ECoG_denoise(:,i)+500*i); hold on;
%     %plot(t_amplifier(t_range(1:end)), data_good(i,1:end)+100*i); hold on;
% end
% title([nickname ' data post-ICA'], 'Interpreter', 'none');
% ylabel('Voltage (\muV)');
% xlabel('Time (s)');


% Compute MUA (500-3000 Hz) and MUA events (passing 5*SD threshold)
d_bpass = designfilt('bandpassiir','FilterOrder',6, ...
'HalfPowerFrequency1',500,'HalfPowerFrequency2',3000, ...
'SampleRate',fs);

% MUA is filtered, full-wave rectified, and LPF for smoothing
data_band_filt = filtfilt(d_bpass, ECoG_denoise);
% data_rect = data_band_filt.^2;
% data_mua = filtfilt(d_low_100, data_rect);

% Plot MUA
figure;
for i = 1:size(data_band_filt,2)
    plot(t_amplifier, data_band_filt(:,i)+500*i); hold on;
end
% for i = 1:length(t_stim)
%     plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
% end
ylim([0 (i+1)*500]);
title('MUA');


% Remove noise if any
%prompt = 'Enter tStart = ';
tStart = 1; %input(prompt); %1; % seconds

%prompt2 = 'Enter tEnd = ';
tEnd = 100; %input(prompt2); %25; % seconds

data_band_filt = data_band_filt(tStart*fs_down:tEnd*fs_down,:);
t_down = t_down(tStart*fs_down:tEnd*fs_down);


% Spike detection thresholds
thresh = 4.5*std(data_band_filt); % modify SD multiplier based on SNR (3, 3.5, 4, or 4.5 usually good)

spikes = zeros(size(data_band_filt));
for i = 1:length(good_channels)
    [pks, loc] = findpeaks(-1*data_band_filt(:,i), 'MinPeakHeight', thresh(i), 'MinPeakDistance', fs*0.001);
    spikes(loc,i) = 1;
end

% Bin spikes and sum all channels
%spikes_sum = sum(spikes,2);

% Plot spike raster and spiking vector
figure;
for i = 1:length(good_channels)
    scatter(t_down, spikes(:,i)*(i+1), '.'); hold on;
end
%plot(t_down, spikes_sum);
xlim([tStart tStart+30]);
xlabel('Time (seconds)');
ylabel('Spike Rate (Spikes/sec) for population - light');
title(nickname, 'Interpreter', 'none');


%% Plot grid of spikes based on each channel's spike times (Supplementary figure 7b)
%channel = 1;
miniTimeStart = -0.001;
miniTimeEnd = 0.002;
time = miniTimeStart*fs_down:1:miniTimeEnd*fs_down;
t_reverse = fliplr(time);
t2 = [time, t_reverse];
c = jet(length(good_channels));

% Plot all event-triggered MUA averages for all channels
%(Supplementary figure 7b)
for k = 1:length(good_channels)
    figure;
    spike_times = find(spikes(:,k) == 1);
    for i = 1:length(good_channels)
        subplot(4,4,find(mapping == good_channels(i)));
        spike_mean = zeros(length(time),1);
        spike_sd = zeros(length(spike_times), length(time));
        for j = 1:length(spike_times)
            %spike_mean = spike_mean + data_band_filt(spike_times(j)+time,i);
            spike_sd(j,:) = data_band_filt(spike_times(j)+time,i);
            %plot(time./fs_down, data_band_filt(spike_times(j)+time,i)); hold on;
        end
        spike_mean = mean(spike_sd);%/length(spike_times);
        spike_sd = std(spike_sd,0,1);
        curve1 = spike_mean + spike_sd;
        curve2 = spike_mean - spike_sd;
        btwn = [curve1, fliplr(curve2)];
%         if find(not_working == i)
%             fill(t2/fs_down*1000, btwn, [1, 0.75, 0.75], 'EdgeColor', 'none'); hold on;
%             plot(time/fs_down*1000, spike_mean, 'LineWidth', 2, 'Color', 'r');
%         else
            fill(t2/fs_down*1000, btwn, [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on;
            plot(time/fs_down*1000, spike_mean, 'LineWidth', 2, 'Color', 'k');
%         end
        
        %plot(time/fs_down*1000, spike_mean, 'Color', c(k,:)); hold on;
        ylim([-30 30]);
        %ylabel('Amplitude (\muV)');
        %xlabel('Time (ms)');
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
    end
    sgtitle(['Ch ' num2str(k) 'nSpikes: ' num2str(length(spike_times))], 'Interpreter', 'none');
    %saveas(gcf, [nickname '_Ch' num2str(good_channels(k))], 'eps');
    %close;%MNW
end
%sgtitle([nickname ' Avg MUA from ' num2str(tStart) ' to ' num2str(tEnd)], 'Interpreter', 'none');
%set(gcf, 'Position', [50 50 1200 600]);

%saveas(gcf, [nickname '_45sd'], 'png');

%% Plot selected channel spike waveforms (Supplementary figure 7c)
channels = [4,6,14,12];

for k = 1:length(channels)
    chK = channels(k);
    figure;
    spike_times = find(spikes(:,chK) == 1);
    for i = 1:length(good_channels)
        subplot(4,4,find(mapping == good_channels(i)));
        spike_mean = zeros(length(time),1);
        spike_sd = zeros(length(spike_times), length(time));
        for j = 1:length(spike_times)
            %spike_mean = spike_mean + data_band_filt(spike_times(j)+time,i);
            spike_sd(j,:) = data_band_filt(spike_times(j)+time,i);
            %plot(time./fs_down, data_band_filt(spike_times(j)+time,i)); hold on;
        end
        spike_mean = mean(spike_sd);%/length(spike_times);
        %spike_sd = std(spike_sd,0,1);
        %curve1 = spike_mean + spike_sd;
        %curve2 = spike_mean - spike_sd;
        %btwn = [curve1, fliplr(curve2)];
%         if find(not_working == i)
%             fill(t2/fs_down*1000, btwn, [1, 0.75, 0.75], 'EdgeColor', 'none'); hold on;
%             plot(time/fs_down*1000, spike_mean, 'LineWidth', 2, 'Color', 'r');
%         else
            %fill(t2/fs_down*1000, btwn, [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on;
            %plot(time/fs_down*1000, spike_mean, 'LineWidth', 2, 'Color', 'k');
%         end
        plot(time/fs_down*1000, spike_sd, 'LineWidth', 0.5); hold on;
        plot(time/fs_down*1000, spike_mean, 'LineWidth', 2, 'Color', 'k');
        
        %plot(time/fs_down*1000, spike_mean, 'Color', c(k,:)); hold on;
        ylim([-30 30]);
        %ylabel('Amplitude (\muV)');
        %xlabel('Time (ms)');
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
    end
    sgtitle(['Ch ' num2str(chK) 'nSpikes: ' num2str(length(spike_times))], 'Interpreter', 'none');
    set(gcf, 'renderer', 'painters');
    %saveas(gcf, [nickname '_Ch' num2str(good_channels(chK)) '_nonaveraged'], 'eps');
    %close;%MNW
end

%% Count number of overlapping events for each channel pairing
overlap_count = zeros(size(spikes,2));
for i = 1:size(spikes,2)
    for j = 1:size(spikes,2)
        overlap_count(i,j) = sum((spikes(:,i) + spikes(:,j)) > 1);
    end
end


% Circularly shift one train and compute mean and sd of spike overlaps
nShifts = 100; % Number of times user wants to shift the data
overlap_mean = zeros(size(spikes,2));
overlap_sd = zeros(size(spikes,2));
overlap_temp_counts = zeros(size(spikes,2),size(spikes,2),nShifts); %zeros(1,nShifts);
n = size(spikes,1);

for i = 1:size(spikes,2)
    %for j = 1:i%size(spikes,2)
        tic;
        for k = 1:nShifts
            temp = circshift(spikes(:,i),randi(n)); %spikes([k+1:end 1:k],i); <- SLOW (~3x slower)%
            overlap_temp_counts(i,:,k) = sum((temp + spikes) > 1,1); %overlap_temp_counts(k) = sum((temp + spikes(:,j)) > 1);
        end
        toc;
        %pd = fitdist(overlap_temp_counts', 'Normal');
        overlap_mean(i,:) = mean(overlap_temp_counts(i,:,:),3); %overlap_mean(i,j) = mean(overlap_temp_counts);
        overlap_sd(i,:) = std(overlap_temp_counts(i,:,:),0,3); %overlap_sd(i,j) = std(overlap_temp_counts);
    %end
    disp(['i = ' num2str(i)]);
end
disp('Done');


% Calculate pval
%y = @(x)exp(-0.5 * ((x - overlap_mean(1,1))./overlap_sd(1,1)).^2) ./ (sqrt(2*pi) .* overlap_sd(1,1));
%figure;
%scatter(overlap_count(1,1),0);

pval = zeros(size(spikes,2));
%pval = integral(y, overlap_mean(1,1)-1*overlap_sd(1,1), overlap_mean(1,1)+overlap_sd(1,1));
for i = 1:size(spikes,2)
    for j = 1:size(spikes,2)
        %y = @(x)normpdf(x, overlap_mean(i,j), overlap_sd(i,j));
        y = @(x)poisspdf(x, overlap_sd(i,j));
        pval(i,j) = integral(y, overlap_count(i,j), Inf);
    end
end

mapped_overlap = zeros(17);
mapped_overlap(good_channels, good_channels) = pval;%overlap_count;% 
figure; pcolor(mapped_overlap');
colorbar;
caxis([0 0.05]);
title(nickname, 'Interpreter', 'none');
%saveas(gcf, [nickname '_pval_10k_45sd_hist'], 'png');

%save('nod13var.mat', 'overlap_mean', 'overlap_sd', 'overlap_count', 'pval', 'overlap_temp_counts', 'good_channels');

%% Plot PDFs
figure;
for i = 1:size(spikes,2)
    for j = 1:i
        mi = good_channels(i);
        mj = good_channels(j);
        subplot(16,16,(mi-1)*16+mj);
        %x = -1:0.01:3;
        %y = @(x)normpdf(x, overlap_mean(i,j), overlap_sd(i,j));
        x = -1:3;
        y = @(x)poisspdf(x, overlap_mean(i,j));
        plot(x,y(x)); hold on;
        if pval(i,j) < 0.05
            scatter(overlap_count(i,j), 0, '*');
        else
            scatter(overlap_count(i,j), 0);
        end
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        ylim([0 2]);
        xlim([-1 3]);
    end
end
%saveas(gcf, [nickname '_pdf_10k_45sd'], 'png');


%% Histogram pvalues
pval = zeros(length(good_channels));
for i = 1:length(good_channels)
    for j = 1:length(good_channels)
        h = histogram(overlap_temp_counts(i,j,:));
        pval(i,j) = sum(h.BinCounts(overlap_count(i,j)+2:end))/nShifts;
    end
end


%% Plot select channel histograms org/cortex (Supplementary figure 7e)
oChannels = [5,6,10];
cChannels = [9,11,12];

figure;
for i = 1:length(cChannels)
    for j = 1:length(oChannels)
        mi = find(good_channels == oChannels(j));
        mj = find(good_channels == cChannels(i));
        subplot(length(cChannels),length(oChannels),(i-1)*length(oChannels)+j);
        h = histogram(overlap_temp_counts(mi,mj,:));
        set(gca, 'xtick', []);
        set(gca, 'ytick', []);
        ylim([0 500]);
        xlim([-1 7]);
        %set(gca, 'yscale', 'log');
        disp([num2str((i-1)*length(oChannels)+j)]);
    end
end
%saveas(gcf, ['NOD18_' num2str(oChannels) ' by ' num2str(cChannels)], 'png');



