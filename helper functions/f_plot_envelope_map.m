function f_plot_envelope_map(data, f_bpass, thresh_std, dt, t, step, fs, t_amplifier, nickname, good_channels, mapping, ica_remove, doSave, t_omit)
    
%     f_bpass = [5,8];
    fs_down = f_bpass(2) * 2;
    
    bpass = designfilt('bandpassiir','FilterOrder',6, ...
    'HalfPowerFrequency1',f_bpass(1),'HalfPowerFrequency2',f_bpass(2), ...
    'SampleRate',fs);

    data_filt = filtfilt(bpass, data);
    data_filt_down = data_filt(1:fs/fs_down:end,:);
    t_down = t_amplifier(1:fs/fs_down:end);
    data_env = envelope(data_filt);
    
    % Get rid of motion epochs
    t_omit_epoch = t_omit.*fs;
    if size(t_omit) > 0
        if t_omit_epoch(1,1) == 0
            t_omit_epoch(1,1) = 1;
        end
    end
    for i = 1:size(t_omit,1)
        epoch = t_omit_epoch(i,1):t_omit_epoch(i,2);
        data_filt(epoch,:) = nan;
        data_env(epoch,:) = nan;
    end
    
    % Compute threshold crossings
    thresh = thresh_std*nanstd(data_env); % modify SD multiplier based on SNR (3, 3.5, 4, or 4.5 usually good)
    
    burst_onset = zeros(size(data_env));
    for i = 1:length(good_channels)
        [pks, loc] = findpeaks(data_env(:,i), 'MinPeakHeight', thresh(i), 'MinPeakDistance', fs);
        burst_onset(loc,i) = 1;
    end
    for i = 1:size(t_omit,1)
        if (t_omit(i,1)+dt(1))*fs < 0
            burst_onset(1:(t_omit(i,2)+dt(2))*fs,:) = 0;
        else
            burst_onset((t_omit(i,1)+dt(1))*fs:(t_omit(i,2)+dt(2))*fs,:) = 0; % Remove peaks around omitted epochs
        end
        
    end
    burst_onset(burst_onset == 0) = nan;
    
    % Plot filtered data and threshold crossings
    f1 = figure;
    set(gcf, 'Position', [6,42,1148,954]);
    hold on;
    for i = 1:size(data,2)
        plot(t_amplifier, data_filt(:,i)+400*i, 'LineWidth', 2);
        plot(t_amplifier, data_env(:,i)+400*i, 'LineWidth', 2);
        scatter(t_amplifier, data_env(:,i).*burst_onset(:,i)+400*i, 'filled'); % Plot threshold crossings
    end
    
    title(['[' num2str(f_bpass) '] Hz filtered signal for ' nickname ' ICA: [' num2str(ica_remove) ']'], 'Interpreter', 'none');
    ylabel('Voltage (\muV)');
    xlabel('Time (s)');
    
    % Compute average Morlet frequency response following each channel's peak
    %dt = [-1, 1]; % seconds
    freq_range = [0,100]; % frequency range of interest
    deltaF = 0.5; % frequency resolution [Hz]
    morletParams = [10,30]; % Range of cycles (low number of cycles => better temporal resolution, high # cycles => better freq resolution)
    
    % Calculate frequency range
    %num_frex = (freq_range(2)-freq_range(1)+1)/deltaF;
    %rangeF = linspace(freq_range(1), freq_range(2), num_frex);

%     time_range = time_range*fs_down_lfp;
%     if time_range(1) == 0
%         time_range(1) = 1;
%     end
%     time_range = time_range(1):time_range(2);

    time_buffer = dt(1)*fs:dt(2)*fs;
    %M = zeros(length(good_channels), length(good_channels), length(rangeF), length(time_buffer));
    env_averages = zeros(length(good_channels), length(good_channels), length(time_buffer));
    
    for i = 1:length(good_channels)
        idx = find(burst_onset(:,i) == 1);
        idx(idx <  fs*dt(1)*-1) = []; % Get rid of indices too close to beginning or end of recording
        idx(idx > length(t_amplifier) - fs*dt(2)) = [];
%         for i2 = 1:size(t_omit,1)
%             idx(bitand(idx > (t_omit(i2,1)+dt(1))*fs,idx < (t_omit(i2,2)+dt(2))*fs)) = []; % Remove peaks too close to removed epochs
%         end
        
        
        for j = 1:length(idx)
            time_range = idx(j) + time_buffer;
            
            env_averages(i,:,:) = squeeze(env_averages(i,:,:)) + data_env(time_range,:)';
            
            % Morlet-based spectrograms
%             for k = 1:length(good_channels)
%                 M(i,k,:,:) = squeeze(M(i,k,:,:)) + my_morlet(data_filt_down(time_range,k), fs_down, freq_range(1), freq_range(2), deltaF, morletParams);
%             end
            
        end
        %M(i,:,:,:) = M(i,:,:,:) ./ length(idx);
        env_averages(i,:,:) = env_averages(i,:,:) ./ length(idx);

        %subtract baseline
%         baseline = mean(M, 3);
%         M = 10*(log10(M) - log10(baseline));
        
        disp(i);
    end
    
    % Plot average envelopes
    f2 = figure;
    f3 = figure;
    amp_max_max = max(env_averages, [], 'all');
    amp_max_max = ceil(amp_max_max/10)*10;
    
    for i = 1:length(good_channels)
%         figure(f1); hold on;
%         for j = 1:length(good_channels)
%             subplot(4,4,find(mapping == good_channels(j)));
%             plot(time_buffer/fs, squeeze(env_averages(i,j,:)),'LineWidth',2); hold on;
%             scatter(delay_max(j),amp_max(j), 'filled');
%             ylim([0 200]);
%         end
%         sgtitle([nickname ' Reference Ch: ' num2str(good_channels(i))], 'Interpreter', 'none');
        
        
        [amp_max, delay_max] = max(squeeze(env_averages(i,:,:)),[],2);
        delay_max = time_buffer(delay_max) ./ fs;
        
        amp_plot = nan(1,16);
        amp_plot(good_channels) = amp_max;
        amp_plot = amp_plot(mapping);
        amp_plot = reshape(amp_plot, [4,4])';
        %amp_plot(isnan(amp_plot)) = interp2(find(~isnan(amp_plot)), amp_plot(~isnan(amp_plot)), find(isnan(amp_plot)),'cubic');
        
        delay_plot = nan(1,16);
        delay_plot(good_channels) = delay_max;
        delay_plot = delay_plot(mapping);
        delay_plot = reshape(delay_plot, [4,4])';
        
        figure(f2);
        subplot(4,4,find(mapping == good_channels(i)));
        imagesc(amp_plot);
        title(['Ref ch = ' num2str(good_channels(i))]);
        colorbar;
        caxis([0, amp_max_max]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        
        figure(f3);
        subplot(4,4,find(mapping == good_channels(i)));
        imagesc(delay_plot);
        title(['Ref ch = ' num2str(good_channels(i))]);
        colorbar;
        caxis([-0.04,0.04]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
    end
    figure(f3);
    sgtitle(['Delay time heat maps freq = [' num2str(f_bpass) ']']);
    set(gcf, 'Position', [3,42,1383,954]);
    figure(f2);
    sgtitle(['Amplitude heat maps freq = [' num2str(f_bpass) ']']);
    set(gcf, 'Position', [3,42,1383,954]);
    
    % Save figures
    if doSave == 1
        saveas(f2, ['_210614_burst_analysis\' nickname '_amp_' num2str(f_bpass(1)) 'to' num2str(f_bpass(2)) 'Hz'], 'jpeg');
        saveas(f3, ['_210614_burst_analysis\' nickname '_delay_' num2str(f_bpass(1)) 'to' num2str(f_bpass(2)) 'Hz'], 'jpeg');
        save(['_210614_burst_analysis\data\' nickname '_envAvg_' num2str(f_bpass(1)) 'to' num2str(f_bpass(2)) 'Hz'], 'env_averages', 'good_channels', 'f_bpass');
    end
    
    
    
    
%     
%     %Plot
%     figure;
%     for i = 1:length(good_channels)
%         subplot(4,4,find(mapping == good_channels(i)));
%         %h = pcolor(t_down_lfp(time_range(1:10:end)),rangeF, M(:,1:10:end));
%         h = pcolor(time_buffer./fs_down,rangeF, squeeze(M(i,1,:,:)));
%         set(h,'EdgeColor','none');
%         %set(gca, 'YScale', 'log');
%         %set(gca, 'FontSize', 25);
%         c = colorbar;
%         c.Label.String = '10*log10(signal/baseline)';
%         colormap jet;
%         caxis([0.5,4]);
%         %set(c, 'visible', 'on');
%         title(['Channel: ' num2str(good_channels(i))]);
%         xlabel('Time (seconds)');
%         ylabel('Frequency (Hz)');
%     end
%     sgtitle([nickname ' Params: [' num2str(morletParams(1)) ',' num2str(morletParams(2)) '] ICA: [' num2str(ica_remove) ']'], 'Interpreter', 'none');
% 
%     %figure; plot(t_down, data_env(:,i)); hold on; plot(get(gca, 'xlim'), [thresh(i),thresh(i)]);
%     figure;
%     hold on;
%     for i = 1:size(burst_onset,2)
%         scatter(t_down, burst_onset(:,i)*200*i, 'filled');
%     end
%     
%     % Make array of min ripple amp latencies based on each spike
%     [idxSpike,chSpike] = find(burst_onset == 1);
%     latencies = zeros(length(idxSpike),length(good_channels));
%     buffer = 0.0625; % seconds
%     for i = 1:length(idxSpike)
%         [~,latencies(i,:)] = min(data_filt_down(idxSpike(i)+(-fs_down*buffer:fs_down*buffer),:));
%     end
%     latencies = latencies - fs_down*buffer - 1;
%     
%     %Average latencies
%     latencies_avg = zeros(length(good_channels),length(good_channels));
%     for i = 1:length(idxSpike)
%         latencies_avg(chSpike(i),:) = latencies_avg(chSpike(i),:) + latencies(i,:) - latencies(i,chSpike(i));
%     end
%     for i = 1:length(good_channels)
%         latencies_avg(i,:) = latencies_avg(i,:) / length(find(chSpike == i));
%     end
%     
%     %Reorder to match physical mapping
%     latency_map = zeros(16,16);
%     latency_map(good_channels, good_channels) = latencies_avg / fs_down; %divide to get units of seconds rather than indices
%     latency_map = latency_map(:,mapping); % remap
%     figure;
%     for i = 1:length(good_channels)
%         subplot(4,4,find(mapping == good_channels(i)));
%         imagesc(reshape(latency_map(good_channels(i),:),[4,4])');
%         colorbar;
%         caxis([-buffer, buffer]);
%     end
%     
%     
%     %Find indices of times of interest
%     idx = zeros(length(t(1):step:t(2)),1);
%     count = 1;
%     for i = t(1):step:t(2)
%         idx(count) = find(t_down >= i,1);
%         count = count+1;
%     end
%     
%     % Make matrix of power values and sort according to physical channel mapping
%     data_plot = zeros(length(idx),16);
%     data_plot(:,good_channels) = data_env(idx, :); %fit to 16 channel array
%     data_plot = data_plot(:,mapping); % Remap according to physical mapping
%     %data_plot = reshape(data_plot,[4,4,:]);
%     
%     f2 = figure;
%     set(gcf, 'Position', [1200,388,560,420]);
%     for i = 1:length(idx)
%         %subplot(5,5,i);
%         figure(f1);
%         if exist('p1')
%             delete(p1);
%         end
%         p1 = plot([t_down(idx(i)),t_down(idx(i))], get(gca,'ylim'), 'k');
%         
%         figure(f2);
%         imagesc(reshape(data_plot(i,:),[4,4])');
%         colorbar;
%         caxis([0,30]);
%         title(['Time = ' num2str(t_down(idx(i)))]);
%         pause(0.1);
%     end
    
end




