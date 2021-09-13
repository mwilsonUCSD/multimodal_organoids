function [trialAverages, t_trial] = f_plotSDTrialAverages2D(smooth, y, light, data, t, fs, stimuli, leadTime, lenTrial, badEpochs, mapping, good_channels, ylab)
% Function to plot the trial average and +-SD for each channel

% Inputs
    % smooth = 0 if no smoothing, else smooth is value to smooth by (seconds)
    % y = limits of y axis
    % light = time limits of light stimuli on (e.g. for 2 Hz, 4 sec stimuli, light = [0,4]
    % data [electrodes x samples]
    % fs = sampling rate [Hz]
    % stimuli = times of stimulation (seconds)
    % leadTime = time prior to stimuli onset (seconds)
    % lenTrial = time posterior to stimuli onset (seconds)
    % badEpochs = vector of trials to disregard/skip
    % impedance mapping = vector of impedances
    % mapping = correct order of electrodes to make 2D plot

% Outputs
    % trialAverages = average trial data
    % t_trial = time/latency of trial average with stimuli onset at t_trial = 0

% convert stimuli from seconds to indicies
for i = 1:length(stimuli)
    stimuli(i) = find(t >= stimuli(i),1);
end

% calculate trial averages for each channel
trialAverages = zeros(size(data,1), round(fs*(lenTrial-leadTime)));
numTrials = length(stimuli);
t_trial = leadTime:1/fs:lenTrial-1/fs;
t_reverse = fliplr(t_trial);
t2 = [t_trial, t_reverse];

% Sum all the trials into a single matrix. The order is the electrode
% mapping order.
figure;
for i = 1:size(data,1)
    subplot(4,4,find(mapping == good_channels(i))); % plot electrode data in correct 2D location
    sdData = zeros(numTrials-length(badEpochs), round(fs*(lenTrial-leadTime)));
    lastUsedJ = 1;
    for j = 1:numTrials
        if find(badEpochs == j)
            continue;
        end
        
        trialAverages(i,:) = trialAverages(i,:) + data(i, stimuli(j)+round(leadTime*fs):(stimuli(j)+round(lenTrial*fs)-1));
        sdData(lastUsedJ, :) = data(i, stimuli(j)+(leadTime*fs):(stimuli(j)+(lenTrial*fs)-1));
        
        lastUsedJ = lastUsedJ+1;
    end
    
    trialAverages(i,:) = trialAverages(i,:)/(numTrials-length(badEpochs));
    %[pk, loc] = findpeaks(-1*trialAverages(i,fs*leadTime*-1:fs*leadTime*-1+fs*0.1),'SortStr','descend');
    sd = std(sdData);
    curve1 = trialAverages(i,:) + sd;
    curve2 = trialAverages(i,:) - sd;
    
    if smooth
        d = smoothdata(trialAverages(i,:), 'gaussian', fs*smooth);
        curve1 = smoothdata(curve1, 'gaussian', fs*smooth);
        curve2 = smoothdata(curve2, 'gaussian', fs*smooth);
    else
        d = trialAverages(i,:);
    end
    
    %plot(t_trial(1:1000:end), curve1(1:1000:end), 'Color', 'w', 'LineWidth', 2); hold on;
    %plot(t_trial(1:1000:end), curve2(1:1000:end), 'Color', 'w', 'LineWidth', 2); hold on;
    if smooth
        btwn = [curve1(1:fs*smooth:end), fliplr(curve2(1:fs*smooth:end))];
        fill(t2(1:fs*smooth:end), btwn, [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on;
        plot(t_trial(1:fs*smooth:end), d(1:fs*smooth:end), 'Color', 'k', 'LineWidth', 2); hold on;
    else
        btwn = [curve1, fliplr(curve2)];
        fill(t2, btwn, [0.8, 0.8, 0.8], 'EdgeColor', 'none'); hold on;
        plot(t_trial, d, 'Color', 'k', 'LineWidth', 2); hold on;
    end

    
    %ylabel(ylab);
    %xlabel('Time (s)');
    %title(['Channel ' num2str(good_channels(i))]);
    
    %text(1, -150, [num2str(t_trial(loc(1)+fs*leadTime*-1)) ' s']); 
    ylim(y);
    xlim([leadTime lenTrial]);
    
    %fill([light(1), light(2), light(2), light(1)], [y(1), y(1), y(1)+10, y(1)+10], [1 0.5 0.25], 'EdgeColor', 'none');
    
    
    %plot([light(1) light(1)], get(gca, 'ylim'), 'Color', 'b'); % Plot stimuli onset line
    %plot([light(2) light(2)], get(gca, 'ylim'), 'Color', 'r'); % Plot stimuli offset line
end

end
