function f_plotPointTrialAverages(spikes, fs, t, leadTime, lenTrial, t_stim, onsets)
% Inputs:
% spikes = binary array of spike times
% fs = frequency
% t = time matching spikes
% leadTime = time leading up to stimulus onset
% lenTrial = time after stimulus onset
% t_stim = times in t where the light onsets happen
% onsets = times that the light flashes

%     leadTime = -1;
%     lenTrial = 9.5;
    t_trial = leadTime:1/fs:lenTrial-1/fs;

    stimuli = zeros(1,length(t_stim));
    for i = 1:length(t_stim)
        stimuli(i) = find(t >= t_stim(i),1);
    end
    
    total = zeros(size(spikes,2), round(fs*(lenTrial-leadTime)));
    for i = 1:16
        for j = 1:20
            total(i,:) = total(i,:) + spikes2(i, (stimuli(j)+leadTime*fs):(stimuli(j)+lenTrial*fs_down-1));
        end
    end

    figure;
    total(total == 0) = nan;
    for i = 1:16
        scatter(t_trial, total(i,:)+i); hold on;
    end

    %onsets = [0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5];
    for i = 1:length(onsets)
        plot([onsets(i) onsets(i)],get(gca,'ylim'),'k','linew',1); hold on;
    end
end