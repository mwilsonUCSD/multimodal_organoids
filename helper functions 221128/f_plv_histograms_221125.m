function [] = f_plv_histograms(phase_specgram_light, phase_specgram_dark, spikes_light, spikes_dark, ind_light_all, ind_dark_all, win, t_amplifier, freq, f, mapping, good_channels)
% Function to compute the PLV of MUA spikes to a given frequency band across all channels, allowing us to plot a 2D PLV colormap.
% The colormap shows which channels the MUA spikes phase locked most to.

% Inputs
    % phase_specgram_light = total spike spectrogram for light spikes [nTotalLightSpikes x nFreq x nTapers x nChannels]
    % phase_specgram_dark = total spike spectrogram for dark spikes [nTotalDarkSpikes x nFreq x nTapers x nChannels]
    % spikes_light = array containing ones when there was a light MUA spike and 0s otherwise [time x channels]
    % spikes_dark = array containing ones when there was a dark MUA spike and 0s otherwise [time x channels]
    % ind_light_all = indices of all the light spikes in relation to t_amplifier [nTotalLightSpikes]
    % ind_dark_all = indices of all the dark spikes in relation to t_amplifier [nTotalDarkSpikes]
    % win = window size [seconds * fs]
    % t_amplifier = aquisition time from amplifier board
    % freq = frequency range of interest [starting frequency, ending frequency]
    % f = frequencies [nFreq]
    % mapping = index order of the channels/how they are positioned in the physical array
    % good_channels = channels with impedances under the threshold
    
% Outputs
    % 

%Calculate indices of frequencies of interest
iFreq = find(f >= freq,1); %index of desired frequency
pval_color = zeros(size(spikes_light,2), size(spikes_light,2));

% Holder variables
% plv_light_final_2D = zeros(size(spikes_light,2), size(spikes_light,2), nFreq_band);
% plv_dark_final_2D = zeros(size(spikes_dark,2), size(spikes_light,2), nFreq_band);

nLight = zeros(size(spikes_light,2),1);
nDark = zeros(size(spikes_light,2),1);
%nSpikes = zeros(size(spikes_light,2),1);
    
f1 = figure;
f2 = figure;
for i = 1:size(spikes_light,2)
    %determine spike times for that channel
    ind_light = find(spikes_light(:,i) == 1);
    ind_dark = find(spikes_dark(:,i) == 1);
    ind_dark(ind_dark < win/2) = []; %get rid of spikes too close to beginning of recording
    ind_dark(ind_dark > length(t_amplifier) - win/2) = []; % get rid of spikes too close to end of recording

    % Equalize number of light and dark spikes for fair bootstrap PLV comparison
    % Note: a random sample of dark spikes is taken before any PLV calculation so the dark plot might change
    % values between runs. But since it's random the results should be similar.
    if length(ind_dark) > length(ind_light)
        temp = randperm(length(ind_dark)); % Select spikes randomly from dark
        ind_dark = ind_dark(temp(1:length(ind_light)));
    end

    % Find spike indices to pull this channel's spikes from the total spike spectrogram.
    % size(i_match_light) should = size(ind_light) and same for dark
    i_match_light = find(ismember(ind_light_all,ind_light));
    i_match_dark = find(ismember(ind_dark_all,ind_dark));
    
    nLight(i) = length(i_match_light);
    nDark(i) = length(i_match_dark);
    %nSpikes(i) = length(i_match_light);
    
    % Generate spectrograms specific to this channel's spikes for PLV calculation
    this_specgram_light = phase_specgram_light(i_match_light, iFreq, :, i);
    this_specgram_light = squeeze(this_specgram_light);
    this_specgram_dark = phase_specgram_dark(i_match_dark, iFreq, :, i);
    this_specgram_dark = squeeze(this_specgram_dark);
    
    % Calculate light angles, itpc, and pval
    phase_light_avg = angle(mean(exp(1i*this_specgram_light), 2));
    [pval_light, ~] = circ_rtest(phase_light_avg);
    itpc1 = abs(squeeze(mean(exp(1i*this_specgram_light),1))); % circular mean across spikes
    itpc_light = squeeze(mean(itpc1)); % mean across tapers
    %itpc_light = abs(mean(exp(1i*phase_light_avg)));
    prefAngle_light = angle(mean(exp(1i*phase_light_avg)));
    
    % Plot light histogram
    figure(f1);
    subplot(4,4,find(mapping == good_channels(i)));
    if pval_light < 0.05
        polarhistogram(phase_light_avg, 36); hold on;
    else
        polarhistogram(phase_light_avg, 36, 'FaceColor', 'r'); hold on;
    end
    polarplot([0, prefAngle_light], [0, itpc_light*50], 'Color', 'm');
    
    rlim([0 6]);
    title(['n = ' num2str(nLight(i)) ' plv: ' num2str(itpc_light)]);
    
    % Calculate dark angles, itpc, and pval
    phase_dark_avg = angle(mean(exp(1i*this_specgram_dark), 2));
    [pval_dark, ~] = circ_rtest(phase_dark_avg);
    itpc2 = abs(squeeze(mean(exp(1i*this_specgram_dark),1))); % circular mean across spikes
    itpc_dark = squeeze(mean(itpc2)); % mean across tapers
    %itpc_dark = abs(mean(exp(1i*phase_dark_avg)));
    prefAngle_dark = angle(mean(exp(1i*phase_dark_avg)));
    
    % Plot dark histogram
    figure(f2);
    subplot(4,4,find(mapping == good_channels(i)));
    if pval_dark < 0.05
        polarhistogram(phase_dark_avg, 36); hold on;
    else
        polarhistogram(phase_dark_avg, 36, 'FaceColor', 'r'); hold on;
    end
    polarplot([0, prefAngle_dark], [0, itpc_dark*50], 'Color', 'm');
    
    rlim([0 6]);
    title(['n = ' num2str(nDark(i)) ' plv: ' num2str(itpc_dark)]);
end

figure(f1);
sgtitle(['Light Histograms for freq = ' num2str(freq) ' Hz']);

figure(f2);
sgtitle(['Dark Histograms for freq = ' num2str(freq) ' Hz']);

end