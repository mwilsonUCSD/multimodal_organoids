function [plv_light_final_2D, plv_dark_final_2D, pval_color] = f_2D_PLV(phase_specgram_light, phase_specgram_dark, spikes_light, spikes_dark, ind_light_all, ind_dark_all, win, t_amplifier, freq, f, nBoots, params)
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
    % nBoots = how many bootstrap samples you want to take
    % params = struct containing tapers, fs, pad, and f_pass
    
% Outputs
    % plv_light_final_2D = 2D phase locking of light MUA spikes to band-filtered LFP of other channels [nChannels x nChannels x nFreq_band]
    % plv_dark_final_2D = 2D phase locking of dark MUA spikes to band-filtered LFP of other channels [nChannels x nChannels x nFreq_band]
    % pval_color = p-value of the light vs dark phase locking 2D arrays [nChannels x nChannels x nFreq_band]

%Calculate indices of frequencies of interest
iFreq = [find(f >= freq(1),1), find(f >= freq(2),1)]; %index of desired frequency band
nFreq_band = iFreq(2)-iFreq(1)+1;
pval_color = zeros(size(spikes_light,2), size(spikes_light,2), nFreq_band);

% Holder variables
plv_light_final_2D = zeros(size(spikes_light,2), size(spikes_light,2), nFreq_band);
plv_dark_final_2D = zeros(size(spikes_dark,2), size(spikes_light,2), nFreq_band);

nSpikes = zeros(size(spikes_light,2),1);
    
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
    elseif length(ind_light) > length(ind_dark)
        temp = randperm(length(ind_light)); % Select spikes randomly from light
        ind_light = ind_light(temp(1:length(ind_dark)));
    end

    %Equalize number of spikes to smallest number in channels
%     temp = randperm(length(ind_dark)); % Select spikes randomly from dark
%     ind_dark = ind_dark(temp(1:44));
%     temp = randperm(length(ind_light)); % Select spikes randomly from light
%     ind_light = ind_light(temp(1:44));
    
    % Find spike indices to pull this channel's spikes from the total spike spectrogram.
    % size(i_match_light) should = size(ind_light) and same for dark
    i_match_light = find(ismember(ind_light_all,ind_light));
    i_match_dark = find(ismember(ind_dark_all,ind_dark));
    
    nSpikes(i) = length(i_match_light);
    
    % Holder variables
    phi_light_2D = zeros(nBoots, nSpikes(i), nFreq_band, params.tapers(2));
    phi_dark_2D = zeros(nBoots, nSpikes(i), nFreq_band, params.tapers(2));
    
    % Create a random bootstrap index matrix to generate samples from
    bootMatrix_light = randi(length(i_match_light), nBoots, nSpikes(i));
    bootMatrix_dark = randi(length(i_match_dark), nBoots, nSpikes(i));
    
    % Generate spectrograms specific to this channel's spikes for PLV calculation
    this_specgram_light = phase_specgram_light(i_match_light, iFreq(1):iFreq(2), :, :);
    this_specgram_light = squeeze(this_specgram_light);
    this_specgram_dark = phase_specgram_dark(i_match_dark, iFreq(1):iFreq(2), :, :);
    this_specgram_dark = squeeze(this_specgram_dark);
    
    tic;
    for k = 1:size(spikes_light,2)
        %Bootstrap samples of spikes for k lfp channels
        for j = 1:nBoots
            phi_light_2D(j,:,:,:) = this_specgram_light(bootMatrix_light(j,:),:,:,k);
            phi_dark_2D(j,:,:,:) = this_specgram_dark(bootMatrix_dark(j,:),:,:,k);
        end
        
        %light PLV
        plv_light = abs(squeeze(mean(exp(1i*phi_light_2D),2))); %mean across spikes
        plv_light_t = squeeze(mean(plv_light, 3)); % mean across tapers
        plv_light_b = mean(plv_light_t, 1);% mean across bootstrapped PLVs
        plv_light_final_2D(i,k,:) = plv_light_b;

        %dark PLV
        plv_dark = abs(squeeze(mean(exp(1i*phi_dark_2D),2))); %mean across spikes
        plv_dark_t = squeeze(mean(plv_dark, 3)); % mean across tapers
        plv_dark_b = mean(plv_dark_t, 1);% mean across bootstrapped PLVs
        plv_dark_final_2D(i,k,:) = plv_dark_b;
        
        % calculate p-value
        temp = plv_light_t > plv_dark_t;
        for j = 1:nFreq_band
            pval_color(i,k,j) = min(length(find(temp(:,j) == 1)), length(find(temp(:,j) == 0)));
        end
    end
    toc;
    disp(i);
end

pval_color = pval_color/nBoots;

end