function [plv_light_final, plv_dark_final, pval, is_significant] = f_plot_PLV_vs_freq(phase_specgram_light, phase_specgram_dark, spikes_light, spikes_dark, ind_light_all, ind_dark_all, t_amplifier, nBoots, nFreq, win, f, mapping, good_channels, params, nickname)
% Function to calculate and plot the phase locking value versus frequency for all channels

% Inputs:
    % phase_specgram_light = total spike spectrogram for light spikes [nTotalLightSpikes x nFreq x nTapers x nChannels]
    % phase_specgram_dark = total spike spectrogram for dark spikes [nTotalDarkSpikes x nFreq x nTapers x nChannels]
    % spikes_light = array containing ones when there was a light MUA spike and 0s otherwise [time x channels]
    % spikes_dark = array containing ones when there was a dark MUA spike and 0s otherwise [time x channels]
    % ind_light_all = indices of all the light spikes in relation to t_amplifier [nTotalLightSpikes]
    % ind_dark_all = indices of all the dark spikes in relation to t_amplifier [nTotalDarkSpikes]
    % t_amplifier = aquisition time from amplifier board
    % nBoots = how many bootstrap samples you want to take
    % nFreq = number of frequencies (resolution of x-axis)
    % win = window size [seconds * fs]
    % f = frequencies [nFreq]
    % mapping = index order of the channels/how they are positioned in the physical array
    % good_channels = channels with impedances under the threshold
    % params = struct containing tapers, fs, pad, and f_pass
    % nickname = nickname of data for figure labeling
    
% Outputs:
    % plv_light = PLV values for all frequencies and channels for light MUA spikes [nFreq x nChannel]
    % plv_dark = PLV values for all frequencies and channels for dark MUA spikes [nFreq x nChannel]
    % pval = p-value between light and dark PLV for all frequencies [nFreq]
    % is_significant = whether the p-value is < 0.05 for each frequency [nFreq]


%declare holder variables
plv_light_final = zeros(nFreq, size(spikes_light,2));
plv_dark_final = zeros(nFreq, size(spikes_light,2));
pval = zeros(nFreq, size(spikes_light,2));

%nLight = zeros(size(spikes_light,2),1);
%nDark = zeros(size(spikes_dark,2),1);
nSpikes = zeros(size(spikes_light,2),1);

for i = 1:size(spikes_light,2) % loop through good channels
%     if i == 4 || i == 7 % || i == 2 || i == 3 || i == 4
%         continue;
%     end
    tic;
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
    
%     %equalize spikes to smallest number of spikes in all channels (hard
%     %coded at the moment)
%     temp = randperm(length(ind_dark)); % Select spikes randomly from dark
%     ind_dark = ind_dark(temp(1:44));
%     temp = randperm(length(ind_light)); % Select spikes randomly from light
%     ind_light = ind_light(temp(1:44));

    % Find spike indices to pull this channel's spikes from the total spike spectrogram.
    % size(i_match_light) should = size(ind_light) and same for dark
    i_match_light = find(ismember(ind_light_all,ind_light));
    i_match_dark = find(ismember(ind_dark_all,ind_dark));

    %Variables to store # of light and dark spikes (should be equal)
    %nLight(i) = length(i_match_light);
    %nDark(i) = length(i_match_dark);
    nSpikes(i) = length(i_match_light);

    % Generate spectrograms specific to this channel's spikes for PLV calculation
    this_specgram_light = phase_specgram_light(i_match_light, :, :, i);%i
    this_specgram_light = squeeze(this_specgram_light);
    this_specgram_dark = phase_specgram_dark(i_match_dark, :, :, i);%i
    this_specgram_dark = squeeze(this_specgram_dark);

    % Generate bootstrap samples
    phi_light = this_specgram_light(randi(nSpikes(i),nBoots*nSpikes(i),1),:,:);
    phi_light = reshape(phi_light, [nBoots, nSpikes(i), nFreq, params.tapers(2)]);

    phi_dark = this_specgram_dark(randi(nSpikes(i),nBoots*nSpikes(i),1),:,:);
    phi_dark = reshape(phi_dark, [nBoots, nSpikes(i), nFreq, params.tapers(2)]);

    %light PLV
    plv_light = abs(squeeze(mean(exp(1i*phi_light),2))); % circular mean across spikes
    plv_light_t = squeeze(mean(plv_light, 3)); % mean across tapers
    plv_light_b = mean(plv_light_t, 1);% mean across bootstrapped PLVs
    plv_light_final(:,i) = plv_light_b;

    %dark PLV
    plv_dark = abs(squeeze(mean(exp(1i*phi_dark),2))); % circular mean across spikes
    plv_dark_t = squeeze(mean(plv_dark, 3)); % mean across tapers
    plv_dark_b = mean(plv_dark_t, 1);% mean across bootstrapped PLVs
    plv_dark_final(:,i) = plv_dark_b;

    % Calculation of p-value
    temp = plv_light_t > plv_dark_t;
    for j = 1:nFreq
        pval(j,i) = min(length(find(temp(:,j) == 1)), length(find(temp(:,j) == 0)));
    end
    toc;
    disp(i);
end

% Calculate significance
pval = pval/1000;
is_significant = double(pval < 0.05);
is_significant(is_significant == 0) = -1;

% Plot PLV vs. frequencies
figure;
for i = 1:size(spikes_light,2)
    subplot(4,4,find(mapping == good_channels(i)));
    %subplot(4,4,i); 
    plot(f, plv_light_final(:,i), 'b'); hold on;
    plot(f, plv_dark_final(:,i), 'r');
    scatter(f, is_significant(:,i) * 0.05, 'k.');
    fill([0,4,4,0], [0,0,0.01,0.01], 'r', 'EdgeColor', 'none');
    fill([4,8,8,4], [0,0,0.01,0.01], [1 0.5 0.25], 'EdgeColor', 'none');
    fill([8,12,12,8], [0,0,0.01,0.01], 'y', 'EdgeColor', 'none');
    fill([12,32,32,12], [0,0,0.01,0.01], 'g', 'EdgeColor', 'none');
    fill([32,80,80,32], [0,0,0.01,0.01], 'c', 'EdgeColor', 'none');
    fill([80,150,150,80], [0,0,0.01,0.01], 'b', 'EdgeColor', 'none');
    ylim([0 0.5]);
    ylabel('PLV');
    xlabel('Frequency (Hz)');
    title(['nSpikes = ' num2str(nSpikes(i))]);
    %set(gca, 'FontSize', 18);
    %title(['F = ' num2str(F(i))]);
    %title(['n spikes = ' num2str(nSpikes(i))]);
end
sgtitle([nickname ' tapers: [' num2str(params.tapers(1)) ', ' num2str(params.tapers(2)) '] Winsize: ' num2str(win/params.Fs) ' sec nBoots: ' num2str(nBoots)], 'Interpreter', 'none');

end