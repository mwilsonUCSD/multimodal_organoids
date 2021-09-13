function [itpc_light, itpc_dark, f] = f_bootstrapPLV(data_lfp, t_amplifier, spikes_light, spikes_dark, nboots, winsize, params)

N=winsize*params.Fs;
nfft=max(2^(nextpow2(N)+params.pad),N);
[f,~]=getfgrid(params.Fs,nfft,params.fpass);
nFreq = length(f);

% Light phi
%phi_light = zeros(length());
itpc_light = zeros(size(spikes_light,2),nboots,nFreq);
itpc_dark = zeros(size(spikes_dark,2),nboots,nFreq);

%allphi_light = zeros(length(ind_light), 26215, 16, 5);
for i = 1:size(spikes_light,2) % num channels
    % Light spikes
    tic;
    ind_light = find(spikes_light(:,i) == 1);
    % calculate spike phases
    allphi_light = zeros(length(ind_light), nFreq, 5);
    for j = 1:length(ind_light)%num spikes
        [allphi_light(j,:,:),f]=mtspectrumc_mnw_allphi(data_lfp(ind_light(j)-params.Fs*winsize/2:ind_light(j)+params.Fs*winsize/2-1,i),params);
        %allphi_light(j,:,i,:) = a;
    end
    allphi_light = exp(1i*allphi_light);
    nTapers = size(allphi_light, 3);
    toc;
    
    %Dark spikes
    tic;
    ind_dark = find(spikes_dark(:,i) == 1);
    ind_dark(ind_dark < params.Fs*winsize/2) = [];
    ind_dark(ind_dark > length(t_amplifier) - params.Fs*winsize/2) = [];
    % calculate spike phases
    allphi_dark = zeros(length(ind_dark), nFreq, 5);
    for j = 1:length(ind_dark)%num spikes
        [allphi_dark(j,:,:),f]=mtspectrumc_mnw_allphi(data_lfp(ind_dark(j)-params.Fs*winsize/2:ind_dark(j)+params.Fs*winsize/2-1,i),params);
        %allphi_light(j,:,i,:) = a;
    end
    allphi_dark = exp(1i*allphi_dark);
    toc;
    
    tic;
    % bootstrap
    temp_itpc_light = zeros(nboots, length(ind_light), length(f), nTapers);
    indBootstrp_light = randi(length(ind_light), nboots, length(ind_light));
    
    temp_itpc_dark = zeros(nboots, length(ind_dark), length(f), nTapers);
    indBootstrp_dark = randi(length(ind_dark), nboots, length(ind_dark));
    
    for j = 1:nboots
        temp_itpc_light(j,:,:,:) = allphi_light(indBootstrp_light(j,:),:,:);
        temp_itpc_dark(j,:,:,:) = allphi_dark(indBootstrp_dark(j,:),:,:);
    end
    itpc_light(i,:,:) = mean(abs(squeeze(mean(temp_itpc_light,2))),3);
    itpc_dark(i,:,:) = mean(abs(squeeze(mean(temp_itpc_dark,2))),3);
    toc;
    
    disp(i);
end

end