function phase = f_plv_XL(signal, spktime, win, nFreq, params)
    % Code written by Xin Liu to help me simulate data for testing PLV
    
    [tapers,pad,Fs,fpass]=getparams(params);
    nTapers = tapers(2);
    nChannels = size(signal, 2);
    phase = zeros(length(spktime),nFreq,nTapers,nChannels);
    signal = change_row_to_column(signal);
    
    N = length(-(win/2)+1:win/2); %M.W. moved this and next three lines out of loop
    nfft = max(2^(nextpow2(N)+pad),N);
    [f,findx] = getfgrid(Fs,nfft,fpass);
    tapers = dpsschk(tapers,N,Fs); % check tapers

    for s = 1:length(spktime)
        data = signal(spktime(s) + (-(win/2)+1:win/2),:);
        % below code is copy-paste from chronux function "mtspectrumc" and
        % modify the final average step
%         N = size(data,1);
%         nfft = max(2^(nextpow2(N)+pad),N);
%         [f,findx] = getfgrid(Fs,nfft,fpass);
%         tapers = dpsschk(tapers,N,Fs); % check tapers
        J = mtfftc(data,tapers,nfft,Fs);
        J = J(findx,:,:);
        phase(s,:,:,:) = angle(J);
    end
    %plv = abs(squeeze(mean(exp(1i*phase),3))); % Circular mean across tapers
end