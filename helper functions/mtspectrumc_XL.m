function [S,f,Serr]=mtspectrumc_XL(data,params)
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
data=change_row_to_column(data);
N=size(data,1);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
J=mtfftc(data,tapers,nfft,Fs);
J=J(findx,:,:);
% S=permute(mean(conj(J).*J,2),[1 3 2]);
S=permute(J,[1 3 2]); % keep the result of individual tapers. X.L.
if trialave; S=squeeze(mean(S,2));else S=squeeze(S);end;
if nargout==3; 
   Serr=specerr(S,J,err,trialave);
end;