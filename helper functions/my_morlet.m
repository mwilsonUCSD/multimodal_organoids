function tf = my_morlet(x,srate,flo,fhi,deltaf, range_cycles)
if ~isrow(x)
    x = x';
end

% frequency parameters
min_freq = flo;
max_freq = fhi;
num_frex = (fhi-flo+1)/deltaf;
frex = linspace(min_freq,max_freq,num_frex);

% other wavelet parameters
%s = logspace(log10(range_cycles(1)),log10(range_cycles(end)),num_frex) ./ (2*pi*frex);
s = linspace(range_cycles(1),range_cycles(end),num_frex) ./ (2*pi*frex);


N_orig=length(x);
% now loop over trials...
N=2^(nextpow2(length(x)));
x=[x,zeros(1,N-length(x))];
% initialize output time-frequency data
tf = zeros(length(frex),N);
dataX = fft(x);
freq_samples=srate*1/N*(-N/2:N/2-1);

% loop over frequencies
for fi=1:length(frex)
    % create wavelet and get its FFT
    % wavelet  = exp(2*1i*pi*frex(fi).*wavtime) .* exp(-wavtime.^2./(2*s(fi)^2));
    % waveletX = fft(wavelet,nConv);
    waveletX = sqrt(pi/(2*s(fi)^2))*exp(-2*pi^2*s(fi)^2*(freq_samples-frex(fi)*ones(1,N)).^2);
    waveletX = waveletX ./ max(waveletX);
    % run convolution
    as = ifft(ifftshift(waveletX) .* dataX);
    %as = as(half_wave+1:end-half_wave);
    
    % compute the amplitude instead of power
    %tf(fi,:) = 2*abs(as);
    
    %compute power instead
    %tf(fi,:) = mean(abs(as).^2,2); % Mean is used for trial averaging
    tf(fi,:) = abs(as).^2;
end

tf=tf(:,1:N_orig);
