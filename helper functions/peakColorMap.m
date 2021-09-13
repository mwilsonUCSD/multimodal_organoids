function peakColorMap(data, t, fs, nPeak, pPeak, leadTime, lenTrial, good_channels, mapping, nickname)
% Function to plot the peak amplitude and time delay

% nPeak = upper bound estimate for time of negative peak
% pPeak = upper bound estimate for time of positive peak (needs to be after nPeak if data has neg peak before pos peak)

negPeak = zeros(size(data,1),1);
posPeak = zeros(size(data,1),1);
iNP = zeros(size(data,1),1);
iPP = zeros(size(data,1),1);

figure;
for i = 1:size(data,1)
    [np, locN] = findpeaks(-1*data(i,fs*(leadTime*-1):fs*leadTime*-1+fs*nPeak),'SortStr','descend');
    %locN = locN+0.02*fs;
    [pp, locP] = findpeaks(data(i,fs*leadTime*-1:fs*leadTime*-1+fs*pPeak),'SortStr','descend');
    
    negPeak(i) = -1*np(1);
    posPeak(i) = pp(1);
    
    iNP(i) = t(locN(1)+fs*leadTime*-1);
    iPP(i) = t(locP(1)+fs*leadTime*-1);
    
    subplot(4,4,find(mapping == good_channels(i)));
    plot(t, data(i,:)); hold on;
    scatter(t(locN(1)+fs*leadTime*-1), -1*np(1)); hold on;
    scatter(t(locP(1)+fs*leadTime*-1), pp(1));
    ylim([-150 150]);
    xlim([leadTime, lenTrial]);
end
sgtitle(nickname, 'Interpreter', 'none');

% reshape the data
peak1D = zeros(16,1);
peak1D(good_channels) = posPeak;
peak1D = peak1D(mapping);
peak2D = reshape(peak1D, [4 4])';
peak2D = flipud(peak2D);
peak2D = interp2(peak2D,2);
%peak2D = [peak2D; peak2D(4,:)];
%peak2D = [peak2D, peak2D(:,4)];

peak1D = zeros(16,1);
peak1D(good_channels) = negPeak;
peak1D = peak1D(mapping);
peak2Dn = reshape(peak1D, [4 4])';
peak2Dn = flipud(peak2Dn);
peak2Dn = interp2(peak2Dn,2);

mappedINP = zeros(16,1);
mappedINP(good_channels) = iNP;
mappedINP = mappedINP(mapping);
np2D = reshape(mappedINP, [4 4])';
np2D = flipud(np2D);
np2D = interp2(np2D,2);

mappedIPP = zeros(16,1);
mappedIPP(good_channels) = iPP;
mappedIPP = mappedIPP(mapping);
pp2D = reshape(mappedIPP, [4 4])';
pp2D = flipud(pp2D);
pp2D = interp2(pp2D,2);

% plot the positive amplitude map
figure;
cmap1 = pcolor(peak2D);
cmap1.FaceColor = 'interp';
c = colorbar;
c.Label.String = '|Amplitude|';
caxis([0,150]);
title(['Positive Peak Amplitudes ' nickname], 'Interpreter', 'none');
set(gcf, 'Position', [200 200 900 750]);

% plot the negative amplitude map
figure;
cmap1 = pcolor(abs(peak2Dn));
cmap1.FaceColor = 'interp';
c = colorbar;
c.Label.String = '|Amplitude|';
caxis([0,150]);
title(['Negative Peak Amplitudes ' nickname], 'Interpreter', 'none');
set(gcf, 'Position', [200 200 900 750]);

% plot the negative delay map
figure;
cmap1 = pcolor(np2D*1000);
cmap1.FaceColor = 'interp';
c = colorbar;
c.Label.String = 'Delay (ms)';
caxis([34 42]);%caxis([min(min(iNP*1000)) max(max(iNP*1000))]);
title(['Negative Peak Delays ' nickname], 'Interpreter', 'none');
set(gcf, 'Position', [200 200 900 750]);

% plot the positive delay map
figure;
cmap1 = pcolor(pp2D*1000);
cmap1.FaceColor = 'interp';
c = colorbar;
c.Label.String = 'Delay (ms)';
caxis([min(min(iPP*1000)) max(max(iPP*1000))]);%35,43]);%
title(['Positive Peak Delays ' nickname], 'Interpreter', 'none');
set(gcf, 'Position', [200 200 900 750]);

end