function itpc = f_PLVfullcolormap(data_hilbert, spikes, nOpp, clim, freq, good_channels, mapping, nickname, nBoot, nSpikes)

itpc = zeros(size(data_hilbert,2),length(good_channels));
angleBoot = zeros(nSpikes,nBoot);
figure;
for j = 1:length(good_channels)
    for i = 1:size(data_hilbert,2)
        hAngles = angle(data_hilbert(spikes(:,j)==1,i));
        iBoot = randi(length(hAngles), nSpikes, nBoot);
        for k = 1:nBoot
            angleBoot(:,k) = hAngles(iBoot(:,k));
        end
        
%         if matchSpikeCount
%             if length(hAngles) > nOpp(j)
%                 randNum = randperm(nOpp(j));
%                 hAngles = hAngles(randNum(1:nOpp(j)));
%             end
%         end
        itpc(i,j) = mean(abs(mean(exp(1i*angleBoot))));%abs(mean(exp(1i*hAngles)));
    end
    % Color map of PLV
    temp = NaN(16,1);
    temp(good_channels) = itpc(:,j);
    itpcMap = temp;
    itpcMap = itpcMap(mapping);
    itpc2D = reshape(itpcMap, [4 4])';
    itpc2D = fillmissing(itpc2D, 'linear',2, 'EndValues','nearest');
    itpc2D = flipud(itpc2D);
    itpc2D = interp2(itpc2D,2);

    subplot(4,4,find(mapping == good_channels(j)));
    cmap1 = pcolor(itpc2D);
    cmap1.FaceColor = 'interp';
    cmap1.EdgeColor = 'none';
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    c = colorbar;
    c.Label.String = 'PLV';
    caxis(clim);
    title(['PLV of ' nickname], 'Interpreter', 'none');
    title(['n = ' num2str(length(hAngles))]);
end
%sgtitle([nickname ' filtered btwn: ' num2str(freq(1)) '-' num2str(freq(2)) ' Hz'], 'Interpreter', 'none');

end