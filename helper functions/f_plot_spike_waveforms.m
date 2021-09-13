function f_plot_spike_waveforms(data, spikeTimes, fs, leadTime, lenSpike, thresh, mapping, good_channels)
    nChannels = size(data,2);
    indAroundSpike = -leadTime*fs:lenSpike*fs;
    tSpike = -leadTime:1/fs:lenSpike;
    tSpike = tSpike*1000;
    
    tempArray = zeros(500,nChannels);
    
    figure;
    for i = 1:nChannels
        %subplot(4,4,find(mapping == good_channels(i)));%(16,nChannels, i);
        spikesHere = find(spikeTimes(:,i) == 1);
        for j = 1:length(spikesHere)
            %bump peak back until lines up with threshold value
            while data(spikesHere(j),i) < -1*thresh(i)
                spikesHere(j) = spikesHere(j)-1;
            end
            
            %plot(tSpike, data(spikesHere(j) + indAroundSpike, i)); hold on;
        end
        tempArray(1:length(spikesHere),i) = spikesHere;
        
        for k = 0:nChannels-1
            subplot(nChannels,nChannels,i+k*nChannels);
            d = data(:,k+1);
            plot(tSpike, mean(d(spikesHere + indAroundSpike)));
            if k == 0
                title(['Spike ch: ' num2str(i) ' Data ch: ' num2str(k+1)]);
            end
            ylim([-50 50]);
        end
        %title(['Channel: ' num2str(i)]);
        %ylim([-120 120]);
        %ylabel('Voltage \muV');
        %xlabel('Time (ms)');
        disp(i);
    end
end