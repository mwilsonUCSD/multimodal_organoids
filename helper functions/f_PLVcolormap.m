function f_PLVcolormap(data_hi, plv, spikes, el, clim, freq, good_channels, mapping, fs, nickname)

temp = zeros(16,1);
temp(good_channels) = plv;
itpc2D = temp(mapping);
itpc2D = reshape(itpc2D, [4 4])';
itpc2D = flipud(itpc2D);
itpc2D = interp2(itpc2D,2);

figure;
cmap1 = pcolor(itpc2D); hold on;
cmap1.FaceColor = 'interp';
cmap1.EdgeColor = 'none';
%set(gca,'xtick',[]);
%set(gca,'ytick',[]);
c = colorbar;
c.Label.String = 'PLV';
caxis([0,clim]);

% Plot average LFP following spikes
range = [-0.25*1/freq(1), 0.75*1/freq(1)];
t_mini = range(1):1/fs:range(2)-1/fs;
t_mini_down = linspace(-1,1,length(t_mini));%t_mini(1:fs/500:end);
%t_mini_norm = normalize(t_mini_down);
yshift = [8,8,5,2,2,5,2,2,5,5,8,11,11,8,11,11];
xshift = [11,8,11,11,8,8,5,2,2,5,2,2,5,5,8,11];

%figure;
for i = 1:length(good_channels)
    %ch = find(mapping == good_channels(i));
    if el
        ind = find(spikes(:,el) == 1);
    else
        ind = find(spikes(:,i) == 1);
    end
    
    %ind_remove = find(spikes_sum_light > 2);
    %ind = setdiff(ind,ind_remove);
    
    % Get rid of indicies too close to beginning or end of trial
    ind(ind <= -1*fs*range(1)) = [];
    ind(ind > length(data_hi)-fs*range(2)+1) = [];
    
    % Create an array of LFPs surrounding spike times
    temp = zeros(length(t_mini),length(ind));
    for j = 1:length(ind)
        temp(:,j) = data_hi(ind(j)+fs*range(1):ind(j)+fs*range(2)-1,i);
        %plot(t_mini, temp(:,j)); hold on;
    end

    plot(t_mini_down + yshift(good_channels(i)), mean(temp,2)/4 + xshift(good_channels(i)), 'LineWidth', 2, 'Color', 'k'); hold on;
    
    %ylabel('Amplitude (\muV)');
    %xlabel('Time (seconds)');
    %title(['Spike Count: ' num2str(j)]);
    %pause(1);
end
title(['PLV of ' nickname ' filtered btwn: ' num2str(freq(1)) '-' num2str(freq(2)) ' Hz'], 'Interpreter', 'none');

end