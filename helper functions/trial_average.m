%% MUA Trial Averages
leadTime = -1;
lenTrial = 7;
badEpochs = [];

trialAverages_mua = zeros((lenTrial - leadTime)*fs, length(good_channels));
j = 5;

temp = zeros((lenTrial - leadTime)*fs, length(t_stim)-length(badEpochs));
lastUsedJ = 1;
for i = 1:length(t_stim)
    if find(badEpochs == i)
        continue;
    end
    ind = find(t_amplifier >= t_stim(i),1);

    temp(:,lastUsedJ) = smoothdata(data_lfp(ind+leadTime*fs:ind+lenTrial*fs-1/fs, j),'gaussian', fs*0.001);
    %sddata(i,:) = std(temp(:,lastUsedJ));
    lastUsedJ = lastUsedJ + 1;
end
trialAverages_mua(:,j) = smoothdata(mean(temp,2), 'gaussian', fs*0.0001);

sd = std(temp,0,2);
curve1 = trialAverages_mua(:,j) + sd;
curve2 = trialAverages_mua(:,j) - sd;
t_trial = leadTime:1/fs:lenTrial-1/fs;
t_reverse = fliplr(t_trial);
t2 = [t_trial, t_reverse];

figure;
% plot(t_trial(1:1000:end), curve1(1:1000:end), 'Color', 'w', 'LineWidth', 2); hold on;
% plot(t_trial(1:1000:end), curve2(1:1000:end), 'Color', 'w', 'LineWidth', 2); hold on;
% btwn = [curve1; fliplr(curve2)];
% fill(t2(1:100:end), btwn(1:100:end), [0.8, 0.8, 0.8], 'EdgeColor', 'none');

plot(leadTime:1/fs:lenTrial-1/fs, smoothdata(mean(temp,2), 'gaussian', 1000), 'k', 'LineWidth', 2); hold on;
ylim([-125 150]);
xlim([leadTime, lenTrial]);
fill([0,5,5,0], -125+[0,0,7,7], [1 0.5 0.25], 'EdgeColor', 'none');
xlabel('Time (s)');
ylabel('Amplitude (\muV)');
set(gca, 'FontSize', 18);
%fill([2,3,3,2], [0,0,1,1], [1 0.5 0.25], 'EdgeColor', 'none');
