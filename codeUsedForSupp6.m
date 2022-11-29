% specify recording file
clc;
clear all;

thresh = 10e6;
% Plot type:
% 's' = scatter
% 'b' = boxplot
% 'l' = lineplot
% 'm' = mean +- sd
% 'e' = errorbar
plotType = 'l';

Root = 'D:\Maddie\Archive_Code Database\';
imTotal = readtable('Electrode Impedances.xlsx');

time = [7,14,21,28,51,71]; % NOD18
ephys_folder = 'Version3\NOD18\';
org = [4,5,6];
im18 = table2array(imTotal(31:36,2:17));
isGood = im18 < thresh;
isGood = isGood([1,2,3,4,5,6],:);
    
cortex = [8,9,11];%setdiff(1:16, org);

data_folder = [Root, ephys_folder];%[Root, exp_date, ephys_folder];
fnamelist = dir(data_folder);
rate_array = zeros(length(fnamelist)-2, 16);
amp_array = zeros(length(fnamelist)-2, 16);

for i = 3:length(fnamelist)
    record_num = i;
    filename = [data_folder, fnamelist(record_num).name];
    load(filename);

    rate_array(i-2,:) = spike_rate;
end
    
% Adjust for high impedance channels
rate_array = rate_array.*isGood;
rate_array(rate_array == 0) = nan;

%time = repmat(time, [1,16]);
t_org = repmat(time, [1,length(org)]);
t_cor = repmat(time, [1,length(cortex)]);
    
figure;
plot(time, abs(rate_array(:,org)), 'r'); hold on;
plot(time, abs(rate_array(:,cortex)), 'b'); hold on;
scatter(t_org, reshape(abs(rate_array(:,org)), [], 1), [], 'r', 'filled'); hold on;
scatter(t_cor, reshape(abs(rate_array(:,cortex)), [], 1), [], 'bs', 'filled'); hold on;
        

xlim([0 80]);
ylim([0 10]);
title([string(ephys_folder(10:14)) ' Spike Rates'], 'fontsize', 15);
xlabel('Days post-implantation', 'fontsize', 15);
ylabel('Spike Rate (Hz)', 'fontsize', 15);
legend('Organoid', 'Cortex');

% meanOrg = nanmedian(abs(rate_array(:,org)),2);
% [rhoO,pO] = corr(meanOrg, time', 'Type', 'Spearman', 'Rows', 'complete');
% meanCortex = nanmedian(abs(rate_array(:,cortex)),2);
% [rhoC,pC] = corr(meanCortex, time', 'Type', 'Spearman', 'Rows', 'complete');

% legend(['Organoid rho= ' num2str(rhoO) ' pval= ' num2str(pO)], ['Cortex rho= ' num2str(rhoC) ' pval= ' num2str(pC)]);


%% Plot spike rasters over time
%clc;
clear all;

ephys_folder = 'Version3\NOD18\';

imTotal = readtable('Electrode Impedances.xlsx');

org = [4,5,6];%[3,4,5,6,7,10]; % NOD18
im = table2array(imTotal(31:36,2:17)); im = im([1,2,3,4,5,6],:);% NOD18

Root = 'D:\Maddie\Archive_Code Database\';
data_folder = [Root, ephys_folder];%[Root, exp_date, ephys_folder];
fnamelist = dir(data_folder);

% Plot spike rasters
figure;
for i = 3:length(fnamelist)
    record_num = i;
    filename = [data_folder, fnamelist(record_num).name];
    load(filename);
    
    bad_channels = find(im(i-2,:) > 10e6);
    
    cor = [8,9,11];%setdiff(1:size(spikes,2), org);
    newOrder = [org, cor];
    isBad = ismember(newOrder, bad_channels);
    spikes_reorg = spikes(:,newOrder);
    spikes_reorg(:,isBad) = [];
    nOrgChannels = length(setdiff(org, bad_channels));
    
    spikes_reorg(spikes_reorg == 0) = nan;
    
    %subplot(2,4,i-2);
    t_down = t_down - t_down(1);
    subset = 10000*[[20 30]; [20 30]; [100 110]; [1 11]; [20 30]; [35 45]];
    for j = 1:size(spikes_reorg,2)
        if find(j <= nOrgChannels)
            scatter(t_down(1:100001), spikes_reorg(subset(i-2,1):subset(i-2,2),j)*(2*j) + (30*(i-3)), 'r.'); hold on;
            %scatter(t_down, spikes_reorg(:,j)*(j+1) + (40*(i-3)), 'r.'); hold on;
        else
            scatter(t_down(1:100001), spikes_reorg(subset(i-2,1):subset(i-2,2),j)*(2*j) + (30*(i-3)), 'b.'); hold on;
            %scatter(t_down, spikes_reorg(:,j)*(j+1) + (40*(i-3)), 'b.'); hold on;
        end
    end
%     spikes_sum_org = sum(spikes(:,org),2);
%     spikes_sum_cor = sum(spikes(:,cor),2);
%     plot(t_down, movsum(spikes_sum_org,10), 'r'); hold on;
%     plot(t_down, movsum(spikes_sum_cor,10), 'b');
%     
    xlim([0 10]);
    ylim([0 180]);
    xlabel('Time (seconds)');
    ylabel('Spike Raster');
    %title(fnamelist(record_num).name, 'Interpreter', 'none');
    disp(['#org: ' num2str(nOrgChannels) ' #cortex: ' num2str(size(spikes_reorg,2)-nOrgChannels)]);
end
sgtitle(string(ephys_folder(10:14)));


