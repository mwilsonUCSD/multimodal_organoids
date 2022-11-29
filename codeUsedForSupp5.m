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
plotType = 'm';


Root = 'D:\Maddie\Archive_Code Database\';
imTotal = readtable('Electrode Impedances.xlsx');

%f1 = figure;
f2 = figure;
for j = 1:4
    switch j
        case 1
            time = [7,14,28,58,78]; % NOD13
            ephys_folder = 'Version3\NOD13\';
            org = [11,12,13,14];
            im13 = table2array(imTotal(1:6,2:17));
            isGood = im13 < thresh;
            isGood = isGood([1,2,4,5,6],:);
        case 2
            time = [13,27,57]; %NOD15
            ephys_folder = 'Version3\NOD15\';
            org = [3,4,5,6,7];
            im15 = table2array(imTotal(16:19,2:17));
            isGood = im15 < thresh;
            isGood = isGood([1,3,4],:);
        case 3
            time = [7,14,21,28,51,71]; % NOD18
            ephys_folder = 'Version3\NOD18\';
            org = [3,4,5,6,7,10];
            im18 = table2array(imTotal(31:36,2:17));
            isGood = im18 < thresh;
            isGood = isGood([1,2,3,4,5,6],:);
        case 4
            time = [6,13,20,27,50,70]; % NOD20
            ephys_folder = 'Version3\NOD20\';
            org = [5,6,7,8,9,10];
            im20 = table2array(imTotal(39:44,2:17));
            isGood = im20 < thresh;
            isGood = isGood([1,2,3,4,5,6],:);
        otherwise
            disp('error');
    end
    
    cortex = setdiff(1:16, org);
    
    data_folder = [Root, ephys_folder];%[Root, exp_date, ephys_folder];
    fnamelist = dir(data_folder);
    rate_array = zeros(length(fnamelist)-2, 16);
    amp_array = zeros(length(fnamelist)-2, 16);

    for i = 3:length(fnamelist)
        record_num = i;
        filename = [data_folder, fnamelist(record_num).name];
        load(filename);

        rate_array(i-2,:) = spike_rate;
        %amp_array(i-2,:) = amp;
    end
    
    % Adjust for high impedance channels
    rate_array = rate_array.*isGood;
    rate_array(rate_array == 0) = nan;
    
    %time = repmat(time, [1,16]);
    t_org = repmat(time, [1,length(org)]);
    t_cor = repmat(time, [1,length(cortex)]);
    
    figure(f2);
    subplot(2,2,j);
    switch plotType
        case 'b'
            boxplot(reshape(abs(rate_array(:,org)), [], 1), t_org, 'Colors', 'r'); hold on;
            boxplot(reshape(abs(rate_array(:,cortex)), [], 1), t_cor, 'Colors', 'b'); hold on;
        case 's'
            scatter(t_org, reshape(abs(rate_array(:,org)), [], 1), [], 'r'); hold on;
            scatter(t_cor, reshape(abs(rate_array(:,cortex)), [], 1), [], 'b'); hold on;
        case 'l'
            plot(time, abs(rate_array(:,org)), 'r'); hold on;
            plot(time, abs(rate_array(:,cortex)), 'b'); hold on;
        case 'm'
            t2 = [time, fliplr(time)];
            mean_org = nanmean(abs(rate_array(:,org)),2);
            sd_org = nanstd(abs(rate_array(:,org)),0,2);
            curve1 = mean_org + sd_org;
            curve2 = mean_org - sd_org;
            btwn = [curve1', fliplr(curve2')];
            fill(t2, btwn, [1, 0.5, 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.3); hold on;
            plot(time, mean_org, '-o', 'MarkerFaceColor', 'r', 'Color', 'r', 'LineWidth', 2); hold on;
        
            mean_cor = nanmean(abs(rate_array(:,cortex)),2);
            sd_cor = nanstd(abs(rate_array(:,cortex)),0,2);
            curve1 = mean_cor + sd_cor;
            curve2 = mean_cor - sd_cor;
            btwn = [curve1', fliplr(curve2')];
            fill(t2, btwn, [0.5, 0.5, 1], 'EdgeColor', 'none', 'FaceAlpha', 0.3); hold on;
            plot(time, mean_cor, '-o', 'MarkerFaceColor', 'b', 'Color', 'b', 'LineWidth', 2); hold on;
        case 'e'
            mean_org = nanmean(abs(rate_array(:,org)),2);
            sd_org = nanstd(abs(rate_array(:,org)),0,2);
            errorbar(time, mean_org, sd_org, 'r', 'LineWidth', 2); hold on;
            mean_cor = nanmean(abs(rate_array(:,cortex)),2);
            sd_cor = nanstd(abs(rate_array(:,cortex)),0,2);
            errorbar(time, mean_cor, sd_cor, 'b', 'LineWidth', 2);
    end
    
    xlim([0 80]);
    ylim([0 10]);
    title(string(ephys_folder(10:14)), 'fontsize', 15);
    xlabel('Days post-implantation', 'fontsize', 15);
    ylabel('Spike Rate (Hz)', 'fontsize', 15);
    
    meanOrg = nanmedian(abs(rate_array(:,org)),2);
    [rhoO,pO] = corr(meanOrg, time', 'Type', 'Spearman', 'Rows', 'complete');
    meanCortex = nanmedian(abs(rate_array(:,cortex)),2);
    [rhoC,pC] = corr(meanCortex, time', 'Type', 'Spearman', 'Rows', 'complete');
    
    if plotType == 'm'
        legend('', ['Organoid rho= ' num2str(rhoO) ' pval= ' num2str(pO)], '', ['Cortex rho= ' num2str(rhoC) ' pval= ' num2str(pC)]);
    else
        legend(['Organoid rho= ' num2str(rhoO) ' pval= ' num2str(pO)], ['Cortex rho= ' num2str(rhoC) ' pval= ' num2str(pC)]);
    end
    
end

figure(f2);
sgtitle('Spike Rates');




