%% Plot impedances vs. time
clc;
clear all;

imTotal = readtable('Electrode Impedances.xlsx');

im13 = table2array(imTotal(1:6,2:17));
im13(im13 > 4e6) = nan;
im14 = table2array(imTotal(9:13,2:17));
im14(im14 > 4e6) = nan;
im15 = table2array(imTotal(16:19,2:17));
im15(im15 > 4e6) = nan;
im16 = table2array(imTotal(22:27,2:17));
im16(im16 > 4e6) = nan;
im18 = table2array(imTotal(31:36,2:17));
im18(im18 > 4e6) = nan;
im20 = table2array(imTotal(39:44,2:17));
im20(im20 > 4e6) = nan;

days13 = [7,14,21,28,58,78];
days14 = [7,14,21,28,78];
days15 = [13,20,27,57];
days16 = [6,13,20,27,34,57];
days18 = [7,14,21,28,51,71];
days20 = [6,13,20,27,50,70];

orgCh13 = [11,12,13,14]; vesCh13 = [1,3,4,15,16]; corCh13 = setdiff(1:16, [orgCh13, vesCh13]);
orgCh14 = [5,7,8,9,11,13,14]; vesCh14 = [1,15,16]; corCh14 = setdiff(1:16, [orgCh14, vesCh14]);
orgCh15 = [3,4,5,6,7]; vesCh15 = [1,13,15,16]; corCh15 = setdiff(1:16, [orgCh15, vesCh15]);
orgCh16 = [11,12,13,14]; vesCh16 = []; corCh16 = setdiff(1:16, [orgCh16, vesCh16]);
orgCh18 = [3,4,5,6,7,10]; vesCh18 = [1,15,16]; corCh18 = setdiff(1:16, [orgCh18, vesCh18]);
orgCh20 = [5,6,7,8,9,10]; vesCh20 = [1,3,15,16]; corCh20 = setdiff(1:16, [orgCh20, vesCh20]);

comb13 = [orgCh13, corCh13];
comb14 = [orgCh14, corCh14];
comb15 = [orgCh15, corCh15];
comb16 = [orgCh16, corCh16];
comb18 = [orgCh18, corCh18];
comb20 = [orgCh20, corCh20];

sd13 = std(im13);
sd14 = std(im14);
sd15 = std(im15);
sd16 = std(im16);
sd18 = std(im18);
sd20 = std(im20);

% Combined, different plots, fill method
figure;
t2 = [days13, fliplr(days13)];
mean_org = nanmean(im13(:,comb13),2)/1e6;
sd_org = nanstd(im13(:,comb13),[],2)/1e6;
curve1 = mean_org + sd_org; curve2 = mean_org - sd_org; btwn = [curve1', fliplr(curve2')];
subplot(2,3,1); fill(t2, btwn, [0.9, 0.9, 0.9], 'EdgeColor', 'none'); hold on; plot(days13, mean_org, '-o', 'MarkerFaceColor', 'k', 'Color', 'k', 'LineWidth', 2);
xlim([1,80]); ylim([0 8]); ylabel('Impedance (MOhm)'); xlabel('DPI'); title('NOD 13'); set(gca, 'FontSize', 15);
t2 = [days14, fliplr(days14)];
mean_org = nanmean(im14(:,comb14),2)/1e6;
sd_org = nanstd(im14(:,comb14),[],2)/1e6;
curve1 = mean_org + sd_org; curve2 = mean_org - sd_org; btwn = [curve1', fliplr(curve2')];
subplot(2,3,2); fill(t2, btwn, [0.9, 0.9, 0.9], 'EdgeColor', 'none'); hold on; plot(days14, mean_org, '-o', 'MarkerFaceColor', 'k', 'Color', 'k', 'LineWidth', 2);
xlim([1,80]); ylim([0 8]); ylabel('Impedance (MOhm)'); xlabel('DPI'); title('NOD 14'); set(gca, 'FontSize', 15);
t2 = [days15, fliplr(days15)];
mean_org = nanmean(im15(:,comb15),2)/1e6;
sd_org = nanstd(im15(:,comb15),[],2)/1e6;
curve1 = mean_org + sd_org; curve2 = mean_org - sd_org; btwn = [curve1', fliplr(curve2')];
subplot(2,3,3); fill(t2, btwn, [0.9, 0.9, 0.9], 'EdgeColor', 'none'); hold on; plot(days15, mean_org, '-o', 'MarkerFaceColor', 'k', 'Color', 'k', 'LineWidth', 2);
xlim([1,80]); ylim([0 8]); ylabel('Impedance (MOhm)'); xlabel('DPI'); title('NOD 15'); set(gca, 'FontSize', 15);
t2 = [days16, fliplr(days16)];
mean_org = nanmean(im16(:,comb16),2)/1e6;
sd_org = nanstd(im16(:,comb16),[],2)/1e6;
curve1 = mean_org + sd_org; curve2 = mean_org - sd_org; btwn = [curve1', fliplr(curve2')];
subplot(2,3,4); fill(t2, btwn, [0.9, 0.9, 0.9], 'EdgeColor', 'none'); hold on; plot(days16, mean_org, '-o', 'MarkerFaceColor', 'k', 'Color', 'k', 'LineWidth', 2);
xlim([1,80]); ylim([0 8]); ylabel('Impedance (MOhm)'); xlabel('DPI'); title('NOD 16'); set(gca, 'FontSize', 15);
t2 = [days18, fliplr(days18)];
mean_org = nanmean(im18(:,comb18),2)/1e6;
sd_org = nanstd(im18(:,comb18),[],2)/1e6;
curve1 = mean_org + sd_org; curve2 = mean_org - sd_org; btwn = [curve1', fliplr(curve2')];
subplot(2,3,5); fill(t2, btwn, [0.9, 0.9, 0.9], 'EdgeColor', 'none'); hold on; plot(days18, mean_org, '-o', 'MarkerFaceColor', 'k', 'Color', 'k', 'LineWidth', 2);
xlim([1,80]); ylim([0 8]); ylabel('Impedance (MOhm)'); xlabel('DPI'); title('NOD 18'); set(gca, 'FontSize', 15);
t2 = [days20, fliplr(days20)];
mean_org = nanmean(im20(:,comb20),2)/1e6;
sd_org = nanstd(im20(:,comb20),[],2)/1e6;
curve1 = mean_org + sd_org; curve2 = mean_org - sd_org; btwn = [curve1', fliplr(curve2')];
subplot(2,3,6); fill(t2, btwn, [0.9, 0.9, 0.9], 'EdgeColor', 'none'); hold on; plot(days20, mean_org, '-o', 'MarkerFaceColor', 'k', 'Color', 'k', 'LineWidth', 2);
xlim([1,80]); ylim([0 8]); ylabel('Impedance (MOhm)'); xlabel('DPI'); title('NOD 20'); set(gca, 'FontSize', 15);





