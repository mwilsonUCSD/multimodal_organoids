%% Code for 1b
% Read in EIS data from Excel file
z_table = readtable('20_10_17 EIS.xlsx', 'Sheet', 'Sheet1');
y_table = readtable('20_10_17 EIS.xlsx', 'Sheet', 'Sheet2');
x = table2array(y_table(2:52,1));

% Get rid of weird bumps in Gamry plots
x(8) = [];
x(37) = [];

% Plot EIS results
z2 = str2double(table2array(z_table(2:52,21:36)));
y2 = str2double(table2array(y_table(2:52,21:36)));
% Get rid of bad trials
z2(:,1) = [];
y2(:,1) = [];
z2(:,4) = [];
y2(:,4) = [];

z2(8,:) = [];
y2(8,:) = [];
z2(37,:) = [];
y2(37,:) = [];


figure;
yyaxis('right');
semilogx(x,y2, 'r-');
ylabel('Y2 - Zphz (deg)');

yyaxis('left');
loglog(x,z2, 'b-');
ylabel('Zmod (Ohm)');
xlabel('Frequency (Hz)');
title('Array 2');

%% Code for 1c
clear all;
clc;

% Plot bar plots of channel impedances
imTotal = readtable('Electrode Impedances.xlsx');
a= [13, 14, 15, 16, 18, 20];

imStart = table2array(imTotal(46:51,2:17));
imStart(imStart > 10e3) = 0;

figure;
bar(1:16,imStart(1,:));
ylim([0 2500]);
xlabel('Electrode');




