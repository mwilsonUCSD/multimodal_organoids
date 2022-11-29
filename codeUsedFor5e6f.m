%% Read in image
close all;
clear all;

n = 1;%Choose what to plot

switch n
    case 1
        % Human puncta density
        a = [15.29377704,16.62813343,15.75566964,25.01878255,22.14547699,12.45307649];
        b = [20.67065657,14.16829771,24.22984247,18.06971303,24.43518011,29.77395896,20.80754833,23.95605893,25.11963894,23.27160011,23.88761305,23.61382952,22.08546496,19.85383261,25.54834274,28.16473929,24.06648398,18.1948418,22.97914284,17.90488416];
        name = categorical({'Boundary','Center'});
        label = '# human puncta/1000 um^3';
        lim = [0,40];
    case 2
        % Mouse puncta density
        a = [12.26581448,10.46956549,8.262745316,20.74386452,20.11313891,11.114036,];
        b = [8.145060038,5.407224731,6.297021206,8.008168273,10.6775577,9.719315339,6.297021206,8.96641063,12.86782594,11.70424594,15.53721537,15.26343184,8.157000996,7.695283958,8.157000996,8.157000996,10.29349616,3.914428117,4.349364574,4.566832803,];
        name = categorical({'Boundary','Center'});
        label = '# mouse puncta/1000 um^3';
        lim = [0,25];
    case 3
        % NM95
        a = [94.382,87.1795,72.0721];
        b = [0,0,0];
        name = categorical({'Organoid','Cortex'});
        name = reordercats(name, {'Organoid','Cortex'});
        label = 'DAB/Hematoxylin (%)';
        lim = [0,100];
    case 4
        % NeuN
        a = [49.23076923,52.77777778,46.34146341];
        b = [92.47311828,86.76470588,80.48780488];
        name = categorical({'Organoid','Cortex'});
        name = reordercats(name, {'Organoid','Cortex'});
        label = 'DAB/Hematoxylin (%)';
        lim = [0,100];
    case 5
        % Ki67
        a = [3.243243243,6.113537118,7.80141844];
        b = [0,0,0];
        name = categorical({'Organoid','Cortex'});
        name = reordercats(name, {'Organoid','Cortex'});
        label = 'DAB/Hematoxylin (%)';
        lim = [0,10];
    case 6
        % Olig2
        a = [6.153846154,11.32075472,5.681818182];
        b = [12.38095238,22.52252252,24.61538462];
        name = categorical({'Organoid','Cortex'});
        name = reordercats(name, {'Organoid','Cortex'});
        label = 'DAB/Hematoxylin (%)';
        lim = [0,30];
    otherwise
        disp('Not valid input')
end

% Student's ttest
[h, p] = ttest2(a,b);

figure;
hold on;
%name = categorical({'Edge','Organoid'});
%name = categorical({'Organoid','Cortex'});
%name = reordercats(name, {'Organoid','Cortex'});
%bar(repmat([1;2],1,3), [a;b]);
bar(name, [mean(a);mean(b)]);
errorbar(name,[mean(a);mean(b)], [std(a);std(b)], 'LineStyle', 'none', 'LineWidth', 1, 'Color', 'k');
scatter(repmat(name(1),length(a),1), a, 'ok');
scatter(repmat(name(2),length(b),1), b, 'ok');
set(gca, 'ylim', lim, 'FontSize', 15);
ylabel(label);




