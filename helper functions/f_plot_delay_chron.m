function f_plot_delay_chron(data_env, good_channels, f_bpass)
    Root = 'C:\Users\Maddie\Documents\Research\Code Database\_210614_burst_analysis\data\';
    ephys_folder = 'NOD14\';
    org_ch = [11,12,13,14];%[5,7,8,9,11,14];
    cor_ch = setdiff(1:16, org_ch);
    fs = 20000;
    t = -1:1/fs:1;
    week = [1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4];
    
    data_folder = [Root, ephys_folder];
    fnamelist = dir(data_folder);
    
    figure; hold on;
    for record_num = 3:length(fnamelist)
        filename = [data_folder, fnamelist(record_num).name];
        load(filename);
        
        [~, delay_max] = max(env_averages,[],3);
        delay_max = t(delay_max);

        [~,org_ch_mapped] = ismember(org_ch, good_channels);
        org_ch_mapped(org_ch_mapped == 0) = [];
        [~,cor_ch_mapped] = ismember(cor_ch, good_channels);
        cor_ch_mapped(cor_ch_mapped == 0) = [];

        subplot(2,3,mod(record_num-3,6)+1); hold on;
        %Organoid -> cortex
        for chO = 1:length(org_ch_mapped)
            for chC = 1:length(cor_ch_mapped)
                %diff = delay_max(chO,chO);
                scatter(week(record_num-2),delay_max(org_ch_mapped(chO),cor_ch_mapped(chC)), 'filled', 'r');
            end
        end

        %Cortex -> organoid
        for chC = 1:length(cor_ch_mapped)
            for chO = 1:length(org_ch_mapped)
                %diff = delay_max(chO,chO);
                scatter(week(record_num-2),delay_max(cor_ch_mapped(chC),org_ch_mapped(chO)), 'filled', 'b');
            end
        end
        
    end
    
    % Add axis labels
    for i = 1:6
        subplot(2,3,i);
        xlim([0,8]);
        ylim([-0.2, 0.2]);
        xlabel('Time (weeks)');
        ylabel('Delay time (s)');
        set(gca, 'FontSize', 15);
    end
    sgtitle(ephys_folder, 'Interpreter', 'none');
    
    
    %load('0217_envAvg_5to7Hz.mat');
%     figure;
%     for i = 1:length(good_channels)
%         subplot(4,4,find(mapping == good_channels(i)));
%         t = dt(1):1/fs:dt(2);
%         plot(t,squeeze(env_averages(1,i,:)));
%     end
    

    
    % Organoid -> Organoid
    
    % Cortex -> Cortex
    
    legend('Organoid -> Cortex', 'Cortex -> Organoid');
end