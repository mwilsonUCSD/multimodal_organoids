function f_analyze_anesthesia(nickname, ica_remove, time_range_orig)
    % Author: Madison Wilson

    %nickname = 'NOD20\0302_run1';
    load(['C:\Users\Maddie\Documents\Martin Exps\Processed Data\' nickname '.mat']);

    % Gather time & frequency info from RHD file, map channels to array, and filter raw data into LFP and MUA
    fs = rhd.frequency_parameters.amplifier_sample_rate;
    impedance = [rhd.amplifier_channels.electrode_impedance_magnitude];

    % mapping of array channels (Gold facing down onto cortex)
    % [4, 3,  1,  16]
    % [5, 6,  2,  15]
    % [7, 10, 14, 13]
    % [8, 9,  11, 12]
    % mapping = [4,3,1,16,5,6,2,15,7,10,14,13,8,9,11,12];
    % 
    impedance_mapped = impedance(9:24); %cut out unused channels
    good_channels = find(impedance_mapped < 10e6); % consider channels below an impedance threshold "good channels"
    % impedance_mapped = impedance_mapped(good_channels);
    % 
    % data_good = rhd.amplifier_data(9:24,:); % cut out unused channels
    % data_good = data_good(good_channels, :); % reorder channels to match array mapping

    % Create time axes for data and stimuli
    t_amplifier = rhd.t_amplifier;

    if isfield(rhd, 'board_adc_data')
        disp("ADC found");
        ind = find(diff(rhd.board_adc_data(1,:) > 0.3) == 1) + 1;
        t_stim = t_amplifier(ind);
    end


    % plot raw data
    % figure;
    % for i = 1:size(data_good,1)
    %     plot(t_amplifier, data_good(i,1:end)+500*i); hold on;
    % end
    % if isfield(rhd, 'board_adc_data')
    %     for i = 1:length(t_stim)
    %         plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
    %     end
    % end
    % 
    % title([nickname ' Raw Data'], 'Interpreter', 'none');
    % ylabel('Voltage (\muV)');
    % xlabel('Time (s)');


    % Evaluate ICA
    % f_plot_ica(Aica,ica_source_activity, good_channels, mapping, t_amplifier);
    % if isfield(rhd, 'board_adc_data')
    %     for i = 1:length(t_stim)
    %         plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1); hold on;
    %     end
    % end


    % Remove noisy ICA components
    %ica_remove = [1];
    good_ind = setdiff(1:length(good_channels), ica_remove); % second input is channel to remove based on imaging artifact
    ECoG_denoise = ica_source_activity(:,good_ind) * Aica(:,good_ind)';


    % Create LPF and MUA data and downsample
    fs_down_lfp = 500;

    d_low = designfilt('lowpassiir','FilterOrder',8, ...
    'PassbandFrequency',250,'PassbandRipple',0.2, ...
    'SampleRate',fs);

    data_lfp = filtfilt(d_low, ECoG_denoise);

    %data_lfp_down_long = data_lfp(1:fs/fs_down:end,:);
    data_lfp_down = data_lfp(1:fs/fs_down_lfp:end,:); % fs_down Hz sampling rate
    t_down_lfp = t_amplifier(1:fs/fs_down_lfp:end);


    % Plot all LFP and MUA (good) channels
    figure;
    hold on;
    for i = 1:size(data_lfp_down,2)
        %plot(t_down, data_lfp_down_long(:,i)+500*i); hold on;
        plot(t_down_lfp, data_lfp_down(:,i)+500*i);
    end
    if isfield(rhd, 'board_adc_data')
        for i = 1:length(t_stim)
            plot([t_stim(i) t_stim(i)],get(gca,'ylim'),'k','linew',1);
        end
    end
    title(['LFP for ' nickname ' ICA: [' num2str(ica_remove) ']'], 'Interpreter', 'none');
    ylabel('Voltage (\muV)');
    xlabel('Time (s)');

    disp(['Calculating Coherence for ' nickname '...']);
    
    % Coherency
    % User params
    %time_range = [0 60] *fs_down_lfp;
    time_range = time_range_orig *fs_down_lfp;

    params.tapers = [10 19];
    params.pad = 0;
    params.Fs = fs_down_lfp;
    params.fpass = [0 150];
    params.err = [1 0.05];
    params.trialave = 0;

    if time_range(1) == 0
        time_range(1) = 1;
    end
    time_range = time_range(1):time_range(2);

    num = length(good_channels);

    % Instantiate heat maps
    delta_heat = zeros(num,num);
    ab_heat = zeros(num,num);
    lowgamma_heat = zeros(num,num);
    highgamma_heat = zeros(num,num);

    % Instantiate frequency bands
    N = length(time_range);
    nfft=max(2^(nextpow2(N)+params.pad),N);
    [f,~]=getfgrid(params.Fs,nfft,params.fpass);

    delta_band = find(f>=2 & f<=4);
    ab_band = find(f>=7 & f<=30);
    lowgamma_band = find(f>=31 & f<=59);
    highgamma_band = find(f>=61 & f<=100);

    figure;
    for i = 1:num
        for j = 1:i
            [C,phi,S12,S1,S2,f,confC,phistd] = coherencyc(data_lfp_down(time_range,j), data_lfp_down(time_range,i), params);

            delta_heat(i,j) = mean(C(delta_band));
            ab_heat(i,j) = mean(C(ab_band));
            lowgamma_heat(i,j) = mean(C(lowgamma_band));
            highgamma_heat(i,j) = mean(C(highgamma_band));

            subplot(num, num, (i-1)*num + j);
            %plot_vector(C,f,'n',confC,colors(mod(i,length(colors))+1));
            plot_vector(C,f,'n',confC,'k');
            set(gca, 'YLabel', []);
            set(gca, 'XLabel', []);
            title(['i = ' num2str(i) ' j = ' num2str(j)]);
            %plot_vector(smoothdata(C,'g',5),f,'n',[],colors(mod(i,length(colors))+1));
        end
        disp(['i = ' num2str(i)]);
    end

    % Heat map of coherency at different frequency bands
    figure;
    temp = zeros(16,16);
    temp(good_channels,good_channels) = delta_heat;
    subplot(1,4,1);
    imagesc(temp);
    title('Delta (2-4 Hz)');
    set(gca, 'FontSize', 14);

    temp = zeros(16,16);
    temp(good_channels,good_channels) = ab_heat;
    subplot(1,4,2);
    imagesc(temp);
    title('Alpha/Beta (7-30 Hz)');
    set(gca, 'FontSize', 14);

    temp = zeros(16,16);
    temp(good_channels,good_channels) = lowgamma_heat;
    subplot(1,4,3);
    imagesc(temp);
    title('Low Gamma (31-59 Hz)');
    set(gca, 'FontSize', 14);

    temp = zeros(16,16);
    temp(good_channels,good_channels) = highgamma_heat;
    subplot(1,4,4);
    imagesc(temp);
    title('High Gamma (61-100 Hz)');
    c = colorbar;
    c.Label.String = 'Coherency';
    set(gca, 'FontSize', 14);

    sgtitle(nickname, 'Interpreter', 'none');
    set(gcf, 'Position', [16,275,1890,400]);
    
    saveas(gcf, ['_210608_coherence_lfp_analysis\' nickname '_heat_' num2str(time_range_orig(1)) 'to' num2str(time_range_orig(2)) 's'], 'jpeg');

end













