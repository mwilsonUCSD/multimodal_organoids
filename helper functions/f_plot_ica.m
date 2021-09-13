function f_plot_ica(Aica,ica_source_activity, good_channels, mapping, t_amplifier)
    ica_templates = nan(16,size(ica_source_activity,2));
    ica_templates(good_channels,:) = Aica;
    ica_templates = ica_templates(mapping,:);
    ica_templates = reshape(ica_templates, 4,4,[]);

    rowcol = [4,4];

    figure;
    for i = 1:size(ica_source_activity,2)
        subplot(rowcol(1),rowcol(2),i);
        h = pcolor([[ica_templates(:,:,i)',zeros(4,1)];zeros(1,5)]);
        set(h, 'EdgeColor','none');
        axis square; axis off;
        set(gca,'Ydir','reverse');
        colorbar;
    end

    figure;
    for i = 1:size(ica_source_activity,2)
        plot(t_amplifier,ica_source_activity(:,i) + 40*i); hold on;
    end
end