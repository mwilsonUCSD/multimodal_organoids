function [A,ica_source_activity] = ica_denoise_fastica(lfp_data,t_low,good_ch,pos)
% lfp_data: N x C matrix
% t_low: N x 1 vector
% good_ch: the good channel numbers
% pos: the position of all the channels

% perform ICA to the LFP data
rng(0);
[icasig, A, W] = fastica(lfp_data(:,good_ch)','approach','symm','epsilon',1e-5);
ica_source_activity = lfp_data(:,good_ch)*W';
% compute and plot the ICA template
ica_templates = nan(16, length(good_ch));
ica_templates(pos(good_ch),:) = A;
ica_templates = reshape(ica_templates, 4, 4,[]);
[rowcol,~] = numSubplots(length(good_ch));
figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:length(good_ch)
    subplot(rowcol(1),rowcol(2),i);
    h = pcolor([[ica_templates(:,:,i)',zeros(4,1)];zeros(1,5)]);
    set(h, 'EdgeColor','none'); colormap jet;
    axis square;
    set(gca,'Ydir','reverse');
    title(sprintf('%2.d',i));
    % caxis([-5,5]);
end

figure;
for i = 1:length(good_ch)
    plot(t_low,ica_source_activity(:,i) + 40*i); hold on;
end


end