function [A,ica_source_activity] = ica_denoise(ripple_env,t_low,mapping)

B = jadeR(ripple_env'); % The inverse mapping matrix
ica_source_activity = (B*ripple_env')'; % The value for demixed sources
A = inv(B);                       % The forward matrix

ica_templates = nan(16,size(ripple_env,2));
ica_templates(mapping,:) = A;
ica_templates = reshape(ica_templates, 4,4,[]);

%[rowcol,~] = numSubplots(size(ripple_env,2));
rowcol = [4,4];

% figure;
% for i = 1:size(ripple_env,2)
%     subplot(rowcol(1),rowcol(2),i);
%     h = pcolor([[ica_templates(:,:,i)',zeros(4,1)];zeros(1,5)]);
%     set(h, 'EdgeColor','none');
%     axis square; axis off;
%     set(gca,'Ydir','reverse');
%     colorbar;
% end
% 
% figure;
% for i = 1:size(ripple_env,2)
%     plot(t_low,ica_source_activity(:,i) + 40*i); hold on;
% end

%ECoG_denoise = ica_source_activity(:,[1,3:11]) * A(:,[1,3:11])';

end