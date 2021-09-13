function [Electrode_CaEvents_Onset, Electrode_CaEvents_Offset, Electrode_CaEvents] = CR_GetCaEvents_WF(Electrode_df_f,frate,Plotting)
if ~iscolumn(Electrode_df_f)
    Electrode_df_f = Electrode_df_f';
end
% Plotting = true or false;
%% Start
% Loess smooth over 1 sec
Electrode_Smooth_loess = smooth(Electrode_df_f,frate/length(Electrode_df_f),'loess');
noise_df_f = abs(Electrode_df_f-Electrode_Smooth_loess);
noise_cutoff = std(Electrode_df_f-Electrode_Smooth_loess); % noise cutoff
Electrode_Smooth_loess_velocity = diff([nan;Electrode_Smooth_loess]);
if Plotting
    figure; plot(Electrode_df_f,'b'); hold on;
    plot(Electrode_Smooth_loess,'k'); hold on; plot(Electrode_Smooth_loess_velocity,'r');
    hold on; line([0 9000],[0 0],'color','k');
end
% Detect events: the first derivative (velocity) of the smoothed fluorescence trace (loess, 1 s) crossed the standard deviation of the inactive velocity trace;
% Inactive velocity trace was derived from periods when the velocity was
% within the standard deviation of the velocity trace
Threshold_inactive = nanstd(Electrode_Smooth_loess_velocity);
Electrode_Smooth_inactive = Electrode_Smooth_loess_velocity <= Threshold_inactive;
if Plotting
    temp = ones(size(Electrode_df_f)); temp(Electrode_Smooth_inactive) = nan;
    hold on; plot(temp.*Electrode_df_f,'c');
end
Threshold_active = nanstd(Electrode_Smooth_loess_velocity(Electrode_Smooth_inactive));
Electrode_Smooth_active = Electrode_Smooth_loess_velocity >= Threshold_active;
if Plotting
    temp = ones(size(Electrode_df_f)); temp(~Electrode_Smooth_active) = nan;
    hold on; plot(temp.*Electrode_df_f,'m');
end
Electrode_Smooth_active_switch = diff([0;Electrode_Smooth_active]);
% Get rough onset: when the velocity above threshold_active
temp_onset = find(Electrode_Smooth_active_switch==1);
if Plotting
    hold on
    for ii = 1:length(temp_onset)
        hold on; line([temp_onset(ii),temp_onset(ii)],[-0.05, 0.1],'color',[0.5 0.5 0.5]);
    end
end
% Get rougn offset (peak time): when the velocoity dropps below zero
for ii = 1:length(temp_onset)
    curr_onset = temp_onset(ii);
    if ~isempty(find(Electrode_Smooth_loess_velocity(curr_onset:end)<0,1,'first'))
        temp_offset(ii,1) = find(Electrode_Smooth_loess_velocity(curr_onset:end)<0,1,'first')+curr_onset-1;
    else
        temp_offset(ii,1) = nan;
    end
end
% Refine
temp_onset(isnan(temp_offset)) = [];
temp_offset(isnan(temp_offset)) = [];
if Plotting
    hold on
    for ii = 1:length(temp_offset)
        hold on; line([temp_offset(ii),temp_offset(ii)],[-0.05, 0.1],'color',[0.5 0.5 1]);
    end
end
% Fill the gap:
temp_offset_diff = diff([0;temp_offset]);
temp_onset(temp_offset_diff == 0) = [];
temp_offset(temp_offset_diff == 0) = [];
if Plotting
    hold on
    for ii = 1:length(temp_offset)
        hold on; line([temp_onset(ii),temp_onset(ii)],[-0.05, 0.1],'color',[0.3 0.3 0.3]);
        hold on; line([temp_offset(ii),temp_offset(ii)],[-0.05, 0.1],'color',[0.3 0.3 1]);
    end
end
% Refine offset (peak time): max df/f +/-5 frames around rougn offset
for ii = 1:length(temp_offset)
    curr_offset = temp_offset(ii);
    if curr_offset-5>0 && curr_offset+5 <= length(Electrode_df_f)
        [~,temp_max] = max(Electrode_df_f([curr_offset-5:curr_offset+5]));
        temp_peak_time(ii,1) = curr_offset+temp_max-6;
    else
        temp_peak_time(ii,1) = nan;
    end
end
if Plotting
    hold on
    for ii = 1:length(temp_peak_time)
        hold on; line([temp_peak_time(ii),temp_peak_time(ii)],[-0.05, 0.1],'color',[0.1 0.1 1]);
    end
end
temp_onset(isnan(temp_peak_time)) = [];
temp_offset(isnan(temp_peak_time)) = [];
temp_peak_time(isnan(temp_peak_time)) = [];

% Refine onset: find baseline: when the velocoity rises above zero before
% peak time, use temp_onffset, more robust
for ii = 1:length(temp_offset)
    curr_offset = temp_offset(ii);
    if ~isempty(find(Electrode_Smooth_loess_velocity(1:curr_offset-1)<0,1,'last'))
        baseline_index(ii,1) = find(Electrode_Smooth_loess_velocity(1:curr_offset-1)<0,1,'last')+1;
        baseline_df_f(ii,1) = Electrode_df_f(baseline_index(ii,1));
    else
        baseline_index(ii,1) = nan;
        baseline_df_f(ii,1) = nan;
    end
end
temp_onset(isnan(baseline_index)) = [];
temp_offset(isnan(baseline_index)) = [];
baseline_index(isnan(baseline_df_f)) = [];
baseline_df_f(isnan(baseline_df_f)) = [];
if Plotting
    hold on
    for ii = 1:length(baseline_index)
        hold on; plot(baseline_index(ii),baseline_df_f(ii),'r*');
    end
end
% Refine onset
for ii = 1:length(baseline_df_f)
    curr_baseline_df_f = baseline_df_f(ii);
    temp_onset(ii,1) = find(abs(Electrode_df_f([baseline_index(ii):temp_offset(ii)])-curr_baseline_df_f)<noise_cutoff,1,'last')+baseline_index(ii)-1;
end
if Plotting
    hold on
    for ii = 1:length(baseline_index)
        hold on; plot(temp_onset(ii),Electrode_df_f(temp_onset(ii)),'g*');
    end
end

Electrode_CaEvents_Onset = temp_onset;
Electrode_CaEvents_Offset = temp_offset;
Electrode_CaEvents = zeros(size(Electrode_df_f));
for ii = 1:length(temp_offset)
    Electrode_CaEvents(temp_onset(ii):temp_offset(ii)) = Electrode_df_f(temp_offset(ii))-baseline_df_f(ii);
end
if Plotting
    figure; plot(Electrode_df_f,'k'); hold on;
    plot(Electrode_CaEvents,'r'); 
end