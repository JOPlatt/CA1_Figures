function handles_outd=decoding_Multi_ROI(handles_choices)
%drgCaImAnInspectDecodingSingleROI
%This function inspects the decoding results of multi ROI decoding
%When run with show_figures=1 it plots results for a couple of ROIs

time_windows=[3.1 4.1];
time_windows_pre=[-1 0];
time_window_lat=[-1.5 10];
% pre_time_window=[-7 -1.5];
MLalgo=6;

%below needs to be edited
ii_out=handles_choices.ii_out;
pre_perFileName=handles_choices.pre_perFileName;
pre_perPathName=handles_choices.pre_perPathName;
% acc_thr=[0.35 0.65];
% dt_span=15;
%

  %This dt_lat is used to find the maximum accuracy attained within dt_lat
%after latency start
sustained_dt=0.2;
lat_fact=1;
mad_pre_accuracy=0.05;

%
% ROI1=1;
% time_to_inspect1=4; %Time to inspect the dFF
%
% remove below once NWB addresses are added
if pre_perPathName(1:5) == "/data"
    pre_perPathName = append('R:/SFTP/Ming Ma/', pre_perPathName(12:end));
end
Val_out = load([pre_perPathName pre_perFileName],'handles_out');



handles_out2= Val_out.handles_out.ii_out(ii_out).handles_out;
no_ROI_draws= handles_out2.no_ROI_draws;
time_span= handles_out2.time_span;
dFF_per_trial_sm= handles_out2.dFF_per_trial_sm;
dFF_per_trial_sp= handles_out2.dFF_per_trial_sp;


%Calculate accuracy in the time window
accuracy_per_ROI = nan(1,no_ROI_draws);
for iiROI=1:no_ROI_draws
    accuracy_per_ROI(1,iiROI) = mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2));
end
handles_outd.accuracy_per_ROI=accuracy_per_ROI;
clear iiROI

accuracy_per_ROI_sh = nan(1,no_ROI_draws);
for iiROI=1:no_ROI_draws
    accuracy_per_ROI_sh(1,iiROI) = mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh(:,(time_span>=time_windows(1))&(time_span<=time_windows(2))),2));
end
handles_outd.accuracy_per_ROI_sh = accuracy_per_ROI_sh;
clear iiROI

%Calculate accuracy in the pre time window
accuracy_per_ROI_pre = nan(1,no_ROI_draws);
for iiROI=1:no_ROI_draws
    accuracy_per_ROI_pre(1,iiROI) = mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_windows_pre(1))&(time_span<=time_windows_pre(2))),2));
end
handles_outd.accuracy_per_ROI_pre = accuracy_per_ROI_pre;
clear iiROI

accuracy_per_ROI_sh_pre = nan(1,no_ROI_draws);
for iiROI=1:no_ROI_draws
    accuracy_per_ROI_sh_pre(1,iiROI) = mean(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh(:,(time_span>=time_windows_pre(1))&(time_span<=time_windows_pre(2))),2));
end
handles_outd.accuracy_per_ROI_sh_pre = accuracy_per_ROI_sh_pre;
clear iiROI

%Quantify outliers in the odor window

ii_t_start=find(time_span>=time_windows(1),1,'first');
ii_t_end=find(time_span<=time_windows(2),1,'last');


dFF_sp_outlier_frac_per_ROI = nan(1,no_ROI_draws);
dFF_sm_outlier_frac_per_ROI = nan(1,no_ROI_draws);

for iiROI=1:no_ROI_draws
    % *need to find the amount diff between ii_t_start and end
    stepNum = ii_t_end - ii_t_start +1;
    dFF_sp_outlier_frac_this_ROI = nan(1,stepNum);
    dFF_sm_outlier_frac_this_ROI = nan(1,stepNum);
    Stepping = 1;
    % * endit above
    for ii_t=ii_t_start:ii_t_end
        
        this_dFF_per_trial_sm=zeros(1,size(dFF_per_trial_sm,1));
        this_dFF_per_trial_sm(1,:)=dFF_per_trial_sm(:,iiROI,ii_t);
        dFF_sm_outlier_frac_this_ROI(1,Stepping) =  sum(isoutlier(this_dFF_per_trial_sm,"median",ThresholdFactor=3))/length(this_dFF_per_trial_sm);
        
        this_dFF_per_trial_sp=zeros(1,size(dFF_per_trial_sp,1));
        this_dFF_per_trial_sp(1,:)=dFF_per_trial_sp(:,iiROI,ii_t);
        dFF_sp_outlier_frac_this_ROI(1,Stepping) = sum(isoutlier(this_dFF_per_trial_sp,"median",ThresholdFactor=3))/length(this_dFF_per_trial_sp);
        Stepping = Stepping + 1;
    end
    dFF_sp_outlier_frac_per_ROI(1,iiROI)=mean(dFF_sp_outlier_frac_this_ROI);
    dFF_sm_outlier_frac_per_ROI(1,iiROI)=mean(dFF_sm_outlier_frac_this_ROI);
end
clear iiROI


%Calculate latency 
% accuracy_per_ROIw = [];
% latency_per_ROI = [];

% *need to look at what this is pulling out
% cropped_time_span_ii=find((time_span>=time_window_lat(1))&(time_span<=time_window_lat(2)));
cropped_time_span = time_span((time_span>=time_window_lat(1))&(time_span<=time_window_lat(2)));

sustained_ii=ceil(sustained_dt/(time_span(2)-time_span(1)));


% accuracy_per_ROIw = nan(1,no_ROI_draws);
latency_per_ROI = nan(1,no_ROI_draws);
for iiROI=1:no_ROI_draws
%     accuracy_per_ROIw(1,iiROI) = prctile(mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_window_lat(1))&(time_span<=time_window_lat(2))),1),95);

    this_acc_timecourse=mean(handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict(:,(time_span>=time_window_lat(1))&(time_span<=time_window_lat(2))),1);

    not_found = 1;
    ii_shift=1;
    while not_found == 1
        this_latency_delta_ii=find(this_acc_timecourse(ii_shift:end)>0.5+lat_fact*mad_pre_accuracy,1,'first');
        checkOne = ~isempty(this_latency_delta_ii) && ((ii_shift+this_latency_delta_ii+sustained_ii) <= length(this_acc_timecourse));
        if checkOne == 1
            checkTwo =  sum(this_acc_timecourse(ii_shift+this_latency_delta_ii-1:ii_shift+this_latency_delta_ii+sustained_ii)<0.5+lat_fact*mad_pre_accuracy);
            if checkTwo == 0
                not_found=0;
                latency_per_ROI(1,iiROI) =cropped_time_span(ii_shift+this_latency_delta_ii-1);
            else
                ii_shift=ii_shift+this_latency_delta_ii;
                if ii_shift>length(this_acc_timecourse)
                    not_found=0;
                    latency_per_ROI(1,iiROI) =NaN;
                end
            end
        else
            not_found = 0;
            latency_per_ROI(1,iiROI) = NaN;
        end
    end
    
      

end

handles_outd.latency_per_ROI=latency_per_ROI;
% handles_outd.accuracy_per_ROIw=accuracy_per_ROIw;
% handles_outd.mad_pre_accuracy=mad_pre_accuracy;

%Save the accuracy timecourse

for iiROI=1:no_ROI_draws
    this_correct_predict=handles_out2.ROI(iiROI).MLalgo(MLalgo).this_correct_predict;
    handles_outd.ROI(iiROI).mean_accuracy=mean(this_correct_predict);
end
handles_outd.time_span= time_span;



