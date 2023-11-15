function [FigEight_choices] = SVZ_ROI_draw(FigEight_choices)
%{
Purpose -
Trains several decoding algorithms with the post odorant and then 
determines what happens throughout the entire timecouse

Inputs - 


Outputs -
Processed ROI with statistics

%}
% below can be removed once NWB addresses are added
pre_perPathName =FigEight_choices.pre_per_PathName;
pre_perFileName =FigEight_choices.pre_per_FileName;
processing_algorithm =FigEight_choices.processing_algorithm;
show_figures =FigEight_choices.show_figures;
%


post_time =FigEight_choices.post_time;
k_fold =FigEight_choices.k_fold;
post_shift =FigEight_choices.post_shift;
MLalgo_to_use =FigEight_choices.MLalgo_to_use;
pre_time =FigEight_choices.pre_time;
p_threshold =FigEight_choices.p_threshold;
dt_p_threshold =FigEight_choices.dt_p_threshold;

ii_cost =FigEight_choices.ii_cost;
no_ROI_draws =FigEight_choices.no_ROI_draws;
no_ROIs =FigEight_choices.no_ROIs;

%
% w = warning ('off','all');
% FigEight_choices.handles_choices=FigEight_choices;

% warning('off')

tic

%Restart random seeds
rng('shuffle');


convert_z=1; %Convert dFF traces to z
dt_span=15; %Seconds shown before and after odor on in the figures, 40 was used before
moving_mean_n=30; %Points used to calculate moving mean for the prediction label figure
no_shuffles=3; %Number of shuffles for per trial shuffling
window_no=2;
time_windows= ...
    [-1 0; 3.1 4.1];

delta_odor=4.127634e+00;
delta_odor_on_reinf_on=4.415787e+00;
delta_reinf=4.078266e-01;

% needs to be removed
load([pre_perPathName pre_perFileName])

if show_figures==1
    fprintf(1, ['\ndrgCaImAn_SVZ_entire_session run for ' pre_perFileName '\n\n']);
    fprintf(1, 'post_time = %d, p_threshold= %d, post_shift= %d, cost %d\n',post_time,p_threshold,post_shift, ii_cost);
end
%
%
if convert_z==1
    for trace_no=1:size(traces,1)
        traces(trace_no,:)=traces(trace_no,:)/std(traces(trace_no,:));
    end
end

% !!needs to be pulled into loaction of images being generated
classifier_names{1}='Linear Discriminant';
classifier_names{2}='Support Vector Machine';
classifier_names{3}='Naive Bayes Classifier';
classifier_names{4}='Neural Network';
classifier_names{5}='Decision tree';
classifier_names{6}='Binomial glm';
%
%

% FigEight_choices.handles_choices=FigEight_choices;
% FigEight_choices.pre_perFileName=pre_perFileName;
% FigEight_choices.pre_perPathName=pre_perPathName;
% FigEight_choices.post_time=post_time;
% FigEight_choices.k_fold=k_fold;
% FigEight_choices.post_shift=post_shift;
% FigEight_choices.pre_time=pre_time;
handles_not_out.MLalgo_to_use=MLalgo_to_use;

figNo=0;



%if the number of ROIs is larger than the number available decrease it to
%the number of ROIs
if no_ROIs>no_traces
    no_ROIs=no_traces;
    no_ROI_draws=1;
end

if no_ROIs==1
    no_ROI_draws=no_traces;
end
%epochs is a vector of the length of time that gives information on
%behavior
% 1=Final Valve
% 6=Hit (on for the duration of odor on)
% 7=Miss
% 8=FA
% 9=CR

%For example Hit||Miss shows S+ odor application times (red)
%and FA||CR gives S- (blue)

%Post points
Nall=size(traces,1);
dt=median(time(2:end)-time(1:end-1));
ii_p_threshold=ceil(dt_p_threshold/dt);
no_points_post=floor(post_time/dt);
no_points_post_shift=floor(post_shift/dt);
no_points_pre=floor(pre_time/dt);


measurements_post=[];
measurements_pre=[];
epochs_sp_post=zeros(1,length(time));
epochs_sp_pre=zeros(1,length(time));
epochs_sm_post=zeros(1,length(time));
epochs_sm_pre=zeros(1,length(time));
which_model_for_traces_loo=no_odor_trials*ones(1,size(traces,2));



%Do both S+ and S-
at_end=0;
this_ii=0;
ii_post=0;
ii_pre=0;
ii=0;
trial_no=0;
ii_sp_post=0;
ii_sm_post=0;
ii_which_model=0;
dt_post_which_model=floor(20/dt); %Points that model will be used beyond the training period
ii_span=ceil(dt_span/dt);
dFF_per_trial_sp=[];
dFF_per_trial_sm=[];
dFFs_sp_per_trial_per_ROI=[];
dFFs_sm_per_trial_per_ROI=[];
hit_per_trial=[];
miss_per_trial=[];
cr_per_trial=[];
fa_per_trial=[];

[hit_per_trial,cr_per_trial,dFFs_sp_per_trial_per_ROI,...
    dFFs_sm_per_trial_per_ROI,dFF_per_trial_sm,dFF_per_trial_sp,training_decisions_post,...
    which_model_for_traces_loo,decisions_per_trial,...
    ii_pointer_to_td,epochs_sp_post,measurements_post,...
    measurements_pre,epochs_sp_pre,ii_post,trial_no...
    ,epochs_sm_post,epochs_sm_pre] = ...
    drgCaImAn_parse_out_trials(dt, dt_span,epochs,no_points_post_shift,no_points_post,no_points_pre,traces,ii_p_threshold,no_odor_trials);


% !!Not needed can be removed Calculate percent correct
handles_out.percent_correct=100*(sum(hit_per_trial)+sum(cr_per_trial))/length(hit_per_trial);
%
%

which_model_for_traces_loo(which_model_for_traces_loo>trial_no)=trial_no;
%

%limiting ROIs to those below p_threshold set above
dFF_size = size(dFFs_sm_per_trial_per_ROI);
p_values=ones(1,dFF_size(2));
for iiROI=1:dFF_size(2)
    dFF_sm=zeros(dFF_size(1),dFF_size(3));
    dFF_sm(:,:)=dFFs_sm_per_trial_per_ROI(:,iiROI,:);
    dFF_sp=zeros(dFF_size(1),dFF_size(3));
    dFF_sp(:,:)=dFFs_sp_per_trial_per_ROI(:,iiROI,:);

    [~,p_values(iiROI)]=ttest2(mean(dFF_sp,2),mean(dFF_sm,2));
end
%
p_value_masks=[];
for iiROI=1:no_ROI_draws

    if no_ROIs~=1
        found_it=0;
        while found_it==0
            these_mask_ROIs=randperm(length(p_values));
            these_mask_ROIs=these_mask_ROIs(1:no_ROIs);
            p_value_mask=false(zeros(1,length(p_values)));
            p_value_mask(these_mask_ROIs)=true;
            found_it=1;
            for jjROI=1:size(p_value_masks,1)
                if sum(p_value_masks(jjROI,:)==p_value_mask)==length(p_values)
                    found_it=0;
                end
            end
        end
    else
        p_value_mask=false(zeros(1,length(p_values)));
        p_value_mask(1,iiROI)= true;
    end

    p_value_masks(iiROI,:)=p_value_mask;
    %     end

    FigEight_choices.p_value_masks=p_value_masks;

    %Trim the number of ROIs in all matrices
    noROIs_before_trimming=size(measurements_post,2);
    %     dFF_per_trial_sp=dFF_per_trial_sp(:,p_value_mask,:);
    %     dFF_per_trial_sm=dFF_per_trial_sm(:,p_value_mask,:);
    measurements_post_trimmed=measurements_post(:,p_value_mask);
    %     measurements_pre_trimmed=measurements_pre(:,p_value_mask);
    traces_trimmed=traces(p_value_mask,:);
    %     no_traces=size(traces,1);

    %Save odor times
    FigEight_choices.sp_times=[];
    FigEight_choices.sp_times_ii=0;
    FigEight_choices.sm_times=[];
    FigEight_choices.sm_times_ii=0;

    %For S+ and S- plot odor on and reinforcement
    for epoch=1:handles.dropcData.epochIndex
        %Epoch 2 is odor on, 3 is odor off
        plot_epoch=(handles.dropcData.epochEvent(epoch)==2)||(handles.dropcData.epochEvent(epoch)==3);
        if plot_epoch
            if handles.dropcData.epochTypeOfOdor(epoch)==handles.dropcProg.splusOdor
                FigEight_choices.sp_times_ii=FigEight_choices.sp_times_ii+1;
                FigEight_choices.sp_times(FigEight_choices.sp_times_ii)=handles.dropcData.epochTime(epoch);
                %             plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces_trimmed+2)*y_shift],...
                %                 '-r','LineWidth',1)
            else
                FigEight_choices.sm_times_ii=FigEight_choices.sm_times_ii+1;
                FigEight_choices.sm_times(FigEight_choices.sm_times_ii)=handles.dropcData.epochTime(epoch);
                %             plot([handles.dropcData.epochTime(epoch) handles.dropcData.epochTime(epoch)], [0 (no_traces_trimmed+2)*y_shift],...
                %                 '-b','LineWidth',1)
            end
        end
    end


    training_decisions_post_sh=zeros(1,length(training_decisions_post));
    training_decisions_post_sh(1,:)=training_decisions_post(randperm(length(training_decisions_post)));

    for ii=1:no_shuffles
        these_shuffled_decisions_per_trial=decisions_per_trial(randperm(length(decisions_per_trial)));
        ww=0;
        for jj=1:length(decisions_per_trial)
            training_decisions_post_sh2(ii,ww+1:ww+no_points_post)=these_shuffled_decisions_per_trial(jj)*ones(1,no_points_post);
            ww=ww+no_points_post;
        end
    end


    FigEight_choices.dt=dt;



    for MLalgo=MLalgo_to_use

        this_cost=[0 ii_cost;ii_cost 0];
        labels=[];
        timepoint_processed=[];
        correct_predict=[];
        correct_predict_shuffled=[];



        Nall_post=size(measurements_post_trimmed,1);


        %leave one trial out
        %Store the training data in a table.


        handles_not_out.MLalgo(MLalgo).models=[];
        handles_not_out.MLalgo(MLalgo).processed_succesfully=1;
        points_masked=floor(no_points_post/k_fold);
        which_model=ones(1,size(measurements_post_trimmed,1));
        no_trials=ii_post/no_points_post;
        for kk=1:no_trials


            training_mask=ones(size(measurements_post_trimmed,1),1);
            training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
            which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;


            these_training_measurements=zeros(sum(training_mask),size(measurements_post_trimmed,2));
            these_training_decisions=zeros(1,sum(training_mask));

            jj=0;
            for ii=1:size(measurements_post_trimmed,1)
                if training_mask(ii)==1
                    jj=jj+1;
                    these_training_measurements(jj,:)=measurements_post_trimmed(ii,:);
                    these_training_decisions(jj)=training_decisions_post(ii);
                end
            end

            tblTrn=[];
            tblTrn = array2table(these_training_measurements);

            %Store the decisions in Y
            Y=these_training_decisions;

            switch MLalgo
                case 1
                    try
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                    catch
                        % Error using ClassificationDiscriminant (line 380)
                        % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                        % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                    end
                case 2
                    handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                case 3
                    %The try catch was entered here because of this
                    %error
                    % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                    % A normal distribution cannot be fit for the combination of class 1 and predictor
                    % these_training_measurements107. The data has zero variance.
                    try
                        handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                    catch
                        handles_not_out.MLalgo(MLalgo).processed_succesfully=0;
                    end
                case 4
                    handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitcnet(tblTrn,Y);
                case 5
                    handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                case 6
                    handles_not_out.MLalgo(MLalgo).models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
            end
        end

        if handles_not_out.MLalgo(MLalgo).processed_succesfully==1
            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
            %Note: for some reason this did not work for net when I did this:
            % [label_traces,score] = predict(Mdl,traces');
            % [label_post,score] = predict(Mdl,measurements_post_trimmed);
            %I had to resort to the for loop:


            if MLalgo==6
                for ii=1:size(traces_trimmed,2)
                    this_time_point=zeros(1,size(traces_trimmed,1));
                    this_time_point(1,:)=traces_trimmed(:,ii);
                    [label,score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                    scores(ii,:)=score;
                    if label>0.5
                        label_traces(ii)=1;
                    else
                        label_traces(ii)=0;
                    end
                end

                for ii=1:size(measurements_post_trimmed,1)
                    this_time_point=zeros(1,size(traces_trimmed,1));
                    this_time_point(1,:)=measurements_post_trimmed(ii,:);
                    [label,score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                    scores_post(ii,:)=score;
                    if label>0.5
                        label_post(ii)=1;
                    else
                        label_post(ii)=0;
                    end
                end
            else
                for ii=1:size(traces_trimmed,2)
                    this_time_point=zeros(1,size(traces_trimmed,1));
                    this_time_point(1,:)=traces_trimmed(:,ii);
                    try
                        [label_traces(ii),score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                        scores(ii,:)=score;
                    catch
                        scores(ii,:)=rand(1,2);
                        if rand(1)>0.5
                            label_traces(ii)=1;
                        else
                            label_traces(ii)=0;
                        end
                    end

                end

                for ii=1:size(measurements_post_trimmed,1)
                    this_time_point=zeros(1,size(traces_trimmed,1));
                    this_time_point(1,:)=measurements_post_trimmed(ii,:);
                    try
                        [label_post(ii),score] = predict(handles_not_out.MLalgo(MLalgo).models(which_model(ii)).Mdl,this_time_point);
                        scores_post(ii,:)=score;
                    catch
                        scores_post(ii,:)=rand(1,2);
                        if rand(1)>0.5
                            label_post(ii)=1;
                        else
                            label_post(ii)=0;
                        end
                    end

                end
            end


            handles_not_out.MLalgo(MLalgo).label_post=label_post;
            handles_out.MLalgo(MLalgo).label_traces=label_traces;
            handles_out.MLalgo(MLalgo).scores=scores;
            handles_not_out.MLalgo(MLalgo).scores_post=scores_post;

            %Now do prediction with shuffled training decisions
            %k-fold cross validation evaluating performance with left-out data
            %Store the training data in a table.

            handles_not_out.MLalgo(MLalgo).sh_models=[];
            points_masked=floor(no_points_post/k_fold);
            which_model=ones(1,size(measurements_post_trimmed,1));
            for kk=1:no_trials


                training_mask=ones(size(measurements_post_trimmed,1),1);
                training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;


                these_training_measurements=zeros(sum(training_mask),size(measurements_post_trimmed,2));
                these_training_decisions=zeros(1,sum(training_mask));


                jj=0;
                for ii=1:size(measurements_post_trimmed,1)
                    if training_mask(ii)==1
                        jj=jj+1;
                        these_training_measurements(jj,:)=measurements_post_trimmed(ii,:);
                        these_training_decisions(jj)=training_decisions_post_sh(ii);
                    end
                end

                tblTrn=[];
                tblTrn = array2table(these_training_measurements);

                %Store the decisions in Y
                Y=these_training_decisions;

                handles_not_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=1;

                switch MLalgo
                    case 1
                        try
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                        catch
                            % Error using ClassificationDiscriminant (line 380)
                            % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                            % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                        end
                    case 2
                        handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                    case 3

                        %The try catch was entered here because of this
                        %error
                        % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                        % A normal distribution cannot be fit for the combination of class 1 and predictor
                        % these_training_measurements107. The data has zero variance.
                        try
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                        catch
                            handles_not_out.MLalgo(MLalgo).sh_models(kk).processed_succesfully=0;
                        end
                    case 4
                        handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                    case 5
                        handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                    case 6
                        handles_not_out.MLalgo(MLalgo).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                end
            end

            %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
            %Note: for some reason this did not work for net when I did this:
            % [label_traces,score] = predict(Mdl,traces');
            % [label_post,score] = predict(Mdl,measurements_post_trimmed);
            %I had to resort to the for loop:



            if MLalgo==6

                for ii=1:size(traces_trimmed,2)
                    this_time_point=zeros(1,size(traces_trimmed,1));
                    this_time_point(1,:)=traces_trimmed(:,ii);
                    [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                    scores_sh(ii,:)=score;
                    if label>0.5
                        label_traces_sh(ii)=1;
                    else
                        label_traces_sh(ii)=0;
                    end
                end

                for ii=1:size(measurements_post_trimmed,1)
                    this_time_point=zeros(1,size(traces_trimmed,1));
                    this_time_point(1,:)=measurements_post_trimmed(ii,:);
                    [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                    scores_post_sh(ii,:)=score;
                    if label>0.5
                        label_post_sh(ii)=1;
                    else
                        label_post_sh(ii)=0;
                    end
                end
            else
                label_traces_sh=zeros(1,size(traces_trimmed,2));
                for ii=1:size(traces_trimmed,2)
                    this_time_point=zeros(1,size(traces_trimmed,1));
                    this_time_point(1,:)=traces_trimmed(:,ii);
                    [label_traces_sh(ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                    scores_sh(ii,:)=score;
                end

                for ii=1:size(measurements_post_trimmed,1)
                    this_time_point=zeros(1,size(traces_trimmed,1));
                    this_time_point(1,:)=measurements_post_trimmed(ii,:);
                    [label_post_sh(ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                    scores_post_sh(ii,:)=score;
                end
            end



            handles_not_out.MLalgo(MLalgo).label_post_sh=label_post_sh;
            handles_out.MLalgo(MLalgo).label_traces_sh=label_traces_sh;
            handles_out.MLalgo(MLalgo).scores_sh=scores_sh;
            handles_not_out.MLalgo(MLalgo).scores_post_sh=scores_post_sh;

            %Now do shuffling on a per trial basis
            for ii_sh=1:no_shuffles


                %Now do prediction with shuffled training decisions
                %k-fold cross validation evaluating performance with left-out data
                %Store the training data in a table.

                handles_not_out.MLalgo(MLalgo).sh2(ii).sh_models=[];
                points_masked=floor(no_points_post/k_fold);
                which_model=ones(1,size(measurements_post_trimmed,1));
                for kk=1:no_trials


                    training_mask=ones(size(measurements_post_trimmed,1),1);
                    training_mask((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=0;
                    which_model((kk-1)*no_points_post+1:(kk-1)*no_points_post+no_points_post)=kk;

                    these_training_measurements=zeros(sum(training_mask),size(measurements_post_trimmed,2));
                    these_training_decisions=zeros(1,sum(training_mask));


                    jj=0;
                    for ii=1:size(measurements_post_trimmed,1)
                        if training_mask(ii)==1
                            jj=jj+1;
                            these_training_measurements(jj,:)=measurements_post_trimmed(ii,:);
                            these_training_decisions(jj)=training_decisions_post_sh2(ii_sh,ii);
                        end
                    end

                    tblTrn=[];
                    tblTrn = array2table(these_training_measurements);

                    %Store the decisions in Y
                    Y=these_training_decisions;

                    handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=1;

                    switch MLalgo
                        case 1
                            try
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcdiscr(tblTrn,Y,'Cost',this_cost);
                            catch
                                % Error using ClassificationDiscriminant (line 380)
                                % Predictor these_training_measurements3 has zero within-class variance. Either exclude this predictor
                                % or set 'discrimType' to 'pseudoLinear' or 'diagLinear'.
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl  = fitcdiscr(tblTrn,Y,'Cost',this_cost,'discrimType','pseudoLinear');
                            end
                        case 2
                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcsvm(tblTrn,Y,'Cost',this_cost);
                        case 3
                            %The try catch was entered here because of this
                            %error
                            % Error using ClassificationNaiveBayes/findNoDataCombos (line 347)
                            % A normal distribution cannot be fit for the combination of class 1 and predictor
                            % these_training_measurements107. The data has zero variance.
                            try
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnb(tblTrn,Y,'Cost',this_cost);
                            catch
                                handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).processed_succesfully=0;
                            end
                        case 4
                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitcnet(tblTrn,Y);
                        case 5
                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitctree(tblTrn,Y,'Cost',this_cost);
                        case 6
                            handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(kk).Mdl = fitglm(these_training_measurements,Y','Distribution','binomial');
                    end
                end

                %Predict labels for the test set. You trained Mdl using a table of data, but you can predict labels using a matrix.
                %Note: for some reason this did not work for net when I did this:
                % [label_traces,score] = predict(Mdl,traces');
                % [label_post,score] = predict(Mdl,measurements_post_trimmed);
                %I had to resort to the for loop:


                if MLalgo==6
                    for ii=1:size(traces_trimmed,2)
                        this_time_point=zeros(1,size(traces_trimmed,1));
                        this_time_point(1,:)=traces_trimmed(:,ii);
                        [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                        scores_sh2(ii_sh,ii,:)=score;
                        if label>0.5
                            label_traces_sh2(ii_sh,ii)=1;
                        else
                            label_traces_sh2(ii_sh,ii)=0;
                        end
                    end


                    for ii=1:size(measurements_post_trimmed,1)
                        this_time_point=zeros(1,size(traces_trimmed,1));
                        this_time_point(1,:)=measurements_post_trimmed(ii,:);
                        [label,score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model(ii)).Mdl,this_time_point);
                        scores_post_sh2(ii_sh,ii,:)=score;
                        if label>0.5
                            label_post_sh2(ii_sh,ii)=1;
                        else
                            label_post_sh2(ii_sh,ii)=0;
                        end
                    end
                else

                    for ii=1:size(traces_trimmed,2)
                        this_time_point=zeros(1,size(traces_trimmed,1));
                        this_time_point(1,:)=traces_trimmed(:,ii);
                        [label_traces_sh2(ii_sh,ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh2(ii_sh).sh_models(which_model_for_traces_loo(ii)).Mdl,this_time_point);
                        scores_sh2(ii_sh,ii,:)=score;
                    end

                    for ii=1:size(measurements_post_trimmed,1)
                        this_time_point=zeros(1,size(traces_trimmed,1));
                        this_time_point(1,:)=measurements_post_trimmed(ii,:);
                        [label_post_sh2(ii_sh,ii),score] = predict(handles_not_out.MLalgo(MLalgo).sh_models(which_model(ii)).Mdl,this_time_point);
                        scores_post_sh2(ii_sh,ii,:)=score;
                    end
                end

                handles_out.MLalgo(MLalgo).label_traces_sh2=label_traces_sh2;
                handles_not_out.MLalgo(MLalgo).label_post_sh2=label_post_sh2;
                handles_out.MLalgo(MLalgo).scores_sh2=scores_sh2;
                handles_not_out.MLalgo(MLalgo).scores_post_sh2=scores_post_sh2;

            end


        end



        if handles_not_out.MLalgo(MLalgo).processed_succesfully==1
            %label is the predicted label, and score is the predicted class
            %posterior probability
            for ii=1:length(training_decisions_post)
                if label_post(ii)==training_decisions_post(ii)
                    correct_predict_tr(ii)=1;
                else
                    correct_predict_tr(ii)=0;
                end
            end


            %Calculate wta for windows of no_points_post
            correct_predict_tr_wta=zeros(1,length(training_decisions_post));
            for ii=1:length(training_decisions_post)-no_points_post
                this_correct=zeros(1,no_points_post);
                for jj=1:no_points_post
                    if label_post(ii+jj-1)==training_decisions_post(ii+jj-1)
                        this_correct(jj)=1;
                    end
                end
                if sum(this_correct)>(no_points_post/2)
                    correct_predict_tr_wta(ii+floor(no_points_post/2))=1;
                else
                    correct_predict_tr_wta(ii+floor(no_points_post/2))=0;
                end
            end

            %Now do shuffled
            %label is the predicted label, and score is the predicted class
            %posterior probability
            for ii=1:length(training_decisions_post_sh)
                if label_post_sh(ii)==training_decisions_post_sh(ii)
                    correct_predict_tr_sh(ii)=1;
                else
                    correct_predict_tr_sh(ii)=0;
                end
            end

            %Calculate wta for windows of no_points_post
            correct_predict_tr_wta_sh=zeros(1,length(training_decisions_post_sh));
            for ii=1:length(training_decisions_post_sh)-no_points_post
                this_correct=zeros(1,no_points_post);
                for jj=1:no_points_post
                    if label_post_sh(ii+jj-1)==training_decisions_post_sh(ii+jj-1)
                        this_correct(jj)=1;
                    end
                end
                if sum(this_correct)>(no_points_post/2)
                    correct_predict_tr_wta_sh(ii+floor(no_points_post/2))=1;
                else
                    correct_predict_tr_wta_sh(ii+floor(no_points_post/2))=0;
                end
            end

            %Now do shuffled per trial
            correct_predict_tr_wta_sh2=zeros(no_shuffles,size(training_decisions_post_sh2,2));
            correct_predict_tr_sh2=zeros(no_shuffles,size(training_decisions_post_sh2,2));
            for ii_sh=1:no_shuffles
                %label is the predicted label, and score is the predicted class
                %posterior probability
                for ii=1:size(training_decisions_post_sh2,2)
                    if label_post_sh2(ii_sh,ii)==training_decisions_post_sh2(ii_sh,ii)
                        correct_predict_tr_sh2(ii_sh,ii)=1;
                    else
                        correct_predict_tr_sh2(ii_sh,ii)=0;
                    end
                end

                %Calculate wta for windows of no_points_post

                for ii=1:size(training_decisions_post_sh2,2)-no_points_post
                    this_correct=zeros(1,no_points_post);
                    for jj=1:no_points_post
                        if label_post_sh2(ii_sh,ii+jj-1)==training_decisions_post_sh2(ii_sh,ii+jj-1)
                            this_correct(jj)=1;
                        end
                    end
                    if sum(this_correct)>(no_points_post/2)
                        correct_predict_tr_wta_sh2(ii_sh,ii+floor(no_points_post/2))=1;
                    else
                        correct_predict_tr_wta_sh2(ii_sh,ii+floor(no_points_post/2))=0;
                    end
                end

            end



            handles_not_out.MLalgo(MLalgo).correct_predict_tr=correct_predict_tr;
            handles_not_out.MLalgo(MLalgo).correct_predict_tr_wta=correct_predict_tr_wta;
            FigEight_choices.ROI(iiROI).MLalgo(MLalgo).accuracy_tr=sum(correct_predict_tr)/length(correct_predict_tr);
            FigEight_choices.ROI(iiROI).MLalgo(MLalgo).accuracy_tr_wta=sum(correct_predict_tr_wta)/length(correct_predict_tr_wta);

            handles_not_out.MLalgo(MLalgo).correct_predict_tr_sh=correct_predict_tr_sh;
            handles_not_out.MLalgo(MLalgo).correct_predict_tr_wta_sh=correct_predict_tr_wta_sh;
            handles_out.MLalgo(MLalgo).accuracy_tr_sh=sum(correct_predict_tr_sh)/length(correct_predict_tr_sh);
            handles_out.MLalgo(MLalgo).accuracy_tr_wta_sh=sum(correct_predict_tr_wta_sh)/length(correct_predict_tr_wta_sh);

            handles_not_out.MLalgo(MLalgo).correct_predict_tr_sh2=correct_predict_tr_sh2;
            handles_not_out.MLalgo(MLalgo).correct_predict_tr_wta_sh2=correct_predict_tr_wta_sh2;
            handles_out.MLalgo(MLalgo).accuracy_tr_sh2=sum(correct_predict_tr_sh2(:))/length(correct_predict_tr_sh2(:));
            handles_out.MLalgo(MLalgo).accuracy_tr_wta_sh2=sum(correct_predict_tr_wta_sh2(:))/length(correct_predict_tr_wta_sh2(:));

            handles_not_out.MLalgo(MLalgo).mean_label_traces=mean(label_traces);
            handles_not_out.MLalgo(MLalgo).var_label_traces=var(label_traces);

            moving_mean_label_traces = movmean(label_traces,moving_mean_n);
            handles_not_out.MLalgo(MLalgo).label_traces=label_traces;

            %Now let's do the carpentry
            %             moving_mean_label_traces_sh = movmean(mean(label_traces_sh2),moving_mean_n);

            moving_mean_label_traces_sh2 = movmean(label_traces_sh2,moving_mean_n);
            handles_not_out.MLalgo(MLalgo).label_traces_sh2=label_traces_sh2;
            moving_mean_label_traces_sh = movmean(label_traces_sh,moving_mean_n);
            handles_not_out.MLalgo(MLalgo).label_traces_sh=label_traces_sh;



            FigEight_choices.time=time;

            %Calculate Shanon enthropy
            pone=sum(label_traces==1)/length(label_traces);
            pzero=sum(label_traces==0)/length(label_traces);
            handles_out.MLalgo(MLalgo).shannon_e=-pzero*log2(pzero) - pone*log2(pone);

            pone=sum(label_traces_sh==1)/length(label_traces_sh);
            pzero=sum(label_traces_sh==0)/length(label_traces_sh);
            handles_out.MLalgo(MLalgo).shannon_e_sh=-pzero*log2(pzero) - pone*log2(pone);

            for ii=1:no_shuffles
                pone=sum(label_traces_sh2(ii,:)==1)/size(label_traces_sh2,2);
                pzero=sum(label_traces_sh2(ii,:)==0)/size(label_traces_sh2,2);
                handles_out.MLalgo(MLalgo).shannon_e_sh2(ii)=-pzero*log2(pzero) - pone*log2(pone);
            end

            %Now let's do accounting and show it in a bar graph

            %post Splus
            post_label_sp=label_traces(logical(epochs_sp_post));
            points_per_cut=no_points_post;
            no_cuts=floor(length(post_label_sp)/points_per_cut);
            mean_post_label_sp=zeros(1,no_cuts);
            for ii=1:no_cuts
                mean_post_label_sp(ii)=mean(post_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            end
            %     sh_post_label_sp=sum(sh_label_traces,1)/size(sh_label_traces,1);
            %     sh_post_label_sp=sh_post_label_sp(logical(epochs_sp_post));
            %     mean_sh_post_label_sp=zeros(1,no_cuts);
            %     for ii=1:no_cuts
            %         mean_sh_post_label_sp(ii)=mean(sh_post_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            %     end

            %post Sminus
            post_label_sm=label_traces(logical(epochs_sm_post));
            points_per_cut=no_points_post;
            no_cuts=floor(length(post_label_sm)/points_per_cut);
            mean_post_label_sm=zeros(1,no_cuts);
            for ii=1:no_cuts
                mean_post_label_sm(ii)=mean(post_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            end
            %     sh_post_label_sm=sum(sh_label_traces,1)/size(sh_label_traces,1);
            %     sh_post_label_sm=sh_post_label_sm(logical(epochs_sm_post));
            %     mean_sh_post_label_sm=zeros(1,no_cuts);
            %     for ii=1:no_cuts
            %         mean_sh_post_label_sm(ii)=mean(sh_post_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            %     end


            %pre Sminus
            pre_label_sm=label_traces(logical(epochs_sm_pre));
            points_per_cut=no_points_pre;
            no_cuts=floor(length(pre_label_sm)/points_per_cut);
            mean_pre_label_sm=zeros(1,no_cuts);
            for ii=1:no_cuts
                mean_pre_label_sm(ii)=mean(pre_label_sm((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            end


            %pre Splus
            pre_label_sp=label_traces(logical(epochs_sp_pre));
            points_per_cut=no_points_pre;
            no_cuts=floor(length(pre_label_sp)/points_per_cut);
            mean_pre_label_sp=zeros(1,no_cuts);
            for ii=1:no_cuts
                mean_pre_label_sp(ii)=mean(pre_label_sp((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            end

            %     sh_pre_label=sum(sh_label_traces,1)/size(sh_label_traces,1);
            %     sh_pre_label=sh_pre_label(logical(epochs_sm_pre+epochs_sp_pre));
            %     mean_sh_pre_label=zeros(1,no_cuts);
            %     for ii=1:no_cuts
            %         mean_sh_pre_label(ii)=mean(sh_pre_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            %     end

            %all
            all_label=label_traces;
            %points_per_cut=no_points_pre;
            points_per_cut=no_points_post;
            no_cuts=floor(length(all_label)/points_per_cut);
            mean_all_label=zeros(1,no_cuts);
            for ii=1:no_cuts
                mean_all_label(ii)=mean(all_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            end
            %     sh_all_label=sum(sh_label_traces,1)/size(sh_label_traces,1);
            %     mean_sh_all_label=zeros(1,no_cuts);
            %     for ii=1:no_cuts
            %         mean_sh_all_label(ii)=mean(sh_all_label((ii-1)*points_per_cut+1:(ii-1)*points_per_cut+points_per_cut));
            %     end



            %Now show average labels for S+, S-, etc for 30 sec
            handles_not_out.MLalgo(MLalgo).mean_all_label=mean_all_label;
            handles_not_out.MLalgo(MLalgo).mean_pre_label_sm=mean_pre_label_sm;
            handles_not_out.MLalgo(MLalgo).mean_pre_label_sp=mean_pre_label_sp;
            handles_not_out.MLalgo(MLalgo).mean_post_label_sm=mean_post_label_sm;
            handles_not_out.MLalgo(MLalgo).mean_post_label_sp=mean_post_label_sp;
            handles_not_out.MLalgo(MLalgo).mean_all_label=mean_all_label;

            at_end=0;
            ii=1;
            sp_ii=0;
            per_trial_sp_timecourse=[];
            epoch_before_sp=[];
            sm_ii=0;
            per_trial_sm_timecourse=[];
            epoch_before_sm=[];
            per_trial_scores_sp=[];
            per_trial_scores_sm=[];

            last_sp_sm=-1;


            while at_end==0
                next_ii=[];
                next_ii_sp=find(epochs_sp_post(ii:end)==1,1,'first');
                next_ii_sm=find(epochs_sm_post(ii:end)==1,1,'first');
                if (~isempty(next_ii_sp))&(~isempty(next_ii_sm))
                    if next_ii_sp<next_ii_sm
                        next_ii=next_ii_sp;
                        this_sp_sm=1;
                    else
                        next_ii=next_ii_sm;
                        this_sp_sm=0;
                    end
                else
                    if ~isempty(next_ii_sp)
                        next_ii=next_ii_sp;
                        this_sp_sm=1;
                    end
                    if ~isempty(next_ii_sm)
                        next_ii=next_ii_sm;
                        this_sp_sm=0;
                    end
                end

                if ~isempty(next_ii)
                    if ((ii+next_ii-ii_span)>0)&((ii+next_ii+ii_span<length(label_traces)))
                        if this_sp_sm==1
                            sp_ii=sp_ii+1;
                            per_trial_sp_timecourse(sp_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                            per_trial_scores_sp(sp_ii,1:2,:)=scores(ii+next_ii-ii_span:ii+next_ii+ii_span,:)';
                            epoch_before_sp(sp_ii)=last_sp_sm;
                            last_sp_sm=1;
                            ii_next_post=find(epochs_sp_post(ii+next_ii:end)==0,1,'first');
                            ii=ii+next_ii+ii_next_post;
                        else
                            sm_ii=sm_ii+1;
                            per_trial_sm_timecourse(sm_ii,:)=label_traces(ii+next_ii-ii_span:ii+next_ii+ii_span);
                            per_trial_scores_sm(sm_ii,1:2,:)=scores(ii+next_ii-ii_span:ii+next_ii+ii_span,:)';
                            epoch_before_sm(sm_ii)=last_sp_sm;
                            last_sp_sm=0;
                            ii_next_post=find(epochs_sm_post(ii+next_ii:end)==0,1,'first');
                            ii=ii+next_ii+ii_next_post;
                        end
                    else
                        if  ((ii+next_ii+ii_span>length(label_traces)))
                            at_end=1;
                        else
                            ii=ii+next_ii;
                        end
                    end
                else
                    at_end=1;
                end

            end

            FigEight_choices.ROI(iiROI).MLalgo(MLalgo).per_trial_sp_timecourse=per_trial_sp_timecourse;
            FigEight_choices.ROI(iiROI).MLalgo(MLalgo).per_trial_sm_timecourse=per_trial_sm_timecourse;



            this_moving_mean_n=10;
            moving_mean_per_trial_sp_timecourse = movmean(per_trial_sp_timecourse',this_moving_mean_n)';
            moving_mean_per_trial_sm_timecourse = movmean(per_trial_sm_timecourse',this_moving_mean_n)';

            FigEight_choices.ROI(iiROI).MLalgo(MLalgo).moving_mean_per_trial_sp_timecourse=moving_mean_per_trial_sp_timecourse;
            FigEight_choices.ROI(iiROI).MLalgo(MLalgo).moving_mean_per_trial_sp_timecourse=moving_mean_per_trial_sp_timecourse;

            time_span=[0:dt:dt*size(per_trial_sp_timecourse,2)]-dt_span+dt;
            time_span=time_span(1:end-1)+post_shift;

            %Calculate correct predict
            for ii_tr=1:size(per_trial_sp_timecourse,1)
                for ii_time=1:size(per_trial_sp_timecourse,2)
                    if per_trial_sp_timecourse(ii_tr,ii_time)==1
                        this_correct_predict(ii_tr,ii_time)=1;
                    else
                        this_correct_predict(ii_tr,ii_time)=0;
                    end
                end
            end

            for ii_tr=1:size(per_trial_sm_timecourse,1)
                for ii_time=1:size(per_trial_sm_timecourse,2)
                    if per_trial_sm_timecourse(ii_tr,ii_time)==0
                        this_correct_predict(ii_tr+sp_ii,ii_time)=1;
                    else
                        this_correct_predict(ii_tr+sp_ii,ii_time)=0;
                    end
                end
            end


            %Calculate correct predict shuffled
            ii_plus=0;
            for ww=1:10
                for ii_tr=1:size(per_trial_sp_timecourse,1)
                    rand_stim=randi([0,1],1,size(per_trial_sp_timecourse,2));
                    for ii_time=1:size(per_trial_sp_timecourse,2)
                        if per_trial_sp_timecourse(ii_tr,ii_time)==rand_stim(ii_time)
                            this_correct_predict_sh(ii_tr+ii_plus,ii_time)=1;
                        else
                            this_correct_predict_sh(ii_tr+ii_plus,ii_time)=0;
                        end
                    end
                end

                ii_plus=ii_plus+sp_ii;

                for ii_tr=1:size(per_trial_sm_timecourse,1)
                    rand_stim=randi([0,1],1,size(per_trial_sp_timecourse,2));
                    for ii_time=1:size(per_trial_sm_timecourse,2)
                        if per_trial_sm_timecourse(ii_tr,ii_time)==rand_stim(ii_time)
                            this_correct_predict_sh(ii_tr+ii_plus,ii_time)=1;
                        else
                            this_correct_predict_sh(ii_tr+ii_plus,ii_time)=0;
                        end
                    end
                end
                ii_plus=ii_plus+sm_ii;
            end

            FigEight_choices.ROI(iiROI).MLalgo(MLalgo).this_correct_predict=this_correct_predict;
            FigEight_choices.ROI(iiROI).MLalgo(MLalgo).this_correct_predict_sh=this_correct_predict_sh;

            %             fprintf(1,'Accuracy for ROI No %d is %d\n', iiROI,mean(this_correct_predict(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2)))))
            %             fprintf(1,'Shuffled trial accuracy for ROI No %d is %d\n', iiROI,mean(this_correct_predict_sh(:,(time_span>=time_windows(window_no,1))&(time_span<=time_windows(window_no,2)))))

            %Calculate correct predict for the odor window



        else
            fprintf(1, [classifier_names{MLalgo} ' was not processed succesfully\n']);
        end




    end


end



%Save other variables to handles_out2
FigEight_choices.time_span=time_span;
FigEight_choices.MLalgo=MLalgo;
FigEight_choices.time_windows=time_windows;
FigEight_choices.wondow_no=window_no;
FigEight_choices.no_ROI_draws=no_ROI_draws;
FigEight_choices.dFF_per_trial_sm=dFF_per_trial_sm;
FigEight_choices.dFF_per_trial_sp=dFF_per_trial_sp;
