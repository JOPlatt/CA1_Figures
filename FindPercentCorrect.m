function [percent_correct] = FindPercentCorrect(trial_amt,epochs)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

idx_max = size(epochs,2);
cK = 1;

hit_per_trial= nan([1,trial_amt]);
miss_per_trial= nan([1,trial_amt]);
cr_per_trial= nan([1,trial_amt]);
fa_per_trial= nan([1,trial_amt]);
at_end=0;
this_ii=0;


%training_decisions is 1 if S+ and 2 if S-
%epochs has masks for the following epochs
% 1 - FV on
% 2 - odor on
% 3 - odor off
% 4 - reinforcement on
% 5 - reinforcement off
% 6 - Hit
% 7 - Miss
% 8 - FA
% 9 - CR
while (at_end==0)
    %6 is hit, 7 is miss, 8 FA and 9 CR
    next_ii_sp=find((epochs(this_ii+1:idx_max)==7)|(epochs(this_ii+1:idx_max)==6),1,'first');
    next_ii_sm=find((epochs(this_ii+1:idx_max)==8)|(epochs(this_ii+1:idx_max)==9),1,'first');
    
    if (isempty(next_ii_sp))&&(isempty(next_ii_sm))
        at_end=1;
    else
        
        if isempty(next_ii_sm)
            %This is S+
            if epochs(this_ii+1+next_ii_sp)==6
                hit_per_trial(cK) = 1;
                miss_per_trial(cK) = 0;
            else
                hit_per_trial(cK) =  0;
                miss_per_trial(cK) = 1;
            end
            cr_per_trial(cK) = 0;
            fa_per_trial(cK) = 0;
            
            delta_next_ii_sp=find(epochs(this_ii+next_ii_sp:idx_max)~=epochs(this_ii+1+next_ii_sp),1,'first');
            this_ii=this_ii+next_ii_sp+delta_next_ii_sp;
            
        end
        
        if isempty(next_ii_sp)
            %This is S-
            if epochs(this_ii+1+next_ii_sm)==9
                cr_per_trial(cK) = 1;
                fa_per_trial(cK) = 0;
            else
                cr_per_trial(cK) = 0;
                fa_per_trial(cK) = 1;
            end
            hit_per_trial(cK) = 0;
            miss_per_trial(cK) = 0;
            
            delta_next_ii_sm=find(epochs(this_ii+next_ii_sm:idx_max)~=epochs(this_ii+1+next_ii_sm),1,'first');
            this_ii=this_ii+next_ii_sm+delta_next_ii_sm;
            
            
        end
        
        if (~isempty(next_ii_sp))&&(~isempty(next_ii_sm))
            if next_ii_sm<next_ii_sp
                %This is S-
%                 next_ii=next_ii_sm;
                
                if epochs(this_ii+1+next_ii_sm)==9
                    cr_per_trial(cK) = 1;
                    fa_per_trial(cK) = 0;
                else
                    cr_per_trial(cK) = 0;
                    fa_per_trial(cK) = 1;
                end
                hit_per_trial(cK) = 0;
                miss_per_trial(cK) = 0;
                
                delta_next_ii_sm=find(epochs(this_ii+next_ii_sm:idx_max)~=epochs(this_ii+1+next_ii_sm),1,'first');
                this_ii=this_ii+next_ii_sm+delta_next_ii_sm;
                
                
            else
                %This is S+
%                 next_ii=next_ii_sp;
                
                if epochs(this_ii+1+next_ii_sp)==6
                    hit_per_trial(cK) = 1;
                    miss_per_trial(cK) = 0;
                else
                    hit_per_trial(cK) = 0;
                    miss_per_trial(cK) = 1;
                end
                cr_per_trial(cK) = 0;
                fa_per_trial(cK) = 0;
                
                delta_next_ii_sp=find(epochs(this_ii+next_ii_sp:idx_max)~=epochs(this_ii+1+next_ii_sp),1,'first');
                this_ii=this_ii+next_ii_sp+delta_next_ii_sp;
                
            end
        end
        
        
    end
    cK = cK +1;
end

%Calculate percent correct
percent_correct=100*(sum(hit_per_trial)+sum(cr_per_trial))/length(hit_per_trial);
end

