function [include] = loadFast (PathName,FileName)
[epochs, epoch_per_trial] = load([PathName FileName],'epochs','epoch_per_trial');
trial_amt = size(epoch_per_trial,2);
%
%determining the percent correct for a given session
[AllSessions.percent_correct.value(fileNo)] = FindPercentCorrect(trial_amt, epochs);

if AllSessions.percent_correct.value(fileNo) > 80
    AllSessions.percent_correct.include(fileNo) = 1;
else
    AllSessions.percent_correct.include(fileNo) = 0;
end
end

