function decode_ROI()
%{

Inputs - 
Table of file locations

Notes - 
fitcnet will not work in Matlab versions earlier than 2021a
%}

[choiceFileName,choiceBatchPathName] = uigetfile({'drgCaImAn_LDAfsdz_choices*.m'},'Select the .m file with all the choices for analysis');


fprintf(1, ['\ndrgCaImAn_batch_pre_per_to_decode_entire_session_multi_ROI_fsdz for ' choiceFileName '\n\n']);

% tempDirName=['temp' choiceFileName(12:end-2)];

addpath(choiceBatchPathName)
BatchFile = eval(choiceFileName(1:end-2) );
BatchFile.choiceFileName=choiceFileName;
BatchFile.choiceBatchPathName=choiceBatchPathName;

all_no_ROIs=[1 2 5 15 2000];
all_no_ROI_draws=[2000 40 40  40 1];


new_no_files=BatchFile.no_files;

if isfield(BatchFile,'processing_algo')
    processing_algo=BatchFile.processing_algo;
else
    processing_algo=1;
end

if isfield(BatchFile,'suffix_out')
    suffix_out=BatchFile.suffix_out;
else
    suffix_out='_dec.mat';
end

if isfield(BatchFile,'first_file')
    first_file=BatchFile.first_file;
else
    first_file=1;
end



%Parallel batch processing for each file
all_files_present=1;
for filNum=first_file:BatchFile.no_files
    
    
    %Make sure that all the files exist
    pre_per_FileName=BatchFile.FileName_pre_per{filNum};
    if iscell(BatchFile.PathName_pre_per)
        pre_per_PathName=BatchFile.PathName_pre_per{filNum};
    else
        pre_per_PathName=BatchFile.PathName_pre_per;
    end
    if pre_per_PathName(1) == 'F'
        pre_per_PathName(1) = 'R';
    end
    
    if exist([pre_per_PathName pre_per_FileName])==0
        fprintf(1, ['Program will be terminated because file No %d, ' pre_per_FileName ' does not exist\n'],filNum);
        all_files_present=0;
    end
    
end


% if exist([handles.PathName_out handles.FileName_out])==0
%     handles_out=[];
%     ii_out=0;
% else
%     load([handles.PathName_out handles.FileName_out])
%     ii_out=handles_out.last_ii_out;
%     first_file=handles_out.last_file_processed+1; 
% end

figNo=0;
show_figures=1;

if all_files_present==1
    
    
    %Process each file separately
    for fileNo=first_file:length(BatchFile.FileName_pre_per)
        tic
        first_toc=toc;
        handles_out=[];
        ii_out=0;
        
        pre_per_PathName=BatchFile.PathName_pre_per{fileNo};
        pre_per_FileName=BatchFile.FileName_pre_per{fileNo};
        if pre_per_PathName(1) == 'F'
            pre_per_PathName(1) = 'R';
        end
        [percent_correct] = drgCaImAnFindPercentCorrect(pre_per_PathName, pre_per_FileName);
        
        if percent_correct>=80
            %Do only for proficient
            for ii_ROI_choices=1:length(all_no_ROIs)
                
                
                
                
                handles_choices.pre_per_PathName=pre_per_PathName;
                handles_choices.pre_per_FileName=pre_per_FileName;
                handles_choices.processing_algorithm=BatchFile.processing_algorithm;
                handles_choices.MLalgo_to_use=BatchFile.MLalgo_to_use;
                handles_choices.dt_p_threshold=BatchFile.dt_p_threshold;
                handles_choices.show_figures=BatchFile.show_figures;
                handles_choices.post_time=BatchFile.post_time;
                handles_choices.k_fold=BatchFile.k_fold;
                handles_choices.post_shift=BatchFile.post_shift;
                handles_choices.pre_time=BatchFile.pre_time;
                handles_choices.ii_cost=BatchFile.ii_cost;
                handles_choices.no_ROI_draws=all_no_ROI_draws(ii_ROI_choices);
                handles_choices.no_ROIs=all_no_ROIs(ii_ROI_choices);
                
                
                
                
                handles_choices.p_threshold=1.1;
                
                ii_out=ii_out+1;
                handles_out.ii_out(ii_out).handles_choices=handles_choices;
                handles_out.ii_out(ii_out).grNo=BatchFile.group(fileNo);
                handles_out.ii_out(ii_out).fileNo=fileNo;
                
                start_toc=toc;
                
                handles_out.ii_out(ii_out).handles_out=drgCaImAn_SVZ_entire_session_randomROIdrawv3(handles_choices);
                
                fprintf(1, 'Data processed for file number %d, number of ROIs= %d\n',fileNo,all_no_ROIs(ii_ROI_choices));
                fprintf(1,'Processing time for number of ROIs= %d is %d hours\n',all_no_ROIs(ii_ROI_choices),(toc-first_toc)/(60*60));
                
            end
            
            %Save output file
            handles_out.last_file_processed=fileNo;
            handles_out.last_ii_out=ii_out;
            handles_out.handles=BatchFile;
            save([pre_per_PathName pre_per_FileName(1:end-4) suffix_out],'handles_out','handles_choices','-v7.3')
            fprintf(1, 'Processing time for file number %d= %d hours\n',toc/(60*60));
            fprintf(1,'\n\n\n')
        end
    end
    
    fprintf(1, 'Total processing time %d hours\n',toc/(60*60));
end






