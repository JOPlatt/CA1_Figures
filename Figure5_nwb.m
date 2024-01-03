%{
creating the NWB file for all inputs needed for figure 5 
%}
%{
path_to_matnwb = '~/Repositories/matnwb'; % change to your own path location
addpath(genpath(pwd));
generateCore();
%}
nwb = NwbFile( ...
    'session_description', 'mouse in open exploration',...
    'identifier', 'Ca1ManuscriptFigures', ...
    'session_start_time', datetime(2018, 4, 25, 2, 30, 3, 'TimeZone', 'local'), ...
    'timestamps_reference_time', datetime(2018, 4, 25, 3, 0, 45, 'TimeZone', 'local'), ...
    'general_institution', 'University of Colorado Anschutz', ... % optional
    'general_related_publications', 'DOI:10.1016/j.neuron.2016.12.011'); % optional


%for figure 5 what is needed to be saved
%{
There are a totle of 5 ROI data sets
each set contains 70 columns
*
time_span, output 1XN, NWB NX70
accuracy_per_ROI, output 1XN, NWB NX70
accuracy_per_ROI_sh, output 1XN, NWB NX70
mean_accuracy, output 1XNXM, NWB NXMX70
accuracy_per_ROI_pre, output 1XN, NWB NX70
accuracy_per_ROI_sh_pre, output 1XN, NWB NX70
latency_per_ROI, output 1XN, NWB NX70

%}
ROI_names = ["ROInum_1", ...
    "ROInum_2", ...
    "ROInum_5", ...
    "ROInum_15", ...
    "ROInum_2000"];
totalSets = 70;
%{
creating process module for time_span
%}
tempData_time = [];
tempData_AccPerRoi = [];
tempData_AccPerRoiSh = [];
tempData_AccPerRoiPre = [];
tempData_AccPerRoiShPre = [];
tempData_LatencyPerRoi = [];
tempData_MeanAccuracy = [];
for rn = 1:size(ROI_names,2)
    index_ticker01 = 0;
    index_ticker02 = 0;
    index_ticker03 = 0;
    index_ticker04 = 0;
    index_ticker05 = 0;
    index_ticker06 = 0;
    index_ticker07 = 0;
    index_ticker08 = 0;
    index_ticker09 = 0;

    for ip = 1:totalSets
        %pulling time data
        index_ticker01 = index_ticker01 + size(FigureEight.session(ip).fig8_ROI(rn).decoding_output.time_span,2);
        tempData_time = [tempData_time; FigureEight.session(ip).fig8_ROI(rn).decoding_output.time_span']; %#ok<AGROW>
        if ip == 1
            index_datatimes01 = index_ticker01;
        else
            index_datatimes01 = [index_datatimes01; index_ticker01]; %#ok<AGROW>
        end
        %pulling accuracy per ROI data
        index_ticker02 = index_ticker02 + size(FigureEight.session(ip).fig8_ROI(rn).decoding_output.accuracy_per_ROI,2);
        tempData_AccPerRoi = [tempData_AccPerRoi; FigureEight.session(ip).fig8_ROI(rn).decoding_output.accuracy_per_ROI']; %#ok<AGROW>
        if ip == 1
            index_datatimes02 = index_ticker02;
        else
            index_datatimes02 = [index_datatimes02; index_ticker02]; %#ok<AGROW>
        end
        %pulling accuracy per ROI sh data
        index_ticker03 = index_ticker03 + size(FigureEight.session(ip).fig8_ROI(rn).decoding_output.accuracy_per_ROI_sh,2);
        tempData_AccPerRoiSh = [tempData_AccPerRoiSh; FigureEight.session(ip).fig8_ROI(rn).decoding_output.accuracy_per_ROI_sh']; %#ok<AGROW>
        if ip == 1
            index_datatimes03 = index_ticker03;
        else
            index_datatimes03 = [index_datatimes03; index_ticker03]; %#ok<AGROW>
        end
        %pulling accuracy per ROI pre data
        index_ticker04 = index_ticker04 + size(FigureEight.session(ip).fig8_ROI(rn).decoding_output.accuracy_per_ROI_pre,2);
        tempData_AccPerRoiPre = [tempData_AccPerRoiPre; FigureEight.session(ip).fig8_ROI(rn).decoding_output.accuracy_per_ROI_pre']; %#ok<AGROW>
        if ip == 1
            index_datatimes04 = index_ticker04;
        else
            index_datatimes04 = [index_datatimes04; index_ticker04]; %#ok<AGROW>
        end
        %pulling accuracy per ROI sh pre data
        index_ticker05 = index_ticker05 + size(FigureEight.session(ip).fig8_ROI(rn).decoding_output.accuracy_per_ROI_sh_pre,2);
        tempData_AccPerRoiShPre = [tempData_AccPerRoiShPre; FigureEight.session(ip).fig8_ROI(rn).decoding_output.accuracy_per_ROI_sh_pre']; %#ok<AGROW>
        if ip == 1
            index_datatimes05 = index_ticker05;
        else
            index_datatimes05 = [index_datatimes05; index_ticker05]; %#ok<AGROW>
        end
        %pulling latency per ROI data
        index_ticker06 = index_ticker06 + size(FigureEight.session(ip).fig8_ROI(rn).decoding_output.latency_per_ROI,2);
        tempData_LatencyPerRoi = [tempData_LatencyPerRoi; FigureEight.session(ip).fig8_ROI(rn).decoding_output.latency_per_ROI']; %#ok<AGROW>
        if ip == 1
            index_datatimes06 = index_ticker06;
        else
            index_datatimes06 = [index_datatimes06; index_ticker06]; %#ok<AGROW>
        end
        %pulling mean accuracy data
        meanSize = size(FigureEight.session(ip).fig8_ROI(rn).decoding_output.ROI);
        
        for ms = 1:meanSize(2)
            index_ticker08 = index_ticker08 + size(FigureEight.session(ip).fig8_ROI(rn).decoding_output.ROI(ms).mean_accuracy,2);
            tempData_MeanAccuracy = [tempData_MeanAccuracy; FigureEight.session(ip).fig8_ROI(rn).decoding_output.ROI(ms).mean_accuracy']; %#ok<AGROW>
            if ip == 1 && ms == 1
                index_datatimes08 = index_ticker08;
            else
                index_datatimes08 = [index_datatimes08; index_ticker08]; %#ok<AGROW>
            end
        end
        index_ticker07 = index_ticker07 + size(tempData_MeanAccuracy,1);
        if ip == 1
            index_datatimes07 = index_ticker07;
        else
            index_datatimes07 = [index_datatimes07; index_ticker07]; %#ok<AGROW>
        end
        index_ticker09 = index_ticker09 + size(index_datatimes08,1);
        if ip == 1
            index_datatimes09 = index_ticker09;
        else
            index_datatimes09 = [index_datatimes09; index_ticker09]; %#ok<AGROW>
        end
    end
    %saving time data
    FullData.(ROI_names(rn)).Time.indexPoints = index_datatimes01;
    FullData.(ROI_names(rn)).Time.timeData = tempData_time;
    %saving accuracy Per ROI data
    FullData.(ROI_names(rn)).AccPerRoi.indexPoints = index_datatimes02;
    FullData.(ROI_names(rn)).AccPerRoi.timeData = tempData_AccPerRoi;
    %saving accuracy Per ROI Sh data
    FullData.(ROI_names(rn)).AccPerRoiSh.indexPoints = index_datatimes03;
    FullData.(ROI_names(rn)).AccPerRoiSh.timeData = tempData_AccPerRoiSh;
    %saving accuracy Per ROI Pre data
    FullData.(ROI_names(rn)).AccPerRoiPre.indexPoints = index_datatimes04;
    FullData.(ROI_names(rn)).AccPerRoiPre.timeData = tempData_AccPerRoiPre;
    %saving accuracy Per ROI Pre data
    FullData.(ROI_names(rn)).AccPerRoiShPre.indexPoints = index_datatimes05;
    FullData.(ROI_names(rn)).AccPerRoiShPre.timeData = tempData_AccPerRoiShPre;
    %saving latency per ROI data
    FullData.(ROI_names(rn)).LatencyPerRoi.indexPoints = index_datatimes06;
    FullData.(ROI_names(rn)).LatencyPerRoi.timeData = tempData_LatencyPerRoi;
    %saving mean accuracy values
    FullData.(ROI_names(rn)).MeanAccuracy.indexPoints = index_datatimes07;
    FullData.(ROI_names(rn)).MeanAccuracy.timeData = tempData_MeanAccuracy;
    %saving mean accuracy value index points for each session
    FullData.(ROI_names(rn)).MeanAccuracyIndexPoints.indexPoints = index_datatimes09;
    FullData.(ROI_names(rn)).MeanAccuracyIndexPoints.timeData = index_datatimes08;

end

%{
Setting up process modules for each ragged data set
%}

desc_vectData = [ ...
    "ragged array: time_span", ...
    "ragged array: accuracy_per_ROI", ...
    "ragged array: accuracy_per_ROI_sh", ...
    "ragged array: accuracy_per_ROI_pre", ...
    "ragged array: accuracy_per_ROI_sh_pre", ...
    "ragged array: latency_per_ROI", ...
    "ragged array: mean_accuracy", ...
    "ragged array: mean_accuracy_indexing", ...
    ];
desc_table = [ ...
    "Figure5 processed data: time_span", ...
    "Figure5 processed data: accuracy_per_ROI", ...
    "Figure5 processed data: accuracy_per_ROI_sh", ...
    "Figure5 processed data: accuracy_per_ROI_pre", ...
    "Figure5 processed data: accuracy_per_ROI_sh_pre", ...
    "Figure5 processed data: latency_per_ROI", ...
    "Figure5 processed data: mean_accuracy", ...
    "Figure5 processed data: mean_accuracy_indexing", ...
    ];
desc_module = [ ...
    "time_span", ...
    "accuracy_per_ROI", ...
    "accuracy_per_ROI_sh", ...
    "accuracy_per_ROI_pre", ...
    "accuracy_per_ROI_sh_pre", ...
    "latency_per_ROI", ...
    "mean_accuracy", ...
    "mean_accuracy_indexing", ...
    ];
desc_type = [ ...
    "ROI_1", ...
    "ROI_2", ...
    "ROI_5", ...
    "ROI_15", ...
    "ROI_2000", ...
    ];
%%

%%
[nwb] = addingprocessModule(nwb,FullData.(ROI_names(1)).Time.timeData,FullData.(ROI_names(1)).Time.indexPoints,desc_vectData(1),desc_table(1),desc_module(1),desc_type(1));

%%
function [file] = addingprocessModule(file,Data,index_datatimes,desc_vectData,desc_table,desc_module,desc_type)
multi_ragged_col = types.hdmf_common.VectorData( ...
    'description', char(desc_vectData),...
    'data', Data ...
);
% Define column with VectorIndex
multi_ragged_index = types.hdmf_common.VectorIndex( ...
    'description', 'index to multi_ragged_col', ...
    'target', types.untyped.ObjectView(multi_ragged_col),'data', index_datatimes ...
);
 
multi_ragged_table = types.hdmf_common.DynamicTable( ...
    'description',char(desc_table), ...
    'colnames', {char(desc_module)}, ...
    char(desc_module), multi_ragged_col, ...
    'col1_index', multi_ragged_index, ...
    'id', types.hdmf_common.ElementIdentifiers('data', (0:69)') ...  % 0-indexed, for compatibility with Python
);


testTwo = types.core.ProcessingModule( ...
    'description',char(desc_module));
testTwo.dynamictable.set(desc_module,multi_ragged_table)

file.processing.set(desc_type,testTwo)


end





















