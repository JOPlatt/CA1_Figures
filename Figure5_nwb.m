%{
creating the NWB file for all inputs needed for figure 5 
%}
%{
path_to_matnwb = '~/Repositories/matnwb'; % change to your own path location
addpath(genpath(pwd));
generateCore();
%}
nwb_data = NwbFile( ...
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
desc_dataTypes = [ ...
    "Time", ...
    "AccPerRoi", ...
    "AccPerRoiSh", ...
    "AccPerRoiPre", ...
    "AccPerRoiShPre", ...
    "LatencyPerRoi", ...
    "MeanAccuracy", ...
    "MeanAccuracyIndexPoints", ...
    ];
ROI_names = [ ...
    "ROInum_1", ...
    "ROInum_2", ...
    "ROInum_5", ...
    "ROInum_15", ...
    "ROInum_2000"];
totalSets = 70;

%{
creating process module for time_span
%}

for rn = 1:size(ROI_names,2)
    tempData_time = [];
    tempData_AccPerRoi = [];
    tempData_AccPerRoiSh = [];
    tempData_AccPerRoiPre = [];
    tempData_AccPerRoiShPre = [];
    tempData_LatencyPerRoi = [];
    tempData_MeanAccuracy = [];
    index_ticker01 = 0; %time
    index_ticker02 = 0; %per ROI
    index_ticker03 = 0; %per ROI sh
    index_ticker04 = 0; %per ROI pre
    index_ticker05 = 0; %per ROI sh pre
    index_ticker06 = 0; %latency
    index_ticker07 = 0; %mean average
    
    index_ticker09 = 0; %session mean average indexs

    for ip = 1:totalSets
        index_ticker08 = 0; %session mean average index values 
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
        index_ticker07 = size(tempData_MeanAccuracy,1);
        if ip == 1
            index_datatimes07 = index_ticker07;
        else
            index_datatimes07 = [index_datatimes07; index_ticker07]; %#ok<AGROW>
        end
        index_ticker09 = size(index_datatimes08,1);
        if ip == 1
            index_datatimes09 = index_ticker09;
        else
            index_datatimes09 = [index_datatimes09; index_ticker09]; %#ok<AGROW>
        end
    end
    %saving time data
    FullData.(ROI_names(rn)).Time.IndexPoints = index_datatimes01;
    FullData.(ROI_names(rn)).Time.Data = tempData_time;
    %saving accuracy Per ROI data
    FullData.(ROI_names(rn)).AccPerRoi.IndexPoints = index_datatimes02;
    FullData.(ROI_names(rn)).AccPerRoi.Data = tempData_AccPerRoi;
    %saving accuracy Per ROI Sh data
    FullData.(ROI_names(rn)).AccPerRoiSh.IndexPoints = index_datatimes03;
    FullData.(ROI_names(rn)).AccPerRoiSh.Data = tempData_AccPerRoiSh;
    %saving accuracy Per ROI Pre data
    FullData.(ROI_names(rn)).AccPerRoiPre.IndexPoints = index_datatimes04;
    FullData.(ROI_names(rn)).AccPerRoiPre.Data = tempData_AccPerRoiPre;
    %saving accuracy Per ROI Pre data
    FullData.(ROI_names(rn)).AccPerRoiShPre.IndexPoints = index_datatimes05;
    FullData.(ROI_names(rn)).AccPerRoiShPre.Data = tempData_AccPerRoiShPre;
    %saving latency per ROI data
    FullData.(ROI_names(rn)).LatencyPerRoi.IndexPoints = index_datatimes06;
    FullData.(ROI_names(rn)).LatencyPerRoi.Data = tempData_LatencyPerRoi;
    %saving mean accuracy values
    FullData.(ROI_names(rn)).MeanAccuracy.IndexPoints = index_datatimes07;
    FullData.(ROI_names(rn)).MeanAccuracy.Data = tempData_MeanAccuracy;
    %saving mean accuracy value index points for each session
    FullData.(ROI_names(rn)).MeanAccuracyIndexPoints.IndexPoints = index_datatimes09;
    FullData.(ROI_names(rn)).MeanAccuracyIndexPoints.Data = index_datatimes08;

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
%% validating that the data is the same for mean average
rn = 1;
Vdata = FullData.(ROI_names(rn)).MeanAccuracy.Data;
Vindex = FullData.(ROI_names(rn)).MeanAccuracy.IndexPoints;
innData = FullData.(ROI_names(rn)).MeanAccuracyIndexPoints.Data;
innindex = FullData.(ROI_names(rn)).MeanAccuracyIndexPoints.IndexPoints;
stPt = 1;
inum = 1;
validData = cell([size(innindex,1),1]);
for iv = 1:size(innindex,2)
    enPt = innindex(iv);
    if enPt <= size(Vdata,2)
        validData{iv,1}{:,:} = Vdata(stPt:enPt);
        if enPt >= Vindex(iNum)
            iNum = iNum +1;
        end

    end
    stPt = enPt+1;
end


%%
for IO = 1:size(desc_type,2)
    for mod = 1:size(desc_module,2)
        [nwb_data] = addingprocessModule(nwb_data,FullData.(ROI_names(IO)).(desc_dataTypes(mod)).Data,FullData.(ROI_names(IO)).(desc_dataTypes(mod)).IndexPoints,desc_vectData(mod),desc_table(mod),desc_module(mod),desc_type(IO));
    end
end

%%

nwbExport(nwb_data,"Figure5_ProcessedData.nwb")
%% pulling the data back out
Figure5_data = nwbRead("Figure5_ProcessedData.nwb",'ignorecache');
%%
name_data = [ ...
    "time_span", ...
    "accuracy_per_ROI", ...
    "accuracy_per_ROI_sh", ...
    "accuracy_per_ROI_pre", ...
    "accuracy_per_ROI_sh_pre", ...
    "latency_per_ROI", ...
    "mean_accuracy", ...
    "mean_accuracy_indexing", ...
    ];
name_ROI = [ ...
    "ROI_1", ...
    "ROI_2", ...
    "ROI_5", ...
    "ROI_15", ...
    "ROI_2000", ...
    ];
FigureEight = struct;
pastend = zeros([70,5]);
%first pull out all the values in there sturcutres. for the last place into
%a group of cells. This will be followed by dividing up into 70 row chuncks





for rnum = 1:size(name_ROI,2) %for all the ROIs
    for ct = 1:70 %for all the sessions
        for pik = 1:size(name_data,2)-1 %for all the datasets within
            %assigning dynamic names
            if ct < 10
                name_session = append("session_0",num2str(ct));
            else
                name_session = append("session_",num2str(ct));
            end
            name_process = append(name_ROI(rnum),'_',name_data(pik));
            %pulling data from NWB file
            temp_table = Figure5_data.processing.get(char(name_process)).dynamictable.get(char(name_data(pik))).getRow(ct);
            %running throuhg a condition that if dataset = 7 the
            %processing is different
            if pik == size(name_data,2)-1
                %name of dataset wihtin the NWB file
                name_index = append(name_ROI(rnum),'_',name_data(pik+1));
                %pulling index points for a single session using NWB files
                temp_indexpoints = Figure5_data.processing.get(char(name_index)).dynamictable.get(char(name_data(pik+1))).getRow(ct);
                tamp_indexpoints = temp_indexpoints{:,:}{1,1};
                %accounting for the index offset after first session
                %pulling offset value
                time_index = Figure5_data.processing.get(char(name_process)).dynamictable.get(char(name_data(pik))).vectordata.get('col1_index').data(:);
                %putting data in correct format
                temp_table = temp_table{1,1}{:,:};
                %pulling the session ROI ending index values
                table_index = Figure5_data.processing.get(char(name_index)).dynamictable.get(char(name_data(pik+1))).getRow(ct); %ending index points for each sessions mean averages
                indexPoints = table_index{1,1}{:,:};
                %creating cell array for data
                average_cell = cell([size(indexPoints,1),1]);
                Round = 1;
                for ti = 1:size(indexPoints,1)
                    RoundEnd = indexPoints(ti);
                    average_cell{ti,1} = temp_table(Round:RoundEnd);
                    Round = RoundEnd + 1;
                end
                FigureFive.(name_session).(name_ROI(rnum)).(name_data(pik)) = cell2table(average_cell);
                FigureFive.(name_session).(name_ROI(rnum)).(name_data(pik)).Properties.VariableNames(1) = "mean_average";
            else
                FigureFive.(name_session).(name_ROI(rnum)).(name_data(pik)) = temp_table{:,:}{1,1}';
            end
            
        end
    end
end


%%
function [file] = addingprocessModule(file,Data,index_datatimes,desc_vectData,desc_table,desc_module,desc_type)
multi_ragged_col = types.hdmf_common.VectorData( ...
    'description', char(desc_vectData),...
    'data', Data ...
);
% Define column with VectorIndex
multi_ragged_index = types.hdmf_common.VectorIndex( ...
    'description', 'index to ragged dataset', ...
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
testTwo.dynamictable.set(desc_module,multi_ragged_table);

processNmae = append(desc_type,'_',desc_module);

file.processing.set(processNmae,testTwo);


end













%                 table_index = nwb.processing.get(char(name_process)).dynamictable.get(char(name_data(pik))+1).getRow(pik); %ending index points for each sessions mean averages
%                 if ct > 1
%                     table_index = table_index - (dataStart-1);
%                 end
%                 average_cell = cell([size(table_index,2),1]);
%                 stpt = pastend(ct,rnum) + 1;
%                 enpt = table_index(pik)
%                 for ti = 1:size(table_index,2)
%                     average_cell{it,1} = temp_table();
%                 end
%                 stpt = enpt +1;
%                 pastend(ct,rnum) = enpt;
%                 temp_table = nwb.processing.get(char(name_process)).dynamictable.get(char(name_data(pik))).vectordata.get('mean_accuracy').data(:);
%                 
%                 fulltemptable = temp_table{:,:}{1,1}';
%                 name_index = append(name_ROI(rnum),'_',name_data(pik+1));
%                 temp_indexpoints = nwb.processing.get(char(name_index)).dynamictable.get(char(name_data(pik+1))).getRow(ct);
%                 tamp_indexpoints = temp_indexpoints{:,:}{1,1};
%                 
%                 if ct > 1
%                     tableindexpoints = temp_indexpoints-pastend(rnum)+1;
%                 else
%                     tableindexpoints = temp_indexpoints;
%                 end
%                 pastend(rnum) = tamp_indexpoints(end);
%                 table_average = cell([size(tableindexpoints,1),1]);
%                 stpt = 1;
%                 try
%                     for ti = 1:size(tableindexpoints,1)
%                         edpt = tableindexpoints(ti);
%                         table_average{ti,1} = fulltemptable(stpt:edpt);
%                         stpt = tableindexpoints(ti)+1;
%                     end
%                 catch
%                     disp("stoppign")
%                 end







