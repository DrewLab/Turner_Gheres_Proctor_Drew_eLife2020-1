function [AnalysisResults] = AnalyzeLaserDoppler_Manuscript2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner%
%
%   Purpose:
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T108','T109','T110','T111','T119','T120','T121','T122','T123'};
modelType = 'Forest';
params.minTime.Rest = 10;   % seconds
params.minTime.Whisk = 7;
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 60;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    % find and load RestData.mat struct
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
    % find and load Manual baseline event information
    manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
    manualBaselineFile = {manualBaselineFileStruct.name}';
    manualBaselineFileID = char(manualBaselineFile);
    load(manualBaselineFileID)
    % find and load RestData.mat struct
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    % find and load SleepData.mat strut
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
    % find and load RestingBaselines.mat struct
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    samplingRate = EventData.flow.data.whisk.samplingRate;
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    
    %% Average behavior-dependent flow
    % Analyze mean laser doppler flow during periods of rest
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    RestPuffCriteria.Fieldname = {'puffDistances'};
    RestPuffCriteria.Comparison = {'gt'};
    RestPuffCriteria.Value = {5};
    [restLogical] = FilterEvents_IOS_Manuscript2020(RestData.flow.data,RestCriteria);
    [puffLogical] = FilterEvents_IOS_Manuscript2020(RestData.flow.data,RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.flow.data.fileIDs(combRestLogical,:);
    restFlowData = RestData.flow.data.NormData(combRestLogical,:);
    restEventTimes = RestData.flow.data.eventTimes(combRestLogical,:);
    restDurations = RestData.flow.data.durations(combRestLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [finalRestFlowData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(restFlowData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    idx = 1;
    for gg = 1:length(finalRestFlowData)
        if isempty(finalRestFlowData{gg,1}) == false
            procRestData{idx,1} = filtfilt(sos,g,finalRestFlowData{gg,1}(1:end)); %#ok<*AGROW>
            idx = idx + 1;
        end
    end
    % analyze correlation coefficient between resting epochs
    for n = 1:length(procRestData)
        restFlowMean(n,1) = mean(procRestData{n,1})*100;   % convert to %
    end
    % save results
    AnalysisResults.(animalID).LDFlow.Rest = restFlowMean;
    
    %% Analyze mean CBV during periods of extended whisking
    % criteria for the FilterEvents data struct
    WhiskCriteria.Fieldname = {'duration','duration','puffDistance'};
    WhiskCriteria.Comparison = {'gt','lt','gt'};
    WhiskCriteria.Value = {2,5,5};
    WhiskPuffCriteria.Fieldname = {'puffDistance'};
    WhiskPuffCriteria.Comparison = {'gt'};
    WhiskPuffCriteria.Value = {5};
    [whiskLogical] = FilterEvents_IOS_Manuscript2020(EventData.flow.data.whisk,WhiskCriteria);
    [puffLogical] = FilterEvents_IOS_Manuscript2020(EventData.flow.data.whisk,WhiskPuffCriteria);
    combWhiskLogical = logical(whiskLogical.*puffLogical);
    whiskFlowData = EventData.flow.data.whisk.NormData(combWhiskLogical,:);
    whiskFileIDs = EventData.flow.data.whisk.fileIDs(combWhiskLogical,:);
    whiskEventTimes = EventData.flow.data.whisk.eventTime(combWhiskLogical,:);
    whiskDurations = EventData.flow.data.whisk.duration(combWhiskLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [finalWhiskData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(whiskFlowData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    for gg = 1:size(finalWhiskData,1)
        procWhiskData(gg,:) = filtfilt(sos,g,finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate));
    end
    % analyze correlation coefficient between resting epochs
    for n = 1:size(procWhiskData,1)
        whiskFlowMean(n,1) = mean(procWhiskData(n,:),2)*100;
    end
    % save results
    AnalysisResults.(animalID).LDFlow.Whisk = whiskFlowMean;
    
    %% Analyze mean CBV during periods of NREM sleep
    % pull data from SleepData.mat structure
    nremData = SleepData.(modelType).NREM.data.DopplerFlow;
    % analyze correlation coefficient between NREM epochs
    idx = 1;
    for n = 1:length(nremData)
        if sum(isnan(nremData{n,1})) == 0
            nremFlowMean(idx,1) = mean(filtfilt(sos,g,nremData{n,1}(1:end)))*100;
            idx = idx + 1;
        end
    end
    % save results
    AnalysisResults.(animalID).LDFlow.NREM = nremFlowMean;
    
    %% Analyze mean CBV during periods of REM sleep
    % pull data from SleepData.mat structure
    remData = SleepData.(modelType).REM.data.DopplerFlow;
    % analyze correlation coefficient between NREM epochs
    idx = 1;
    for n = 1:length(remData)
        if sum(isnan(remData{n,1})) == 0
            remFlowMean(idx,1) = mean(filtfilt(sos,g,remData{n,1}(1:end)))*100;
            idx = idx + 1;
        end
    end
    % save results
    AnalysisResults.(animalID).LDFlow.REM = remFlowMean;
    
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
