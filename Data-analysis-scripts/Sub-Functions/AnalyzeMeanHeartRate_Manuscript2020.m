function [AnalysisResults] = AnalyzeMeanHeartRate_Manuscript2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Determine the average heart rate during different behaviors.
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
modelType = 'SVM';
params.minTime.Rest = 10;   % seconds
params.minTime.Whisk = 5;
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 60;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    % Character list of all ProcData files
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    % find and load RestData.mat struct
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    % find and load Manual baseline event information
    manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
    manualBaselineFile = {manualBaselineFileStruct.name}';
    manualBaselineFileID = char(manualBaselineFile);
    load(manualBaselineFileID)
    % find and load EventData.mat struct
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
    % find and load RestingBaselines.mat strut
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % find and load SleepData.mat strut
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
    % identify animal's ID and pull important infortmat
    fileBreaks = strfind(restDataFileID, '_');
    animalID = restDataFileID(1:fileBreaks(1)-1);
    WhiskCriteria.Fieldname = {'duration','puffDistance'};
    WhiskCriteria.Comparison = {'gt','gt'};
    WhiskCriteria.Value = {5,5};
    WhiskPuffCriteria.Fieldname = {'puffDistance'};
    WhiskPuffCriteria.Comparison = {'gt'};
    WhiskPuffCriteria.Value = {5};
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    RestPuffCriteria.Fieldname = {'puffDistances'};
    RestPuffCriteria.Comparison = {'gt'};
    RestPuffCriteria.Value = {5};
    
    %% Analyze heart rate during long whisking events
    [whiskLogical] = FilterEvents_IOS_Manuscript2020(EventData.CBV.LH.whisk,WhiskCriteria);
    [puffLogical] = FilterEvents_IOS_Manuscript2020(EventData.CBV.LH.whisk,WhiskPuffCriteria);
    combWhiskLogical = logical(whiskLogical.*puffLogical);
    [allWhiskFileIDs] = EventData.CBV.LH.whisk.fileIDs(combWhiskLogical,:);
    [allWhiskEventTimes] = EventData.CBV.LH.whisk.eventTime(combWhiskLogical,:);
    [allWhiskDurations] = EventData.CBV.LH.whisk.duration(combWhiskLogical,:);
    [allWhiskCBVData] = EventData.CBV.LH.whisk.data(combWhiskLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [~,finalWhiskFileIDs,finalWhiskDurations,finalWhiskEventTimes] = DecimateRestData_Manuscript2020(allWhiskCBVData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
    clear whiskingHeartRate
    for a = 1:length(finalWhiskFileIDs)
        whiskFileID = [animalID '_' finalWhiskFileIDs{a,1} '_ProcData.mat'];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            if strcmp(whiskFileID,procDataFileID) == true
                load(whiskFileID)
                heartRate = ProcData.data.heartRate;
                eventTime = round(finalWhiskEventTimes(a,1));
                duration = round(finalWhiskDurations(a,1));
                try
                    whiskingHeartRate(a,1) = mean(heartRate(eventTime:eventTime + duration)); %#ok<*AGROW>
                catch
                    whiskingHeartRate(a,1) = mean(heartRate(1:eventTime + duration));
                end
                break
            end
        end
    end
    % save results
    AnalysisResults.(animalID).MeanHR.Whisk = whiskingHeartRate;
    
    %% Analyze heart rate during rest data
    % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
    [restLogical] = FilterEvents_IOS_Manuscript2020(RestData.CBV.LH,RestCriteria);
    [puffLogical] = FilterEvents_IOS_Manuscript2020(RestData.CBV.LH,RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.CBV.LH.fileIDs(combRestLogical,:);
    restEventTimes = RestData.CBV.LH.eventTimes(combRestLogical,:);
    restDurations = RestData.CBV.LH.durations(combRestLogical,:);
    restCBVData = RestData.CBV.LH.data(combRestLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [~,finalRestFileList,finalRestDurations,finalRestEventTimes] = DecimateRestData_Manuscript2020(restCBVData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    clear restingHeartRate
    for a = 1:length(finalRestFileList)
        restFileID = [animalID '_' finalRestFileList{a,1} '_ProcData.mat'];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            if strcmp(restFileID,procDataFileID) == true
                load(restFileID)
                heartRate = ProcData.data.heartRate;
                eventTime = floor(finalRestEventTimes(a,1));
                duration = floor(finalRestDurations(a,1));
                try
                    restingHeartRate(a,1) = mean(heartRate(eventTime:eventTime + duratiob));
                catch
                    restingHeartRate(a,1) = mean(heartRate(1:eventTime + duration));
                end
                break
            end
        end
    end
    % save results
    AnalysisResults.(animalID).MeanHR.Rest = restingHeartRate;
    
    %% Analyze heart rate during periods of NREM sleep
    % pull data from SleepData.mat structure
    nremData = SleepData.(modelType).NREM.data.HeartRate;
    % analyze correlation coefficient between NREM epochs
    for n = 1:length(nremData)
        nremHRMean(n,1) = mean(nremData{n,1}(1:end));
    end
    % save results
    AnalysisResults.(animalID).MeanHR.NREM = nremHRMean;
    
    %% Analyze heart rate during periods of REM sleep
    % pull data from SleepData.mat structure
    remData = SleepData.(modelType).REM.data.HeartRate;
    % analyze correlation coefficient between REM epochs
    for n = 1:length(remData)
        remHRMean(n,1) = mean(remData{n,1}(1:end));
    end
    % save results
    AnalysisResults.(animalID).MeanHR.REM = remHRMean;
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
