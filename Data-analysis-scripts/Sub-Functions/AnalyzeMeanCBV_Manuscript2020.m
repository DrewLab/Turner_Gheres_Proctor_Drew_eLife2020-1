function [AnalysisResults] = AnalyzeMeanCBV_Manuscript2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the mean hemodynamic signal during each behavioral state
%________________________________________________________________________________________________________________________

%% function parameters
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
Iso_animalIDs = {'T108','T109','T110','T111','T119','T120','T121','T122','T123'};
modelTypes = {'SVM','Ensemble','Forest','Manual'};
params.minTime.Rest = 10;   % seconds
params.Offset = 2;
params.minTime.Whisk = params.Offset + 5;
params.minTime.Stim = params.Offset + 2;
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 60;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(IOS_animalIDs,animalID))
    dataLocation = [rootFolder '\' animalID '\Bilateral Imaging\'];
    cd(dataLocation)
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
    fileBreaks = strfind(restDataFileID,'_');
    animalID = restDataFileID(1:fileBreaks(1)-1);
    samplingRate = RestData.CBV.adjLH.CBVCamSamplingRate;
    % lowpass filter and detrend each segment
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    WhiskCriteria.Fieldname = {'duration','duration','puffDistance'};
    WhiskCriteria.Comparison = {'gt','lt','gt'};
    WhiskCriteria.Value = {2,5,5};
    WhiskPuffCriteria.Fieldname = {'puffDistance'};
    WhiskPuffCriteria.Comparison = {'gt'};
    WhiskPuffCriteria.Value = {5};
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    RestPuffCriteria.Fieldname = {'puffDistances'};
    RestPuffCriteria.Comparison = {'gt'};
    RestPuffCriteria.Value = {5};
    
    %% Analyze mean CBV during periods of rest
    [restLogical] = FilterEvents_IOS_Manuscript2020(RestData.CBV_HbT.adjLH,RestCriteria);
    [puffLogical] = FilterEvents_IOS_Manuscript2020(RestData.CBV_HbT.adjLH,RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.CBV_HbT.adjLH.fileIDs(combRestLogical,:);
    restEventTimes = RestData.CBV_HbT.adjLH.eventTimes(combRestLogical,:);
    restDurations = RestData.CBV_HbT.adjLH.durations(combRestLogical,:);
    LH_RestingData = RestData.CBV_HbT.adjLH.data(combRestLogical,:);
    RH_RestingData = RestData.CBV_HbT.adjRH.data(combRestLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [LH_finalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS_Manuscript2020(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    clear LH_ProcRestData LH_restCBVMean
    clear RH_ProcRestData RH_restCBVMean
    for gg = 1:length(LH_finalRestData)
        LH_ProcRestData{gg,1} = filtfilt(sos,g,LH_finalRestData{gg,1}); %#ok<*AGROW>
        RH_ProcRestData{gg,1} = filtfilt(sos,g,RH_finalRestData{gg,1});
    end
    % analyze mean HbT during resting epochs
    for n = 1:length(LH_ProcRestData)
        LH_restCBVMean(n,1) = mean(LH_ProcRestData{n,1}(1:end));
        RH_restCBVMean(n,1) = mean(RH_ProcRestData{n,1}(1:end));
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.MeanAdjLH = LH_restCBVMean;
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.MeanAdjRH = RH_restCBVMean;
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.IndAdjLH = LH_ProcRestData;
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.IndAdjRH = RH_ProcRestData;
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.FileIDs = finalRestFileIDs;
    
    %% Analyze mean CBV during periods of extended whisking
    [whiskLogical] = FilterEvents_IOS_Manuscript2020(EventData.CBV_HbT.adjLH.whisk,WhiskCriteria);
    [puffLogical] = FilterEvents_IOS_Manuscript2020(EventData.CBV_HbT.adjLH.whisk,WhiskPuffCriteria);
    combWhiskLogical = logical(whiskLogical.*puffLogical);
    whiskFileIDs = EventData.CBV_HbT.adjLH.whisk.fileIDs(combWhiskLogical,:);
    whiskEventTimes = EventData.CBV_HbT.adjLH.whisk.eventTime(combWhiskLogical,:);
    whiskDurations = EventData.CBV_HbT.adjLH.whisk.duration(combWhiskLogical,:);
    LH_whiskData = EventData.CBV_HbT.adjLH.whisk.data(combWhiskLogical,:);
    RH_whiskData = EventData.CBV_HbT.adjRH.whisk.data(combWhiskLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [LH_finalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS_Manuscript2020(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    [RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    clear LH_ProcWhiskData LH_whiskCBVMean LH_whiskCBV
    clear RH_ProcWhiskData RH_whiskCBVMean RH_whiskCBV
    for gg = 1:size(LH_finalWhiskData,1)
        LH_ProcWhiskData_temp = filtfilt(sos,g,LH_finalWhiskData(gg,:));
        LH_ProcWhiskData(gg,:) = LH_ProcWhiskData_temp - mean(LH_ProcWhiskData_temp(1:params.Offset*samplingRate));
        RH_ProcWhiskData_temp = filtfilt(sos,g,RH_finalWhiskData(gg,:));
        RH_ProcWhiskData(gg,:) = RH_ProcWhiskData_temp - mean(RH_ProcWhiskData_temp(1:params.Offset*samplingRate));
    end
    % analyze mean HbT during whisking epochs
    for n = 1:size(LH_ProcWhiskData,1)
        LH_whiskCBVMean{n,1} = mean(LH_ProcWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        RH_whiskCBVMean{n,1} = mean(RH_ProcWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        LH_whiskCBV{n,1} = LH_ProcWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
        RH_whiskCBV{n,1} = RH_ProcWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.MeanAdjLH = cell2mat(LH_whiskCBVMean);
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.MeanAdjRH = cell2mat(RH_whiskCBVMean);
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.IndAdjLH = LH_whiskCBV;
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.IndAdjRH = RH_whiskCBV;
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.FileIDs = finalWhiskFileIDs;
    
    %% Analyze mean CBV during periods of stimulations
    % Criteria for the FilterEvents data struct
    stimCriteriaA.Value = {'RPadSol'};
    stimCriteriaA.Fieldname = {'solenoidName'};
    stimCriteriaA.Comparison = {'equal'};
    stimCriteriaB.Value = {'LPadSol'};
    stimCriteriaB.Fieldname = {'solenoidName'};
    stimCriteriaB.Comparison = {'equal'};
    % filter the EventData.mat structure for stimulus events that meet the desired criteria
    LH_stimFilter = FilterEvents_IOS_Manuscript2020(EventData.CBV_HbT.adjLH.stim,stimCriteriaA);
    RH_stimFilter = FilterEvents_IOS_Manuscript2020(EventData.CBV_HbT.adjRH.stim,stimCriteriaB);
    [LH_stimHbTData] = EventData.CBV_HbT.adjLH.stim.data(LH_stimFilter,:);
    [RH_stimHbTData] = EventData.CBV_HbT.adjRH.stim.data(RH_stimFilter,:);
    [LH_stimFileIDs] = EventData.CBV_HbT.adjLH.stim.fileIDs(LH_stimFilter,:);
    [RH_stimFileIDs] = EventData.CBV_HbT.adjRH.stim.fileIDs(RH_stimFilter,:);
    [LH_stimEventTimes] = EventData.CBV_HbT.adjLH.stim.eventTime(LH_stimFilter,:);
    [RH_stimEventTimes] = EventData.CBV_HbT.adjRH.stim.eventTime(RH_stimFilter,:);
    LH_stimDurations = zeros(length(LH_stimEventTimes),1);
    RH_stimDurations = zeros(length(RH_stimEventTimes),1);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [LH_finalStimData,LH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS_Manuscript2020(LH_stimHbTData,LH_stimFileIDs,LH_stimDurations,LH_stimEventTimes,ManualDecisions);
    [RH_finalStimData,RH_finalStimFileIDs,~,~] = RemoveInvalidData_IOS_Manuscript2020(RH_stimHbTData,RH_stimFileIDs,RH_stimDurations,RH_stimEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    clear LH_ProcStimData LH_stimCBVMean LH_stimCBV
    clear RH_ProcStimData RH_stimCBVMean RH_stimCBV
    for gg = 1:size(LH_finalStimData,1)
        LH_ProcStimData_temp = filtfilt(sos,g,LH_finalStimData(gg,:));
        LH_ProcStimData(gg,:) = LH_ProcStimData_temp - mean(LH_ProcStimData_temp(1:params.Offset*samplingRate));
    end
    %
    for gg = 1:size(RH_finalStimData,1)
        RH_ProcStimData_temp = filtfilt(sos,g,RH_finalStimData(gg,:));
        RH_ProcStimData(gg,:) = RH_ProcStimData_temp - mean(RH_ProcStimData_temp(1:params.Offset*samplingRate));
    end
    % analyze mean HbT during whisking epochs
    for n = 1:size(LH_ProcStimData,1)
        LH_stimCBVMean{n,1} = mean(LH_ProcStimData(n,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
        LH_stimCBV{n,1} = LH_ProcStimData(n,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
    end
    %
    for n = 1:size(RH_ProcStimData,1)
        RH_stimCBVMean{n,1} = mean(RH_ProcStimData(n,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate),2);
        RH_stimCBV{n,1} = RH_ProcStimData(n,(params.Offset + 1)*samplingRate:params.minTime.Stim*samplingRate);
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.MeanAdjLH = cell2mat(LH_stimCBVMean);
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.MeanAdjRH = cell2mat(RH_stimCBVMean);
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.IndAdjLH = LH_stimCBV;
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.IndAdjRH = RH_stimCBV;
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.LH_FileIDs = LH_finalStimFileIDs;
    AnalysisResults.(animalID).MeanCBV.Stim.CBV_HbT.RH_FileIDs = RH_finalStimFileIDs;
    
    %% Analyze mean CBV during periods of NREM sleep
    for xx = 1:length(modelTypes)
        modelType = modelTypes{1,xx};
        % pull data from SleepData.mat structure
        LH_nremData = SleepData.(modelType).NREM.data.CBV_HbT.LH;
        RH_nremData = SleepData.(modelType).NREM.data.CBV_HbT.RH;
        nremFileIDs = SleepData.(modelType).NREM.FileIDs;
        clear LH_nremCBVMean  
        clear RH_nremCBVMean 
        % analyze mean HbT during NREM epochs
        for n = 1:length(LH_nremData)
            LH_nremCBVMean(n,1) = mean(filtfilt(sos,g,LH_nremData{n,1}(1:end)));
            RH_nremCBVMean(n,1) = mean(filtfilt(sos,g,RH_nremData{n,1}(1:end)));
        end
        % save results
        AnalysisResults.(animalID).MeanCBV.NREM.(modelType).CBV_HbT.MeanAdjLH = LH_nremCBVMean;
        AnalysisResults.(animalID).MeanCBV.NREM.(modelType).CBV_HbT.MeanAdjRH = RH_nremCBVMean;
        AnalysisResults.(animalID).MeanCBV.NREM.(modelType).CBV_HbT.IndAdjLH = LH_nremData;
        AnalysisResults.(animalID).MeanCBV.NREM.(modelType).CBV_HbT.IndAdjRH = RH_nremData;
        AnalysisResults.(animalID).MeanCBV.NREM.(modelType).CBV_HbT.FileIDs = nremFileIDs;
        
        %% Analyze mean CBV during periods of REM sleep
        % pull data from SleepData.mat structure
        LH_remData = SleepData.(modelType).REM.data.CBV_HbT.LH;
        RH_remData = SleepData.(modelType).REM.data.CBV_HbT.RH;
        remFileIDs = SleepData.(modelType).REM.FileIDs;
        clear LH_remCBVMean 
        clear RH_remCBVMean 
        % analyze mean HbT during REM epochs
        for n = 1:length(LH_remData)
            LH_remCBVMean(n,1) = mean(filtfilt(sos,g,LH_remData{n,1}(1:end)));
            RH_remCBVMean(n,1) = mean(filtfilt(sos,g,RH_remData{n,1}(1:end)));
        end
        % save results
        AnalysisResults.(animalID).MeanCBV.REM.(modelType).CBV_HbT.MeanAdjLH = LH_remCBVMean;
        AnalysisResults.(animalID).MeanCBV.REM.(modelType).CBV_HbT.MeanAdjRH = RH_remCBVMean;
        AnalysisResults.(animalID).MeanCBV.REM.(modelType).CBV_HbT.IndAdjLH = LH_remData;
        AnalysisResults.(animalID).MeanCBV.REM.(modelType).CBV_HbT.IndAdjRH = RH_remData;
        AnalysisResults.(animalID).MeanCBV.REM.(modelType).CBV_HbT.FileIDs = remFileIDs;
    end
    
    %% Analyze mean CBV during periods of Isolfurane
    if any(strcmp(Iso_animalIDs,animalID))
        dataLocation = [rootFolder '\' animalID '\Isoflurane Trials\'];
        cd(dataLocation)
        % pull ProcData file
        procDataFileStruct = dir('*_ProcData.mat');
        procDataFile = {procDataFileStruct.name}';
        procDataFileID = char(procDataFile);
        load(procDataFileID)
        % extract left and right CBV changes during the last 100 seconds of data
        isoLH_HbT = ProcData.data.CBV_HbT.adjLH((end - samplingRate*100):end);
        filtIsoLH_HbT = filtfilt(sos,g,isoLH_HbT);
        isoRH_HbT = ProcData.data.CBV_HbT.adjRH((end - samplingRate*100):end);
        filtIsoRH_HbT = filtfilt(sos,g,isoRH_HbT);
        AnalysisResults.(animalID).MeanCBV.Iso.CBV_HbT.adjLH = mean(filtIsoLH_HbT);
        AnalysisResults.(animalID).MeanCBV.Iso.CBV_HbT.adjRH = mean(filtIsoRH_HbT);
        AnalysisResults.(animalID).MeanCBV.Iso.CBV_HbT.FileIDs = procDataFileID;
    end
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
