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
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120'};
Iso_animalIDs = {'T108','T109','T110','T111','T119','T120'};
baselineType = 'manualSelection';
modelTypes = {'SVM','Ensemble','Forest','Manual'};
params.minTime.Rest = 10;   % seconds

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
    
    %% Analyze mean CBV during periods of rest
    [restLogical] = FilterEvents_IOS(RestData.CBV_HbT.adjLH,RestCriteria);
    [puffLogical] = FilterEvents_IOS(RestData.CBV_HbT.adjLH,RestPuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.CBV_HbT.adjLH.fileIDs(combRestLogical,:);
    restEventTimes = RestData.CBV_HbT.adjLH.eventTimes(combRestLogical,:);
    restDurations = RestData.CBV_HbT.adjLH.durations(combRestLogical,:);
    LH_RestingData = RestData.CBV_HbT.adjLH.data(combRestLogical,:);
    RH_RestingData = RestData.CBV_HbT.adjRH.data(combRestLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [LH_finalRestData,~,~,~] = DecimateRestData_Manuscript2020(LH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [RH_finalRestData,~,~,~] = DecimateRestData_Manuscript2020(RH_RestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    % lowpass filter and detrend each segment
    [B,A] = butter(3,1/(samplingRate/2),'low');
    clear LH_ProcRestData
    clear RH_ProcRestData
    for g = 1:length(LH_finalRestData)
        LH_ProcRestData{g,1} = filtfilt(B,A,LH_finalRestData{g,1}); %#ok<*AGROW>
        RH_ProcRestData{g,1} = filtfilt(B,A,RH_finalRestData{g,1});
    end
    % analyze correlation coefficient between resting epochs
    for n = 1:length(LH_ProcRestData)
        LH_restCBVMean(n,1) = mean(LH_ProcRestData{n,1});
        RH_restCBVMean(n,1) = mean(RH_ProcRestData{n,1});
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.adjLH = LH_restCBVMean;
    AnalysisResults.(animalID).MeanCBV.Rest.CBV_HbT.adjRH = RH_restCBVMean;
    
    %% Analyze mean CBV during periods of extended whisking
    [whiskLogical] = FilterEvents_IOS(EventData.CBV_HbT.adjLH.whisk,WhiskCriteria);
    [puffLogical] = FilterEvents_IOS(EventData.CBV_HbT.adjLH.whisk,WhiskPuffCriteria);
    combWhiskLogical = logical(whiskLogical.*puffLogical);
    whiskFileIDs = EventData.CBV_HbT.adjLH.whisk.fileIDs(combWhiskLogical,:);
    whiskEventTimes = EventData.CBV_HbT.adjLH.whisk.eventTime(combWhiskLogical,:);
    whiskDurations = EventData.CBV_HbT.adjLH.whisk.duration(combWhiskLogical,:);
    LH_whiskData = EventData.CBV_HbT.adjLH.whisk.data(combWhiskLogical,:);
    RH_whiskData = EventData.CBV_HbT.adjRH.whisk.data(combWhiskLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    LH_finalWhiskData = DecimateRestData_Manuscript2020(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    RH_finalWhiskData = DecimateRestData_Manuscript2020(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    % lowpass filter and detrend each segment
    [B,A] = butter(3,1/(samplingRate/2),'low');
    clear LH_ProcWhiskData
    clear RH_ProcWhiskData
    for g = 1:size(LH_finalWhiskData,1)
        LH_ProcWhiskData(g,:) = filtfilt(B,A,LH_finalWhiskData(g,:));
        RH_ProcWhiskData(g,:) = filtfilt(B,A,RH_finalWhiskData(g,:));
    end
    % analyze correlation coefficient between resting epochs
    for n = 1:size(LH_ProcWhiskData,1)
        LH_whiskCBVMean{n,1} = mean(LH_ProcWhiskData(n,samplingRate*2:end),2);
        RH_whiskCBVMean{n,1} = mean(RH_ProcWhiskData(n,samplingRate*2:end),2);
    end
    % save results
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.adjLH = cell2mat(LH_whiskCBVMean);
    AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.adjRH = cell2mat(RH_whiskCBVMean);
    
    %% Analyze mean CBV during periods of NREM sleep
    for xx = 1:length(modelTypes)
        modelType = modelTypes{1,xx};
        % pull data from SleepData.mat structure
        LH_nremData = SleepData.(modelType).NREM.data.CBV_HbT.LH;
        RH_nremData = SleepData.(modelType).NREM.data.CBV_HbT.RH;
        % analyze correlation coefficient between NREM epochs
        clear LH_nremCBVMean RH_nremCBVMean
        for n = 1:length(LH_nremData)
            LH_nremCBVMean(n,1) = mean(LH_nremData{n,1});
            RH_nremCBVMean(n,1) = mean(RH_nremData{n,1});
        end
        % save results
        AnalysisResults.(animalID).MeanCBV.NREM.(modelType).CBV_HbT.adjLH = LH_nremCBVMean;
        AnalysisResults.(animalID).MeanCBV.NREM.(modelType).CBV_HbT.adjRH = RH_nremCBVMean;
        
        %% Analyze mean CBV during periods of REM sleep
        % pull data from SleepData.mat structure
        LH_remData = SleepData.(modelType).REM.data.CBV_HbT.LH;
        RH_remData = SleepData.(modelType).REM.data.CBV_HbT.RH;
        % analyze correlation coefficient between NREM epochs
        clear LH_remCBVMean RH_remCBVMean
        for n = 1:length(LH_remData)
            LH_remCBVMean(n,1) = mean(LH_remData{n,1});
            RH_remCBVMean(n,1) = mean(RH_remData{n,1});
        end
        % save results
        AnalysisResults.(animalID).MeanCBV.REM.(modelType).CBV_HbT.adjLH = LH_remCBVMean;
        AnalysisResults.(animalID).MeanCBV.REM.(modelType).CBV_HbT.adjRH = RH_remCBVMean;
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
        filtIsoLH_HbT = filtfilt(B,A,isoLH_HbT);
        isoRH_HbT = ProcData.data.CBV_HbT.adjRH((end - samplingRate*100):end);
        filtIsoRH_HbT = filtfilt(B,A,isoRH_HbT);
        AnalysisResults.(animalID).MeanCBV.Iso.CBV_HbT.adjLH = mean(filtIsoLH_HbT);
        AnalysisResults.(animalID).MeanCBV.Iso.CBV_HbT.adjRH = mean(filtIsoRH_HbT);
    end
end
cd(rootFolder)

end
