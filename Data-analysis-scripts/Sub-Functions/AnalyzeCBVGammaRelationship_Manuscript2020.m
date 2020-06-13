function [AnalysisResults] = AnalyzeCBVGammaRelationship_Manuscript2020(animalID,rootFolder,AnalysisResults)
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
modelType = 'Forest';
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
    LH_HbTrestingData = RestData.CBV_HbT.adjLH.data(combRestLogical,:);
    RH_HbTrestingData = RestData.CBV_HbT.adjRH.data(combRestLogical,:);
    LH_gammaRestingData = RestData.cortical_LH.gammaBandPower.data(combRestLogical,:);
    RH_gammaRestingData = RestData.cortical_RH.gammaBandPower.data(combRestLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [LH_HbTfinalRestData,finalRestFileIDs,~,~] = RemoveInvalidData_IOS_Manuscript2020(LH_HbTrestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [RH_HbTfinalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(RH_HbTrestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [LH_gammaFinalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(LH_gammaRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    [RH_gammaFinalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(RH_gammaRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    clear LH_HbTprocRestData LH_HbTrestCBVMean LH_gammaProcRestData LH_gammaRestCBVMean
    clear RH_HbTprocRestData RH_HbTrestCBVMean RH_gammaProcRestData RH_gammaRestCBVMean
    for gg = 1:length(LH_HbTfinalRestData)
        LH_HbTprocRestData{gg,1} = filtfilt(sos,g,LH_HbTfinalRestData{gg,1}); %#ok<*AGROW>
        RH_HbTprocRestData{gg,1} = filtfilt(sos,g,RH_HbTfinalRestData{gg,1});
        LH_gammaProcRestData{gg,1} = filtfilt(sos,g,LH_gammaFinalRestData{gg,1});
        RH_gammaProcRestData{gg,1} = filtfilt(sos,g,RH_gammaFinalRestData{gg,1});
    end
    % analyze mean HbT during resting epochs
    for n = 1:length(LH_HbTprocRestData)
        LH_HbTrestCBVMean(n,1) = mean(LH_HbTprocRestData{n,1}(1:end));
        RH_HbTrestCBVMean(n,1) = mean(RH_HbTprocRestData{n,1}(1:end));
        LH_gammaRestCBVMean(n,1) = mean(LH_gammaProcRestData{n,1}(1:end));
        RH_gammaRestCBVMean(n,1) = mean(RH_gammaProcRestData{n,1}(1:end));
    end
    % save results
    AnalysisResults.(animalID).HbTvsGamma.Rest.HbT.MeanAdjLH = LH_HbTrestCBVMean;
    AnalysisResults.(animalID).HbTvsGamma.Rest.HbT.MeanAdjRH = RH_HbTrestCBVMean;
    AnalysisResults.(animalID).HbTvsGamma.Rest.HbT.IndAdjLH = LH_HbTprocRestData;
    AnalysisResults.(animalID).HbTvsGamma.Rest.HbT.IndAdjRH = RH_HbTprocRestData;
    AnalysisResults.(animalID).HbTvsGamma.Rest.Gamma.MeanAdjLH = LH_gammaRestCBVMean;
    AnalysisResults.(animalID).HbTvsGamma.Rest.Gamma.MeanAdjRH = RH_gammaRestCBVMean;
    AnalysisResults.(animalID).HbTvsGamma.Rest.Gamma.IndAdjLH = LH_gammaProcRestData;
    AnalysisResults.(animalID).HbTvsGamma.Rest.Gamma.IndAdjRH = RH_gammaProcRestData;
    AnalysisResults.(animalID).HbTvsGamma.Rest.fileIDs = finalRestFileIDs;
    
    %% Analyze mean HbT during periods of extended whisking
    [whiskLogical] = FilterEvents_IOS_Manuscript2020(EventData.CBV_HbT.adjLH.whisk,WhiskCriteria);
    [puffLogical] = FilterEvents_IOS_Manuscript2020(EventData.CBV_HbT.adjLH.whisk,WhiskPuffCriteria);
    combWhiskLogical = logical(whiskLogical.*puffLogical);
    whiskFileIDs = EventData.CBV_HbT.adjLH.whisk.fileIDs(combWhiskLogical,:);
    whiskEventTimes = EventData.CBV_HbT.adjLH.whisk.eventTime(combWhiskLogical,:);
    whiskDurations = EventData.CBV_HbT.adjLH.whisk.duration(combWhiskLogical,:);
    LH_HbTwhiskData = EventData.CBV_HbT.adjLH.whisk.data(combWhiskLogical,:);
    RH_HbTwhiskData = EventData.CBV_HbT.adjRH.whisk.data(combWhiskLogical,:);
    LH_gammaWhiskData = EventData.cortical_LH.gammaBandPower.whisk.NormData(combWhiskLogical,:);
    RH_gammaWhiskData = EventData.cortical_RH.gammaBandPower.whisk.NormData(combWhiskLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [LH_HbTfinalWhiskData,finalWhiskFileIDs,~,~] = RemoveInvalidData_IOS_Manuscript2020(LH_HbTwhiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    [RH_HbTfinalWhiskData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(RH_HbTwhiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    [LH_gammaFinalWhiskData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(LH_gammaWhiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    [RH_gammaFinalWhiskData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(RH_gammaWhiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    clear LH_HbTprocWhiskData LH_HbTwhiskCBVMean LH_HbTwhiskCBV LH_gammaProcWhiskData LH_gammaWhiskCBVMean LH_gammaWhiskCBV
    clear RH_HbTprocWhiskData RH_HbTwhiskCBVMean RH_HbTwhiskCBV RH_gammaProcWhiskData RH_gammaWhiskCBVMean RH_gammaWhiskCBV
    for gg = 1:size(LH_HbTfinalWhiskData,1)
        LH_HbTprocWhiskData_temp = filtfilt(sos,g,LH_HbTfinalWhiskData(gg,:));
        LH_HbTprocWhiskData(gg,:) = LH_HbTprocWhiskData_temp - mean(LH_HbTprocWhiskData_temp(1:params.Offset*samplingRate));
        RH_HbTprocWhiskData_temp = filtfilt(sos,g,RH_HbTfinalWhiskData(gg,:));
        RH_HbTprocWhiskData(gg,:) = RH_HbTprocWhiskData_temp - mean(RH_HbTprocWhiskData_temp(1:params.Offset*samplingRate));
        LH_gammaProcWhiskData_temp = filtfilt(sos,g,LH_gammaFinalWhiskData(gg,:));
        LH_gammaProcWhiskData(gg,:) = LH_gammaProcWhiskData_temp;%LH_gammaProcWhiskData_temp - mean(LH_gammaProcWhiskData_temp(1:params.Offset*samplingRate));
        RH_gammaProcWhiskData_temp = filtfilt(sos,g,RH_gammaFinalWhiskData(gg,:));
        RH_gammaProcWhiskData(gg,:) = RH_gammaProcWhiskData_temp;%RH_gammaProcWhiskData_temp - mean(RH_gammaProcWhiskData_temp(1:params.Offset*samplingRate));
    end
    % analyze mean HbT during whisking epochs
    for n = 1:size(LH_HbTprocWhiskData,1)
        LH_HbTwhiskCBVMean{n,1} = mean(LH_HbTprocWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        RH_HbTwhiskCBVMean{n,1} = mean(RH_HbTprocWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        LH_HbTwhiskCBV{n,1} = LH_HbTprocWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
        RH_HbTwhiskCBV{n,1} = RH_HbTprocWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
        LH_gammaWhiskCBVMean{n,1} = mean(LH_gammaProcWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        RH_gammaWhiskCBVMean{n,1} = mean(RH_gammaProcWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate),2);
        LH_gammaWhiskCBV{n,1} = LH_gammaProcWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
        RH_gammaWhiskCBV{n,1} = RH_gammaProcWhiskData(n,params.Offset*samplingRate:params.minTime.Whisk*samplingRate);
    end
    % save results
    AnalysisResults.(animalID).HbTvsGamma.Whisk.HbT.MeanAdjLH = cell2mat(LH_HbTwhiskCBVMean);
    AnalysisResults.(animalID).HbTvsGamma.Whisk.HbT.MeanAdjRH = cell2mat(RH_HbTwhiskCBVMean);
    AnalysisResults.(animalID).HbTvsGamma.Whisk.HbT.IndAdjLH = LH_HbTwhiskCBV;
    AnalysisResults.(animalID).HbTvsGamma.Whisk.HbT.IndAdjRH = RH_HbTwhiskCBV;
    AnalysisResults.(animalID).HbTvsGamma.Whisk.Gamma.MeanAdjLH = cell2mat(LH_gammaWhiskCBVMean);
    AnalysisResults.(animalID).HbTvsGamma.Whisk.Gamma.MeanAdjRH = cell2mat(RH_gammaWhiskCBVMean);
    AnalysisResults.(animalID).HbTvsGamma.Whisk.Gamma.IndAdjLH = LH_gammaWhiskCBV;
    AnalysisResults.(animalID).HbTvsGamma.Whisk.Gamma.IndAdjRH = RH_gammaWhiskCBV;
    AnalysisResults.(animalID).HbTvsGamma.Whisk.fileIDs = finalWhiskFileIDs;
    
    %% Analyze mean CBV during periods of NREM sleep
    % pull data from SleepData.mat structure
    LH_HbTnremData = SleepData.(modelType).NREM.data.CBV_HbT.LH;
    RH_HbTnremData = SleepData.(modelType).NREM.data.CBV_HbT.RH;
    LH_gammaNremData = SleepData.(modelType).NREM.data.cortical_LH.gammaBandPower;
    RH_gammaNremData = SleepData.(modelType).NREM.data.cortical_RH.gammaBandPower;
    nremFileIDs = SleepData.(modelType).NREM.FileIDs;
    clear LH_HbTnremCBVMean LH_gammaNremCBVMean
    clear RH_HbTnremCBVMean RH_gammaNremCBVMean
    % analyze mean HbT during NREM epochs
    for n = 1:length(LH_HbTnremData)
        LH_HbTnremMean(n,1) = mean(filtfilt(sos,g,LH_HbTnremData{n,1}(1:end)));
        RH_HbTnremMean(n,1) = mean(filtfilt(sos,g,RH_HbTnremData{n,1}(1:end)));
        LH_gammaNremMean(n,1) = mean(filtfilt(sos,g,LH_gammaNremData{n,1}(1:end)));
        RH_gammaNremMean(n,1) = mean(filtfilt(sos,g,RH_gammaNremData{n,1}(1:end)));
    end
    % save results
    AnalysisResults.(animalID).HbTvsGamma.NREM.HbT.MeanAdjLH = LH_HbTnremMean;
    AnalysisResults.(animalID).HbTvsGamma.NREM.HbT.MeanAdjRH = RH_HbTnremMean;
    AnalysisResults.(animalID).HbTvsGamma.NREM.HbT.IndAdjLH = LH_HbTnremData;
    AnalysisResults.(animalID).HbTvsGamma.NREM.HbT.IndAdjRH = RH_HbTnremData;
    AnalysisResults.(animalID).HbTvsGamma.NREM.Gamma.MeanAdjLH = LH_gammaNremMean;
    AnalysisResults.(animalID).HbTvsGamma.NREM.Gamma.MeanAdjRH = RH_gammaNremMean;
    AnalysisResults.(animalID).HbTvsGamma.NREM.Gamma.IndAdjLH = LH_gammaNremData;
    AnalysisResults.(animalID).HbTvsGamma.NREM.Gamma.IndAdjRH = RH_gammaNremData;
    AnalysisResults.(animalID).HbTvsGamma.NREM.fileIDs = nremFileIDs;
    
    %% Analyze mean CBV during periods of REM sleep
    % pull data from SleepData.mat structure
    LH_HbTremData = SleepData.(modelType).REM.data.CBV_HbT.LH;
    RH_HbTremData = SleepData.(modelType).REM.data.CBV_HbT.RH;
    LH_gammaRemData = SleepData.(modelType).REM.data.cortical_LH.gammaBandPower;
    RH_gammaRemData = SleepData.(modelType).REM.data.cortical_RH.gammaBandPower;
    remFileIDs = SleepData.(modelType).REM.FileIDs;
    clear LH_HbTremCBVMean LH_gammaRemCBVMean
    clear RH_HbTremCBVMean RH_gammaRemCBVMean
    % analyze mean HbT during REM epochs
    for n = 1:length(LH_HbTremData)
        LH_HbTremMean(n,1) = mean(filtfilt(sos,g,LH_HbTremData{n,1}(1:end)));
        RH_HbTremMean(n,1) = mean(filtfilt(sos,g,RH_HbTremData{n,1}(1:end)));
        LH_gammaRemMean(n,1) = mean(filtfilt(sos,g,LH_gammaRemData{n,1}(1:end)));
        RH_gammaRemMean(n,1) = mean(filtfilt(sos,g,RH_gammaRemData{n,1}(1:end)));
    end
    % save results
    AnalysisResults.(animalID).HbTvsGamma.REM.HbT.MeanAdjLH = LH_HbTremMean;
    AnalysisResults.(animalID).HbTvsGamma.REM.HbT.MeanAdjRH = RH_HbTremMean;
    AnalysisResults.(animalID).HbTvsGamma.REM.HbT.IndAdjLH = LH_HbTremData;
    AnalysisResults.(animalID).HbTvsGamma.REM.HbT.IndAdjRH = RH_HbTremData;
    AnalysisResults.(animalID).HbTvsGamma.REM.Gamma.MeanAdjLH = LH_gammaRemMean;
    AnalysisResults.(animalID).HbTvsGamma.REM.Gamma.MeanAdjRH = RH_gammaRemMean;
    AnalysisResults.(animalID).HbTvsGamma.REM.Gamma.IndAdjLH = LH_gammaRemData;
    AnalysisResults.(animalID).HbTvsGamma.REM.Gamma.IndAdjRH = RH_gammaRemData;
    AnalysisResults.(animalID).HbTvsGamma.REM.fileIDs = remFileIDs;
    
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
