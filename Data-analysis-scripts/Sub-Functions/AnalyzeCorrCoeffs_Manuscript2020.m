function [AnalysisResults] = AnalyzeCorrCoeffs_Manuscript2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Pearson's correlation coefficient between hemodynamic or neural envelope signals
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
modelType = 'Forest';
params.minTime.Rest = 10;   % seconds
params.minTime.Whisk = 7;   % 5 seconds after epoch
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 60;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
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
    samplingRate = RestData.CBV_HbT.adjLH.CBVCamSamplingRate;
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
    % lowpass filter
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % go through each valid data type for behavior-based correlation analysis
    for a = 1:length(dataTypes)
        dataType = dataTypes{1,a};
        
        %% Analyze Pearson's correlation coefficient during periods of rest
        % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
        if strcmp(dataType,'CBV_HbT') == true
            [restLogical] = FilterEvents_IOS_Manuscript2020(RestData.(dataType).adjLH,RestCriteria);
            [puffLogical] = FilterEvents_IOS_Manuscript2020(RestData.(dataType).adjLH,RestPuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            restFileIDs = RestData.(dataType).adjLH.fileIDs(combRestLogical,:);
            restEventTimes = RestData.(dataType).adjLH.eventTimes(combRestLogical,:);
            restDurations = RestData.(dataType).adjLH.durations(combRestLogical,:);
            LH_unstimRestingData = RestData.(dataType).adjLH.data(combRestLogical,:);
            RH_unstimRestingData = RestData.(dataType).adjRH.data(combRestLogical,:);
        else
            [restLogical] = FilterEvents_IOS_Manuscript2020(RestData.cortical_LH.(dataType),RestCriteria);
            [puffLogical] = FilterEvents_IOS_Manuscript2020(RestData.cortical_LH.(dataType),RestPuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
            restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
            restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
            LH_unstimRestingData = RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
            RH_unstimRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
        end
        % decimate the file list to only include those files that occur within the desired number of target minutes
        [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(LH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(RH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
        % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
        clear LH_ProcRestData
        clear RH_ProcRestData
        for gg = 1:length(LH_finalRestData)
            LH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,LH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant'); %#ok<*AGROW>
            RH_ProcRestData{gg,1} = detrend(filtfilt(sos,g,RH_finalRestData{gg,1}(1:params.minTime.Rest*samplingRate)),'constant');
        end
        % analyze correlation coefficient between resting epochs
        for n = 1:length(LH_ProcRestData)
            rest_CC = corrcoef(LH_ProcRestData{n,1},RH_ProcRestData{n,1});
            rest_R(n,1) = rest_CC(2,1);
        end
        meanRest_R = mean(rest_R);
        stdRest_R = std(rest_R,0,1);
        % save results
        AnalysisResults.(animalID).CorrCoeff.Rest.(dataType).R = rest_R;
        AnalysisResults.(animalID).CorrCoeff.Rest.(dataType).meanR = meanRest_R;
        AnalysisResults.(animalID).CorrCoeff.Rest.(dataType).stdR = stdRest_R;
        
        %% Analyze Pearson's correlation coefficient during periods of extended whisking
        % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
        if strcmp(dataType,'CBV_HbT') == true
            [whiskLogical] = FilterEvents_IOS_Manuscript2020(EventData.(dataType).adjLH.whisk,WhiskCriteria);
            [puffLogical] = FilterEvents_IOS_Manuscript2020(EventData.(dataType).adjLH.whisk,WhiskPuffCriteria);
            combWhiskLogical = logical(whiskLogical.*puffLogical);
            whiskFileIDs = EventData.(dataType).adjLH.whisk.fileIDs(combWhiskLogical,:);
            whiskEventTimes = EventData.(dataType).adjLH.whisk.eventTime(combWhiskLogical,:);
            whiskDurations = EventData.(dataType).adjLH.whisk.duration(combWhiskLogical,:);
            LH_whiskData = EventData.(dataType).adjLH.whisk.data(combWhiskLogical,:);
            RH_whiskData = EventData.(dataType).adjRH.whisk.data(combWhiskLogical,:);
        else
            [whiskLogical] = FilterEvents_IOS_Manuscript2020(EventData.cortical_LH.(dataType).whisk,WhiskCriteria);
            [puffLogical] = FilterEvents_IOS_Manuscript2020(EventData.cortical_LH.(dataType).whisk,WhiskPuffCriteria);
            combWhiskLogical = logical(whiskLogical.*puffLogical);
            whiskFileIDs = EventData.cortical_LH.(dataType).whisk.fileIDs(combWhiskLogical,:);
            whiskEventTimes = EventData.cortical_LH.(dataType).whisk.eventTime(combWhiskLogical,:);
            whiskDurations = EventData.cortical_LH.(dataType).whisk.duration(combWhiskLogical,:);
            LH_whiskData = EventData.cortical_LH.(dataType).whisk.NormData(combWhiskLogical,:);
            RH_whiskData = EventData.cortical_RH.(dataType).whisk.NormData(combWhiskLogical,:);
        end
        % decimate the file list to only include those files that occur within the desired number of target minutes
        [LH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(LH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        [RH_finalWhiskData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(RH_whiskData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
        % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
        % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
        % lowpass filter and detrend each segment
        clear LH_ProcWhiskData
        clear RH_ProcWhiskData
        for gg = 1:size(LH_finalWhiskData,1)
            LH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,LH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
            RH_ProcWhiskData(gg,:) = detrend(filtfilt(sos,g,RH_finalWhiskData(gg,2*samplingRate:params.minTime.Whisk*samplingRate)),'constant');
        end
        % analyze correlation coefficient between resting epochs
        for n = 1:size(LH_ProcWhiskData,1)
            whisk_CC = corrcoef(LH_ProcWhiskData(n,:),RH_ProcWhiskData(n,:));
            whisk_R(n,1) = whisk_CC(2,1);
        end
        meanWhisk_R = mean(whisk_R);
        stdWhisk_R = std(whisk_R,0,1);
        % save results
        AnalysisResults.(animalID).CorrCoeff.Whisk.(dataType).R = whisk_R;
        AnalysisResults.(animalID).CorrCoeff.Whisk.(dataType).meanR = meanWhisk_R;
        AnalysisResults.(animalID).CorrCoeff.Whisk.(dataType).stdR = stdWhisk_R;
        
        %% Analyze Pearson's correlation coefficient during periods of awake data
        zz = 1;
        clear LH_AwakeData RH_AwakeData LH_ProcAwakeData RH_ProcAwakeData
        LH_AwakeData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS_Manuscript2020(procDataFileID);
            strDay = ConvertDate_IOS_Manuscript2020(allDataFileDate);
            scoringLabels = [];
            for cc = 1:length(ScoringResults.fileIDs)
                if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                    scoringLabels = ScoringResults.labels{cc,1};
                end
            end
            % check labels for sleep
            if sum(strcmp(scoringLabels,'Not Sleep')) > 168   % 6 bins (180 total) or 30 seconds of sleep
                load(procDataFileID)
                puffs = ProcData.data.solenoids.LPadSol;
                if isempty(puffs) == true
                    if strcmp(dataType,'CBV_HbT') == true
                        LH_AwakeData{zz,1} = ProcData.data.(dataType).adjLH;
                        RH_AwakeData{zz,1} = ProcData.data.(dataType).adjRH;
                    else
                        LH_AwakeData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                        RH_AwakeData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                    end
                    zz = zz + 1;
                end
            end
        end
        for gg = 1:length(LH_AwakeData,1)
            LH_ProcAwakeData{gg,1} = detrend(filtfilt(sos,g,LH_AwakeData{gg,1}),'constant');
            RH_ProcAwakeData{gg,1} = detrend(filtfilt(sos,g,RH_AwakeData{gg,1}),'constant');
        end
        % analyze correlation coefficient between resting epochs
        for n = 1:length(LH_ProcAwakeData)
            awake_CC = corrcoef(LH_ProcAwakeData{n,1},RH_ProcAwakeData{n,1});
            awake_R(n,1) = awake_CC(2,1);
        end
        meanAwake_R = mean(awake_R);
        stdAwake_R = std(awake_R,0,1);
        % save results
        AnalysisResults.(animalID).CorrCoeff.Awake.(dataType).R = awake_R;
        AnalysisResults.(animalID).CorrCoeff.Awake.(dataType).meanR = meanAwake_R;
        AnalysisResults.(animalID).CorrCoeff.Awake.(dataType).stdR = stdAwake_R;
        
        %% Analyze Pearson's correlation coefficient during periods of Sleep data
        zz = 1;
        clear LH_SleepData RH_SleepData LH_ProcSleepData RH_ProcSleepData
        LH_SleepData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,allDataFileDate,allDataFileID] = GetFileInfo_IOS_Manuscript2020(procDataFileID);
            strDay = ConvertDate_IOS_Manuscript2020(allDataFileDate);
            scoringLabels = [];
            for cc = 1:length(ScoringResults.fileIDs)
                if strcmp(allDataFileID,ScoringResults.fileIDs{cc,1}) == true
                    scoringLabels = ScoringResults.labels{cc,1};
                end
            end
            % check labels for sleep
            if sum(strcmp(scoringLabels,'Not Sleep')) > 168   % 6 bins (180 total) or 30 seconds of sleep
                load(procDataFileID)
                puffs = ProcData.data.solenoids.LPadSol;
                if isempty(puffs) == true
                    if strcmp(dataType,'CBV_HbT') == true
                        LH_SleepData{zz,1} = ProcData.data.(dataType).adjLH;
                        RH_SleepData{zz,1} = ProcData.data.(dataType).adjRH;
                    else
                        LH_SleepData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                        RH_SleepData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                    end
                    zz = zz + 1;
                end
            end
        end
        for gg = 1:length(LH_SleepData)
            LH_ProcSleepData{gg,1} = detrend(filtfilt(sos,g,LH_SleepData{gg,1}),'constant');
            RH_ProcSleepData{gg,1} = detrend(filtfilt(sos,g,RH_SleepData{gg,1}),'constant');
        end
        % analyze correlation coefficient between resting epochs
        for n = 1:length(LH_ProcSleepData)
            sleep_CC = corrcoef(LH_ProcSleepData{n,1},RH_ProcSleepData{n,1});
            sleep_R(n,1) = sleep_CC(2,1);
        end
        meanSleep_R = mean(sleep_R);
        stdSleep_R = std(sleep_R,0,1);
        % save results
        AnalysisResults.(animalID).CorrCoeff.Sleep.(dataType).R = sleep_R;
        AnalysisResults.(animalID).CorrCoeff.Sleep.(dataType).meanR = meanSleep_R;
        AnalysisResults.(animalID).CorrCoeff.Sleep.(dataType).stdR = stdSleep_R;
        
        %% Analyze Pearson's correlation coefficient during periods of all data
        zz = 1;
        clear LH_AllData RH_AllData LH_ProcAllData RH_ProcAllData
        LH_AllData = [];
        for bb = 1:size(procDataFileIDs,1)
            procDataFileID = procDataFileIDs(bb,:);
            [~,allDataFileDate,~] = GetFileInfo_IOS_Manuscript2020(procDataFileID);
            strDay = ConvertDate_IOS_Manuscript2020(allDataFileDate);
            load(procDataFileID)
            puffs = ProcData.data.solenoids.LPadSol;
            if isempty(puffs) == true
                if strcmp(dataType,'CBV_HbT') == true
                    LH_AllData{zz,1} = ProcData.data.(dataType).adjLH;
                    RH_AllData{zz,1} = ProcData.data.(dataType).adjRH;
                else
                    LH_AllData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                    RH_AllData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                end
                zz = zz + 1;
            end
        end
        for gg = 1:length(LH_AllData)
            LH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,LH_AllData{gg,1}),'constant');
            RH_ProcAllData{gg,1} = detrend(filtfilt(sos,g,RH_AllData{gg,1}),'constant');
        end
        % analyze correlation coefficient between resting epochs
        for n = 1:length(LH_ProcAllData)
            all_CC = corrcoef(LH_ProcAllData{n,1},RH_ProcAllData{n,1});
            all_R(n,1) = all_CC(2,1);
        end
        meanAll_R = mean(all_R);
        stdAll_R = std(all_R,0,1);
        % save results
        AnalysisResults.(animalID).CorrCoeff.All.(dataType).R = all_R;
        AnalysisResults.(animalID).CorrCoeff.All.(dataType).meanR = meanAll_R;
        AnalysisResults.(animalID).CorrCoeff.All.(dataType).stdR = stdAll_R;
        
        %% Analyze Pearson's correlation coefficient during periods of NREM sleep
        % pull data from SleepData.mat structure
        if strcmp(dataType,'CBV_HbT') == true
            LH_nremData = SleepData.(modelType).NREM.data.(dataType).LH;
            RH_nremData = SleepData.(modelType).NREM.data.(dataType).RH;
        else
            LH_nremData = SleepData.(modelType).NREM.data.cortical_LH.(dataType);
            RH_nremData = SleepData.(modelType).NREM.data.cortical_RH.(dataType);
        end
        % detrend - data is already lowpass filtered
        for j = 1:length(LH_nremData)
            LH_nremData{j,1} = detrend(filtfilt(sos,g,LH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
            RH_nremData{j,1} = detrend(filtfilt(sos,g,RH_nremData{j,1}(1:params.minTime.NREM*samplingRate)),'constant');
        end
        % analyze correlation coefficient between NREM epochs
        for n = 1:length(LH_nremData)
            nrem_CC = corrcoef(LH_nremData{n,1},RH_nremData{n,1});
            nrem_R(n,1) = nrem_CC(2,1);
        end
        meanNREM_R = mean(nrem_R);
        stdNREM_R = std(nrem_R,0,1);
        % save results
        AnalysisResults.(animalID).CorrCoeff.NREM.(dataType).R = nrem_R;
        AnalysisResults.(animalID).CorrCoeff.NREM.(dataType).meanR = meanNREM_R;
        AnalysisResults.(animalID).CorrCoeff.NREM.(dataType).stdR = stdNREM_R;
        
        %% Analyze Pearson's correlation coefficient during periods of REM sleep
        % pull data from SleepData.mat structure
        if strcmp(dataType,'CBV_HbT') == true
            LH_remData = SleepData.(modelType).REM.data.(dataType).LH;
            RH_remData = SleepData.(modelType).REM.data.(dataType).RH;
        else
            LH_remData = SleepData.(modelType).REM.data.cortical_LH.(dataType);
            RH_remData = SleepData.(modelType).REM.data.cortical_RH.(dataType);
        end
        % detrend - data is already lowpass filtered
        for m = 1:length(LH_remData)
            LH_remData{m,1} = detrend(filtfilt(sos,g,LH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
            RH_remData{m,1} = detrend(filtfilt(sos,g,RH_remData{m,1}(1:params.minTime.REM*samplingRate)),'constant');
        end
        % analyze correlation coefficient between NREM epochs
        for n = 1:length(LH_remData)
            rem_CC = corrcoef(LH_remData{n,1},RH_remData{n,1});
            rem_R(n,1) = rem_CC(2,1);
        end
        meanREM_R = mean(rem_R);
        stdREM_R = std(rem_R,0,1);
        % save results
        AnalysisResults.(animalID).CorrCoeff.REM.(dataType).R = rem_R;
        AnalysisResults.(animalID).CorrCoeff.REM.(dataType).meanR = meanREM_R;
        AnalysisResults.(animalID).CorrCoeff.REM.(dataType).stdR = stdREM_R;
    end
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
