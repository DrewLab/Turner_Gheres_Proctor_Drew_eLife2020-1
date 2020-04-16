function [AnalysisResults] = AnalyzeMeanVesselDiameter_Manuscript2020(animalID,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyze the mean hemodynamic signal during each behavioral state
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T115','T116','T117','T118','T125','T126'};
modelType = 'Manual';
params.minTime.Rest = 10;   % seconds
params.minTime.Whisk = 7;   % seconds
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 60;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '\' animalID '\2P Data\'];
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
    samplingRate = RestData.vesselDiameter.data.samplingRate;
    % lowpass filter and detrend each segment
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    WhiskCriteria.Fieldname = {'duration','duration'};
    WhiskCriteria.Comparison = {'gt','lt'};
    WhiskCriteria.Value = {2,5};
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    
    %% Analyze mean CBV during periods of rest
    [restLogical] = FilterEvents_2P_Manuscript2020(RestData.vesselDiameter.data,RestCriteria);
    combRestLogical = logical(restLogical);
    restVesselData = RestData.vesselDiameter.data.data(combRestLogical,:);
    restFileIDs = RestData.vesselDiameter.data.fileIDs(combRestLogical,:);
    restVesselIDs = RestData.vesselDiameter.data.vesselIDs(combRestLogical,:);
    restDurations = RestData.vesselDiameter.data.durations(combRestLogical,:);
    restEventTimes = RestData.vesselDiameter.data.eventTimes(combRestLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [finalRestVesselData,finalRestFileIDs,finalRestVesselIDs,~,~] = RemoveInvalidData_2P_Manuscript2020(restVesselData,restFileIDs,restVesselIDs,restDurations,restEventTimes,ManualDecisions);
    % go through the data and normalize + filter each rest epoch based on individual vessels
    uniqueRestVesselIDs = unique(finalRestVesselIDs);
    for aa = 1:length(uniqueRestVesselIDs)
        cc = 1;
        for bb = 1:length(finalRestVesselIDs)
            if strcmp(uniqueRestVesselIDs{aa,1},finalRestVesselIDs{bb,1})
                strDay = ConvertDate_2P_Manuscript2020(finalRestFileIDs{bb,1}(1:6));
                tempRestData.(uniqueRestVesselIDs{aa,1}){cc,1} = filtfilt(sos,g,((finalRestVesselData{bb,1} - RestingBaselines.manualSelection.vesselDiameter.data.(uniqueRestVesselIDs{aa,1}).(strDay))/RestingBaselines.manualSelection.vesselDiameter.data.(uniqueRestVesselIDs{aa,1}).(strDay)));
                cc = cc + 1;
            end
        end
    end
    % take the average of each vessel's individual resting event
    for dd = 1:length(uniqueRestVesselIDs)
        for ee = 1:length(tempRestData.(uniqueRestVesselIDs{dd,1}))
            tempRestDataMeans.(uniqueRestVesselIDs{dd,1})(ee,1) = max(tempRestData.(uniqueRestVesselIDs{dd,1}){ee,1});
        end
    end
    % take the average of each vessel's total resting events
    for ff = 1:length(uniqueRestVesselIDs)
        AnalysisResults.(animalID).MeanVesselDiameter.Rest.(uniqueRestVesselIDs{ff,1}) = mean(tempRestDataMeans.(uniqueRestVesselIDs{ff,1}))*100;
    end
    
    %% Analyze mean CBV during periods of extended whisking
    [whiskLogical] = FilterEvents_2P_Manuscript2020(EventData.vesselDiameter.data.whisk,WhiskCriteria);
    combWhiskLogical = logical(whiskLogical);
    whiskVesselData = EventData.vesselDiameter.data.whisk.data(combWhiskLogical,:);
    whiskFileIDs = EventData.vesselDiameter.data.whisk.fileIDs(combWhiskLogical,:);
    whiskVesselIDs = EventData.vesselDiameter.data.whisk.vesselIDs(combWhiskLogical,:);
    whiskDurations = EventData.vesselDiameter.data.whisk.duration(combWhiskLogical,:);
    whiskEventTimes = EventData.vesselDiameter.data.whisk.eventTime(combWhiskLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [finalWhiskVesselData,finalWhiskFileIDs,finalWhiskVesselIDs,~,~] = RemoveInvalidData_2P_Manuscript2020(whiskVesselData,whiskFileIDs,whiskVesselIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    uniqueWhiskVesselIDs = unique(finalWhiskVesselIDs);
    for aa = 1:length(uniqueWhiskVesselIDs)
        cc = 1;
        for bb = 1:length(finalWhiskVesselIDs)
            if strcmp(uniqueWhiskVesselIDs{aa,1},finalWhiskVesselIDs{bb,1})
                strDay = ConvertDate_2P_Manuscript2020(finalWhiskFileIDs{bb,1}(1:6));
                tempWhiskData.(uniqueWhiskVesselIDs{aa,1}){cc,1} = filtfilt(sos,g,((finalWhiskVesselData(bb,:) - RestingBaselines.manualSelection.vesselDiameter.data.(uniqueWhiskVesselIDs{aa,1}).(strDay))/RestingBaselines.manualSelection.vesselDiameter.data.(uniqueWhiskVesselIDs{aa,1}).(strDay)));
                cc = cc + 1;
            end
        end
    end
    % take the average of each vessel's individual resting event
    for dd = 1:length(uniqueWhiskVesselIDs)
        for ee = 1:length(tempWhiskData.(uniqueWhiskVesselIDs{dd,1}))
            tempWhiskDataMeans.(uniqueWhiskVesselIDs{dd,1})(ee,1) = max(tempWhiskData.(uniqueWhiskVesselIDs{dd,1}){ee,1}(2*samplingRate:params.minTime.Whisk*samplingRate));
        end
    end
    % take the average of each vessel's total resting events
    for ff = 1:length(uniqueWhiskVesselIDs)
        AnalysisResults.(animalID).MeanVesselDiameter.Whisk.(uniqueWhiskVesselIDs{ff,1}) = mean(tempWhiskDataMeans.(uniqueWhiskVesselIDs{ff,1}))*100;
    end
    
    %% Analyze mean CBV during periods of NREM sleep
    if isfield(SleepData.(modelType),'NREM') == true
        % pull data from SleepData.mat structure
        nremVesselData = SleepData.(modelType).NREM.data.vesselDiameter.data;
        nremVesselIDs = SleepData.(modelType).NREM.VesselIDs;
        %
        uniqueNREMVesselIDs = unique(nremVesselIDs);
        for aa = 1:length(uniqueNREMVesselIDs)
            cc = 1;
            for bb = 1:length(nremVesselIDs)
                if strcmp(uniqueNREMVesselIDs{aa,1},nremVesselIDs{bb,1})
                    tempNREMData.(uniqueNREMVesselIDs{aa,1}){cc,1} = filtfilt(sos,g,nremVesselData{bb,1});
                    cc = cc + 1;
                end
            end
        end
        % take the average of each vessel's individual resting event
        for dd = 1:length(uniqueNREMVesselIDs)
            for ee = 1:length(tempNREMData.(uniqueNREMVesselIDs{dd,1}))
                tempNREMDataMeans.(uniqueNREMVesselIDs{dd,1})(ee,1) = max(tempNREMData.(uniqueNREMVesselIDs{dd,1}){ee,1});
            end
        end
        % take the average of each vessel's total resting events
        for ff = 1:length(uniqueNREMVesselIDs)
            AnalysisResults.(animalID).MeanVesselDiameter.NREM.(uniqueNREMVesselIDs{ff,1}) = mean(tempNREMDataMeans.(uniqueNREMVesselIDs{ff,1}))*100;
        end
    end
    
    %% Analyze mean CBV during periods of REM sleep
    if isfield(SleepData.(modelType),'REM') == true
        % pull data from SleepData.mat structure
        remVesselData = SleepData.(modelType).REM.data.vesselDiameter.data;
        remVesselIDs = SleepData.(modelType).REM.VesselIDs;
        %
        uniqueREMVesselIDs = unique(remVesselIDs);
        for aa = 1:length(uniqueREMVesselIDs)
            cc = 1;
            for bb = 1:length(remVesselIDs)
                if strcmp(uniqueREMVesselIDs{aa,1},remVesselIDs{bb,1})
                    tempREMData.(uniqueREMVesselIDs{aa,1}){cc,1} = filtfilt(sos,g,remVesselData{bb,1});
                    cc = cc + 1;
                end
            end
        end
        % take the average of each vessel's individual resting event
        for dd = 1:length(uniqueREMVesselIDs)
            for ee = 1:length(tempREMData.(uniqueREMVesselIDs{dd,1}))
                tempREMDataMeans.(uniqueREMVesselIDs{dd,1})(ee,1) = max(tempREMData.(uniqueREMVesselIDs{dd,1}){ee,1});
            end
        end
        % take the average of each vessel's total resting events
        for ff = 1:length(uniqueREMVesselIDs)
            AnalysisResults.(animalID).MeanVesselDiameter.REM.(uniqueREMVesselIDs{ff,1}) = mean(tempREMDataMeans.(uniqueREMVesselIDs{ff,1}))*100;
        end
    end
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end


end
