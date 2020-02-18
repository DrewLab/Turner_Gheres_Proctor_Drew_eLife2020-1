function [AnalysisResults] = AnalyzeLaserDoppler_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner%
%
%   Purpose:
%________________________________________________________________________________________________________________________

%% function parameters
IOS_animalIDs = {'T108','T109','T110','T111','T119','T120'};
modelType = 'SVM';
params.minTime.Rest = 10;   % seconds
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 30;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(IOS_animalIDs,animalID))
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
    % determine the animal's ID use the EventData.mat file's name for the current folder
    fileBreaks = strfind(eventDataFileID,'_');
    animalID = eventDataFileID(1:fileBreaks(1)-1);
    
    %% Whisking-evoked responses
    % criteria for the FilterEvents data struct
    whiskCriteriaA.Fieldname = {'duration','duration','puffDistance'};
    whiskCriteriaA.Comparison = {'gt','lt','gt'};
    whiskCriteriaA.Value = {0.5,2,5};
    whiskCriteriaB.Fieldname = {'duration','duration','puffDistance'};
    whiskCriteriaB.Comparison = {'gt','lt','gt'};
    whiskCriteriaB.Value = {2,5,5};
    whiskCriteriaC.Fieldname = {'duration','puffDistance'};
    whiskCriteriaC.Comparison = {'gt','gt'};
    whiskCriteriaC.Value = {5,5};
    PuffCriteria.Fieldname = {'puffDistance'};
    PuffCriteria.Comparison = {'gt'};
    PuffCriteria.Value = {5};
    whiskCriteriaNames = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
    % filter the EventData.mat structure for whisking events that meet the desired criteria
    % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
    samplingRate = EventData.flow.data.whisk.samplingRate;
    timeVector = (0:(EventData.flow.data.whisk.epoch.duration*samplingRate))/samplingRate - EventData.flow.data.whisk.epoch.offset;
    offset = EventData.flow.data.whisk.epoch.offset;
    for b = 1:length(whiskCriteriaNames)
        whiskCriteriaName = whiskCriteriaNames{1,b};
        if strcmp(whiskCriteriaName,'ShortWhisks') == true
            WhiskCriteria = whiskCriteriaA;
        elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
            WhiskCriteria = whiskCriteriaB;
        elseif strcmp(whiskCriteriaName,'LongWhisks') == true
            WhiskCriteria = whiskCriteriaC;
        end
        [whiskLogical] = FilterEvents_IOS(EventData.flow.data.whisk,WhiskCriteria);
        [puffLogical] = FilterEvents_IOS(EventData.flow.data.whisk,PuffCriteria);
        combWhiskLogical = logical(whiskLogical.*puffLogical);
        [allWhiskFlowData] = EventData.flow.data.whisk.NormData(combWhiskLogical,:);
        [allWhiskFileIDs] = EventData.flow.data.whisk.fileIDs(combWhiskLogical,:);
        [allWhiskEventTimes] = EventData.flow.data.whisk.eventTime(combWhiskLogical,:);
        allWhiskDurations = EventData.flow.data.whisk.eventTime(combWhiskLogical,:);
        % decimate the file list to only include those files that occur within the desired number of target minutes
        [finalWhiskFlowData,~,~,~] = DecimateRestData_Manuscript2020(allWhiskFlowData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        clear procWhiskFlowData
        for f = 1:size(finalWhiskFlowData,1)
            whiskFlowArray = finalWhiskFlowData(f,:);
            filtWhiskFlowArray = sgolayfilt(whiskFlowArray,3,17);
            procWhiskFlowData(f,:) = filtWhiskFlowArray - mean(filtWhiskFlowArray(1:(offset*samplingRate))); %#ok<*AGROW>
        end
        meanWhiskFlowData = mean(procWhiskFlowData,1)*100;
        stdWhiskFlowData = std(procWhiskFlowData,0,1)*100;
        % save results
        AnalysisResults.(animalID).EvokedAvgs.Whisk.DopplerFlow.(whiskCriteriaName).flowMean = meanWhiskFlowData;
        AnalysisResults.(animalID).EvokedAvgs.Whisk.DopplerFlow.(whiskCriteriaName).flowStD = stdWhiskFlowData;
        AnalysisResults.(animalID).EvokedAvgs.Whisk.DopplerFlow.(whiskCriteriaName).timeVector = timeVector;
        % Save figures if desired
        if strcmp(saveFigs,'y') == true
            whiskEvoked = figure;
            plot(timeVector,meanWhiskFlowData,'k')
            hold on
            plot(timeVector,meanWhiskFlowData + stdWhiskFlowData,'color',colors_Manuscript2020('battleship grey'))
            plot(timeVector,meanWhiskFlowData - stdWhiskFlowData,'color',colors_Manuscript2020('battleship grey'))
            title([animalID ' whisking-evoked flow averages - ' whiskCriteriaName])
            xlabel('Time (sec)')
            ylabel('Perfusion Units (%)')
            axis tight
            axis square
            % save figure
            [pathstr,~,~] = fileparts(cd);
            dirpath = [pathstr '/Figures/Evoked Responses/'];
            if ~exist(dirpath, 'dir')
                mkdir(dirpath);
            end
            savefig(whiskEvoked, [dirpath animalID '_DopplerFlow_' whiskCriteriaName '_WhiskEvokedAverages']);
            close(whiskEvoked)
        end
    end
    
    %% Stimulus-evoked responses
    % Criteria for the FilterEvents data struct
    stimCriteriaA.Value = {'LPadSol'};
    stimCriteriaA.Fieldname = {'solenoidName'};
    stimCriteriaA.Comparison = {'equal'};
    stimCriteriaB.Value = {'RPadSol'};
    stimCriteriaB.Fieldname = {'solenoidName'};
    stimCriteriaB.Comparison = {'equal'};
    stimCriteriaC.Value = {'AudSol'};
    stimCriteriaC.Fieldname = {'solenoidName'};
    stimCriteriaC.Comparison = {'equal'};
    stimCriteriaNames = {'stimCriteriaA','stimCriteriaB','stimCriteriaC'};
    % filter the EventData.mat structure for stimulus events that meet the desired criteria
    for j = 1:length(stimCriteriaNames)
        stimCriteriaName = stimCriteriaNames{1,j};
        if strcmp(stimCriteriaName,'stimCriteriaA') == true
            stimCriteria = stimCriteriaA;
            solenoid = 'LPadSol';
        elseif strcmp(stimCriteriaName,'stimCriteriaB') == true
            stimCriteria = stimCriteriaB;
            solenoid = 'RPadSol';
        elseif strcmp(stimCriteriaName,'stimCriteriaC') == true
            stimCriteria = stimCriteriaC;
            solenoid = 'AudSol';
        end
        allStimFilter = FilterEvents_IOS(EventData.flow.data.stim,stimCriteria);
        [allStimFlowData] = EventData.flow.data.stim.data(allStimFilter,:);
        [allStimFileIDs] = EventData.flow.data.stim.fileIDs(allStimFilter,:);
        [allStimEventTimes] = EventData.flow.data.stim.eventTime(allStimFilter,:);
        allStimDurations = zeros(length(allStimEventTimes),1);
        % decimate the file list to only include those files that occur within the desired number of target minutes
        [finalStimFlowData,~,~,~] = DecimateRestData_Manuscript2020(allStimFlowData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
        % lowpass filter each whisking event and mean-subtract by the first 2 seconds
        clear procStimFlowData
        for p = 1:size(finalStimFlowData,1)
            stimFlowArray = finalStimFlowData(p,:);
            filtStimFlowArray = sgolayfilt(stimFlowArray,3,17);
            procStimFlowData(p,:) = filtStimFlowArray - mean(filtStimFlowArray(1:(offset*samplingRate)));
        end
        meanStimFlowData = mean(procStimFlowData,1)*100;
        stdStimFlowData = std(procStimFlowData,0,1)*100;
        % save results
        AnalysisResults.(animalID).EvokedAvgs.Stim.DopperFlow.(solenoid).flowMean = meanStimFlowData;
        AnalysisResults.(animalID).EvokedAvgs.Stim.DopperFlow.(solenoid).flowStD = stdStimFlowData;
        AnalysisResults.(animalID).EvokedAvgs.Stim.DopperFlow.(solenoid).timeVector = timeVector;
        % Save figures if desired
        if strcmp(saveFigs,'y') == true
            % summary figure
            stimEvoked = figure;
            plot(timeVector,meanStimFlowData,'k')
            hold on
            plot(timeVector,meanStimFlowData + stdStimFlowData,'color',colors_Manuscript2020('battleship grey'))
            plot(timeVector,meanStimFlowData - stdStimFlowData,'color',colors_Manuscript2020('battleship grey'))
            title([animalID ' stimulus-evoked flow averages - ' solenoid])
            xlabel('Time (sec)')
            ylabel('Perfusion Units (%)')
            axis tight
            axis square
            % save figure
            [pathstr,~,~] = fileparts(cd);
            dirpath = [pathstr '/Figures/Evoked Responses/'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(stimEvoked,[dirpath animalID '_DopperFlow_' solenoid '_StimEvokedAverages']);
            close(stimEvoked)
        end
    end
    
    %% Average behavior-dependent flow
    % Analyze mean laser doppler flow during periods of rest
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    [restLogical] = FilterEvents_IOS(RestData.flow.data,RestCriteria);
    [puffLogical] = FilterEvents_IOS(RestData.flow.data,PuffCriteria);
    combRestLogical = logical(restLogical.*puffLogical);
    restFileIDs = RestData.flow.data.fileIDs(combRestLogical,:);
    restFlowData = RestData.flow.data.NormData(combRestLogical,:);
    restEventTimes = RestData.flow.data.eventTimes(combRestLogical,:);
    restDurations = RestData.flow.data.durations(combRestLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [finalRestFlowData,~,~,~] = DecimateRestData_Manuscript2020(restFlowData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    % lowpass filter and detrend each segment
    [B,A] = butter(3,1/(samplingRate/2),'low');
    clear procRestData
    for g = 1:length(finalRestFlowData)
        procRestData{g,1} = filtfilt(B,A,finalRestFlowData{g,1});
    end
    % analyze correlation coefficient between resting epochs
    for n = 1:length(procRestData)
        restFlowMean(n,1) = mean(procRestData{n,1});
    end
    % save results
    AnalysisResults.(animalID).LDFlow.Rest = restFlowMean;
    
    %% Analyze mean CBV during periods of extended whisking
    % criteria for the FilterEvents data struct
    [whiskLogical] = FilterEvents_IOS(EventData.flow.data.whisk,whiskCriteriaC);
    [puffLogical] = FilterEvents_IOS(EventData.flow.data.whisk,PuffCriteria);
    combWhiskLogical = logical(whiskLogical.*puffLogical);
    whiskFlowData = EventData.flow.data.whisk.NormData(combWhiskLogical,:);
    whiskFileIDs = EventData.flow.data.whisk.fileIDs(combWhiskLogical,:);
    whiskEventTimes = EventData.flow.data.whisk.eventTimes(combWhiskLogical,:);
    whiskDurations = EventData.flow.data.whisk.durations(combWhiskLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [finalWhiskData,~,~,~] = DecimateRestData_Manuscript2020(whiskFlowData,whiskFileIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    % lowpass filter and detrend each segment
    [B,A] = butter(3,1/(samplingRate/2),'low');
    clear procWhiskData
    for g = 1:size(finalWhiskData,1)
        procWhiskData(g,:) = filtfilt(B,A,finalWhiskData(g,:));
    end
    % analyze correlation coefficient between resting epochs
    for n = 1:size(procWhiskData,1)
        whiskFlowMean{n,1} = mean(procWhiskData(n,samplingRate*2:end),2);
    end
    % save results
    AnalysisResults.(animalID).LDFlow.Whisk = cell2mat(whiskFlowMean);
    
    %% Analyze mean CBV during periods of NREM sleep
    % pull data from SleepData.mat structure
    nremData = SleepData.(modelType).NREM.data.DopplerFlow;
    % analyze correlation coefficient between NREM epochs
    for n = 1:length(nremData)
        nremFlowMean(n,1) = mean(nremData{n,1});
    end
    % save results
    AnalysisResults.(animalID).LDFlow.NREM = nremFlowMean;
    
    %% Analyze mean CBV during periods of REM sleep
    % pull data from SleepData.mat structure
    remData = SleepData.(modelType).REM.data.DopplerFlow;
    % analyze correlation coefficient between NREM epochs
    for n = 1:length(remData)
        remFlowMean(n,1) = mean(remData{n,1});
    end
    % save results
    AnalysisResults.(animalID).LDFlow.REM = remFlowMean;
    
    %% Analyze mean CBV during periods of Isolfurane
    dataLocation = [rootFolder '/' animalID '/Isoflurane Trials/'];
    cd(dataLocation)
    % pull ProcData file
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFile = {procDataFileStruct.name}';
    procDataFileID = char(procDataFile);
    load(procDataFileID)
    isoFlow = ProcData.data.flow((end - samplingRate*100):end);
    normIsoFlow = (isoFlow - RestingBaselines.(baselineType).flow.data.(strDay))./(RestingBaselines.(baselineType).flow.data.(strDay));
    filtIsoFlow = filtfilt(D,C,normIsoFlow)*100;
    filtIsoFlow = filtIsoFlow - mean(filtIsoFlow(1:300*samplingRate));
    AnalysisResults.(animalID).LDFlow.REM = mean(filtIsoFlow);
end
cd(rootFolder)

end
