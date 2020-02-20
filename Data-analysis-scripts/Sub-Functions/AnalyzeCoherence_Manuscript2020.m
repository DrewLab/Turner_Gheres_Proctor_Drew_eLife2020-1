function [AnalysisResults] = AnalyzeCoherence_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Analyze the spectral coherence between bilateral hemodynamic and neural signals.
%________________________________________________________________________________________________________________________

%% function parameters
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120'};
dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
modelTypes = {'SVM','Ensemble','Forest','Manual'};
params.minTime.Rest = 10;   % seconds
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 30;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(IOS_animalIDs,animalID))
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
    samplingRate = RestData.CBV_HbT.adjLH.CBVCamSamplingRate;
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    PuffCriteria.Fieldname = {'puffDistances'};
    PuffCriteria.Comparison = {'gt'};
    PuffCriteria.Value = {5};
    % go through each valid data type for behavior-based coherence analysis
    for a = 1:length(dataTypes)
        dataType = dataTypes{1,a};
        
        %% Analyze coherence during periods of rest
        % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
        if strcmp(dataType,'CBV_HbT') == true
            [restLogical] = FilterEvents_IOS_Manuscript2020(RestData.(dataType).adjLH,RestCriteria);
            [puffLogical] = FilterEvents_IOS_Manuscript2020(RestData.(dataType).adjLH,PuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            restFileIDs = RestData.(dataType).adjLH.fileIDs(combRestLogical,:);
            restEventTimes = RestData.(dataType).adjLH.eventTimes(combRestLogical,:);
            restDurations = RestData.(dataType).adjLH.durations(combRestLogical,:);
            LH_unstimRestingData = RestData.(dataType).adjLH.data(combRestLogical,:);
            RH_unstimRestingData = RestData.(dataType).adjRH.data(combRestLogical,:);
        else
            [restLogical] = FilterEvents_IOS_Manuscript2020(RestData.cortical_LH.(dataType),RestCriteria);
            [puffLogical] = FilterEvents_IOS_Manuscript2020(RestData.cortical_LH.(dataType),PuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            restFileIDs = RestData.cortical_LH.(dataType).fileIDs(combRestLogical,:);
            restEventTimes = RestData.cortical_LH.(dataType).eventTimes(combRestLogical,:);
            restDurations = RestData.cortical_LH.(dataType).durations(combRestLogical,:);
            LH_unstimRestingData = RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
            RH_unstimRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
        end
        % decimate the file list to only include those files that occur within the desired number of target minutes        
        [LH_finalRestData,~,~,~] = DecimateRestData_Manuscript2020(LH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [RH_finalRestData,~,~,~] = DecimateRestData_Manuscript2020(RH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
        % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
        % lowpass filter and detrend each segment
        [B,A] = butter(3,1/(samplingRate/2),'low');
        clear LH_ProcRestData RH_ProcRestData
        for g = 1:length(LH_finalRestData)
            if length(LH_finalRestData{g,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(LH_finalRestData{g,1});
                LH_restPad = (ones(1,restChunkSampleDiff))*LH_finalRestData{g,1}(end);
                RH_restPad = (ones(1,restChunkSampleDiff))*RH_finalRestData{g,1}(end);
                LH_ProcRestData{g,1} = horzcat(LH_finalRestData{g,1},LH_restPad); %#ok<*AGROW>
                RH_ProcRestData{g,1} = horzcat(RH_finalRestData{g,1},RH_restPad);
                LH_ProcRestData{g,1} = detrend(filtfilt(B,A,LH_ProcRestData{g,1}),'constant');
                RH_ProcRestData{g,1} = detrend(filtfilt(B,A,RH_ProcRestData{g,1}),'constant');
            else
                LH_ProcRestData{g,1} = detrend(filtfilt(B,A,LH_finalRestData{g,1}(1:(params.minTime.Rest*samplingRate))),'constant');
                RH_ProcRestData{g,1} = detrend(filtfilt(B,A,RH_finalRestData{g,1}(1:(params.minTime.Rest*samplingRate))),'constant');
            end
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        LH_restData = zeros(length(LH_ProcRestData{1,1}),length(LH_ProcRestData));
        RH_restData = zeros(length(RH_ProcRestData{1,1}),length(RH_ProcRestData));
        for h = 1:length(LH_ProcRestData)
            LH_restData(:,h) = LH_ProcRestData{h,1};
            RH_restData(:,h) = RH_ProcRestData{h,1};
        end
        % parameters for coherencyc - information available in function
        params.tapers = [3,5];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;   % Sampling Rate
        params.fpass = [0,1];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the coherence between desired signals
        [C_RestData,~,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc_Manuscript2020(LH_restData,RH_restData,params);
        % save data and figures
        AnalysisResults.(animalID).Coherence.Rest.(dataType).C = C_RestData;
        AnalysisResults.(animalID).Coherence.Rest.(dataType).f = f_RestData;
        AnalysisResults.(animalID).Coherence.Rest.(dataType).confC = confC_RestData;
        AnalysisResults.(animalID).Coherence.Rest.(dataType).cErr = cErr_RestData;
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            restCoherence = figure;
            plot(f_RestData,C_RestData,'k')
            hold on;
            plot(f_RestData,cErr_RestData,'color',colors_Manuscript2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Coherence');
            title([animalID  ' ' dataType ' coherence for resting data']);
            set(gca,'Ticklength',[0,0]);
            legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
            set(legend,'FontSize',6);
            ylim([0,1])
            xlim([0,1])
            axis square
            [pathstr,~,~] = fileparts(cd);
            dirpath = [pathstr '/Figures/Coherence/'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(restCoherence,[dirpath animalID '_Rest_' dataType '_Coherence']);
            close(restCoherence)
        end
        
        %% Analyze coherence during periods of NREM sleep
        % pull data from SleepData.mat structure
        for xx = 1:length(modelTypes)
            modelType = modelTypes{1,xx};
            if strcmp(dataType,'CBV_HbT') == true
                LH_nremData = SleepData.(modelType).NREM.data.(dataType).LH;
                RH_nremData = SleepData.(modelType).NREM.data.(dataType).RH;
            else
                LH_nremData = SleepData.(modelType).NREM.data.cortical_LH.(dataType);
                RH_nremData = SleepData.(modelType).NREM.data.cortical_RH.(dataType);
            end
            % detrend - data is already lowpass filtered
            for j = 1:length(LH_nremData)
                LH_nremData{j,1} = detrend(LH_nremData{j,1}(1:(params.minTime.NREM*samplingRate)),'constant');
                RH_nremData{j,1} = detrend(RH_nremData{j,1}(1:(params.minTime.NREM*samplingRate)),'constant');
            end
            % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            LH_nrem = zeros(length(LH_nremData{1,1}),length(LH_nremData));
            RH_nrem = zeros(length(RH_nremData{1,1}),length(RH_nremData));
            for k = 1:length(LH_nremData)
                LH_nrem(:,k) = LH_nremData{k,1};
                RH_nrem(:,k) = RH_nremData{k,1};
            end
            % calculate the coherence between desired signals
            [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc_Manuscript2020(LH_nrem,RH_nrem,params);
            % save data and figures
            AnalysisResults.(animalID).Coherence.NREM.(modelType).(dataType).C = C_nrem;
            AnalysisResults.(animalID).Coherence.NREM.(modelType).(dataType).f = f_nrem;
            AnalysisResults.(animalID).Coherence.NREM.(modelType).(dataType).confC = confC_nrem;
            AnalysisResults.(animalID).Coherence.NREM.(modelType).(dataType).cErr = cErr_nrem;
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                nremCoherence = figure;
                plot(f_nrem,C_nrem,'k')
                hold on;
                plot(f_nrem,cErr_nrem,'color',colors_Manuscript2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Coherence');
                title([animalID  ' ' dataType ' coherence for ' modelType ' NREM data']);
                set(gca,'Ticklength',[0,0]);
                legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
                set(legend,'FontSize',6);
                ylim([0,1])
                xlim([0,1])
                axis square
                savefig(nremCoherence,[dirpath animalID '_' modelType '_NREM_' dataType '_Coherence']);
                close(nremCoherence)
            end
            
            %% Analyze coherence during periods of REM sleep
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
                LH_remData{m,1} = detrend(LH_remData{m,1}(1:(params.minTime.REM*samplingRate)),'constant');
                RH_remData{m,1} = detrend(RH_remData{m,1}(1:(params.minTime.REM*samplingRate)),'constant');
            end
            % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            LH_rem = zeros(length(LH_remData{1,1}),length(LH_remData));
            RH_rem = zeros(length(RH_remData{1,1}),length(RH_remData));
            for n = 1:length(LH_remData)
                LH_rem(:,n) = LH_remData{n,1};
                RH_rem(:,n) = RH_remData{n,1};
            end
            % calculate the coherence between desired signals
            [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc_Manuscript2020(LH_rem,RH_rem,params);
            % save data and figures
            AnalysisResults.(animalID).Coherence.REM.(modelType).(dataType).C = C_rem;
            AnalysisResults.(animalID).Coherence.REM.(modelType).(dataType).f = f_rem;
            AnalysisResults.(animalID).Coherence.REM.(modelType).(dataType).confC = confC_rem;
            AnalysisResults.(animalID).Coherence.REM.(modelType).(dataType).cErr = cErr_rem;
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                remCoherence = figure;
                plot(f_rem,C_rem,'k')
                hold on;
                plot(f_rem,cErr_rem,'color',colors_Manuscript2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Coherence');
                title([animalID  ' ' dataType ' coherence for ' modelType 'REM data']);
                set(gca,'Ticklength',[0,0]);
                legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
                set(legend,'FontSize',6);
                ylim([0,1])
                xlim([0,1])
                axis square
                savefig(remCoherence,[dirpath animalID '_' modelType '_REM_' dataType '_Coherence']);
                close(remCoherence)
            end
        end
    end
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
