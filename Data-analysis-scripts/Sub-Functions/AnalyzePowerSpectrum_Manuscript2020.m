function [AnalysisResults] = AnalyzePowerSpectrum_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Analyze the spectral power of hemodynamic and neural signals.
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
modelType = 'Forest';
params.minTime.Rest = 10;   % seconds
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 60;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    % character list of all ProcData file IDs
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
    % find and load Forest_ScoringResults.mat struct
    forestScoringResultsFileID = 'Forest_ScoringResults.mat';
    load(forestScoringResultsFileID,'-mat')
    % identify animal's ID and pull important infortmat
    fileBreaks = strfind(restDataFileID, '_');
    animalID = restDataFileID(1:fileBreaks(1)-1);
    samplingRate = RestData.CBV_HbT.adjLH.CBVCamSamplingRate;
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    PuffCriteria.Fieldname = {'puffDistances'};
    PuffCriteria.Comparison = {'gt'};
    PuffCriteria.Value = {5};
    % lowpass filter and detrend each segment
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % go through each valid data type for behavior-based power spectrum analysis
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        
        %% Analyze power spectra during periods of rest
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
            LH_unstimRestingData =RestData.cortical_LH.(dataType).NormData(combRestLogical,:);
            RH_unstimRestingData = RestData.cortical_RH.(dataType).NormData(combRestLogical,:);
            Hip_unstimRestingData = RestData.hippocampus.(dataType).NormData(combRestLogical,:);
        end
        % decimate the file list to only include those files that occur within the desired number of target minutes
        [LH_finalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(LH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [RH_finalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(RH_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        if strcmp(dataType,'CBV_HbT') == false
            [Hip_finalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(Hip_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        end
        % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
        % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
        clear LH_ProcRestData
        clear RH_ProcRestData
        clear Hip_ProcRestData
        for bb = 1:length(LH_finalRestData)
            if length(LH_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                restChunkSampleDiff = params.minTime.Rest*samplingRate - length(LH_finalRestData{bb,1});
                LH_restPad = (ones(1,restChunkSampleDiff))*LH_finalRestData{bb,1}(end);
                RH_restPad = (ones(1,restChunkSampleDiff))*RH_finalRestData{bb,1}(end);
                LH_ProcRestData{bb,1} = horzcat(LH_finalRestData{bb,1},LH_restPad); %#ok<*AGROW>
                RH_ProcRestData{bb,1} = horzcat(RH_finalRestData{bb,1},RH_restPad);
                LH_ProcRestData{bb,1} = filtfilt(sos,g,detrend(LH_ProcRestData{bb,1},'constant'));
                RH_ProcRestData{bb,1} = filtfilt(sos,g,detrend(RH_ProcRestData{bb,1},'constant'));
                if strcmp(dataType,'CBV_HbT') == false
                    Hip_restPad = (ones(1,restChunkSampleDiff))*Hip_finalRestData{bb,1}(end);
                    Hip_ProcRestData{bb,1} = horzcat(Hip_finalRestData{bb,1},Hip_restPad);
                    Hip_ProcRestData{bb,1} = filtfilt(sos,g,detrend(Hip_ProcRestData{bb,1},'constant'));
                end
            else
                LH_ProcRestData{bb,1} = filtfilt(sos,g,detrend(LH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                RH_ProcRestData{bb,1} = filtfilt(sos,g,detrend(RH_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                if strcmp(dataType,'CBV_HbT') == false
                    Hip_ProcRestData{bb,1} = filtfilt(sos,g,detrend(Hip_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                end
            end
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        LH_restData = zeros(length(LH_ProcRestData{1,1}),length(LH_ProcRestData));
        RH_restData = zeros(length(RH_ProcRestData{1,1}),length(RH_ProcRestData));
        if strcmp(dataType,'CBV_HbT') == false
            Hip_restData = zeros(length(Hip_ProcRestData{1,1}),length(Hip_ProcRestData));
        end
        for cc = 1:length(LH_ProcRestData)
            LH_restData(:,cc) = LH_ProcRestData{cc,1};
            RH_restData(:,cc) = RH_ProcRestData{cc,1};
            if strcmp(dataType,'CBV_HbT') == false
                Hip_restData(:,cc) = Hip_ProcRestData{cc,1};
            end
        end
        % parameters for mtspectrumc - information available in function
        params.tapers = [1,1];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;   % Sampling Rate
        params.fpass = [0,0.5];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the power spectra of the desired signals
        [LH_rest_S,LH_rest_f,LH_rest_sErr] = mtspectrumc_Manuscript2020(LH_restData,params);
        [RH_rest_S,RH_rest_f,RH_rest_sErr] = mtspectrumc_Manuscript2020(RH_restData,params);
        if strcmp(dataType,'CBV_HbT') == false
            [Hip_rest_S,Hip_rest_f,Hip_rest_sErr] = mtspectrumc_Manuscript2020(Hip_restData,params);
        end
        % save data and figures
        AnalysisResults.(animalID).PowerSpectra.Rest.(dataType).adjLH.S = LH_rest_S;
        AnalysisResults.(animalID).PowerSpectra.Rest.(dataType).adjLH.f = LH_rest_f;
        AnalysisResults.(animalID).PowerSpectra.Rest.(dataType).adjLH.sErr = LH_rest_sErr;
        AnalysisResults.(animalID).PowerSpectra.Rest.(dataType).adjRH.S = RH_rest_S;
        AnalysisResults.(animalID).PowerSpectra.Rest.(dataType).adjRH.f = RH_rest_f;
        AnalysisResults.(animalID).PowerSpectra.Rest.(dataType).adjRH.sErr = RH_rest_sErr;
        if strcmp(dataType,'CBV_HbT') == false
            AnalysisResults.(animalID).PowerSpectra.Rest.(dataType).Hip.S = Hip_rest_S;
            AnalysisResults.(animalID).PowerSpectra.Rest.(dataType).Hip.f = Hip_rest_f;
            AnalysisResults.(animalID).PowerSpectra.Rest.(dataType).Hip.sErr = Hip_rest_sErr;
        end
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            % awake rest summary figures
            LH_RestPower = figure;
            loglog(LH_rest_f,LH_rest_S,'k')
            hold on;
            loglog(LH_rest_f,LH_rest_sErr,'color',colors_Manuscript2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Power');
            title([animalID  ' adjLH ' dataType ' Power during awake rest']);
            set(gca,'Ticklength',[0,0]);
            legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
            set(legend,'FontSize',6);
            xlim([0.1,0.5])
            axis square
            set(gca,'box','off')
            RH_RestPower = figure;
            loglog(RH_rest_f,RH_rest_S,'k')
            hold on;
            loglog(RH_rest_f,RH_rest_sErr,'color',colors_Manuscript2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Power');
            title([animalID  ' adjRH ' dataType ' Power during awake rest']);
            set(gca,'Ticklength',[0,0]);
            legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
            set(legend,'FontSize',6);
            xlim([0.1,0.5])
            axis square
            set(gca,'box','off')
            if strcmp(dataType,'CBV_HbT') == false
                Hip_RestPower = figure;
                loglog(Hip_rest_f,Hip_rest_S,'k')
                hold on;
                loglog(Hip_rest_f,Hip_rest_sErr,'color',colors_Manuscript2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' Hippocampal ' dataType ' Power during awake rest']);
                set(gca,'Ticklength',[0,0]);
                legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
                set(legend,'FontSize',6);
                xlim([0.1,0.5])
                axis square
                set(gca,'box','off')
            end
            [pathstr, ~, ~] = fileparts(cd);
            dirpath = [pathstr '/Figures/Power Spectrum/'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(LH_RestPower,[dirpath animalID '_Rest_LH_' dataType '_PowerSpectra']);
            close(LH_RestPower)
            savefig(RH_RestPower,[dirpath animalID '_Rest_RH_' dataType '_PowerSpectra']);
            close(RH_RestPower)
            if strcmp(dataType,'CBV_HbT') == false
                savefig(Hip_RestPower,[dirpath animalID '_Rest_Hippocampal_' dataType '_PowerSpectra']);
                close(Hip_RestPower)
            end
        end
        
        %% Analyze coherence during awake periods with no sleep scores
        zz = 1;
        clear LH_AwakeData RH_AwakeData Hip_AwakeData LH_ProcAwakeData RH_ProcAwakeData Hip_ProcAwakeData
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
            if sum(strcmp(scoringLabels,'Not Sleep')) > 170   % 6 bins (180 total) or 30 seconds of sleep
                load(procDataFileID)
                if strcmp(dataType,'CBV_HbT') == true
                    LH_AwakeData{zz,1} = ProcData.data.(dataType).adjLH;
                    RH_AwakeData{zz,1} = ProcData.data.(dataType).adjRH;
                else
                    LH_AwakeData{zz,1} = (ProcData.data.cortical_LH.(dataType) - RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_LH.(dataType).(strDay);
                    RH_AwakeData{zz,1} = (ProcData.data.cortical_RH.(dataType) - RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay))./RestingBaselines.manualSelection.cortical_RH.(dataType).(strDay);
                    Hip_AwakeData{zz,1} = (ProcData.data.hippocampus.(dataType) - RestingBaselines.manualSelection.hippocampus.(dataType).(strDay))./RestingBaselines.manualSelection.hippocampus.(dataType).(strDay);
                end
                zz = zz + 1;
            end
        end
        if isempty(LH_AwakeData) == false
            % process
            for bb = 1:length(LH_AwakeData)
                LH_ProcAwakeData{bb,1} = filtfilt(sos,g,detrend(LH_AwakeData{bb,1},'constant'));
                RH_ProcAwakeData{bb,1} = filtfilt(sos,g,detrend(RH_AwakeData{bb,1},'constant'));
                if strcmp(dataType,'CBV_HbT') == false
                    Hip_ProcAwakeData{bb,1} = filtfilt(sos,g,detrend(Hip_AwakeData{bb,1},'constant'));
                end
            end
            % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            LH_awakeData = zeros(length(LH_ProcAwakeData{1,1}),length(LH_ProcAwakeData));
            RH_awakeData = zeros(length(RH_ProcAwakeData{1,1}),length(RH_ProcAwakeData));
            if strcmp(dataType,'CBV_HbT') == false
                Hip_awakeData = zeros(length(Hip_ProcAwakeData{1,1}),length(Hip_ProcAwakeData));
            end
            for cc = 1:length(LH_ProcAwakeData)
                LH_awakeData(:,cc) = LH_ProcAwakeData{cc,1};
                RH_awakeData(:,cc) = RH_ProcAwakeData{cc,1};
                if strcmp(dataType,'CBV_HbT') == false
                    Hip_awakeData(:,cc) = Hip_ProcAwakeData{cc,1};
                end
            end
            % parameters for mtspectrumc - information available in function
            params.tapers = [5,9];   % Tapers [n, 2n - 1]
            params.pad = 1;
            params.Fs = samplingRate;   % Sampling Rate
            params.fpass = [0,0.5];   % Pass band [0, nyquist]
            params.trialave = 1;
            params.err = [2,0.05];
            % calculate the power spectra of the desired signals
            [LH_awake_S,LH_awake_f,LH_awake_sErr] = mtspectrumc_Manuscript2020(LH_awakeData,params);
            [RH_awake_S,RH_awake_f,RH_awake_sErr] = mtspectrumc_Manuscript2020(RH_awakeData,params);
            if strcmp(dataType,'CBV_HbT') == false
                [Hip_awake_S,Hip_awake_f,Hip_awake_sErr] = mtspectrumc_Manuscript2020(Hip_awakeData,params);
            end
            % save data and figures
            AnalysisResults.(animalID).PowerSpectra.Awake.(dataType).adjLH.S = LH_awake_S;
            AnalysisResults.(animalID).PowerSpectra.Awake.(dataType).adjLH.f = LH_awake_f;
            AnalysisResults.(animalID).PowerSpectra.Awake.(dataType).adjLH.sErr = LH_awake_sErr;
            AnalysisResults.(animalID).PowerSpectra.Awake.(dataType).adjRH.S = RH_awake_S;
            AnalysisResults.(animalID).PowerSpectra.Awake.(dataType).adjRH.f = RH_awake_f;
            AnalysisResults.(animalID).PowerSpectra.Awake.(dataType).adjRH.sErr = RH_awake_sErr;
            if strcmp(dataType,'CBV_HbT') == false
                AnalysisResults.(animalID).PowerSpectra.Awake.(dataType).Hip.S = Hip_awake_S;
                AnalysisResults.(animalID).PowerSpectra.Awake.(dataType).Hip.f = Hip_awake_f;
                AnalysisResults.(animalID).PowerSpectra.Awake.(dataType).Hip.sErr = Hip_awake_sErr;
            end
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                % awake awake summary figures
                LH_AwakePower = figure;
                loglog(LH_awake_f,LH_awake_S,'k')
                hold on;
                loglog(LH_awake_f,LH_awake_sErr,'color',colors_Manuscript2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' adjLH ' dataType ' Power during awake awake']);
                set(gca,'Ticklength',[0,0]);
                legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
                set(legend,'FontSize',6);
                xlim([0.1,0.5])
                axis square
                set(gca,'box','off')
                RH_AwakePower = figure;
                loglog(RH_awake_f,RH_awake_S,'k')
                hold on;
                loglog(RH_awake_f,RH_awake_sErr,'color',colors_Manuscript2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' adjRH ' dataType ' Power during awake awake']);
                set(gca,'Ticklength',[0,0]);
                legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
                set(legend,'FontSize',6);
                xlim([0.1,0.5])
                axis square
                set(gca,'box','off')
                if strcmp(dataType,'CBV_HbT') == false
                    Hip_AwakePower = figure;
                    loglog(Hip_awake_f,Hip_awake_S,'k')
                    hold on;
                    loglog(Hip_awake_f,Hip_awake_sErr,'color',colors_Manuscript2020('battleship grey'))
                    xlabel('Freq (Hz)');
                    ylabel('Power');
                    title([animalID  ' Hippocampal ' dataType ' Power during awake awake']);
                    set(gca,'Ticklength',[0,0]);
                    legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
                    set(legend,'FontSize',6);
                    xlim([0.1,0.5])
                    axis square
                    set(gca,'box','off')
                end
                [pathstr, ~, ~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Power Spectrum/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(LH_AwakePower,[dirpath animalID '_Awake_LH_' dataType '_PowerSpectra']);
                close(LH_AwakePower)
                savefig(RH_AwakePower,[dirpath animalID '_Awake_RH_' dataType '_PowerSpectra']);
                close(RH_AwakePower)
                if strcmp(dataType,'CBV_HbT') == false
                    savefig(Hip_AwakePower,[dirpath animalID '_Awake_Hippocampal_' dataType '_PowerSpectra']);
                    close(Hip_AwakePower)
                end
            end
        end
        
        %% Analyze power spectra during periods of NREM sleep
        % pull data from SleepData.mat structure
        if strcmp(dataType,'CBV_HbT') == true
            LH_nremData = SleepData.(modelType).NREM.data.(dataType).LH;
            RH_nremData = SleepData.(modelType).NREM.data.(dataType).RH;
        else
            LH_nremData = SleepData.(modelType).NREM.data.cortical_LH.(dataType);
            RH_nremData = SleepData.(modelType).NREM.data.cortical_RH.(dataType);
            Hip_nremData = SleepData.(modelType).NREM.data.hippocampus.(dataType);
        end
        % detrend - data is already lowpass filtered
        for dd = 1:length(LH_nremData)
            LH_nremData{dd,1} = filtfilt(sos,g,detrend(LH_nremData{dd,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
            RH_nremData{dd,1} = filtfilt(sos,g,detrend(RH_nremData{dd,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
            if strcmp(dataType,'CBV_HbT') == false
                Hip_nremData{dd,1} = filtfilt(sos,g,detrend(Hip_nremData{dd,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
            end
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        LH_nrem = zeros(length(LH_nremData{1,1}),length(LH_nremData));
        RH_nrem = zeros(length(RH_nremData{1,1}),length(RH_nremData));
        if strcmp(dataType,'CBV_HbT') == false
            Hip_nrem = zeros(length(Hip_nremData{1,1}),length(Hip_nremData));
        end
        for ee = 1:length(LH_nremData)
            LH_nrem(:,ee) = LH_nremData{ee,1};
            RH_nrem(:,ee) = RH_nremData{ee,1};
            if strcmp(dataType,'CBV_HbT') == false
                Hip_nrem(:,ee) = Hip_nremData{ee,1};
            end
        end
        % parameters for mtspectrumc - information available in function
        params.tapers = [1,1];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;   % Sampling Rate
        params.fpass = [0,0.5];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the power spectra of the desired signals
        [LH_nrem_S,LH_nrem_f,LH_nrem_sErr] = mtspectrumc_Manuscript2020(LH_nrem,params);
        [RH_nrem_S,RH_nrem_f,RH_nrem_sErr] = mtspectrumc_Manuscript2020(RH_nrem,params);
        if strcmp(dataType,'CBV_HbT') == false
            [Hip_nrem_S,Hip_nrem_f,Hip_nrem_sErr] = mtspectrumc_Manuscript2020(Hip_nrem,params);
        end
        % save data and figures
        AnalysisResults.(animalID).PowerSpectra.NREM.(dataType).adjLH.S = LH_nrem_S;
        AnalysisResults.(animalID).PowerSpectra.NREM.(dataType).adjLH.f = LH_nrem_f;
        AnalysisResults.(animalID).PowerSpectra.NREM.(dataType).adjLH.sErr = LH_nrem_sErr;
        AnalysisResults.(animalID).PowerSpectra.NREM.(dataType).adjRH.S = RH_nrem_S;
        AnalysisResults.(animalID).PowerSpectra.NREM.(dataType).adjRH.f = RH_nrem_f;
        AnalysisResults.(animalID).PowerSpectra.NREM.(dataType).adjRH.sErr = RH_nrem_sErr;
        if strcmp(dataType,'CBV_HbT') == false
            AnalysisResults.(animalID).PowerSpectra.NREM.(dataType).Hip.S = Hip_nrem_S;
            AnalysisResults.(animalID).PowerSpectra.NREM.(dataType).Hip.f = Hip_nrem_f;
            AnalysisResults.(animalID).PowerSpectra.NREM.(dataType).Hip.sErr = Hip_nrem_sErr;
        end
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            LH_nremPower = figure;
            loglog(LH_nrem_f,LH_nrem_S,'k')
            hold on;
            loglog(LH_nrem_f,LH_nrem_sErr,'color',colors_Manuscript2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Power');
            title([animalID  ' adjLH ' dataType ' Power during NREM']);
            set(gca,'Ticklength',[0,0]);
            legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
            set(legend,'FontSize',6);
            xlim([0.1,0.5])
            axis square
            set(gca,'box','off')
            RH_nremPower = figure;
            loglog(RH_nrem_f,RH_nrem_S,'k')
            hold on;
            loglog(RH_nrem_f,RH_nrem_sErr,'color',colors_Manuscript2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Power');
            title([animalID  ' adjRH ' dataType ' Power during NREM']);
            set(gca,'Ticklength',[0,0]);
            legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
            set(legend,'FontSize',6);
            xlim([0.1,0.5])
            axis square
            set(gca,'box','off')
            if strcmp(dataType,'CBV_HbT') == false
                Hip_nremPower = figure;
                loglog(Hip_nrem_f,Hip_nrem_S,'k')
                hold on;
                loglog(Hip_nrem_f,Hip_nrem_sErr,'color',colors_Manuscript2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' Hippocampal ' dataType ' Power during NREM']);
                set(gca,'Ticklength',[0,0]);
                legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
                set(legend,'FontSize',6);
                xlim([0.1,0.5])
                axis square
                set(gca,'box','off')
            end
            savefig(LH_nremPower,[dirpath animalID '_NREM_LH_' dataType '_PowerSpectra']);
            close(LH_nremPower)
            savefig(RH_nremPower,[dirpath animalID '_NREM_RH_' dataType '_PowerSpectra']);
            close(RH_nremPower)
            if strcmp(dataType,'CBV_HbT') == false
                savefig(Hip_nremPower,[dirpath animalID '_NREM_Hippocampal_' dataType '_PowerSpectra']);
                close(Hip_nremPower)
            end
        end
        
        %% Analyze power spectra during periods of REM sleep
        % pull data from SleepData.mat structure
        if strcmp(dataType,'CBV_HbT') == true
            LH_remData = SleepData.(modelType).REM.data.(dataType).LH;
            RH_remData = SleepData.(modelType).REM.data.(dataType).RH;
        else
            LH_remData = SleepData.(modelType).REM.data.cortical_LH.(dataType);
            RH_remData = SleepData.(modelType).REM.data.cortical_RH.(dataType);
            Hip_remData = SleepData.(modelType).REM.data.hippocampus.(dataType);
        end
        % detrend - data is already lowpass filtered
        for ff = 1:length(LH_remData)
            LH_remData{ff,1} = filtfilt(sos,g,detrend(LH_remData{ff,1}(1:(params.minTime.REM*samplingRate)),'constant'));
            RH_remData{ff,1} = filtfilt(sos,g,detrend(RH_remData{ff,1}(1:(params.minTime.REM*samplingRate)),'constant'));
            if strcmp(dataType,'CBV_HbT') == false
                Hip_remData{ff,1} = filtfilt(sos,g,detrend(Hip_remData{ff,1}(1:(params.minTime.REM*samplingRate)),'constant'));
            end
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        LH_rem = zeros(length(LH_remData{1,1}),length(LH_remData));
        RH_rem = zeros(length(RH_remData{1,1}),length(RH_remData));
        if strcmp(dataType,'CBV_HbT') == false
            Hip_rem = zeros(length(Hip_remData{1,1}),length(Hip_remData));
        end
        for gg = 1:length(LH_remData)
            LH_rem(:,gg) = LH_remData{gg,1};
            RH_rem(:,gg) = RH_remData{gg,1};
            if strcmp(dataType,'CBV_HbT') == false
                Hip_rem(:,gg) = Hip_remData{gg,1};
            end
        end
        % parameters for mtspectrumc - information available in function
        params.tapers = [1,1];   % Tapers [n, 2n - 1]
        params.pad = 1;
        params.Fs = samplingRate;   % Sampling Rate
        params.fpass = [0,0.5];   % Pass band [0, nyquist]
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the power spectra of the desired signals
        [LH_rem_S,LH_rem_f,LH_rem_sErr] = mtspectrumc_Manuscript2020(LH_rem,params);
        [RH_rem_S,RH_rem_f,RH_rem_sErr] = mtspectrumc_Manuscript2020(RH_rem,params);
        if strcmp(dataType,'CBV_HbT') == false
            [Hip_rem_S,Hip_rem_f,Hip_rem_sErr] = mtspectrumc_Manuscript2020(Hip_rem,params);
        end
        %save data and figures
        AnalysisResults.(animalID).PowerSpectra.REM.(dataType).adjLH.S = LH_rem_S;
        AnalysisResults.(animalID).PowerSpectra.REM.(dataType).adjLH.f = LH_rem_f;
        AnalysisResults.(animalID).PowerSpectra.REM.(dataType).adjLH.sErr = LH_rem_sErr;
        AnalysisResults.(animalID).PowerSpectra.REM.(dataType).adjRH.S = RH_rem_S;
        AnalysisResults.(animalID).PowerSpectra.REM.(dataType).adjRH.f = RH_rem_f;
        AnalysisResults.(animalID).PowerSpectra.REM.(dataType).adjRH.sErr = RH_rem_sErr;
        if strcmp(dataType,'CBV_HbT') == false
            AnalysisResults.(animalID).PowerSpectra.REM.(dataType).Hip.S = Hip_rem_S;
            AnalysisResults.(animalID).PowerSpectra.REM.(dataType).Hip.f = Hip_rem_f;
            AnalysisResults.(animalID).PowerSpectra.REM.(dataType).Hip.sErr = Hip_rem_sErr;
        end
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            LH_remPower = figure;
            loglog(LH_rem_f,LH_rem_S,'k')
            hold on;
            loglog(LH_rem_f,LH_rem_sErr,'color',colors_Manuscript2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Power');
            title([animalID  ' adjLH ' dataType ' Power during REM']);
            set(gca,'Ticklength',[0,0]);
            legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
            set(legend,'FontSize',6);
            xlim([0.1,0.5])
            axis square
            set(gca,'box','off')
            RH_remPower = figure;
            loglog(RH_rem_f,RH_rem_S,'k')
            hold on;
            loglog(RH_rem_f,RH_rem_sErr,'color',colors_Manuscript2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Power');
            title([animalID  ' adjRH ' dataType ' Power during REM']);
            set(gca,'Ticklength',[0,0]);
            legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
            set(legend,'FontSize',6);
            xlim([0.1,0.5])
            axis square
            set(gca,'box','off')
            if strcmp(dataType,'CBV_HbT') == false
                Hip_remPower = figure;
                loglog(Hip_rem_f,Hip_rem_S,'k')
                hold on;
                loglog(Hip_rem_f,Hip_rem_sErr,'color',colors_Manuscript2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' Hippocampal ' dataType ' Power during REM']);
                set(gca,'Ticklength',[0,0]);
                legend('Coherence','Jackknife Lower','JackknifeUpper','Location','Southeast');
                set(legend,'FontSize',6);
                xlim([0.1,0.5])
                axis square
                set(gca,'box','off')
            end
            savefig(LH_remPower,[dirpath animalID '_REM_LH_' dataType '_PowerSpectra']);
            close(LH_remPower)
            savefig(RH_remPower,[dirpath animalID '_REM_RH_' dataType '_PowerSpectra']);
            close(RH_remPower)
            if strcmp(dataType,'CBV_HbT') == false
                savefig(Hip_remPower,[dirpath animalID '_REM_Hippocampal_' dataType '_PowerSpectra']);
                close(Hip_remPower)
            end
        end
    end
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end