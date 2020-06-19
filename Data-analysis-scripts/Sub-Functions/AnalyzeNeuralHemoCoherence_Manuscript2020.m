function [AnalysisResults] = AnalyzeNeuralHemoCoherence_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Analyze the spectral coherence between bilateral hemodynamic and neural signals.
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
dataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
hemDataTypes = {'adjLH','adjRH'};
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
    % find and load RestingBaselines.mat struct
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % find and load SleepData.mat struct
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
    % find and load Forest_ScoringResults.mat struct
    forestScoringResultsFileID = 'Forest_ScoringResults.mat';
    load(forestScoringResultsFileID,'-mat')
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
    % lowpass filter and detrend each segment
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % go through each valid data type for behavior-based coherence analysis
    for zzz = 1:length(hemDataTypes)
        hemDataType = hemDataTypes{1,zzz};
        for aa = 1:length(dataTypes)
            dataType = dataTypes{1,aa};
            %% Analyze coherence during periods of rest
            % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
            [restLogical] = FilterEvents_IOS_Manuscript2020(RestData.CBV_HbT.(hemDataType),RestCriteria);
            [puffLogical] = FilterEvents_IOS_Manuscript2020(RestData.CBV_HbT.(hemDataType),PuffCriteria);
            combRestLogical = logical(restLogical.*puffLogical);
            restFileIDs = RestData.CBV_HbT.(hemDataType).fileIDs(combRestLogical,:);
            restEventTimes = RestData.CBV_HbT.(hemDataType).eventTimes(combRestLogical,:);
            restDurations = RestData.CBV_HbT.(hemDataType).durations(combRestLogical,:);
            HbT_unstimRestingData = RestData.CBV_HbT.(hemDataType).data(combRestLogical,:);
            Gamma_unstimRestingData = RestData.(['cortical_' hemDataType(4:5)]).(dataType).NormData(combRestLogical,:);
            % decimate the file list to only include those files that occur within the desired number of target minutes
            [HbT_finalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(HbT_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
            [Gamma_finalRestData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(Gamma_unstimRestingData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
            % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
            % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
            clear HbT_ProcRestData Gamma_ProcRestData
            for bb = 1:length(HbT_finalRestData)
                if length(HbT_finalRestData{bb,1}) < params.minTime.Rest*samplingRate
                    restChunkSampleDiff = params.minTime.Rest*samplingRate - length(HbT_finalRestData{bb,1});
                    HbT_restPad = (ones(1,restChunkSampleDiff))*HbT_finalRestData{bb,1}(end);
                    Gamma_restPad = (ones(1,restChunkSampleDiff))*Gamma_finalRestData{bb,1}(end);
                    HbT_ProcRestData{bb,1} = horzcat(HbT_finalRestData{bb,1},HbT_restPad); %#ok<*AGROW>
                    Gamma_ProcRestData{bb,1} = horzcat(Gamma_finalRestData{bb,1},Gamma_restPad);
                    HbT_ProcRestData{bb,1} = filtfilt(sos,g,detrend(HbT_ProcRestData{bb,1},'constant'));
                    Gamma_ProcRestData{bb,1} = filtfilt(sos,g,detrend(Gamma_ProcRestData{bb,1},'constant'));
                else
                    HbT_ProcRestData{bb,1} = filtfilt(sos,g,detrend(HbT_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                    Gamma_ProcRestData{bb,1} = filtfilt(sos,g,detrend(Gamma_finalRestData{bb,1}(1:(params.minTime.Rest*samplingRate)),'constant'));
                end
            end
            % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            HbT_restData = zeros(length(HbT_ProcRestData{1,1}),length(HbT_ProcRestData));
            Gamma_restData = zeros(length(Gamma_ProcRestData{1,1}),length(Gamma_ProcRestData));
            for cc = 1:length(HbT_ProcRestData)
                HbT_restData(:,cc) = HbT_ProcRestData{cc,1};
                Gamma_restData(:,cc) = Gamma_ProcRestData{cc,1};
            end
            % parameters for coherencyc - information available in function
            params.tapers = [1,1];   % Tapers [n, 2n - 1]
            params.pad = 1;
            params.Fs = samplingRate;   % Sampling Rate
            params.fpass = [0,0.5];   % Pass band [0, nyquist]
            params.trialave = 1;
            params.err = [2,0.05];
            % calculate the coherence between desired signals
            [C_RestData,~,~,~,~,f_RestData,confC_RestData,~,cErr_RestData] = coherencyc_Manuscript2020(HbT_restData,Gamma_restData,params);
            % save data and figures
            AnalysisResults.(animalID).NeuralHemoCoherence.Rest.(dataType).(hemDataType).C = C_RestData;
            AnalysisResults.(animalID).NeuralHemoCoherence.Rest.(dataType).(hemDataType).f = f_RestData;
            AnalysisResults.(animalID).NeuralHemoCoherence.Rest.(dataType).(hemDataType).confC = confC_RestData;
            AnalysisResults.(animalID).NeuralHemoCoherence.Rest.(dataType).(hemDataType).cErr = cErr_RestData;
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
                xlim([0,0.5])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Coherence/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(restCoherence,[dirpath animalID '_Rest_' dataType '_' hemDataType '_Coherence']);
                close(restCoherence)
            end
            
            %% Analyze coherence during awake periods with no sleep scores
            zz = 1;
            clear HbT_AwakeData Gamma_AwakeData HbT_ProcAwakeData Gamma_ProcAwakeData
            HbT_AwakeData = [];
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
                        HbT_AwakeData{zz,1} = ProcData.data.CBV_HbT.(hemDataType);
                        Gamma_AwakeData{zz,1} = (ProcData.data.(['cortical_' hemDataType(4:5)]).(dataType) - RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay))./RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay);
                        zz = zz + 1;
                    end
                end
            end
            % process
            if isempty(HbT_AwakeData) == false
                for bb = 1:length(HbT_AwakeData)
                    HbT_ProcAwakeData{bb,1} = filtfilt(sos,g,detrend(HbT_AwakeData{bb,1},'constant'));
                    Gamma_ProcAwakeData{bb,1} = filtfilt(sos,g,detrend(Gamma_AwakeData{bb,1},'constant'));
                end
                % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
                HbT_awakeData = zeros(length(HbT_ProcAwakeData{1,1}),length(HbT_ProcAwakeData));
                Gamma_awakeData = zeros(length(Gamma_ProcAwakeData{1,1}),length(Gamma_ProcAwakeData));
                for cc = 1:length(HbT_ProcAwakeData)
                    HbT_awakeData(:,cc) = HbT_ProcAwakeData{cc,1};
                    Gamma_awakeData(:,cc) = Gamma_ProcAwakeData{cc,1};
                end
                % parameters for coherencyc - information available in function
                params.tapers = [10,19];   % Tapers [n, 2n - 1]
                params.pad = 1;
                params.Fs = samplingRate;   % Sampling Rate
                params.fpass = [0,0.5];   % Pass band [0, nyquist]
                params.trialave = 1;
                params.err = [2,0.05];
                % calculate the coherence between desired signals
                [C_AwakeData,~,~,~,~,f_AwakeData,confC_AwakeData,~,cErr_AwakeData] = coherencyc_Manuscript2020(HbT_awakeData,Gamma_awakeData,params);
                % save data and figures
                AnalysisResults.(animalID).NeuralHemoCoherence.Awake.(dataType).(hemDataType).C = C_AwakeData;
                AnalysisResults.(animalID).NeuralHemoCoherence.Awake.(dataType).(hemDataType).f = f_AwakeData;
                AnalysisResults.(animalID).NeuralHemoCoherence.Awake.(dataType).(hemDataType).confC = confC_AwakeData;
                AnalysisResults.(animalID).NeuralHemoCoherence.Awake.(dataType).(hemDataType).cErr = cErr_AwakeData;
                % save figures if desired
                if strcmp(saveFigs,'y') == true
                    awakeCoherence = figure;
                    plot(f_AwakeData,C_AwakeData,'k')
                    hold on;
                    plot(f_AwakeData,cErr_AwakeData,'color',colors_Manuscript2020('battleship grey'))
                    xlabel('Freq (Hz)');
                    ylabel('Coherence');
                    title([animalID  ' ' dataType ' coherence for awake data']);
                    set(gca,'Ticklength',[0,0]);
                    legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
                    set(legend,'FontSize',6);
                    ylim([0,1])
                    xlim([0,0.5])
                    axis square
                    set(gca,'box','off')
                    [pathstr,~,~] = fileparts(cd);
                    dirpath = [pathstr '/Figures/Coherence/'];
                    if ~exist(dirpath,'dir')
                        mkdir(dirpath);
                    end
                    savefig(awakeCoherence,[dirpath animalID '_Awake_' dataType '_' hemDataType '_Coherence']);
                    close(awakeCoherence)
                end
            end
            
            %% Analyze coherence during sleep periods with no sleep scores
            zz = 1;
            clear HbT_SleepData Gamma_SleepData HbT_ProcSleepData Gamma_ProcSleepData
            HbT_SleepData = [];
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
                if sum(strcmp(scoringLabels,'Not Sleep')) < 12   % 6 bins (180 total) or 30 seconds of sleep
                    load(procDataFileID)
                    puffs = ProcData.data.solenoids.LPadSol;
                    if isempty(puffs) == true
                        HbT_SleepData{zz,1} = ProcData.data.CBV_HbT.(hemDataType);
                        Gamma_SleepData{zz,1} = (ProcData.data.(['cortical_' hemDataType(4:5)]).(dataType) - RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay))./RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay);
                        zz = zz + 1;
                    end
                end
            end
            % process
            if isempty(HbT_SleepData) == false
                for bb = 1:length(HbT_SleepData)
                    HbT_ProcSleepData{bb,1} = filtfilt(sos,g,detrend(HbT_SleepData{bb,1},'constant'));
                    Gamma_ProcSleepData{bb,1} = filtfilt(sos,g,detrend(Gamma_SleepData{bb,1},'constant'));
                end
                % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
                HbT_sleepData = zeros(length(HbT_ProcSleepData{1,1}),length(HbT_ProcSleepData));
                Gamma_sleepData = zeros(length(Gamma_ProcSleepData{1,1}),length(Gamma_ProcSleepData));
                for cc = 1:length(HbT_ProcSleepData)
                    HbT_sleepData(:,cc) = HbT_ProcSleepData{cc,1};
                    Gamma_sleepData(:,cc) = Gamma_ProcSleepData{cc,1};
                end
                % parameters for coherencyc - information available in function
                params.tapers = [10,19];   % Tapers [n, 2n - 1]
                params.pad = 1;
                params.Fs = samplingRate;   % Sampling Rate
                params.fpass = [0,0.5];   % Pass band [0, nyquist]
                params.trialave = 1;
                params.err = [2,0.05];
                % calculate the coherence between desired signals
                [C_SleepData,~,~,~,~,f_SleepData,confC_SleepData,~,cErr_SleepData] = coherencyc_Manuscript2020(HbT_sleepData,Gamma_sleepData,params);
                % save data and figures
                AnalysisResults.(animalID).NeuralHemoCoherence.Sleep.(dataType).(hemDataType).C = C_SleepData;
                AnalysisResults.(animalID).NeuralHemoCoherence.Sleep.(dataType).(hemDataType).f = f_SleepData;
                AnalysisResults.(animalID).NeuralHemoCoherence.Sleep.(dataType).(hemDataType).confC = confC_SleepData;
                AnalysisResults.(animalID).NeuralHemoCoherence.Sleep.(dataType).(hemDataType).cErr = cErr_SleepData;
                % save figures if desired
                if strcmp(saveFigs,'y') == true
                    sleepCoherence = figure;
                    plot(f_SleepData,C_SleepData,'k')
                    hold on;
                    plot(f_SleepData,cErr_SleepData,'color',colors_Manuscript2020('battleship grey'))
                    xlabel('Freq (Hz)');
                    ylabel('Coherence');
                    title([animalID  ' ' dataType ' coherence for sleep data']);
                    set(gca,'Ticklength',[0,0]);
                    legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
                    set(legend,'FontSize',6);
                    ylim([0,1])
                    xlim([0,0.5])
                    axis square
                    set(gca,'box','off')
                    [pathstr,~,~] = fileparts(cd);
                    dirpath = [pathstr '/Figures/Coherence/'];
                    if ~exist(dirpath,'dir')
                        mkdir(dirpath);
                    end
                    savefig(sleepCoherence,[dirpath animalID '_Sleep_' dataType '_' hemDataType '_Coherence']);
                    close(sleepCoherence)
                end
            end
            
            %% Analyze coherence during awake periods with no sleep scores
            zz = 1;
            clear HbT_AllUnstimData Gamma_AllUnstimData HbT_ProcAllUnstimData Gamma_ProcAllUnstimData
            HbT_AllUnstimData = [];
            for bb = 1:size(procDataFileIDs,1)
                procDataFileID = procDataFileIDs(bb,:);
                [~,allUnstimDataFileDate,~] = GetFileInfo_IOS_Manuscript2020(procDataFileID);
                strDay = ConvertDate_IOS_Manuscript2020(allUnstimDataFileDate);
                load(procDataFileID)
                puffs = ProcData.data.solenoids.LPadSol;
                if isempty(puffs) == true
                    HbT_AllUnstimData{zz,1} = ProcData.data.CBV_HbT.(hemDataType);
                    Gamma_AllUnstimData{zz,1} = (ProcData.data.(['cortical_' hemDataType(4:5)]).(dataType) - RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay))./RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay);
                    zz = zz + 1;
                end
            end
            % process
            if isempty(HbT_AllUnstimData) == false
                for bb = 1:length(HbT_AllUnstimData)
                    HbT_ProcAllUnstimData{bb,1} = filtfilt(sos,g,detrend(HbT_AllUnstimData{bb,1},'constant'));
                    Gamma_ProcAllUnstimData{bb,1} = filtfilt(sos,g,detrend(Gamma_AllUnstimData{bb,1},'constant'));
                end
                % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
                HbT_allUnstimData = zeros(length(HbT_ProcAllUnstimData{1,1}),length(HbT_ProcAllUnstimData));
                Gamma_allUnstimData = zeros(length(Gamma_ProcAllUnstimData{1,1}),length(Gamma_ProcAllUnstimData));
                for cc = 1:length(HbT_ProcAllUnstimData)
                    HbT_allUnstimData(:,cc) = HbT_ProcAllUnstimData{cc,1};
                    Gamma_allUnstimData(:,cc) = Gamma_ProcAllUnstimData{cc,1};
                end
                % parameters for coherencyc - information available in function
                params.tapers = [10,19];   % Tapers [n, 2n - 1]
                params.pad = 1;
                params.Fs = samplingRate;   % Sampling Rate
                params.fpass = [0,0.5];   % Pass band [0, nyquist]
                params.trialave = 1;
                params.err = [2,0.05];
                % calculate the coherence between desired signals
                [C_AllUnstimData,~,~,~,~,f_AllUnstimData,confC_AllUnstimData,~,cErr_AllUnstimData] = coherencyc_Manuscript2020(HbT_allUnstimData,Gamma_allUnstimData,params);
                % save data and figures
                AnalysisResults.(animalID).NeuralHemoCoherence.All.(dataType).(hemDataType).C = C_AllUnstimData;
                AnalysisResults.(animalID).NeuralHemoCoherence.All.(dataType).(hemDataType).f = f_AllUnstimData;
                AnalysisResults.(animalID).NeuralHemoCoherence.All.(dataType).(hemDataType).confC = confC_AllUnstimData;
                AnalysisResults.(animalID).NeuralHemoCoherence.All.(dataType).(hemDataType).cErr = cErr_AllUnstimData;
                % save figures if desired
                if strcmp(saveFigs,'y') == true
                    allUnstimCoherence = figure;
                    plot(f_AllUnstimData,C_AllUnstimData,'k')
                    hold on;
                    plot(f_AllUnstimData,cErr_AllUnstimData,'color',colors_Manuscript2020('battleship grey'))
                    xlabel('Freq (Hz)');
                    ylabel('Coherence');
                    title([animalID  ' ' dataType ' coherence for all unstimulated data']);
                    set(gca,'Ticklength',[0,0]);
                    legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
                    set(legend,'FontSize',6);
                    ylim([0,1])
                    xlim([0,0.5])
                    axis square
                    set(gca,'box','off')
                    [pathstr,~,~] = fileparts(cd);
                    dirpath = [pathstr '/Figures/Coherence/'];
                    if ~exist(dirpath,'dir')
                        mkdir(dirpath);
                    end
                    savefig(allUnstimCoherence,[dirpath animalID '_AllUnstim_' dataType '_' hemDataType '_Coherence']);
                    close(allUnstimCoherence)
                end
            end
            
            %% Analyze coherence during periods of NREM sleep
            % pull data from SleepData.mat structure
            
            HbT_nremData = SleepData.(modelType).NREM.data.CBV_HbT.(hemDataType(4:5));
            Gamma_nremData = SleepData.(modelType).NREM.data.(['cortical_' hemDataType(4:5)]).(dataType);
            % detrend - data is already lowpass filtered
            for ee = 1:length(HbT_nremData)
                HbT_nremData{ee,1} = filtfilt(sos,g,detrend(HbT_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
                Gamma_nremData{ee,1} = filtfilt(sos,g,detrend(Gamma_nremData{ee,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
            end
            % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            HbT_nrem = zeros(length(HbT_nremData{1,1}),length(HbT_nremData));
            Gamma_nrem = zeros(length(Gamma_nremData{1,1}),length(Gamma_nremData));
            for ff = 1:length(HbT_nremData)
                HbT_nrem(:,ff) = HbT_nremData{ff,1};
                Gamma_nrem(:,ff) = Gamma_nremData{ff,1};
            end
            % parameters for coherencyc - information available in function
            params.tapers = [3,5];   % Tapers [n, 2n - 1]
            params.pad = 1;
            params.Fs = samplingRate;   % Sampling Rate
            params.fpass = [0,0.5];   % Pass band [0, nyquist]
            params.trialave = 1;
            params.err = [2,0.05];
            % calculate the coherence between desired signals
            [C_nrem,~,~,~,~,f_nrem,confC_nrem,~,cErr_nrem] = coherencyc_Manuscript2020(HbT_nrem,Gamma_nrem,params);
            % save data and figures
            AnalysisResults.(animalID).NeuralHemoCoherence.NREM.(dataType).(hemDataType).C = C_nrem;
            AnalysisResults.(animalID).NeuralHemoCoherence.NREM.(dataType).(hemDataType).f = f_nrem;
            AnalysisResults.(animalID).NeuralHemoCoherence.NREM.(dataType).(hemDataType).confC = confC_nrem;
            AnalysisResults.(animalID).NeuralHemoCoherence.NREM.(dataType).(hemDataType).cErr = cErr_nrem;
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
                ylim([0.1,0.5])
                xlim([0,0.5])
                axis square
                set(gca,'box','off')
                savefig(nremCoherence,[dirpath animalID '_' modelType '_NREM_' dataType '_' hemDataType '_Coherence']);
                close(nremCoherence)
            end
            
            %% Analyze coherence during periods of REM sleep
            % pull data from SleepData.mat structure
            HbT_remData = SleepData.(modelType).REM.data.CBV_HbT.(hemDataType(4:5));
            Gamma_remData = SleepData.(modelType).REM.data.(['cortical_' hemDataType(4:5)]).(dataType);
            % detrend - data is already lowpass filtered
            for gg = 1:length(HbT_remData)
                HbT_remData{gg,1} = filtfilt(sos,g,detrend(HbT_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant'));
                Gamma_remData{gg,1} = filtfilt(sos,g,detrend(Gamma_remData{gg,1}(1:(params.minTime.REM*samplingRate)),'constant'));
            end
            % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
            HbT_rem = zeros(length(HbT_remData{1,1}),length(HbT_remData));
            Gamma_rem = zeros(length(Gamma_remData{1,1}),length(Gamma_remData));
            for hh = 1:length(HbT_remData)
                HbT_rem(:,hh) = HbT_remData{hh,1};
                Gamma_rem(:,hh) = Gamma_remData{hh,1};
            end
            % parameters for coherencyc - information available in function
            params.tapers = [5,9];   % Tapers [n, 2n - 1]
            params.pad = 1;
            params.Fs = samplingRate;   % Sampling Rate
            params.fpass = [0,0.5];   % Pass band [0, nyquist]
            params.trialave = 1;
            params.err = [2,0.05];
            % calculate the coherence between desired signals
            [C_rem,~,~,~,~,f_rem,confC_rem,~,cErr_rem] = coherencyc_Manuscript2020(HbT_rem,Gamma_rem,params);
            % save data and figures
            AnalysisResults.(animalID).NeuralHemoCoherence.REM.(dataType).(hemDataType).C = C_rem;
            AnalysisResults.(animalID).NeuralHemoCoherence.REM.(dataType).(hemDataType).f = f_rem;
            AnalysisResults.(animalID).NeuralHemoCoherence.REM.(dataType).(hemDataType).confC = confC_rem;
            AnalysisResults.(animalID).NeuralHemoCoherence.REM.(dataType).(hemDataType).cErr = cErr_rem;
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
                ylim([0.1,0.5])
                xlim([0,0.5])
                axis square
                set(gca,'box','off')
                savefig(remCoherence,[dirpath animalID '_' modelType '_REM_' dataType '_' hemDataType '_Coherence']);
                close(remCoherence)
            end
        end
    end
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
