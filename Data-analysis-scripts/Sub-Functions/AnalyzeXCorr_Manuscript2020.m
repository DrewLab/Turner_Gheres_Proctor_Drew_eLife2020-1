function [AnalysisResults] = AnalyzeXCorr_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Analyze the cross-correlation between a hemodynamic signal and a spectrogram during different behaviors.
%________________________________________________________________________________________________________________________

%% function parameters
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120'};
dataTypes = {'adjLH','adjRH'};
modelType = 'SVM';
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
    % find and load AllSpecStruct.mat struct
    allSpecStructFileStruct = dir('*_AllSpecStruct.mat');
    allSpecStructFile = {allSpecStructFileStruct.name}';
    allSpecStructFileID = char(allSpecStructFile);
    load(allSpecStructFileID)
    % determine the animal's ID use the RestData.mat file's name for the current folder
    fileBreaks = strfind(restDataFileID,'_');
    animalID = restDataFileID(1:fileBreaks(1)-1);
    % go through each valid data type for behavior-based cross-correlation analysis
    for z = 1:length(dataTypes)
        dataType = dataTypes{1,z};
        neuralDataType = ['cortical_' dataType(4:end)];
        % pull a few necessary numbers from the RestData.mat struct such as trial duration and sampling rate
        samplingRate = RestData.CBV_HbT.LH.CBVCamSamplingRate;
        trialDuration_sec = RestData.CBV_HbT.LH.trialDuration_sec;   % sec
        sleepBinWidth = 5;   % sec
        oneSecSpecFs = 10;   % sec   5 for fiveSec, 10 for oneSec
        frequencyDiff = 3;   % Hz    6 for fiveSec, 3 for oneSec
        
        %% Cross-correlation analysis for resting data
        % set criteria for rest event filter
        RestCriteria.Fieldname = {'durations'};
        RestCriteria.Comparison = {'gt'};
        RestCriteria.Value = {params.minTime.Rest};
        PuffCriteria.Fieldname = {'puffDistances'};
        PuffCriteria.Comparison = {'gt'};
        PuffCriteria.Value = {5};
        % filter the RestData structure for events that meet the desired criteria
        [restLogical] = FilterEvents_IOS(RestData.CBV_HbT.(dataType),RestCriteria);
        [puffLogical] = FilterEvents_IOS(RestData.CBV_HbT.(dataType),PuffCriteria);
        combRestLogical = logical(restLogical.*puffLogical);
        restFileIDs = RestData.CBV_HbT.(dataType).fileIDs(combRestLogical,:);
        restDurations = RestData.CBV_HbT.(dataType).durations(combRestLogical,:);
        restEventTimes = RestData.CBV_HbT.(dataType).eventTimes(combRestLogical,:);
        restingHbTData = RestData.CBV_HbT.(dataType).data(combRestLogical,:);
        restingMUAData = RestData.(neuralDataType).muaPower.NormData(combRestLogical,:);
        % decimate the file list to only include those files that occur within the desired number of target minutes
        [restFinalRestHbTData,restFinalFileIDs,restFinalDurations,restFinalEventTimes] = DecimateRestData_Manuscript2020(restingHbTData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        [restFinalRestMUAData,~,~,~] = DecimateRestData_Manuscript2020(restingMUAData,restFileIDs,restDurations,restEventTimes,ManualDecisions);
        zz = 1;
        for e = 1:length(restFinalFileIDs)
            restFileID = restFinalFileIDs{e,1};
            % check whether the event occurs in the appropriate time frame
            restStartTime = ceil(restFinalEventTimes(e,1)*10)/10; % *10/10 used to round to first decimal place in a floor/ceil fashion.
            restDuration = floor(restFinalDurations(e,1)*10)/10;
            if restStartTime >= 0.5 && (restStartTime + restDuration) <= (trialDuration_sec - 0.5)
                % remove the number of samples due to rounding up to start and rounding down to end. This is done to keep the HbT/MUA vectores aligned positionally with the upcoming
                % spectral analysis which is at 10 Hz
                leadSamples = round((restStartTime - restFinalEventTimes(e,1))*samplingRate);
                lagSamples = round((restFinalDurations(e,1) - restDuration)*samplingRate);
                % Load in CBV_HbT from rest period
                restHbT = restFinalRestHbTData{e,1};
                restMUA = restFinalRestMUAData{e,1};
                % remove leading/lag samples due to rounding to nearest 0.1 up/0.1 down
                restSnipHbT = restHbT(1 + leadSamples:end - lagSamples);
                restSnipMUA = restMUA(1 + leadSamples:end - lagSamples);
                % low pass filter the epoch below 1 Hz
                [B,A] = butter(3,1/(samplingRate/2),'low');
                restFiltHbT = filtfilt(B,A,restSnipHbT);
                restFiltMUA = filtfilt(B,A,restSnipMUA);
                % only take the first 10 seconds of the epoch. occassionally a sample gets lost from rounding during the
                % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
                if length(restFiltHbT) < params.minTime.Rest*samplingRate
                    restChunkSampleDiff = params.minTime.Rest*samplingRate - length(restFiltHbT);
                    restPadHbT = (ones(1,restChunkSampleDiff))*restFiltHbT(end);
                    restPadMUA = (ones(1,restChunkSampleDiff))*restFiltMUA(end);
                    restShortHbT = horzcat(restFiltHbT,restPadHbT);
                    restShortMUA = horzcat(restFiltMUA,restPadMUA);
                else
                    restShortHbT = restFiltHbT(1:params.minTime.Rest*samplingRate);
                    restShortMUA = restFiltMUA(1:params.minTime.Rest*samplingRate);
                end
                % downsample the 10 second epoch to 5 Hz
                restDsHbT = downsample(restShortHbT,frequencyDiff);
                restDsMUA = downsample(restShortMUA,frequencyDiff);
                % mean subtract the downsampled epoch
                restProcData.HbT{zz,1} = detrend(restDsHbT,'constant');
                restProcData.MUA{zz,1} = detrend(restDsMUA,'constant');
                % extract LFP from spectrograms associated with the whisking indecies
                specDataFileID = [animalID '_' restFileID '_SpecData.mat'];
                clear S_data
                for g = 1:length(AllSpecData.(neuralDataType).fileIDs)
                    if strcmp(AllSpecData.(neuralDataType).fileIDs{g,1},specDataFileID) == true
                        rest_S = AllSpecData.(neuralDataType).oneSec.normS{g,1};
                        rest_T = round(AllSpecData.(neuralDataType).oneSec.T{g,1},1);
                        rest_F = AllSpecData.(neuralDataType).oneSec.F{g,1};
                    end
                end
                restStartTimeIndex = find(rest_T == restStartTime);
                restDurationIndex = find(rest_T == round((restStartTime + restDuration),1));
                restS_Vals = rest_S(:,restStartTimeIndex:restDurationIndex);
                % only take the first min rest time seconds
                shortRestS_Vals = restS_Vals(:,1:params.minTime.Rest*oneSecSpecFs);
                % mean subtract with detrend and lowpass filter each column
                restProcData.S{zz,1} = detrend(shortRestS_Vals','constant')';
                zz = zz + 1;
            end
            % set parameters for cross-correlation analysis
            restHbTvLFPzhold = [];
            restLagTime = 5;   % seconds
            restFrequency = oneSecSpecFs;   % Hz
            restMaxLag = restLagTime*restFrequency;
            restHbTvLFPxcVals = ones(length(rest_F),2*restMaxLag + 1);
            % run cross-correlation analysis - average through time
            for xx = 1:length(restProcData.HbT)
                for b = 1:size(restProcData.S{xx, 1}, 1)
                    restHbTarray = restProcData.HbT{xx,1};
                    restMUAarray = restProcData.MUA{xx,1};
                    restNeuralArray = restProcData.S{xx,1}(b,:);
                    [restHbTvLFPxcVals(b,:),restLFP_lags] = xcorr(restHbTarray,restNeuralArray,restMaxLag,'coeff');
                end
                [restHbTvMUAxcVals(xx,:),restMUA_lags] = xcorr(restHbTarray,restMUAarray,restMaxLag,'coeff'); %#ok<*AGROW>
                restHbTvLFPzhold = cat(3,restHbTvLFPzhold,restHbTvLFPxcVals);
            end
            restMeanHbTvLFPxcVals = mean(restHbTvLFPzhold,3);
            restMeanHbTvMUAxcVals = mean(restHbTvMUAxcVals,1);
            restStdHbTvMUAxcVals = std(restHbTvMUAxcVals,0,1);
        end
        % save data and figures
        AnalysisResults.(animalID).XCorr.Rest.(dataType).LFP_lags = restLFP_lags;
        AnalysisResults.(animalID).XCorr.Rest.(dataType).MUA_lags = restMUA_lags;
        AnalysisResults.(animalID).XCorr.Rest.(dataType).F = rest_F;
        AnalysisResults.(animalID).XCorr.Rest.(dataType).HbTvLFPxcVals = restMeanHbTvLFPxcVals;
        AnalysisResults.(animalID).XCorr.Rest.(dataType).HbTvMUAxcVals = restMeanHbTvMUAxcVals;
        AnalysisResults.(animalID).XCorr.Rest.(dataType).HbTvMUAxcVals_std = restStdHbTvMUAxcVals;
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            titleID = strrep(dataType,'_',' ');
            RestingXCorr = figure;
            sgtitle([animalID ' ' titleID ' resting cross-correlation'])
            subplot(2,1,1)
            plot(restMUA_lags,restMeanHbTvMUAxcVals,'k')
            hold on
            plot(restMUA_lags,restMeanHbTvMUAxcVals + restStdHbTvMUAxcVals,'color',colors_Manuscript2020('battleship grey'))
            plot(restMUA_lags,restMeanHbTvMUAxcVals - restStdHbTvMUAxcVals,'color',colors_Manuscript2020('battleship grey'))
            title('MUA XCorr')
            xticks([-restMaxLag -restMaxLag/2 0 restMaxLag/2 restMaxLag])
            xticklabels({'-5','-2.5','0','2.5','5'})
            xlim([-restLagTime*restFrequency restLagTime*restFrequency])
            xlabel('Lags (sec)')
            ylabel('Cross-correlation')
            axis xy
            axis square
            subplot(2,1,2)
            imagesc(restLFP_lags,rest_F,restMeanHbTvLFPxcVals)
            title('LFP XCorr')
            xticks([-restMaxLag -restMaxLag/2 0 restMaxLag/2 restMaxLag])
            xticklabels({'-5','-2.5','0','2.5','5'})
            xlim([-restLagTime*restFrequency restLagTime*restFrequency])
            xlabel('Lags (sec)')
            ylabel('Freq (Hz)')
            ylim([1,100])
            colorbar
            axis xy
            axis square
            [pathstr,~,~] = fileparts(cd);
            dirpath = [pathstr '/Figures/XCorr/'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(RestingXCorr,[dirpath animalID '_' dataType '_RestingXCorr']);
            close(RestingXCorr)
        end
        
        %% Cross-correlation analysis for NREM sleep data
        NREM_sleepTime = params.minTime.NREM;   % seconds
        NREM_allSleepFileIDs = SleepData.(modelType).NREM.FileIDs;
        NREM_uniqueSleepFileIDs = unique(SleepData.(modelType).NREM.FileIDs);
        k = 1;
        for m = 1:length(NREM_uniqueSleepFileIDs)
            % pull out the bin times (there may be multiple events) in each unique NREM sleep file
            NREM_uniqueSleepFileID = char(NREM_uniqueSleepFileIDs(m));
            n = 1;
            clear NREM_binTimes
            for ee = 1:length(NREM_allSleepFileIDs)
                NREM_sleepFileID = char(NREM_allSleepFileIDs(ee));
                if strcmp(NREM_uniqueSleepFileID,NREM_sleepFileID)
                    NREM_binTimes{n,1} = SleepData.(modelType).NREM.BinTimes{ee,1};
                    n = n + 1;
                end
            end
            % pull out the Spectrogram data that matches the unique NREM sleep file
            NREM_specDataFileID = [animalID '_' NREM_uniqueSleepFileID '_SpecData.mat'];
            load(NREM_specDataFileID)
            NREM_S_Data = SpecData.(neuralDataType).oneSec.normS;
            for q = 1:length(NREM_binTimes)
                NREM_Bins = NREM_binTimes{q,1};
                NREM_startTime = NREM_Bins(1) - sleepBinWidth;
                NREM_endTime = NREM_Bins(end);
                if NREM_startTime > 5 && NREM_endTime < trialDuration_sec
                    NREM_startTimeIndex = find(rest_T == NREM_startTime);
                    NREM_durationIndex = find(rest_T == NREM_endTime);
                    NREM_sleepNeuralVals{k,1} = NREM_S_Data(:,NREM_startTimeIndex:NREM_durationIndex);
                    editIndex{k,1} = {'none'};
                elseif NREM_startTime == 5 && length(NREM_Bins) >= 7
                    NREM_startTime = NREM_Bins(2) - sleepBinWidth;
                    NREM_endTime = NREM_Bins(end);
                    NREM_startTimeIndex = find(rest_T == NREM_startTime);
                    NREM_durationIndex = find(rest_T == NREM_endTime);
                    NREM_sleepNeuralVals{k,1} = NREM_S_Data(:,NREM_startTimeIndex:NREM_durationIndex);
                    editIndex{k,1} = {'leading'};
                elseif NREM_endTime == 900 && length(NREM_Bins) >= 7
                    NREM_startTime = NREM_Bins(1) - sleepBinWidth;
                    NREM_endTime = NREM_Bins(end - 1);
                    NREM_startTimeIndex = find(rest_T == NREM_startTime);
                    NREM_durationIndex = find(rest_T == NREM_endTime);
                    NREM_sleepNeuralVals{k,1} = NREM_S_Data(:,NREM_startTimeIndex:NREM_durationIndex);
                    editIndex{k,1} = {'lagging'};
                else
                    NREM_sleepNeuralVals{k,1} = [];
                    editIndex{k,1} = {'delete'};
                end
                k = k + 1;
            end
        end
        % detrend spectrogram neural values
        for r = 1:length(NREM_sleepNeuralVals)
            NREM_indSleepNeuralVals = NREM_sleepNeuralVals{r,1};
            if isempty(NREM_indSleepNeuralVals) == false
                NREM_indSleepNeuralVals = NREM_indSleepNeuralVals(:,1:NREM_sleepTime*oneSecSpecFs)';
                NREM_dtSleepNeuralVals{r,1} = (detrend(NREM_indSleepNeuralVals,'constant'))';
            else
                NREM_dtSleepNeuralVals{r,1} = [];
            end
        end
        % adjust HbT and MUA events to match the edits made to the length of each spectrogram
        qx = 1;
        for s = 1:length(NREM_dtSleepNeuralVals)
            if isempty(NREM_dtSleepNeuralVals) == false
                NREM_finalSleepNeuralVals{qx,1} = NREM_dtSleepNeuralVals{s,1};
                if strcmp(editIndex{s,1},'none') == true
                    NREM_HbTVals = SleepData.(modelType).NREM.data.CBV_HbT.(dataType(4:end)){s,1}(1:NREM_sleepTime*samplingRate);
                    NREM_MUAVals = SleepData.(modelType).NREM.data.(neuralDataType).muaPower{s,1}(1:NREM_sleepTime*samplingRate);
                    NREM_finalHbTVals{qx,1} = detrend(downsample(NREM_HbTVals,frequencyDiff),'constant');
                    NREM_finalMUAVals{qx,1} = detrend(downsample(NREM_MUAVals,frequencyDiff),'constant');
                elseif strcmp(editIndex{s,1},'leading') == true
                    NREM_HbTVals = SleepData.(modelType).NREM.data.CBV_HbT.(dataType(4:end)){s,1}((samplingRate*sleepBinWidth) + 1:(NREM_sleepTime*samplingRate + samplingRate*sleepBinWidth));
                    NREM_MUAVals = SleepData.(modelType).NREM.data.(neuralDataType).muaPower{s,1}((samplingRate*sleepBinWidth) + 1:(NREM_sleepTime*samplingRate + samplingRate*sleepBinWidth));
                    NREM_finalHbTVals{qx,1} = detrend(downsample(NREM_HbTVals,frequencyDiff),'constant');
                    NREM_finalMUAVals{qx,1} = detrend(downsample(NREM_MUAVals,frequencyDiff),'constant');
                elseif strcmp(editIndex{s,1},'lagging') == true
                    NREM_HbTVals = SleepData.(modelType).NREM.data.CBV_HbT.(dataType(4:end)){s,1}(1:NREM_sleepTime*samplingRate);
                    NREM_MUAVals = SleepData.(modelType).NREM.data.(neuralDataType).muaPower{s,1}(1:NREM_sleepTime*samplingRate);
                    NREM_finalHbTVals{qx,1} = detrend(downsample(NREM_HbTVals,frequencyDiff),'constant');
                    NREM_finalMUAVals{qx,1} = detrend(downsample(NREM_MUAVals,frequencyDiff),'constant');
                elseif strcmp(editIndex{s,1},'delete') == true
                    % remove HbT/MUA from final file
                end
                qx = qx + 1;
            end
        end
        % run cross-correlation analysis - average through time
        NREM_F = SpecData.(neuralDataType).oneSec.F;
        NREM_HbTvLFPzHold = [];
        NREM_lagTime = 15;   % Seconds
        NREM_frequency = oneSecSpecFs;   % Hz
        NREM_maxLag = NREM_lagTime*NREM_frequency;
        NREM_HbTvLFPxcVals = ones(size(NREM_indSleepNeuralVals,2),2*NREM_maxLag + 1);
        for t = 1:length(NREM_finalSleepNeuralVals)
            for u = 1:size(NREM_finalSleepNeuralVals{t,1},1)
                NREM_HbT_array = NREM_finalHbTVals{t,1};
                NREM_MUA_array = NREM_finalMUAVals{t,1};
                NREM_Neural_array = NREM_finalSleepNeuralVals{t,1}(u,:);
                [NREM_HbTvLFPxcVals(u,:),NREM_LFP_lags] = xcorr(NREM_HbT_array,NREM_Neural_array,NREM_maxLag,'coeff');
            end
            [NREM_HbTvMUAxcVals(t,:),NREM_MUA_lags] = xcorr(NREM_HbT_array,NREM_MUA_array,NREM_maxLag,'coeff');
            NREM_HbTvLFPzHold = cat(3,NREM_HbTvLFPzHold,NREM_HbTvLFPxcVals);
        end
        NREM_meanHbTvLFPxcVals = mean(NREM_HbTvLFPzHold,3);
        NREM_meanHbTvMUAxcVals = mean(NREM_HbTvMUAxcVals,1);
        NREM_stdHbTvMUAxcVals = std(NREM_HbTvMUAxcVals,0,1);
        % save data and figures
        AnalysisResults.(animalID).XCorr.NREM.(dataType).LFP_lags = NREM_LFP_lags;
        AnalysisResults.(animalID).XCorr.NREM.(dataType).MUA_lags = NREM_MUA_lags;
        AnalysisResults.(animalID).XCorr.NREM.(dataType).F = NREM_F;
        AnalysisResults.(animalID).XCorr.NREM.(dataType).HbTvLFPxcVals = NREM_meanHbTvLFPxcVals;
        AnalysisResults.(animalID).XCorr.NREM.(dataType).HbTvMUAxcVals = NREM_meanHbTvMUAxcVals;
        AnalysisResults.(animalID).XCorr.NREM.(dataType).HbTvMUAxcVals_std = NREM_stdHbTvMUAxcVals;
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            NREMXCorr = figure;
            subplot(2,1,1)
            sgtitle([animalID ' ' titleID ' NREM cross-correlation'])
            plot(NREM_MUA_lags,NREM_meanHbTvMUAxcVals,'k')
            hold on
            plot(NREM_MUA_lags,NREM_meanHbTvMUAxcVals + NREM_stdHbTvMUAxcVals,'color',colors_Manuscript2020('battleship grey'))
            plot(NREM_MUA_lags,NREM_meanHbTvMUAxcVals - NREM_stdHbTvMUAxcVals,'color',colors_Manuscript2020('battleship grey'))
            title('MUA XCorr')
            xticks([-NREM_maxLag -NREM_maxLag/2 0 NREM_maxLag/2 NREM_maxLag])
            xticklabels({'-15','-7.5','0','7.5','15'})
            xlim([-NREM_lagTime*NREM_frequency NREM_lagTime*NREM_frequency])
            xlabel('Lags (sec)')
            ylabel('Cross-correlation')
            axis xy
            axis square
            subplot(2,1,2)
            imagesc(NREM_LFP_lags,NREM_F,NREM_meanHbTvLFPxcVals)
            title('LFP XCorr')
            xticks([-NREM_maxLag -NREM_maxLag/2 0 NREM_maxLag/2 NREM_maxLag])
            xticklabels({'-15','-7.5','0','7.5','15'})
            xlim([-NREM_lagTime*NREM_frequency NREM_lagTime*NREM_frequency])
            xlabel('Lags (sec)')
            ylabel('Freq (Hz)')
            ylim([1,100])
            colorbar
            axis xy
            axis square
            savefig(NREMXCorr,[dirpath animalID '_' dataType '_NREMXCorr']);
            close(NREMXCorr)
        end
        
        %% Cross-correlation analysis for REM sleep data
        REM_sleepTime = params.minTime.REM;   % seconds
        REM_allSleepFileIDs = SleepData.(modelType).REM.FileIDs;
        REM_uniqueSleepFileIDs = unique(SleepData.(modelType).REM.FileIDs);
        k = 1;
        clear editIndex
        for m = 1:length(REM_uniqueSleepFileIDs)
            % pull out the bin times (there may be multiple events) in each unique NREM sleep file
            REM_uniqueSleepFileID = char(REM_uniqueSleepFileIDs(m));
            n = 1;
            clear REM_binTimes
            for ee = 1:length(REM_allSleepFileIDs)
                REM_sleepFileID = char(REM_allSleepFileIDs(ee));
                if strcmp(REM_uniqueSleepFileID,REM_sleepFileID)
                    REM_binTimes{n,1} = SleepData.(modelType).REM.BinTimes{ee,1};
                    n = n + 1;
                end
            end
            % pull out the Spectrogram data that matches the unique NREM sleep file
            REM_specDataFileID = [animalID '_' REM_uniqueSleepFileID '_SpecData.mat'];
            load(REM_specDataFileID)
            REM_S_Data = SpecData.(neuralDataType).oneSec.normS;
            for q = 1:length(REM_binTimes)
                REM_Bins = REM_binTimes{q,1};
                REM_startTime = REM_Bins(1) - sleepBinWidth;
                REM_endTime = REM_Bins(end);
                if REM_startTime > 5 && REM_endTime < trialDuration_sec
                    REM_startTimeIndex = find(rest_T == REM_startTime);
                    REM_durationIndex = find(rest_T == REM_endTime);
                    REM_sleepNeuralVals{k,1} = REM_S_Data(:,REM_startTimeIndex:REM_durationIndex);
                    editIndex{k,1} = {'none'};
                elseif REM_startTime == 5 && length(REM_Bins) >= 7
                    REM_startTime = REM_Bins(2) - sleepBinWidth;
                    REM_endTime = REM_Bins(end);
                    REM_startTimeIndex = find(rest_T == REM_startTime);
                    REM_durationIndex = find(rest_T == REM_endTime);
                    REM_sleepNeuralVals{k,1} = REM_S_Data(:,REM_startTimeIndex:REM_durationIndex);
                    editIndex{k,1} = {'leading'};
                elseif REM_endTime == 900 && length(REM_Bins) >= 7
                    REM_startTime = REM_Bins(1) - sleepBinWidth;
                    REM_endTime = REM_Bins(end - 1);
                    REM_startTimeIndex = find(rest_T == REM_startTime);
                    REM_durationIndex = find(rest_T == REM_endTime);
                    REM_sleepNeuralVals{k,1} = REM_S_Data(:,REM_startTimeIndex:REM_durationIndex);
                    editIndex{k,1} = {'lagging'};
                else
                    REM_sleepNeuralVals{k,1} = [];
                    editIndex{k,1} = {'delete'};
                end
                k = k + 1;
            end
        end
        % detrend spectrogram neural values
        for r = 1:length(REM_sleepNeuralVals)
            REM_indSleepNeuralVals = REM_sleepNeuralVals{r,1};
            if isempty(REM_indSleepNeuralVals) == false
                REM_indSleepNeuralVals = REM_indSleepNeuralVals(:,1:REM_sleepTime*oneSecSpecFs)';
                REM_dtSleepNeuralVals{r,1} = (detrend(REM_indSleepNeuralVals,'constant'))';
            else
                REM_dtSleepNeuralVals{r,1} = [];
            end
        end
        % adjust HbT and MUA events to match the edits made to the length of each spectrogram
        qx = 1;
        for s = 1:length(REM_dtSleepNeuralVals)
            if isempty(REM_dtSleepNeuralVals) == false
                REM_finalSleepNeuralVals{qx,1} = REM_dtSleepNeuralVals{s,1};
                if strcmp(editIndex{s,1},'none') == true
                    REM_HbTVals = SleepData.(modelType).REM.data.CBV_HbT.(dataType(4:end)){s,1}(1:REM_sleepTime*samplingRate);
                    REM_MUAVals = SleepData.(modelType).REM.data.(neuralDataType).muaPower{s,1}(1:REM_sleepTime*samplingRate);
                    REM_finalHbTVals{qx,1} = detrend(downsample(REM_HbTVals,frequencyDiff),'constant');
                    REM_finalMUAVals{qx,1} = detrend(downsample(REM_MUAVals,frequencyDiff),'constant');
                elseif strcmp(editIndex{s,1},'leading') == true
                    REM_HbTVals = SleepData.(modelType).REM.data.CBV_HbT.(dataType(4:end)){s,1}((samplingRate*sleepBinWidth) + 1:(REM_sleepTime*samplingRate + samplingRate*sleepBinWidth));
                    REM_MUAVals = SleepData.(modelType).REM.data.(neuralDataType).muaPower{s,1}((samplingRate*sleepBinWidth) + 1:(REM_sleepTime*samplingRate + samplingRate*sleepBinWidth));
                    REM_finalHbTVals{qx,1} = detrend(downsample(REM_HbTVals,frequencyDiff),'constant');
                    REM_finalMUAVals{qx,1} = detrend(downsample(REM_MUAVals,frequencyDiff),'constant');
                elseif strcmp(editIndex{s,1},'lagging') == true
                    REM_HbTVals = SleepData.(modelType).REM.data.CBV_HbT.(dataType(4:end)){s,1}(1:REM_sleepTime*samplingRate);
                    REM_MUAVals = SleepData.(modelType).REM.data.(neuralDataType).muaPower{s,1}(1:REM_sleepTime*samplingRate);
                    REM_finalHbTVals{qx,1} = detrend(downsample(REM_HbTVals,frequencyDiff),'constant');
                    REM_finalMUAVals{qx,1} = detrend(downsample(REM_MUAVals,frequencyDiff),'constant');
                elseif strcmp(editIndex{s,1},'delete') == true
                    % remove HbT/MUA from final file
                end
                qx = qx + 1;
            end
        end
        % run cross-correlation analysis - average through time
        REM_F = SpecData.(neuralDataType).oneSec.F;
        REM_HbTvLFPzHold = [];
        REM_lagTime = 15;   % Seconds
        REM_frequency = oneSecSpecFs;   % Hz
        REM_maxLag = REM_lagTime*REM_frequency;
        REM_HbTvLFPxcVals = ones(size(REM_indSleepNeuralVals,2),2*REM_maxLag + 1);
        for t = 1:length(REM_finalSleepNeuralVals)
            for u = 1:size(REM_finalSleepNeuralVals{t,1},1)
                REM_HbT_array = REM_finalHbTVals{t,1};
                REM_MUA_array = REM_finalMUAVals{t,1};
                REM_Neural_array = REM_finalSleepNeuralVals{t,1}(u,:);
                [REM_HbTvLFPxcVals(u,:),REM_LFP_lags] = xcorr(REM_HbT_array,REM_Neural_array,REM_maxLag,'coeff');
            end
            [REM_HbTvMUAxcVals(t,:),REM_MUA_lags] = xcorr(REM_HbT_array,REM_MUA_array,REM_maxLag,'coeff');
            REM_HbTvLFPzHold = cat(3,REM_HbTvLFPzHold,REM_HbTvLFPxcVals);
        end
        REM_meanHbTvLFPxcVals = mean(REM_HbTvLFPzHold,3);
        REM_meanHbTvMUAxcVals = mean(REM_HbTvMUAxcVals,1);
        REM_stdHbTvMUAxcVals = std(REM_HbTvMUAxcVals,0,1);
        % save data and figures
        AnalysisResults.(animalID).XCorr.REM.(dataType).LFP_lags = REM_LFP_lags;
        AnalysisResults.(animalID).XCorr.REM.(dataType).MUA_lags = REM_MUA_lags;
        AnalysisResults.(animalID).XCorr.REM.(dataType).F = REM_F;
        AnalysisResults.(animalID).XCorr.REM.(dataType).HbTvLFPxcVals = REM_meanHbTvLFPxcVals;
        AnalysisResults.(animalID).XCorr.REM.(dataType).HbTvMUAxcVals = REM_meanHbTvMUAxcVals;
        AnalysisResults.(animalID).XCorr.REM.(dataType).HbTvMUAxcVals_std = REM_stdHbTvMUAxcVals;
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            REMXCorr = figure;
            subplot(2,1,1)
            sgtitle([animalID ' ' titleID ' REM cross-correlation'])
            plot(REM_MUA_lags,REM_meanHbTvMUAxcVals,'k')
            hold on
            plot(REM_MUA_lags,REM_meanHbTvMUAxcVals + REM_stdHbTvMUAxcVals,'color',colors_Manuscript2020('battleship grey'))
            plot(REM_MUA_lags,REM_meanHbTvMUAxcVals - REM_stdHbTvMUAxcVals,'color',colors_Manuscript2020('battleship grey'))
            title('MUA XCorr')
            xticks([-REM_maxLag -REM_maxLag/2 0 REM_maxLag/2 REM_maxLag])
            xticklabels({'-15','-7.5','0','7.5','15'})
            xlim([-REM_lagTime*REM_frequency REM_lagTime*REM_frequency])
            xlabel('Lags (sec)')
            ylabel('Cross-correlation')
            axis xy
            axis square
            subplot(2,1,2)
            imagesc(REM_LFP_lags,REM_F,REM_meanHbTvLFPxcVals)
            title('LFP XCorr')
            xticks([-REM_maxLag -REM_maxLag/2 0 REM_maxLag/2 REM_maxLag])
            xticklabels({'-15','-7.5','0','7.5','15'})
            xlim([-REM_lagTime*REM_frequency REM_lagTime*REM_frequency])
            xlabel('Lags (sec)')
            ylabel('Freq (Hz)')
            ylim([1,100])
            colorbar
            axis xy
            axis square
            savefig(REMXCorr,[dirpath animalID '_' dataType '_REMXCorr']);
            close(REMXCorr)
        end
    end
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
