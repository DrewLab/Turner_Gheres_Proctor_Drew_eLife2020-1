function [AnalysisResults] = AnalyzeEvokedResponses_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Use epochs from the EventData.mat struct to determine the average hemodynamic and neural responses to
%            both volitional whisking and whisker stimuli
%________________________________________________________________________________________________________________________

%% function parameters
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111'};
dataTypes = {'adjLH','adjRH'};

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
    % find and load RestingBaselines.mat struct
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % find and load AllSpecStruct.mat struct
    allSpecStructFileStruct = dir('*_AllSpecStruct.mat');
    allSpecStructFile = {allSpecStructFileStruct.name}';
    allSpecStructFileID = char(allSpecStructFile);
    load(allSpecStructFileID)
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
    PuffCriteria.Fieldname = {'puffDistances'};
    PuffCriteria.Comparison = {'gt'};
    PuffCriteria.Value = {5};
    whiskCriteriaNames = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
    % filter the EventData.mat structure for whisking events that meet the desired criteria
    for a = 1:length(dataTypes)
        dataType = dataTypes{1,a};
        neuralDataType = ['cortical_' dataType(4:end)];
        % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
        samplingRate = EventData.CBV_HbT.(dataType).whisk.samplingRate;
        trialDuration_sec = EventData.CBV_HbT.(dataType).whisk.trialDuration_sec;
        timeVector = (0:(EventData.CBV_HbT.(dataType).whisk.epoch.duration*samplingRate))/samplingRate - EventData.CBV_HbT.(dataType).whisk.epoch.offset;
        offset = EventData.CBV_HbT.(dataType).whisk.epoch.offset;
        for b = 1:length(whiskCriteriaNames)
            whiskCriteriaName = whiskCriteriaNames{1,b};
            if strcmp(whiskCriteriaName,'ShortWhisks') == true
                WhiskCriteria = whiskCriteriaA;
            elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
                WhiskCriteria = whiskCriteriaB;
            elseif strcmp(whiskCriteriaName,'LongWhisks') == true
                WhiskCriteria = whiskCriteriaC;
            end
            [whiskLogical] = FilterEvents_IOS(EventData.CBV_HbT.(dataType).whisk,WhiskCriteria);
            [puffLogical] = FilterEvents_IOS(EventData.CBV_HbT.(dataType).whisk,PuffCriteria);
            combWhiskLogical = logical(whiskLogical.*puffLogical);
            [allWhiskHbTData] = EventData.CBV_HbT.(dataType).whisk.data(combWhiskLogical,:);
            [allWhiskCBVData] = EventData.CBV.(dataType).whisk.NormData(combWhiskLogical,:);
            [allWhiskCorticalMUAData] = EventData.(neuralDataType).muaPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskHippocampalMUAData] = EventData.hippocampus.muaPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskFileIDs] = EventData.CBV_HbT.(dataType).whisk.fileIDs(combWhiskLogical,:);
            [allWhiskEventTimes] = EventData.CBV_HbT.(dataType).whisk.eventTime(combWhiskLogical,:);
            allWhiskDurations = EventData.CBV_HbT.(dataType).whisk.duration(combWhiskLogical,:);
            % decimate the file list to only include those files that occur within the desired number of target minutes
            [finalWhiskHbTData,finalWhiskFileIDs,~,finalWhiskFileEventTimes] = DecimateRestData_Manuscript2020(allWhiskHbTData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCBVData,~,~,~] = DecimateRestData_Manuscript2020(allWhiskCBVData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCorticalMUAData,~,~,~] = DecimateRestData_Manuscript2020(allWhiskCorticalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskHippocampalMUAData,~,~,~] = DecimateRestData_Manuscript2020(allWhiskHippocampalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            % lowpass filter each whisking event and mean-subtract by the first 2 seconds
            clear procWhiskHbTData procWhiskCBVData procWhiskCorticalMUAData procWhiskHippocampalMUAData finalWhiskStartTimes finalWhiskEndTimes finalWFileIDs
            qx = 1;
            for f = 1:size(finalWhiskHbTData,1)
                whiskStartTime = round(finalWhiskFileEventTimes(f,1),1) - 2;
                whiskEndTime = whiskStartTime + 10;
                finalWhiskFileID = finalWhiskFileIDs{f,1};
                if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                    whiskHbTarray = finalWhiskHbTData(f,:);
                    whiskCBVarray = finalWhiskCBVData(f,:);
                    whiskCorticalMUAarray = finalWhiskCorticalMUAData(f,:);
                    whiskHippocampalMUAarray = finalWhiskHippocampalMUAData(f,:);
                    filtWhiskHbTarray = sgolayfilt(whiskHbTarray,3,17);
                    filtWhiskCBVarray = sgolayfilt(whiskCBVarray,3,17);
                    filtWhiskCorticalMUAarray = sgolayfilt(whiskCorticalMUAarray,3,17);
                    filtWhiskHippocampalMUAarray = sgolayfilt(whiskHippocampalMUAarray,3,17);
                    procWhiskHbTData(qx,:) = filtWhiskHbTarray - mean(filtWhiskHbTarray(1:(offset*samplingRate))); %#ok<*AGROW>
                    procWhiskCBVData(qx,:) = filtWhiskCBVarray - mean(filtWhiskCBVarray(1:(offset*samplingRate)));
                    procWhiskCorticalMUAData(qx,:) = filtWhiskCorticalMUAarray - mean(filtWhiskCorticalMUAarray(1:(offset*samplingRate)));
                    procWhiskHippocampalMUAData(qx,:) = filtWhiskHippocampalMUAarray - mean(filtWhiskHippocampalMUAarray(1:(offset*samplingRate)));
                    finalWhiskStartTimes(qx,1) = whiskStartTime;
                    finalWhiskEndTimes(qx,1) = whiskEndTime;
                    finalWFileIDs{qx,1} = finalWhiskFileID;
                    qx = qx + 1;
                end
            end
            meanWhiskHbTData = mean(procWhiskHbTData,1);
            stdWhiskHbTData = std(procWhiskHbTData,0,1);
            meanWhiskCBVData = mean(procWhiskCBVData,1)*100;
            stdWhiskCBVData = std(procWhiskCBVData,0,1)*100;
            meanWhiskCorticalMUAData = mean(procWhiskCorticalMUAData,1)*100;
            stdWhiskCorticalMUAData = std(procWhiskCorticalMUAData,0,1)*100;
            meanWhiskHippocampalMUAData = mean(procWhiskHippocampalMUAData,1)*100;
            stdWhiskHippocampalMUAData = std(procWhiskHippocampalMUAData,0,1)*100;
            % extract LFP from spectrograms associated with the whisking indecies
            whiskCorticalZhold = [];
            whiskHippocampalZhold = [];
            for g = 1:length(finalWFileIDs)
                % load normalized one-second bin data from each file
                whiskFileID = finalWFileIDs{g,1};
                whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecData.mat'];
                whiskSpecField = neuralDataType;
                for h = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                    if strcmp(AllSpecData.(whiskSpecField).fileIDs{h,1},whiskSpecDataFileID) == true
                        whiskCorticalS_Data = AllSpecData.(whiskSpecField).oneSec.normS{h,1};
                        whiskHippocampalS_Data = AllSpecData.hippocampus.oneSec.normS{h,1};
                        F = AllSpecData.(whiskSpecField).oneSec.F{h,1};
                        T = round(AllSpecData.(whiskSpecField).oneSec.T{h,1},1);
                    end
                end
                whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(g,1),1));
                whiskDurationIndex = find(T == round(finalWhiskEndTimes(g,1),1));
                whiskCorticalS_Vals = whiskCorticalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
                whiskHippocampalS_Vals = whiskHippocampalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
                % mean subtract each row with detrend
                transpWhiskCorticalS_Vals = whiskCorticalS_Vals';   % Transpose since detrend goes down columns
                transpWhiskHippocampalS_Vals = whiskHippocampalS_Vals';
                dTWhiskCorticalS_Vals = detrend(transpWhiskCorticalS_Vals,'constant');
                dTWhiskHippocampalS_Vals = detrend(transpWhiskHippocampalS_Vals,'constant');
                whiskCorticalZhold = cat(3,whiskCorticalZhold,dTWhiskCorticalS_Vals');   % transpose back to original orientation
                whiskHippocampalZhold = cat(3,whiskHippocampalZhold,dTWhiskHippocampalS_Vals');
            end
            % figure time/frequency axis and average each S data matrix through time
            meanWhiskCorticalS = mean(whiskCorticalZhold,3);
            meanWhiskHippocampalS = mean(whiskHippocampalZhold,3);
            % Save figures if desired
            if strcmp(saveFigs,'y') == true
                whiskEvoked = figure;
                sgtitle([animalID ' ' dataType ' whisking-evoked averages - ' whiskCriteriaName])
                subplot(2,3,1);
                plot(timeVector,meanWhiskCorticalMUAData,'k')
                hold on
                plot(timeVector,meanWhiskCorticalMUAData + stdWhiskCorticalMUAData,'color',colors_Manuscript2020('battleship grey'))
                plot(timeVector,meanWhiskCorticalMUAData - stdWhiskCorticalMUAData,'color',colors_Manuscript2020('battleship grey'))
                title('Cortical MUA')
                xlabel('Time (sec)')
                ylabel('Fold-change (Norm Power)')
                axis tight
                axis square
                subplot(2,3,4);
                plot(timeVector,meanWhiskHippocampalMUAData,'k')
                hold on
                plot(timeVector,meanWhiskHippocampalMUAData + stdWhiskHippocampalMUAData,'color',colors_Manuscript2020('battleship grey'))
                plot(timeVector,meanWhiskHippocampalMUAData - stdWhiskHippocampalMUAData,'color',colors_Manuscript2020('battleship grey'))
                title([animalID ' ' dataType ' ' whiskCriteriaName ' whisking-evoked averages'])
                title('Hippocampal MUA')
                xlabel('Time (sec)')
                ylabel('Fold-change (Norm Power)')
                axis tight
                axis square
                subplot(2,3,2);
                imagesc(T,F,meanWhiskCorticalS)
                title('Cortical LFP')
                xlabel('Time (sec)')
                ylabel('Freq (Hz)')
                ylim([1 100])
                caxis([-0.5 1])
                set(gca,'Ticklength',[0,0])
                axis xy
                axis square
                subplot(2,3,5);
                imagesc(T,F,meanWhiskHippocampalS)
                title('Hippocampal LFP')
                xlabel('Time (sec)')
                ylabel('Freq (Hz)')
                ylim([1,100])
                caxis([-0.5,1])
                set(gca,'Ticklength',[0 0])
                axis xy
                axis square
                subplot(2,3,[3,6]);
                plot(timeVector,meanWhiskHbTData,'k')
                hold on
                plot(timeVector,meanWhiskHbTData + stdWhiskHbTData,'color',colors_Manuscript2020('battleship grey'))
                plot(timeVector,meanWhiskHbTData - stdWhiskHbTData,'color',colors_Manuscript2020('battleship grey'))
                title('Hemodynamic response')
                xlabel('Time (sec)')
                ylabel('\DeltaHbT (\muM)')
                axis tight
                axis square
                % save figure
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Evoked Responses/'];
                if ~exist(dirpath, 'dir')
                    mkdir(dirpath);
                end
                savefig(whiskEvoked,[dirpath animalID '_' dataType '_' whiskCriteriaName '_WhiskEvokedAverages']);
                close(whiskEvoked)
            end
            % save results
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV_HbT.HbT = meanWhiskHbTData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV_HbT.HbTStD = stdWhiskHbTData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.CBV = meanWhiskCBVData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.CBVStD = stdWhiskCBVData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.hippocampalData = meanWhiskHippocampalMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.hippocampalStD = stdWhiskHippocampalMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).timeVector = timeVector;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.corticalS = meanWhiskCorticalS;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.hippocampalS = meanWhiskHippocampalS;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.T = T;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.F = F;
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
            allStimFilter = FilterEvents_IOS(EventData.CBV_HbT.(dataType).stim,stimCriteria);
            [allStimHbTData] = EventData.CBV_HbT.(dataType).stim.data(allStimFilter,:);
            [allStimCBVData] = EventData.CBV.(dataType).stim.NormData(allStimFilter,:);
            [allStimCortMUAData] = EventData.(neuralDataType).muaPower.stim.NormData(allStimFilter,:);
            [allStimHipMUAData] = EventData.hippocampus.muaPower.stim.NormData(allStimFilter,:);
            [allStimFileIDs] = EventData.CBV_HbT.(dataType).stim.fileIDs(allStimFilter,:);
            [allStimEventTimes] = EventData.CBV_HbT.(dataType).stim.eventTime(allStimFilter,:);
            allStimDurations = zeros(length(allStimEventTimes),1);
            % decimate the file list to only include those files that occur within the desired number of target minutes
            [finalStimHbTData,finalStimFileIDs,~,finalStimFileEventTimes] = DecimateRestData_Manuscript2020(allStimHbTData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCBVData,~,~,~] = DecimateRestData_Manuscript2020(allStimCBVData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortMUAData,~,~,~] = DecimateRestData_Manuscript2020(allStimCortMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimHipMUAData,~,~,~] = DecimateRestData_Manuscript2020(allStimHipMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);       
            % lowpass filter each whisking event and mean-subtract by the first 2 seconds
            clear procStimHbTData procStimCBVData procStimCortMUAData procStimHipMUAData
            qx = 1;
            for p = 1:size(finalStimHbTData,1)
                stimStartTime = round(finalStimFileEventTimes(p,1),1) - 2;
                stimEndTime = stimStartTime + 10;
                finalStimFileID = finalStimFileIDs{p,1};
                if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)  
                    stimHbTarray = finalStimHbTData(p,:);
                    stimCBVarray = finalStimCBVData(p,:);
                    stimCortMUAarray = finalStimCortMUAData(p,:);
                    stimHipMUAarray = finalStimHipMUAData(p,:);
                    filtStimHbTarray = sgolayfilt(stimHbTarray,3,17);
                    filtStimCBVarray = sgolayfilt(stimCBVarray,3,17)*100;
                    filtStimCortMUAarray = sgolayfilt(stimCortMUAarray,3,17);
                    filtStimHipMUAarray = sgolayfilt(stimHipMUAarray,3,17);
                    procStimHbTData(p,:) = filtStimHbTarray - mean(filtStimHbTarray(1:(offset*samplingRate)));
                    procStimCBVData(p,:) = filtStimCBVarray - mean(filtStimCBVarray(1:(offset*samplingRate)));
                    procStimCortMUAData(p,:) = filtStimCortMUAarray - mean(filtStimCortMUAarray(1:(offset*samplingRate)));
                    procStimHipMUAData(p,:) = filtStimHipMUAarray - mean(filtStimHipMUAarray(1:(offset*samplingRate)));
                    finalStimStartTimes(qx,1) = stimStartTime;
                    finalStimEndTimes(qx,1) = stimEndTime;
                    finalSFileIDs{qx,1} = finalStimFileID;
                    qx = qx + 1;
                end
            end
            meanStimHbTData = mean(procStimHbTData,1);
            stdStimHbTData = std(procStimHbTData,0,1);
            meanStimCBVData = mean(procStimCBVData,1)*100;
            stdStimCBVData = std(procStimCBVData,0,1)*100;
            meanStimCortMUAData = mean(procStimCortMUAData,1)*100;
            stdStimCortMUAData = std(procStimCortMUAData,0,1)*100;
            meanStimHipMUAData = mean(procStimHipMUAData,1)*100;
            stdStimHipMUAData = std(procStimHipMUAData,0,1)*100;
            % extract LFP from spectrograms associated with the stimuli indecies
            stimCortZhold = [];
            stimHipZhold = [];
            for q = 1:length(finalSFileIDs)
                % load normalized one-second bin data from each file
                stimFileID = finalSFileIDs{q,1};
                stimSpecDataFileID = [animalID '_' stimFileID '_SpecData.mat'];
                stimSpecField = neuralDataType;
                for r = 1:length(AllSpecData.(stimSpecField).fileIDs)
                    if strcmp(AllSpecData.(stimSpecField).fileIDs{r,1},stimSpecDataFileID) == true
                        stimCorticalS_Data = AllSpecData.(stimSpecField).oneSec.normS{r,1};
                        stimHippocampalS_Data = AllSpecData.hippocampus.oneSec.normS{r,1};
                    end
                end
                stimStartTimeIndex = find(T == round(finalStimStartTimes(q,1),1));
                stimDurationIndex = find(T == round(finalStimEndTimes(q,1),1));
                stimCorticalS_Vals = stimCorticalS_Data(:,stimStartTimeIndex:stimDurationIndex);
                stimHippocampalS_Vals = stimHippocampalS_Data(:,stimStartTimeIndex:stimDurationIndex);
                % mean subtract each row with detrend
                transpStimCorticalS_Vals = stimCorticalS_Vals';   % Transpose since detrend goes down columns
                transpStimHippocampalS_Vals = stimHippocampalS_Vals';   % Transpose since detrend goes down columns
                dTStimCortS_Vals = detrend(transpStimCorticalS_Vals,'constant');
                dTStimHipS_Vals = detrend(transpStimHippocampalS_Vals,'constant');
                stimCortZhold = cat(3,stimCortZhold,dTStimCortS_Vals');   % transpose back to original orientation
                stimHipZhold = cat(3,stimHipZhold,dTStimHipS_Vals');   % transpose back to original orientation
            end
            % figure time/frequency axis and average each S data matrix through time
            meanStimCortS = mean(stimCortZhold,3);
            meanStimHipS = mean(stimHipZhold,3);
            % Save figures if desired
            if strcmp(saveFigs,'y') == true
                stimEvoked = figure;
                sgtitle([animalID ' ' dataType ' ' solenoid ' stimulus-evoked averages'])
                subplot(2,3,1);
                plot(timeVector,meanStimCortMUAData,'k')
                hold on
                plot(timeVector,meanStimCortMUAData + stdStimCortMUAData,'color',colors_Manuscript2020('battleship grey'))
                plot(timeVector,meanStimCortMUAData - stdStimCortMUAData,'color',colors_Manuscript2020('battleship grey'))
                title('Cortical MUA')
                xlabel('Time (sec)')
                ylabel('Fold-change (Norm Power)')
                axis tight
                axis square
                subplot(2,3,2);
                imagesc(T,F,meanStimCortS)
                title('Cortical MUA')
                xlabel('Time (sec)')
                ylabel('Freq (Hz)')
                ylim([1,100])
                caxis([-0.5,1])
                set(gca,'Ticklength',[0,0])
                axis xy
                axis square
                subplot(2,3,4);
                plot(timeVector,meanStimHipMUAData,'k')
                hold on
                plot(timeVector,meanStimHipMUAData + stdStimHipMUAData,'color',colors_Manuscript2020('battleship grey'))
                plot(timeVector,meanStimHipMUAData - stdStimHipMUAData,'color',colors_Manuscript2020('battleship grey'))
                title('Hippocampal MUA')
                xlabel('Time (sec)')
                ylabel('Fold-change (Norm Power)')
                axis tight
                axis square
                subplot(2,3,5);
                imagesc(T,F,meanStimHipS)
                title('Hippocampal MUA')
                xlabel('Time (sec)')
                ylabel('Freq (Hz)')
                ylim([1 100])
                caxis([-0.5 1])
                set(gca,'Ticklength',[0,0])
                axis xy
                axis square
                subplot(2,3,[3,6]);
                plot(timeVector,meanStimHbTData,'k')
                hold on
                plot(timeVector,meanStimHbTData + stdStimHbTData,'color',colors_Manuscript2020('battleship grey'))
                plot(timeVector,meanStimHbTData - stdStimHbTData,'color',colors_Manuscript2020('battleship grey'))
                title('Hemodynamics')
                xlabel('Time (sec)')
                ylabel('\DeltaHbT (\muM)')
                axis tight
                axis square
                savefig(stimEvoked,[dirpath animalID '_' dataType '_' solenoid '_StimEvokedAverages']);
                close(stimEvoked)
            end
            % save results
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).CBV_HbT.HbT = meanStimHbTData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).CBV_HbT.HbTStD = stdStimHbTData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).CBV.CBV = meanStimCBVData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).CBV.CBVStD = stdStimCBVData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).MUA.corticalData = meanStimCortMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).MUA.corticalStD = stdStimCortMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).MUA.hippocampalData = meanStimHipMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).MUA.hippocampalStD = stdStimHipMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).timeVector = timeVector;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).LFP.corticalS = meanStimCortS;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).LFP.hippocampalS = meanStimHipS;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).LFP.T = T;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).LFP.F = F;
        end
    end
    cd(rootFolder)
end

end
