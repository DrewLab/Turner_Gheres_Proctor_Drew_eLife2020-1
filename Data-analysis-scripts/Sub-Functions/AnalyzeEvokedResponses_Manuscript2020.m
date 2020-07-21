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
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
dataTypes = {'adjLH','adjRH'};

%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    % find and load EventData.mat struct
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
    allSpecStructFileStruct = dir('*_AllSpecStructB.mat');
    allSpecStructFile = {allSpecStructFileStruct.name}';
    allSpecStructFileID = char(allSpecStructFile);
    load(allSpecStructFileID)
    % forest ID sctruct
    forestScoringResultsID = 'Forest_ScoringResults.mat';
    load(forestScoringResultsID,'-mat')
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
    whiskCriteriaNames = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
    % filter the EventData.mat structure for whisking events that meet the desired criteria
    for aa = 1:length(dataTypes)
        dataType = dataTypes{1,aa};
        neuralDataType = ['cortical_' dataType(4:end)];
        % pull a few necessary numbers from the EventData.mat struct such as trial duration and sampling rate
        samplingRate = EventData.CBV_HbT.(dataType).whisk.samplingRate;
        specSamplingRate = 10;
        trialDuration_sec = EventData.CBV_HbT.(dataType).whisk.trialDuration_sec;
        timeVector = (0:(EventData.CBV_HbT.(dataType).whisk.epoch.duration*samplingRate))/samplingRate - EventData.CBV_HbT.(dataType).whisk.epoch.offset;
        offset = EventData.CBV_HbT.(dataType).whisk.epoch.offset;
        for bb = 1:length(whiskCriteriaNames)
            whiskCriteriaName = whiskCriteriaNames{1,bb};
            if strcmp(whiskCriteriaName,'ShortWhisks') == true
                WhiskCriteria = whiskCriteriaA;
            elseif strcmp(whiskCriteriaName,'IntermediateWhisks') == true
                WhiskCriteria = whiskCriteriaB;
            elseif strcmp(whiskCriteriaName,'LongWhisks') == true
                WhiskCriteria = whiskCriteriaC;
            end
            [whiskLogical] = FilterEvents_IOS_Manuscript2020(EventData.CBV_HbT.(dataType).whisk,WhiskCriteria);
            combWhiskLogical = logical(whiskLogical);
            [allWhiskHbTData] = EventData.CBV_HbT.(dataType).whisk.data(combWhiskLogical,:);
            [allWhiskCBVData] = EventData.CBV.(dataType).whisk.NormData(combWhiskLogical,:);
            [allWhiskCorticalMUAData] = EventData.(neuralDataType).muaPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskHippocampalMUAData] = EventData.hippocampus.muaPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskCorticalGamData] = EventData.(neuralDataType).gammaBandPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskHippocampalGamData] = EventData.hippocampus.gammaBandPower.whisk.NormData(combWhiskLogical,:);
            [allWhiskFileIDs] = EventData.CBV_HbT.(dataType).whisk.fileIDs(combWhiskLogical,:);
            [allWhiskEventTimes] = EventData.CBV_HbT.(dataType).whisk.eventTime(combWhiskLogical,:);
            allWhiskDurations = EventData.CBV_HbT.(dataType).whisk.duration(combWhiskLogical,:);
            % decimate the file list to only include those files that occur within the desired number of target minutes
            [finalWhiskHbTData,finalWhiskFileIDs,~,finalWhiskFileEventTimes] = RemoveInvalidData_IOS_Manuscript2020(allWhiskHbTData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCBVData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(allWhiskCBVData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCorticalMUAData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(allWhiskCorticalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskHippocampalMUAData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(allWhiskHippocampalMUAData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskCorticalGamData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(allWhiskCorticalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            [finalWhiskHippocampalGamData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(allWhiskHippocampalGamData,allWhiskFileIDs,allWhiskDurations,allWhiskEventTimes,ManualDecisions);
            % lowpass filter each whisking event and mean-subtract by the first 2 seconds
            clear procWhiskHbTData procWhiskCBVData procWhiskCorticalMUAData procWhiskHippocampalMUAData procWhiskCorticalGamData procWhiskHippocampalGamData finalWhiskStartTimes finalWhiskEndTimes finalWhiskFiles
            dd = 1;
            for cc = 1:size(finalWhiskHbTData,1)
                whiskStartTime = round(finalWhiskFileEventTimes(cc,1),1) - 2;
                whiskEndTime = whiskStartTime + 12;
                finalWhiskFileID = finalWhiskFileIDs{cc,1};
                if whiskStartTime >= 0.5 && whiskEndTime <= (trialDuration_sec - 0.5)
                    whiskHbTarray = finalWhiskHbTData(cc,:);
                    whiskCBVarray = finalWhiskCBVData(cc,:);
                    whiskCorticalMUAarray = finalWhiskCorticalMUAData(cc,:);
                    whiskHippocampalMUAarray = finalWhiskHippocampalMUAData(cc,:);
                    whiskCorticalGamArray = finalWhiskCorticalGamData(cc,:);
                    whiskHippocampalGamArray = finalWhiskHippocampalGamData(cc,:);
                    filtWhiskHbTarray = sgolayfilt(whiskHbTarray,3,17);
                    filtWhiskCBVarray = sgolayfilt(whiskCBVarray,3,17);
                    filtWhiskCorticalMUAarray = sgolayfilt(whiskCorticalMUAarray,3,17);
                    filtWhiskHippocampalMUAarray = sgolayfilt(whiskHippocampalMUAarray,3,17);
                    filtWhiskCorticalGamArray = sgolayfilt(whiskCorticalGamArray,3,17);
                    filtWhiskHippocampalGamArray = sgolayfilt(whiskHippocampalGamArray,3,17);
                    procWhiskHbTData(dd,:) = filtWhiskHbTarray - mean(filtWhiskHbTarray(1:(offset*samplingRate))); %#ok<*AGROW>
                    procWhiskCBVData(dd,:) = filtWhiskCBVarray - mean(filtWhiskCBVarray(1:(offset*samplingRate)));
                    procWhiskCorticalMUAData(dd,:) = filtWhiskCorticalMUAarray - mean(filtWhiskCorticalMUAarray(1:(offset*samplingRate)));
                    procWhiskHippocampalMUAData(dd,:) = filtWhiskHippocampalMUAarray - mean(filtWhiskHippocampalMUAarray(1:(offset*samplingRate)));
                    procWhiskCorticalGamData(dd,:) = filtWhiskCorticalGamArray - mean(filtWhiskCorticalGamArray(1:(offset*samplingRate)));
                    procWhiskHippocampalGamData(dd,:) = filtWhiskHippocampalGamArray - mean(filtWhiskHippocampalGamArray(1:(offset*samplingRate)));
                    finalWhiskStartTimes(dd,1) = whiskStartTime;
                    finalWhiskEndTimes(dd,1) = whiskEndTime;
                    finalWhiskFiles{dd,1} = finalWhiskFileID;
                    dd = dd + 1;
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
            meanWhiskCorticalGamData = mean(procWhiskCorticalGamData,1)*100;
            stdWhiskCorticalGamData = std(procWhiskCorticalGamData,0,1)*100;
            meanWhiskHippocampalGamData = mean(procWhiskHippocampalGamData,1)*100;
            stdWhiskHippocampalGamData = std(procWhiskHippocampalGamData,0,1)*100;
            % extract LFP from spectrograms associated with the whisking indecies
            whiskCorticalZhold = [];
            whiskHippocampalZhold = [];
            for ee = 1:length(finalWhiskFiles)
                % load normalized one-second bin data from each file
                whiskFileID = finalWhiskFiles{ee,1};
                whiskSpecDataFileID = [animalID '_' whiskFileID '_SpecDataB.mat'];
                whiskSpecField = neuralDataType;
                for ff = 1:length(AllSpecData.(whiskSpecField).fileIDs)
                    if strcmp(AllSpecData.(whiskSpecField).fileIDs{ff,1},whiskSpecDataFileID) == true
                        whiskCorticalS_Data = AllSpecData.(whiskSpecField).normS{ff,1};
                        whiskHippocampalS_Data = AllSpecData.hippocampus.normS{ff,1};
                        F = AllSpecData.(whiskSpecField).F{ff,1};
                        T = round(AllSpecData.(whiskSpecField).T{ff,1},1);
                    end
                end
                whiskStartTimeIndex = find(T == round(finalWhiskStartTimes(ee,1),1));
                whiskStartTimeIndex = whiskStartTimeIndex(1);
                whiskDurationIndex = find(T == round(finalWhiskEndTimes(ee,1),1));
                whiskDurationIndex = whiskDurationIndex(end);
                whiskCorticalS_Vals = whiskCorticalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
                whiskHippocampalS_Vals = whiskHippocampalS_Data(:,whiskStartTimeIndex:whiskDurationIndex);
                % mean subtract each row with detrend
                transpWhiskCorticalS_Vals = whiskCorticalS_Vals';   % Transpose since detrend goes down columns
                transpWhiskHippocampalS_Vals = whiskHippocampalS_Vals';
                dTWhiskCorticalS_Vals = transpWhiskCorticalS_Vals;
                dTWhiskCorticalS_Vals = dTWhiskCorticalS_Vals(1:12*specSamplingRate + 1,:);
                dTWhiskHippocampalS_Vals = transpWhiskHippocampalS_Vals;
                dTWhiskHippocampalS_Vals = dTWhiskHippocampalS_Vals(1:12*specSamplingRate + 1,:);
                whiskCorticalZhold = cat(3,whiskCorticalZhold,dTWhiskCorticalS_Vals');   % transpose back to original orientation
                whiskHippocampalZhold = cat(3,whiskHippocampalZhold,dTWhiskHippocampalS_Vals');
            end
            % figure time/frequency axis and average each S data matrix through time
            meanWhiskCorticalS = mean(whiskCorticalZhold,3);
            meanWhiskHippocampalS = mean(whiskHippocampalZhold,3);
%             % Save figures if desired
%             if strcmp(saveFigs,'y') == true
%                 whiskEvoked = figure;
%                 sgtitle([animalID ' ' dataType ' whisking-evoked averages - ' whiskCriteriaName])
                T2 = -2:(1/specSamplingRate):10;
%                 subplot(2,3,1);
%                 p1 = plot(timeVector,meanWhiskCorticalMUAData,'k');
%                 hold on
%                 plot(timeVector,meanWhiskCorticalMUAData + stdWhiskCorticalMUAData,'color',colors_Manuscript2020('battleship grey'))
%                 plot(timeVector,meanWhiskCorticalMUAData - stdWhiskCorticalMUAData,'color',colors_Manuscript2020('battleship grey'))
%                 p2 = plot(timeVector,meanWhiskCorticalGamData,'r');
%                 hold on
%                 plot(timeVector,meanWhiskCorticalGamData + stdWhiskCorticalGamData,'color',colors_Manuscript2020('deep carrot orange'))
%                 plot(timeVector,meanWhiskCorticalGamData - stdWhiskCorticalGamData,'color',colors_Manuscript2020('deep carrot orange'))
%                 title('Cortical MUA/Gam')
%                 xlabel('Time (sec)')
%                 ylabel('Fold-change (Norm Power)')
%                 legend([p1,p2],'MUA','Gam')
%                 axis tight
%                 axis square
%                 set(gca,'box','off')
%                 subplot(2,3,4);
%                 plot(timeVector,meanWhiskHippocampalMUAData,'k')
%                 hold on
%                 plot(timeVector,meanWhiskHippocampalMUAData + stdWhiskHippocampalMUAData,'color',colors_Manuscript2020('battleship grey'))
%                 plot(timeVector,meanWhiskHippocampalMUAData - stdWhiskHippocampalMUAData,'color',colors_Manuscript2020('battleship grey'))
%                 plot(timeVector,meanWhiskHippocampalGamData,'r')
%                 hold on
%                 plot(timeVector,meanWhiskHippocampalGamData + stdWhiskHippocampalGamData,'color',colors_Manuscript2020('deep carrot orange'))
%                 plot(timeVector,meanWhiskHippocampalGamData - stdWhiskHippocampalGamData,'color',colors_Manuscript2020('deep carrot orange'))
%                 title([animalID ' ' dataType ' ' whiskCriteriaName ' whisking-evoked averages'])
%                 title('Hippocampal MUA/Gam')
%                 xlabel('Time (sec)')
%                 ylabel('Fold-change (Norm Power)')
%                 axis tight
%                 axis square
%                 set(gca,'box','off')
%                 subplot(2,3,2);
%                 imagesc(T2,F,meanWhiskCorticalS)
%                 title('Cortical LFP')
%                 xlabel('Time (sec)')
%                 ylabel('Freq (Hz)')
%                 ylim([1 100])
%                 caxis([-0.5 1])
%                 set(gca,'Ticklength',[0,0])
%                 axis xy
%                 axis square
%                 subplot(2,3,5);
%                 imagesc(T2,F,meanWhiskHippocampalS)
%                 title('Hippocampal LFP')
%                 xlabel('Time (sec)')
%                 ylabel('Freq (Hz)')
%                 ylim([1,100])
%                 caxis([-0.5,1])
%                 set(gca,'Ticklength',[0 0])
%                 axis xy
%                 axis square
%                 set(gca,'box','off')
%                 subplot(2,3,[3,6]);
%                 plot(timeVector,meanWhiskHbTData,'k')
%                 hold on
%                 plot(timeVector,meanWhiskHbTData + stdWhiskHbTData,'color',colors_Manuscript2020('battleship grey'))
%                 plot(timeVector,meanWhiskHbTData - stdWhiskHbTData,'color',colors_Manuscript2020('battleship grey'))
%                 title('Hemodynamic response')
%                 xlabel('Time (sec)')
%                 ylabel('\DeltaHbT (\muM)')
%                 axis tight
%                 axis square
%                 set(gca,'box','off')
%                 % save figure
%                 [pathstr,~,~] = fileparts(cd);
%                 dirpath = [pathstr '/Figures/Evoked Responses/'];
%                 if ~exist(dirpath, 'dir')
%                     mkdir(dirpath);
%                 end
%                 savefig(whiskEvoked,[dirpath animalID '_' dataType '_' whiskCriteriaName '_WhiskEvokedAverages']);
%                 close(whiskEvoked)
%             end
            % save results
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV_HbT.HbT = meanWhiskHbTData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV_HbT.HbTStD = stdWhiskHbTData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.CBV = meanWhiskCBVData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).CBV.CBVStD = stdWhiskCBVData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.corticalData = meanWhiskCorticalMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.corticalStD = stdWhiskCorticalMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.hippocampalData = meanWhiskHippocampalMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).MUA.hippocampalStD = stdWhiskHippocampalMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).Gam.corticalData = meanWhiskCorticalGamData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).Gam.corticalStD = stdWhiskCorticalGamData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).Gam.hippocampalData = meanWhiskHippocampalGamData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).Gam.hippocampalStD = stdWhiskHippocampalGamData;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).timeVector = timeVector;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.corticalS = meanWhiskCorticalS;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.hippocampalS = meanWhiskHippocampalS;
            AnalysisResults.(animalID).EvokedAvgs.Whisk.(dataType).(whiskCriteriaName).LFP.T = T2;
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
        for gg = 1:length(stimCriteriaNames)
            stimCriteriaName = stimCriteriaNames{1,gg};
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
            allStimFilter = FilterEvents_IOS_Manuscript2020(EventData.CBV_HbT.(dataType).stim,stimCriteria);
            [allStimHbTData] = EventData.CBV_HbT.(dataType).stim.data(allStimFilter,:);
            [allStimCBVData] = EventData.CBV.(dataType).stim.NormData(allStimFilter,:);
            [allStimCortMUAData] = EventData.(neuralDataType).muaPower.stim.NormData(allStimFilter,:);
            [allStimHipMUAData] = EventData.hippocampus.muaPower.stim.NormData(allStimFilter,:);
            [allStimCortGamData] = EventData.(neuralDataType).gammaBandPower.stim.NormData(allStimFilter,:);
            [allStimHipGamData] = EventData.hippocampus.gammaBandPower.stim.NormData(allStimFilter,:);
            [allStimFileIDs] = EventData.CBV_HbT.(dataType).stim.fileIDs(allStimFilter,:);
            [allStimEventTimes] = EventData.CBV_HbT.(dataType).stim.eventTime(allStimFilter,:);
            allStimDurations = zeros(length(allStimEventTimes),1);
            % decimate the file list to only include those files that occur within the desired number of target minutes
            [finalStimHbTData,finalStimFileIDs,~,finalStimFileEventTimes] = RemoveInvalidData_IOS_Manuscript2020(allStimHbTData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCBVData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(allStimCBVData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortMUAData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(allStimCortMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimHipMUAData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(allStimHipMUAData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimCortGamData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(allStimCortGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            [finalStimHipGamData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(allStimHipGamData,allStimFileIDs,allStimDurations,allStimEventTimes,ManualDecisions);
            % lowpass filter each whisking event and mean-subtract by the first 2 seconds
            clear procStimHbTData procStimCBVData procStimCortMUAData procStimHipMUAData procStimCortGamData procStimHipGamData finalStimStartTimes finalStimEndTimes finalStimFiles
            ii = 1;
            for hh = 1:size(finalStimHbTData,1)
                stimStartTime = round(finalStimFileEventTimes(hh,1),1) - 2;
                stimEndTime = stimStartTime + 12;
                finalStimFileID = finalStimFileIDs{hh,1};
                if stimStartTime >= 0.5 && stimEndTime <= (trialDuration_sec - 0.5)
                    stimHbTarray = finalStimHbTData(hh,:);
                    stimCBVarray = finalStimCBVData(hh,:);
                    stimCortMUAarray = finalStimCortMUAData(hh,:);
                    stimHipMUAarray = finalStimHipMUAData(hh,:);
                    stimCortGamArray = finalStimCortGamData(hh,:);
                    stimHipGamArray = finalStimHipGamData(hh,:);
                    filtStimHbTarray = sgolayfilt(stimHbTarray,3,17);
                    filtStimCBVarray = sgolayfilt(stimCBVarray,3,17);
                    filtStimCortMUAarray = sgolayfilt(stimCortMUAarray,3,17);
                    filtStimHipMUAarray = sgolayfilt(stimHipMUAarray,3,17);
                    filtStimCortGamArray = sgolayfilt(stimCortGamArray,3,17);
                    filtStimHipGamArray = sgolayfilt(stimHipGamArray,3,17);
                    procStimHbTData(hh,:) = filtStimHbTarray - mean(filtStimHbTarray(1:(offset*samplingRate)));
                    procStimCBVData(hh,:) = filtStimCBVarray - mean(filtStimCBVarray(1:(offset*samplingRate)));
                    procStimCortMUAData(hh,:) = filtStimCortMUAarray - mean(filtStimCortMUAarray(1:(offset*samplingRate)));
                    procStimHipMUAData(hh,:) = filtStimHipMUAarray - mean(filtStimHipMUAarray(1:(offset*samplingRate)));
                    procStimCortGamData(hh,:) = filtStimCortGamArray - mean(filtStimCortGamArray(1:(offset*samplingRate)));
                    procStimHipGamData(hh,:) = filtStimHipGamArray - mean(filtStimHipGamArray(1:(offset*samplingRate)));
                    finalStimStartTimes(ii,1) = stimStartTime;
                    finalStimEndTimes(ii,1) = stimEndTime;
                    finalStimFiles{ii,1} = finalStimFileID;
                    ii = ii + 1;
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
            meanStimCortGamData = mean(procStimCortGamData,1)*100;
            stdStimCortGamData = std(procStimCortGamData,0,1)*100;
            meanStimHipGamData = mean(procStimHipGamData,1)*100;
            stdStimHipGamData = std(procStimHipGamData,0,1)*100;
            % extract LFP from spectrograms associated with the stimuli indecies
            stimCortZhold = [];
            stimHipZhold = [];
            for jj = 1:length(finalStimFiles)
                % load normalized one-second bin data from each file
                stimFileID = finalStimFiles{jj,1};
                stimSpecDataFileID = [animalID '_' stimFileID '_SpecDataB.mat'];
                stimSpecField = neuralDataType;
                for kk = 1:length(AllSpecData.(stimSpecField).fileIDs)
                    if strcmp(AllSpecData.(stimSpecField).fileIDs{kk,1},stimSpecDataFileID) == true
                        stimCorticalS_Data = AllSpecData.(stimSpecField).normS{kk,1};
                        stimHippocampalS_Data = AllSpecData.hippocampus.normS{kk,1};
                    end
                end
                stimStartTimeIndex = find(T == round(finalStimStartTimes(jj,1),1));
                stimStartTimeIndex = stimStartTimeIndex(1);
                stimDurationIndex = find(T == round(finalStimEndTimes(jj,1),1));
                stimDurationIndex = stimDurationIndex(end);
                stimCorticalS_Vals = stimCorticalS_Data(:,stimStartTimeIndex:stimDurationIndex);
                stimHippocampalS_Vals = stimHippocampalS_Data(:,stimStartTimeIndex:stimDurationIndex);
                % mean subtract each row with detrend
                transpStimCorticalS_Vals = stimCorticalS_Vals';
                transpStimHippocampalS_Vals = stimHippocampalS_Vals';
                dTStimCortS_Vals = transpStimCorticalS_Vals;
                dTStimCortS_Vals = dTStimCortS_Vals(1:12*specSamplingRate + 1,:);
                dTStimHipS_Vals = transpStimHippocampalS_Vals;
                dTStimHipS_Vals = dTStimHipS_Vals(1:12*specSamplingRate + 1,:);
                stimCortZhold = cat(3,stimCortZhold,dTStimCortS_Vals');
                stimHipZhold = cat(3,stimHipZhold,dTStimHipS_Vals');
            end
            % figure time/frequency axis and average each S data matrix through time
            meanStimCortS = mean(stimCortZhold,3);
            meanStimHipS = mean(stimHipZhold,3);
            % Save figures if desired
%             if strcmp(saveFigs,'y') == true
%                 stimEvoked = figure;
%                 sgtitle([animalID ' ' dataType ' ' solenoid ' stimulus-evoked averages'])
%                 subplot(2,3,1);
%                 p1 = plot(timeVector,meanStimCortMUAData,'k');
%                 hold on
%                 plot(timeVector,meanStimCortMUAData + stdStimCortMUAData,'color',colors_Manuscript2020('battleship grey'))
%                 plot(timeVector,meanStimCortMUAData - stdStimCortMUAData,'color',colors_Manuscript2020('battleship grey'))
%                 p2 = plot(timeVector,meanStimCortGamData,'r');
%                 hold on
%                 plot(timeVector,meanStimCortGamData + stdStimCortGamData,'color',colors_Manuscript2020('deep carrot orange'))
%                 plot(timeVector,meanStimCortGamData - stdStimCortGamData,'color',colors_Manuscript2020('deep carrot orange'))
%                 title('Cortical MUA/Gam')
%                 xlabel('Time (sec)')
%                 ylabel('Fold-change (Norm Power)')
%                 legend([p1,p2],'MUA','Gam')
%                 axis tight
%                 axis square
%                 set(gca,'box','off')
%                 subplot(2,3,2);
%                 imagesc(T2,F,(meanStimCortS))
%                 title('Cortical MUA')
%                 xlabel('Time (sec)')
%                 ylabel('Freq (Hz)')
%                 ylim([1,100])
%                 caxis([-0.5,1])
%                 set(gca,'Ticklength',[0,0])
%                 axis xy
%                 axis square
%                 set(gca,'box','off')
%                 subplot(2,3,4);
%                 plot(timeVector,meanStimHipMUAData,'k')
%                 hold on
%                 plot(timeVector,meanStimHipMUAData + stdStimHipMUAData,'color',colors_Manuscript2020('battleship grey'))
%                 plot(timeVector,meanStimHipMUAData - stdStimHipMUAData,'color',colors_Manuscript2020('battleship grey'))
%                 plot(timeVector,meanStimHipGamData,'r')
%                 hold on
%                 plot(timeVector,meanStimHipGamData + stdStimHipGamData,'color',colors_Manuscript2020('deep carrot orange'))
%                 plot(timeVector,meanStimHipGamData - stdStimHipGamData,'color',colors_Manuscript2020('deep carrot orange'))
%                 title('Hippocampal MUA/Gam')
%                 xlabel('Time (sec)')
%                 ylabel('Fold-change (Norm Power)')
%                 axis tight
%                 axis square
%                 set(gca,'box','off')
%                 subplot(2,3,5);
%                 imagesc(T2,F,meanStimHipS)
%                 title('Hippocampal MUA')
%                 xlabel('Time (sec)')
%                 ylabel('Freq (Hz)')
%                 ylim([1,100])
%                 caxis([-0.5,1])
%                 set(gca,'Ticklength',[0,0])
%                 axis xy
%                 axis square
%                 set(gca,'box','off')
%                 subplot(2,3,[3,6]);
%                 plot(timeVector,meanStimHbTData,'k')
%                 hold on
%                 plot(timeVector,meanStimHbTData + stdStimHbTData,'color',colors_Manuscript2020('battleship grey'))
%                 plot(timeVector,meanStimHbTData - stdStimHbTData,'color',colors_Manuscript2020('battleship grey'))
%                 title('Hemodynamics')
%                 xlabel('Time (sec)')
%                 ylabel('\DeltaHbT (\muM)')
%                 axis tight
%                 axis square
%                 set(gca,'box','off')
%                 savefig(stimEvoked,[dirpath animalID '_' dataType '_' solenoid '_StimEvokedAverages']);
%                 close(stimEvoked)
%             end
            % save results
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).count = size(procStimHipMUAData,1); 
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).CBV_HbT.HbT = meanStimHbTData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).CBV_HbT.HbTStD = stdStimHbTData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).CBV.CBV = meanStimCBVData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).CBV.CBVStD = stdStimCBVData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).MUA.corticalData = meanStimCortMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).MUA.corticalStD = stdStimCortMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).MUA.hippocampalData = meanStimHipMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).MUA.hippocampalStD = stdStimHipMUAData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).Gam.corticalData = meanStimCortGamData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).Gam.corticalStD = stdStimCortGamData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).Gam.hippocampalData = meanStimHipGamData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).Gam.hippocampalStD = stdStimHipGamData;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).timeVector = timeVector;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).LFP.corticalS = meanStimCortS;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).LFP.hippocampalS = meanStimHipS;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).LFP.T = T2;
            AnalysisResults.(animalID).EvokedAvgs.Stim.(dataType).(solenoid).LFP.F = F;
        end
    end
end
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')
end

