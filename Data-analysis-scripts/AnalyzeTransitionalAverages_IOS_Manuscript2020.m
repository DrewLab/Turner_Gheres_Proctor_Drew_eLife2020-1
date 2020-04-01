function [AnalysisResults] = AnalyzeTransitionalAverages_IOS_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose:
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};

%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    % load model
    modelDirectory = [rootFolder '\Support Files\'];
    cd(modelDirectory)
    modelName = 'IOS_SVM_SleepScoringModel.mat';
    load(modelName)
    % go to data and load the model files
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    modelDataFileStruct = dir('*_ModelData.mat');
    modelDataFile = {modelDataFileStruct.name}';
    modelDataFileIDs = char(modelDataFile);
    baselineFileStruct = dir('*_RestingBaselines.mat');
    baselineFile = {baselineFileStruct.name}';
    baselineFileID = char(baselineFile);
    load(baselineFileID)
    samplingRate = 30;
    % go through each file and sleep score the data
    for a = 1:size(modelDataFileIDs,1)
        modelDataFileID = modelDataFileIDs(a,:);
        if a == 1
            load(modelDataFileID)
            dataLength = size(paramsTable,1);
            joinedTable = paramsTable;
            joinedFileList = cell(size(paramsTable,1),1);
            joinedFileList(:) = {modelDataFileID};
        else
            load(modelDataFileID)
            fileIDCells = cell(size(paramsTable,1),1);
            fileIDCells(:) = {modelDataFileID};
            joinedTable = vertcat(joinedTable,paramsTable); %#ok<*AGROW>
            joinedFileList = vertcat(joinedFileList,fileIDCells);
        end
    end
    scoringTable = joinedTable;
    [labels,~] = predict(SVM_MDL,scoringTable);
    % apply a logical patch on the REM events
    REMindex = strcmp(labels,'REM Sleep');
    numFiles = length(labels)/dataLength;
    reshapedREMindex = reshape(REMindex,dataLength,numFiles);
    patchedREMindex = [];
    % patch missing REM indeces due to theta band falling off
    for b = 1:size(reshapedREMindex,2)
        remArray = reshapedREMindex(:,b);
        patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
        patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
    end
    % change labels for each event
    for c = 1:length(labels)
        if patchedREMindex(c,1) == 1
            labels{c,1} = 'REM Sleep';
        end
    end
    % convert strings to numbers for easier comparisons
    labelNumbers = zeros(length(labels),1);
    for d = 1:length(labels)
        if strcmp(labels{d,1},'Not Sleep') == true
            labelNumbers(d,1) = 1;
        elseif strcmp(labels{d,1},'NREM Sleep') == true
            labelNumbers(d,1) = 2;
        elseif strcmp(labels{d,1},'REM Sleep') == true
            labelNumbers(d,1) = 3;
        end
    end
    % reshape
    fileIDs = unique(joinedFileList);
    fileLabels = reshape(labelNumbers,dataLength,size(modelDataFileIDs,1))';
    stringArray = zeros(1,12);
    for d = 1:length(transitions)
        transition = transitions{1,d};
        if strcmp(transition,'AWAKEtoNREM') == true
            for e = 1:12
                if e <= 6
                    stringArray(1,e) = 1;
                else
                    stringArray(1,e) = 2;
                end
            end
        elseif strcmp(transition,'NREMtoAWAKE') == true
            for e = 1:12
                if e <= 6
                    stringArray(1,e) = 2;
                else
                    stringArray(1,e) = 1;
                end
            end
        elseif strcmp(transition,'NREMtoREM') == true
            for e = 1:12
                if e <= 6
                    stringArray(1,e) = 2;
                else
                    stringArray(1,e) = 3;
                end
            end
        elseif strcmp(transition,'REMtoAWAKE') == true
            for e = 1:12
                if e <= 6
                    stringArray(1,e) = 3;
                else
                    stringArray(1,e) = 1;
                end
            end
        end
        % go through and pull out all indeces of matching strings
        idx = 1;
        for f = 1:length(fileIDs)
            fileID = fileIDs{f,1};
            labelArray = fileLabels(f,:);
            indeces = strfind(labelArray,stringArray);
            if isempty(indeces) == false
                for g = 1:length(indeces)
                    data.(transition).files{idx,1} = fileID;
                    data.(transition).startInd(idx,1) = indeces(1,g);
                    idx = idx + 1;
                end
            end
        end
    end
    % extract data
    for h = 1:length(transitions)
        transition = transitions{1,h};
        iqx = 1;
        for i = 1:length(data.(transition).files)
            file = data.(transition).files{i,1};
            startBin = data.(transition).startInd(i,1);
            if startBin > 1 && startBin < (180 - 12)
                [animalID,fileDate,fileID] = GetFileInfo_IOS_Manuscript2020(file);
                strDay = ConvertDate_IOS_Manuscript2020(fileDate);
                procDataFileID = [animalID '_' fileID '_ProcData.mat'];
                load(procDataFileID)
                specDataFileID = [animalID '_' fileID '_SpecData.mat'];
                load(specDataFileID)
                startTime = (startBin - 1)*5;   % sec
                endTime = startTime + (12*5);   % sec
                % whisking data
                [B,A] = butter(3,10/(samplingRate/2),'low');
                filtWhiskAngle = filtfilt(B,A,ProcData.data.whiskerAngle(startTime*samplingRate + 1:endTime*samplingRate));
                % heart rate data
                heartRate = ProcData.data.heartRate(startTime + 1:endTime);
                % EMG
                EMG = (ProcData.data.EMG.emg(startTime*samplingRate + 1:endTime*samplingRate) - RestingBaselines.manualSelection.EMG.emg.(strDay));
                % spectrogram data
                cortical_LHnormS = SpecData.cortical_LH.oneSec.normS;
                cortical_RHnormS = SpecData.cortical_RH.oneSec.normS;
                hippocampusNormS = SpecData.hippocampus.oneSec.normS;
                T = round(SpecData.cortical_LH.oneSec.T,1);
                F = SpecData.cortical_LH.oneSec.F;
                specStartIndex = find(T == startTime);
                specEndIndex = find(T == endTime);
                LH_cortSpec = cortical_LHnormS(:,specStartIndex + 1:specEndIndex);
                RH_cortSpec = cortical_RHnormS(:,specStartIndex + 1:specEndIndex);
                Hip_spec = hippocampusNormS(:,specStartIndex + 1:specEndIndex);
                T_short = T(1:size(LH_cortSpec,2));
                % HbT data
                [D,C] = butter(3,1/(samplingRate/2),'low');
                LH_HbT = ProcData.data.CBV_HbT.adjLH;
                RH_HbT = ProcData.data.CBV_HbT.adjRH;
                filtLH_HbT = filtfilt(D,C,LH_HbT(startTime*samplingRate + 1:endTime*samplingRate));
                filtRH_HbT = filtfilt(D,C,RH_HbT(startTime*samplingRate + 1:endTime*samplingRate));
                data.(transition).whisk(iqx,:) = filtWhiskAngle;
                data.(transition).HR(iqx,:) = heartRate;
                data.(transition).EMG(iqx,:) = EMG;
                data.(transition).LH_cort(:,:,iqx) = LH_cortSpec;
                data.(transition).RH_cort(:,:,iqx) = RH_cortSpec;
                data.(transition).Hip(:,:,iqx) = Hip_spec;
                data.(transition).T_short = T_short;
                data.(transition).F = F;
                data.(transition).LH_HbT(iqx,:) = filtLH_HbT;
                data.(transition).RH_HbT(iqx,:) = filtRH_HbT;
                iqx = iqx + 1;
            end
        end
    end
    % take averages of each behavior
    for d = 1:length(transitions)
        transition = transitions{1,d};
        AnalysisResults.(animalID).Transitions.(transition).whisk = mean(data.(transition).whisk,1);
        AnalysisResults.(animalID).Transitions.(transition).HR = mean(data.(transition).HR,1);
        AnalysisResults.(animalID).Transitions.(transition).EMG = mean(data.(transition).EMG,1);
        AnalysisResults.(animalID).Transitions.(transition).Hip = mean(data.(transition).Hip,3);
        AnalysisResults.(animalID).Transitions.(transition).T = data.(transition).T_short;
        AnalysisResults.(animalID).Transitions.(transition).F = data.(transition).F;
        allCort = cat(3,data.(transition).LH_cort,data.(transition).RH_cort);
        allHbT = cat(1,data.(transition).LH_HbT,data.(transition).RH_HbT);
        AnalysisResults.(animalID).Transitions.(transition).Cort = mean(allCort,3);
        AnalysisResults.(animalID).Transitions.(transition).HbT = mean(allHbT,1);
    end
    % save figures if desired
    if strcmp(saveFigs,'y') == true
        for e = 1:length(transitions)
            transition = transitions{1,e};
            transitionFig = figure;
            sgtitle([animalID ' average ' transition])
            % behavior
            subplot(3,1,1)
            p1 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).HbT))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).HbT,'color',colors_Manuscript2020('rich black'),'LineWidth',2);
            ylabel('\DeltaHbT (%)')
            xlim([0,samplingRate*60])
            yyaxis right
            p2 = plot((1:length(AnalysisResults.(animalID).Transitions.(transition).EMG))/samplingRate,AnalysisResults.(animalID).Transitions.(transition).EMG,'color',colors_Manuscript2020('deep carrot orange'),'LineWidth',2);
            ylabel('EMG (Volts^2)')
            set(gca,'TickLength',[0,0])
            set(gca,'Xticklabel',[])
            set(gca,'box','off')
            axis tight
            legend([p1,p2],'HbT','EMG')
            % LH cort neural
            subplot(3,1,2)
            semilog_imagesc_Manuscript2020(AnalysisResults.(animalID).Transitions.(transition).T,AnalysisResults.(animalID).Transitions.(transition).F,AnalysisResults.(animalID).Transitions.(transition).Cort,'y')
            axis xy
            caxis([-1,2])
            ylabel('Frequency (Hz)')
            set(gca,'Yticklabel','10^1')
            xlim([0,60])
            set(gca,'TickLength',[0,0])
            set(gca,'Xticklabel',[])
            set(gca,'box','off')
            yyaxis right
            ylabel('Cortical LFP')
            % hip
            subplot(3,1,3)
            semilog_imagesc_Manuscript2020(AnalysisResults.(animalID).Transitions.(transition).T,AnalysisResults.(animalID).Transitions.(transition).F,AnalysisResults.(animalID).Transitions.(transition).Hip,'y')
            caxis([-1,2])
            xlabel('Time (sec)')
            ylabel('Frequency (Hz)')
            set(gca,'Yticklabel','10^1')
            xlim([0,60])
            set(gca,'TickLength',[0,0])
            set(gca,'box','off')
            yyaxis right
            ylabel('Hippocampal LFP')
            
            [pathstr,~,~] = fileparts(cd);
            dirpath = [pathstr '/Figures/Transitions/'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(transitionFig,[dirpath animalID '_' transition '_Transition']);
            close(transitionFig)
        end
    end
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
