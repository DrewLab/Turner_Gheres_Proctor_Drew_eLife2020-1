function [AnalysisResults] = AnalyzePowerSpectrum_2P_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Analyze the spectral power of hemodynamic and neural signals.
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T115','T116','T117','T118','T125','T126'};
modelType = 'Manual';
params.minTime.Rest = 10;   % seconds
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 70;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/2P Data/'];
    cd(dataLocation)
    % Character list of all MergedData files
    mergedDirectory = dir('*_MergedData.mat');
    mergedDataFiles = {mergedDirectory.name}';
    mergedDataFileIDs = char(mergedDataFiles);
    % character list of all TrainingData files
    trainingDirectory = dir('*_TrainingData.mat');
    trainingDataFiles = {trainingDirectory.name}';
    trainingDataFileIDs = char(trainingDataFiles);
    % find and load EventData.mat struct
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
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
    fileBreaks = strfind(restDataFileID, '_');
    animalID = restDataFileID(1:fileBreaks(1) - 1);
    samplingRate = RestData.vesselDiameter.data.samplingRate;
    RestCriteria.Fieldname = {'durations'};
    RestCriteria.Comparison = {'gt'};
    RestCriteria.Value = {params.minTime.Rest};
    WhiskCriteria.Fieldname = {'duration','duration'};
    WhiskCriteria.Comparison = {'gt','lt'};
    WhiskCriteria.Value = {1,5};
    % lowpass filter and detrend each segment
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    
    %% Analyze power spectra during periods of rest
    % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
    [restLogical] = FilterEvents_2P_Manuscript2020(RestData.vesselDiameter.data,RestCriteria);
    restLogical = logical(restLogical);
    restingData = RestData.vesselDiameter.data.data(restLogical,:);
    restFileIDs = RestData.vesselDiameter.data.fileIDs(restLogical,:);
    restVesselIDs = RestData.vesselDiameter.data.vesselIDs(restLogical,:);
    restEventTimes = RestData.vesselDiameter.data.eventTimes(restLogical,:);
    restDurations = RestData.vesselDiameter.data.durations(restLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [finalRestData,finalRestFileIDs,finalRestVesselIDs,~,~] = RemoveInvalidData_2P_Manuscript2020(restingData,restFileIDs,restVesselIDs,restDurations,restEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    for aa = 1:length(finalRestData)
        restStrDay = ConvertDate_2P_Manuscript2020(finalRestFileIDs{aa,1}(1:6));
        if length(finalRestData{aa,1}) < params.minTime.Rest*samplingRate
            restChunkSampleDiff = params.minTime.Rest*samplingRate - length(finalRestData{aa,1});
            restPad = (ones(1,restChunkSampleDiff))*finalRestData{aa,1}(end);
            arrayRestData = horzcat(finalRestData{aa,1},restPad); %#ok<*AGROW>
            normRestData = (arrayRestData - RestingBaselines.vesselDiameter.data.(finalRestVesselIDs{aa,1}).(restStrDay))./RestingBaselines.vesselDiameter.data.(finalRestVesselIDs{aa,1}).(restStrDay);
            procRestData{aa,1} = filtfilt(sos,g,detrend(normRestData,'constant'));
        else
            normRestData = (finalRestData{aa,1}(1:(params.minTime.Rest*samplingRate)) - RestingBaselines.manualSelection.vesselDiameter.data.(finalRestVesselIDs{aa,1}).(restStrDay))./RestingBaselines.manualSelection.vesselDiameter.data.(finalRestVesselIDs{aa,1}).(restStrDay);
            procRestData{aa,1} = filtfilt(sos,g,detrend(normRestData,'constant'));
        end
    end
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontal)
    reshapedRestData = zeros(length(procRestData{1,1}),length(procRestData));
    for bb = 1:length(procRestData)
        reshapedRestData(:,bb) = procRestData{bb,1};
    end
    % remove veins from the artery list
    uniqueRestVesselIDs = unique(finalRestVesselIDs);
    dd = 1;
    for cc = 1:length(uniqueRestVesselIDs)
        if strcmp(uniqueRestVesselIDs{cc,1}(1),'V') == false
            restArterioleIDs{dd,1} = uniqueRestVesselIDs{cc,1};
            dd = dd + 1;
        end
    end
    % split the data based on different arteries
    for ee = 1:length(restArterioleIDs)
        restArterioleID = restArterioleIDs{ee,1};
        gg = 1;
        for ff = 1:length(finalRestVesselIDs)
            if strcmp(restArterioleID,finalRestVesselIDs{ff,1}) == true
                restArterioleData.(restArterioleID)(:,gg) = reshapedRestData(:,ff);
                gg = gg + 1;
            end
        end
    end
    % parameters for mtspectrumc - information available in function
    params.tapers = [3,5];
    params.pad = 1;
    params.Fs = samplingRate;
    params.fpass = [0,1];
    params.trialave = 1;
    params.err = [2,0.05];
    % calculate the power spectra of the desired signals
    for hh = 1:length(restArterioleIDs)
        restVID = restArterioleIDs{hh,1};
        [rest_S{hh,1},rest_f{hh,1},rest_sErr{hh,1}] = mtspectrumc_Manuscript2020(restArterioleData.(restVID),params);
        % save data and figures
        AnalysisResults.(animalID).PowerSpectra.Rest.(restVID).S = rest_S{hh,1};
        AnalysisResults.(animalID).PowerSpectra.Rest.(restVID).f = rest_f{hh,1};
        AnalysisResults.(animalID).PowerSpectra.Rest.(restVID).sErr = rest_sErr{hh,1};
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            % awake rest summary figures
            RestPower = figure;
            loglog(rest_f{hh,1},rest_S{hh,1},'k')
            hold on;
            loglog(rest_f{hh,1},rest_sErr{hh,1},'color',colors_Manuscript2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Power');
            title([animalID  ' ' restVID ' vessel diameter power during awake rest']);
            set(gca,'Ticklength',[0,0]);
            xlim([0,1])
            axis square
            set(gca,'box','off')
            [pathstr,~,~] = fileparts(cd);
            dirpath = [pathstr '/Figures/Vessel Power Spectrum/'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(RestPower,[dirpath animalID '_' restVID '_RestPowerSpectra']);
            close(RestPower)
        end
    end
    
    %% Analyze power spectra during periods of whisking
    % use the RestCriteria we specified earlier to find unstim resting events that are greater than the criteria
    [whiskLogical] = FilterEvents_2P_Manuscript2020(EventData.vesselDiameter.data.whisk,WhiskCriteria);
    whiskLogical = logical(whiskLogical);
    whiskingData = EventData.vesselDiameter.data.whisk.data(whiskLogical,:);
    whiskFileIDs = EventData.vesselDiameter.data.whisk.fileIDs(whiskLogical,:);
    whiskVesselIDs = EventData.vesselDiameter.data.whisk.vesselIDs(whiskLogical,:);
    whiskEventTimes = EventData.vesselDiameter.data.whisk.eventTime(whiskLogical,:);
    whiskDurations = EventData.vesselDiameter.data.whisk.duration(whiskLogical,:);
    % decimate the file list to only include those files that occur within the desired number of target minutes
    [finalWhiskData,finalWhiskFileIDs,finalWhiskVesselIDs,~,~] = RemoveInvalidData_2P_Manuscript2020(whiskingData,whiskFileIDs,whiskVesselIDs,whiskDurations,whiskEventTimes,ManualDecisions);
    % only take the first 10 seconds of the epoch. occassionunstimy a sample gets lost from rounding during the
    % original epoch create so we can add a sample of two back to the end for those just under 10 seconds
    for aa = 1:length(finalWhiskData)
        whiskStrDay = ConvertDate_2P_Manuscript2020(finalWhiskFileIDs{aa,1}(1:6));
        normWhiskData = (finalWhiskData(aa,:) - RestingBaselines.manualSelection.vesselDiameter.data.(finalWhiskVesselIDs{aa,1}).(whiskStrDay))./RestingBaselines.manualSelection.vesselDiameter.data.(finalWhiskVesselIDs{aa,1}).(whiskStrDay);
        procWhiskData{aa,1} = sgolayfilt(detrend(normWhiskData,'constant'),3,17);
    end
    % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontal)
    reshapedWhiskData = zeros(length(procWhiskData{1,1}),length(procWhiskData));
    for bb = 1:length(procWhiskData)
        reshapedWhiskData(:,bb) = procWhiskData{bb,1};
    end
    % remove veins from the artery list
    uniqueWhiskVesselIDs = unique(finalWhiskVesselIDs);
    dd = 1;
    for cc = 1:length(uniqueWhiskVesselIDs)
        if strcmp(uniqueWhiskVesselIDs{cc,1}(1),'V') == false
            whiskArterioleIDs{dd,1} = uniqueWhiskVesselIDs{cc,1};
            dd = dd + 1;
        end
    end
    % split the data based on different arteries
    for ee = 1:length(whiskArterioleIDs)
        whiskArterioleID = whiskArterioleIDs{ee,1};
        gg = 1;
        for ff = 1:length(finalWhiskVesselIDs)
            if strcmp(whiskArterioleID,finalWhiskVesselIDs{ff,1}) == true
                whiskArterioleData.(whiskArterioleID)(:,gg) = reshapedWhiskData(1:params.minTime.Rest*samplingRate,ff);
                gg = gg + 1;
            end
        end
    end
    % parameters for mtspectrumc - information available in function
    params.tapers = [3,5];
    params.pad = 1;
    params.Fs = samplingRate;
    params.fpass = [0,1];
    params.trialave = 1;
    params.err = [2,0.05];
    % calculate the power spectra of the desired signals
    for hh = 1:length(whiskArterioleIDs)
        whiskVID = whiskArterioleIDs{hh,1};
        [whisk_S{hh,1},whisk_f{hh,1},whisk_sErr{hh,1}] = mtspectrumc_Manuscript2020(whiskArterioleData.(whiskVID),params);
        % save data and figures
        AnalysisResults.(animalID).PowerSpectra.Whisk.(whiskVID).S = whisk_S{hh,1};
        AnalysisResults.(animalID).PowerSpectra.Whisk.(whiskVID).f = whisk_f{hh,1};
        AnalysisResults.(animalID).PowerSpectra.Whisk.(whiskVID).sErr = whisk_sErr{hh,1};
        % save figures if desired
        if strcmp(saveFigs,'y') == true
            % awake rest summary figures
            WhiskPower = figure;
            loglog(whisk_f{hh,1},whisk_S{hh,1},'k')
            hold on;
            loglog(whisk_f{hh,1},whisk_sErr{hh,1},'color',colors_Manuscript2020('battleship grey'))
            xlabel('Freq (Hz)');
            ylabel('Power');
            title([animalID  ' ' whiskVID ' vessel diameter power during awake whisking']);
            set(gca,'Ticklength',[0,0]);
            xlim([0,1])
            axis square
            set(gca,'box','off')
            [pathstr,~,~] = fileparts(cd);
            dirpath = [pathstr '/Figures/Vessel Power Spectrum/'];
            if ~exist(dirpath,'dir')
                mkdir(dirpath);
            end
            savefig(WhiskPower,[dirpath animalID '_' whiskVID '_WhiskPowerSpectra']);
            close(WhiskPower)
        end
    end
    
    %% All data from trials with no sleep
    allData = [];
    for aa = 1:size(mergedDataFileIDs,1)
        mergedDataFileID = mergedDataFileIDs(aa,:);
        [~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P_Manuscript2020(mergedDataFileID);
        if strcmp(vesselID(1),'V') == false
            load(mergedDataFileID,'-mat')
            trainingDataFileID = trainingDataFileIDs(aa,:);
            load(trainingDataFileID,'-mat')
            strDay = ConvertDate_2P_Manuscript2020(fileDate);
            if sum(strcmp(trainingTable.behavState,'Not Sleep')) > 170
                if isfield(allData,vesselID) == false
                    allData.(vesselID) = [];
                    binWhisking.(vesselID) = [];
                end
                vesselDiam = MergedData.data.vesselDiameter.data;
                normVesselDiam = filtfilt(sos,g,detrend((vesselDiam - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay))./ RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay),'constant'));
                allData.(vesselID) = horzcat(allData.(vesselID),normVesselDiam');
                binWhisk = MergedData.data.binWhiskerAngle;
                [linkedBinarizedWhiskers] = LinkBinaryEvents_2P_Manuscript2020(gt(binWhisk,0),[round(30/3),0]);
                binWhiskingPercent = sum(linkedBinarizedWhiskers)/length(linkedBinarizedWhiskers)*100 ;
                binWhisking.(vesselID) = horzcat(binWhisking.(vesselID),binWhiskingPercent);
            end
        end
    end
    % parameters for mtspectrumc - information available in function
    params.tapers = [3,5];
    params.pad = 1;
    params.Fs = samplingRate;
    params.fpass = [0,1];
    params.trialave = 1;
    params.err = [2,0.05];
    if isempty(allData) == false
        allDataVesselIDs = fieldnames(allData);
        for bb = 1:length(allDataVesselIDs)
            allDataVID = allDataVesselIDs{bb,1};
            [allData_S{bb,1},allData_f{bb,1},allData_sErr{bb,1}] = mtspectrumc_Manuscript2020(allData.(allDataVID),params);
            % save data and figures
            AnalysisResults.(animalID).PowerSpectra.AllData.(allDataVID).whiskingPerc = mean(binWhisking.(allDataVID));
            AnalysisResults.(animalID).PowerSpectra.AllData.(allDataVID).S = allData_S{bb,1};
            AnalysisResults.(animalID).PowerSpectra.AllData.(allDataVID).f = allData_f{bb,1};
            AnalysisResults.(animalID).PowerSpectra.AllData.(allDataVID).sErr = allData_sErr{bb,1};
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                % awake rest summary figures
                allDataPower = figure;
                loglog(allData_f{bb,1},allData_S{bb,1},'k')
                hold on;
                loglog(allData_f{bb,1},allData_sErr{bb,1},'color',colors_Manuscript2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' ' allDataVID ' vessel diameter power during all awake data']);
                set(gca,'Ticklength',[0,0]);
                xlim([0,1])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Vessel Power Spectrum/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(allDataPower,[dirpath animalID '_' allDataVID '_AllDataPowerSpectra']);
                close(allDataPower)
            end
        end
    end
    
    %% Analyze power spectra during periods of NREM sleep
    if isfield(SleepData.(modelType),'NREM') == true
        % pull data from SleepData.mat structure
        nremData = SleepData.(modelType).NREM.data.vesselDiameter.data;
        nremVesselIDs = SleepData.(modelType).NREM.VesselIDs;
        % detrend - data is already lowpass filtered
        for aa = 1:length(nremData)
            nremData{aa,1} = filtfilt(sos,g,detrend(nremData{aa,1}(1:(params.minTime.NREM*samplingRate)),'constant'));
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        nremReshape = zeros(length(nremData{1,1}),length(nremData));
        for bb = 1:length(nremData)
            nremReshape(:,bb) = nremData{bb,1};
        end
        % remove veins from the artery list
        uniqueNREMVesselIDs = unique(nremVesselIDs);
        dd = 1;
        for cc = 1:length(uniqueNREMVesselIDs)
            if strcmp(uniqueNREMVesselIDs{cc,1}(1),'V') == false
                nremArterioleIDs{dd,1} = uniqueNREMVesselIDs{cc,1};
                dd = dd + 1;
            end
        end
        % split the data based on different arteries
        for ee = 1:length(nremArterioleIDs)
            nremArterioleID = nremArterioleIDs{ee,1};
            gg = 1;
            for ff = 1:length(nremVesselIDs)
                if strcmp(nremArterioleID,nremVesselIDs{ff,1}) == true
                    nremArterioleData.(nremArterioleID)(:,gg) = nremReshape(:,ff);
                    gg = gg + 1;
                end
            end
        end
        % parameters for mtspectrumc - information available in function
        params.tapers = [3,5];
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,1];
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the power spectra of the desired signals
        for bb = 1:length(nremArterioleIDs)
            nremVID = nremArterioleIDs{bb,1};
            [nrem_S{bb,1},nrem_f{bb,1},nrem_sErr{bb,1}] = mtspectrumc_Manuscript2020(nremArterioleData.(nremVID),params);
            % save data and figures
            AnalysisResults.(animalID).PowerSpectra.NREM.(nremVID).S = nrem_S{bb,1};
            AnalysisResults.(animalID).PowerSpectra.NREM.(nremVID).f = nrem_f{bb,1};
            AnalysisResults.(animalID).PowerSpectra.NREM.(nremVID).sErr = nrem_sErr{bb,1};
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                % awake rest summary figures
                nremPower = figure;
                loglog(nrem_f{bb,1},nrem_S{bb,1},'k')
                hold on;
                loglog(nrem_f{bb,1},nrem_sErr{bb,1},'color',colors_Manuscript2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' ' nremVID ' vessel diameter power during NREM']);
                set(gca,'Ticklength',[0,0]);
                xlim([0,1])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Vessel Power Spectrum/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(nremPower,[dirpath animalID '_' nremVID '_NREMPowerSpectra']);
                close(nremPower)
            end
        end
    end
    
    %% Analyze power spectra during periods of REM sleep
    if isfield(SleepData.(modelType),'REM') == true
        % pull data from SleepData.mat structure
        remData = SleepData.(modelType).REM.data.vesselDiameter.data;
        remVesselIDs = SleepData.(modelType).REM.VesselIDs;
        % detrend - data is already lowpass filtered
        for aa = 1:length(remData)
            remDataA{aa,1} = filtfilt(sos,g,detrend(remData{aa,1}(1:(params.minTime.REM*samplingRate)),'constant'));
        end
        % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
        remReshape = zeros(length(remDataA{1,1}),length(remDataA));
        for bb = 1:length(remDataA)
            remReshape(:,bb) = remDataA{bb,1};
        end
        % remove veins from the artery list
        uniqueREMVesselIDs = unique(remVesselIDs);
        dd = 1;
        for cc = 1:length(uniqueREMVesselIDs)
            if strcmp(uniqueREMVesselIDs{cc,1}(1),'V') == false
                remArterioleIDs{dd,1} = uniqueREMVesselIDs{cc,1};
                dd = dd + 1;
            end
        end
        % split the data based on different arteries
        for ee = 1:length(remArterioleIDs)
            remArterioleID = remArterioleIDs{ee,1};
            gg = 1;
            for ff = 1:length(remVesselIDs)
                if strcmp(remArterioleID,remVesselIDs{ff,1}) == true
                    remArterioleData.(remArterioleID)(:,gg) = remReshape(:,ff);
                    gg = gg + 1;
                end
            end
        end
        % parameters for mtspectrumc - information available in function
        params.tapers = [3,5];
        params.pad = 1;
        params.Fs = samplingRate;
        params.fpass = [0,1];
        params.trialave = 1;
        params.err = [2,0.05];
        % calculate the power spectra of the desired signals
        for bb = 1:length(remArterioleIDs)
            remVID = remArterioleIDs{bb,1};
            [rem_S{bb,1},rem_f{bb,1},rem_sErr{bb,1}] = mtspectrumc_Manuscript2020(remArterioleData.(remVID),params);
            % save data and figures
            AnalysisResults.(animalID).PowerSpectra.REM.(remVID).S = rem_S{bb,1};
            AnalysisResults.(animalID).PowerSpectra.REM.(remVID).f = rem_f{bb,1};
            AnalysisResults.(animalID).PowerSpectra.REM.(remVID).sErr = rem_sErr{bb,1};
            % save figures if desired
            if strcmp(saveFigs,'y') == true
                % awake rest summary figures
                remPower = figure;
                loglog(rem_f{bb,1},rem_S{bb,1},'k')
                hold on;
                loglog(rem_f{bb,1},rem_sErr{bb,1},'color',colors_Manuscript2020('battleship grey'))
                xlabel('Freq (Hz)');
                ylabel('Power');
                title([animalID  ' ' remVID ' vessel diameter power during REM']);
                set(gca,'Ticklength',[0,0]);
                xlim([0,1])
                axis square
                set(gca,'box','off')
                [pathstr,~,~] = fileparts(cd);
                dirpath = [pathstr '/Figures/Vessel Power Spectrum/'];
                if ~exist(dirpath,'dir')
                    mkdir(dirpath);
                end
                savefig(remPower,[dirpath animalID '_' remVID '_REMPowerSpectra']);
                close(remPower)
            end
        end
    end
end
cd(rootFolder)
save('AnalysisResults.mat','AnalysisResults')
end
