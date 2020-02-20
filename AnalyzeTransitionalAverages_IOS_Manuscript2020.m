function [AnalysisResults] = AnalyzeTransitionalAverages_IOS_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: 
%________________________________________________________________________________________________________________________

%% function parameters
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120'};
modelTypes = {'SVM'};
params.minTime.Rest = 30;   % seconds
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 30;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(IOS_animalIDs,animalID))
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
    [labels,~] = predict(MDL,scoringTable);
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
    for d = 1:length(labels)
        if patchedREMindex(d,1) == 1
            labels{d,1} = 'REM Sleep';
        end
    end
    
    
    
    
    
    
        save('AnalysisResults.mat','AnalysisResults')

    % EMG
    % HR
    % Whisker Angle
    % Spectrogram
    % HbT
    
    
    
    
    
    
    
    
end