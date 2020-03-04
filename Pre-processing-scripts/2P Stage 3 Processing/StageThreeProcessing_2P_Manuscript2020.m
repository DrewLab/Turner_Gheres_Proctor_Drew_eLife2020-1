%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: 1) Categorize behavioral (rest,whisk,stim) data using previously processed data structures, add 'flags'  
%            2) Create a temporary RestData structure that contains periods of rest - use this for initial figures
%            3) Analyze neural data and create different spectrograms for each file's electrodes
%            4) Uses periods when animal is not being stimulated or moving to establish an initial baseline
%            5) Manually select awake files for a slightly different baseline not based on hard time vals
%            6) Use the best baseline to convert reflectance changes to total hemoglobin
%            7) Re-create the RestData structure now that we can deltaHbT
%            8) Create an EventData structure looking at the different data types after whisking or stimulation
%            9) Apply the resting baseline to each data type to create a percentage change 
%            10) Use the time indeces of the resting baseline file to apply a percentage change to the spectrograms
%            11) Use the time indeces of the resting baseline file to create a reflectance pixel-based baseline
%            12) Generate a summary figure for all of the analyzed and processed data
%________________________________________________________________________________________________________________________

%% BLOCK PURPOSE: [0] Load the script's necessary variables and data structures.
% Clear the workspace variables and command window.
clc;
clear;
disp('Analyzing Block [0] Preparing the workspace and loading variables.'); disp(' ')
% Character list of all MergedData files
mergedDirectory = dir('*_MergedData.mat');
mergedDataFiles = {mergedDirectory.name}';
mergedDataFileIDs = char(mergedDataFiles);
[animalID,~,~,~,~,~] = GetFileInfo2_2P_Manuscript2020(mergedDataFileIDs(1,:));
dataTypes = {'vesselDiameter','corticalNeural','hippocampalNeural','EMG'};
neuralDataTypes = {'corticalNeural','hippocampalNeural'};
specNeuralDataTypes = {'rawCorticalNeural','rawHippocampalNeural'};

%% BLOCK PURPOSE: [1] Categorize data 
disp('Analyzing Block [1] Categorizing data.'); disp(' ')
for a = 1:size(mergedDataFileIDs,1)
    mergedDataFileID = mergedDataFileIDs(a,:);
    disp(['Analyzing file ' num2str(a) ' of ' num2str(size(mergedDataFileIDs,1)) '...']); disp(' ')
    CategorizeData_2P_Manuscript2020(mergedDataFileID)
end

%% BLOCK PURPOSE: [2] Create RestData data structure.
disp('Analyzing Block [2] Creating RestData struct for vessels and neural data.'); disp(' ')
[RestData] = ExtractRestingData_2P_Manuscript2020(mergedDataFileIDs,dataTypes);
    
%% BLOCK PURPOSE: [3] Analyze the spectrogram for each session.
disp('Analyzing Block [3] Analyzing the spectrogram for each file and normalizing by the resting baseline.'); disp(' ')
CreateTrialSpectrograms_2P_Manuscript2020(mergedDataFileIDs,specNeuralDataTypes);

%% BLOCK PURPOSE: [4] Create EventData data structure.
disp('Analyzing Block [3] Creating EventData struct for vessels and neural data.'); disp(' ')
[EventData] = ExtractEventTriggeredData_2P_Manuscript2020(mergedDataFileIDs,dataTypes);

%% BLOCK PURPOSE: [5] Create Baselines data structure
disp('Analyzing Block [4] Create Baselines struct for CBV and neural data.'); disp(' ')
baselineType = 'setDuration';
trialDuration_sec = 900;
targetMinutes = 30;
[RestingBaselines] = CalculateRestingBaselines_2P_Manuscript2020(animalID,targetMinutes,trialDuration_sec,RestData);
% Find spectrogram baselines for each day
specDirectory = dir('*_SpecData.mat');
specDataFiles = {specDirectory.name}';
specDataFileIDs = char(specDataFiles);
[RestingBaselines] = CalculateSpectrogramBaselines_2P_Manuscript2020(animalID,neuralDataTypes,trialDuration_sec,specDataFileIDs,RestingBaselines,baselineType);
% Normalize spectrogram by baseline
NormalizeSpectrograms_2P_Manuscript2020(specDataFileIDs,neuralDataTypes,RestingBaselines);

% 
% 
% 
% %% BLOCK PURPOSE [5] Analyze the spectrogram for each session.
% disp('Analyzing Block [5] Analyzing the spectrogram for each file and normalizing by the resting baseline.'); disp(' ')
% neuralDataTypes = {'corticalNeural','hippocampalNeural'};
% CreateTrialSpectrograms_2P_Manuscript2020(mergedDataFileIDs,neuralDataTypes);
% 
% % Find spectrogram baselines for each day
% specDirectory = dir('*_SpecData.mat');
% specDataFiles = {specDirectory.name}';
% specDataFiles = char(specDataFiles);
% [RestingBaselines] = CalculateSpectrogramBaselines_2P(animalID,trialDuration_Sec,specDataFiles,RestingBaselines,neuralDataTypes);
% 
% % Normalize spectrogram by baseline
% NormalizeSpectrograms_2P(specDataFiles,RestingBaselines,neuralDataTypes);
% 
% %% BLOCK PURPOSE [6]
% saveFigs = 'y';
% GenerateSingleFigures_2P(mergedDataFileIDs,RestingBaselines,saveFigs)
% 
% disp('Two Photon Stage Three Processing - Complete.'); disp(' ')
