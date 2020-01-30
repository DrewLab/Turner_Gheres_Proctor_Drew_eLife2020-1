function [] = MainScript_Manuscript2020()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Generates KLT's main and supplemental figs for the 2020 sleep paper. 
%
%            Scripts used to pre-process the original data are located in the folder "Pre-processing-scripts".
%            Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________

clear; clc; close all
%% Make sure the current directory is 'TurnerFigs-Manuscript2020' and that the code repository is present.
currentFolder = pwd;
addpath(genpath(currentFolder));
fileparts = strsplit(currentFolder,filesep);
if ismac
    rootFolder = fullfile(filesep,fileparts{1:end});
else
    rootFolder = fullfile(fileparts{1:end});
end
% Add root folder to Matlab's working directory.
addpath(genpath(rootFolder))

%% Run the data analysis. The progress bars will show the analysis progress.
dataSummary = dir('ComparisonData.mat');
% If the analysis structure has already been created, load it and skip the analysis.
if ~isempty(dataSummary)
    load(dataSummary.name);
    disp('Loading analysis results and generating figures...'); disp(' ')
else
    multiWaitbar_Manuscript2020('Analyzing coherence',0,'Color',[0.720000 0.530000 0.040000]); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing power spectra',0,'Color',[0.720000 0.530000 0.040000]); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing evoked responses',0,'Color',[0.720000 0.530000 0.040000]); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing cross correlation',0,'Color',[0.720000 0.530000 0.040000]); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing Pearson''s correlation coefficients',0,'Color',[0.720000 0.530000 0.040000]); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing behavioral hemodynamics',0,'Color',[0.720000 0.530000 0.040000]); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing behavioral heart rate',0,'Color',[0.720000 0.530000 0.040000]); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing hemodynamic response functions',0,'Color',[0.720000 0.530000 0.040000]); pause(0.25);
    % Run analysis and output a structure containing all the analyzed data.
    [ComparisonData] = AnalyzeData_Manuscript2020(rootFolder);
    multiWaitbar_Manuscript2020('CloseAll');
end

%% Informational figures with function dependencies for the various analysis and the time per vessel.
% To view individual summary figures, change the value of line 72 to false. You will then be prompted to manually select
% any number of figures (CTL-A for all) inside any of the five folders. You can only do one animal at a time.
% functionNames = {'MainScript_Manuscript2020','StageOneProcessing_Manuscript2020','StageTwoProcessing_Manuscript2020','StageThreeProcessing_Manuscript2020'};
% functionList = {};
% for a = 1:length(functionNames)
%     [functionList] = GetFuncDependencies_Manuscript2020(a,functionNames{1,a},functionNames,functionList);
% end

%% Individual figures can be re-run after the analysis has completed.
AvgCoherence_Manuscript2020
AvgPowerSpectra_Manuscript2020
AvgXCorr_Manuscript2020
AvgStim_Manuscript2020
AvgWhisk_Manuscript2020
AvgCorrCoeff_Manuscript2020
AvgCBVandHeartRate_Manuscript2020
AvgResponseFunctionPredictions_Manuscript2020
disp('MainScript Analysis - Complete'); disp(' ')
end

function [AnalysisResults] = AnalyzeData_Manuscript2020(rootFolder)
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111'};
AnalysisResults = [];   % pre-allocate the results structure as empty

%% BLOCK PURPOSE: [1] Analyze the coherence between bilateral hemispheres (IOS)
for a = 1:length(animalIDs)
    [AnalysisResults] = AnalyzeCoherence_Manuscript2020(animalIDs{1,a},rootFolder,AnalysisResults);
    multiWaitbar_Manuscript2020('Analyzing coherence','Value',a/length(animalIDs));
end

%% BLOCK PURPOSE: [2] Analyze the power spectra of each single hemisphere (IOS)
for b = 1:length(animalIDs)
    [AnalysisResults] = AnalyzePowerSpectrum_Manuscript2020(animalIDs{1,b},rootFolder,AnalysisResults);
    multiWaitbar_Manuscript2020('Analyzing power spectra','Value',b/length(animalIDs));
end

%% BLOCK PURPOSE: [3] Analyze the cross-correlation between local neural activity and hemodynamics (IOS)
for b = 1:length(animalIDs)
    [AnalysisResults] = AnalyzeXCorr_Manuscript2020(animalIDs{1,c},AnalysisResults);
    multiWaitbar_Manuscript2020('Analyzing cross correlation','Value',c/length(animalIDs));
end

% %% BLOCK PURPOSE: [4] 
% for b = 1:length(animalIDs)
%     [AnalysisResults] = AnalyzeXCorr_Manuscript2020(animalIDs{1,b},AnalysisResults);
%     multiWaitbar_Manuscript2020('Analyzing cross correlation','Value',b/length(animalIDs));
% end

% %% BLOCK PURPOSE: [5] 
% for b = 1:length(animalIDs)
%     [AnalysisResults] = AnalyzeXCorr_Manuscript2020(animalIDs{1,b},AnalysisResults);
%     multiWaitbar_Manuscript2020('Analyzing cross correlation','Value',b/length(animalIDs));
% end

% %% BLOCK PURPOSE: [6] 
% for b = 1:length(animalIDs)
%     [AnalysisResults] = AnalyzeXCorr_Manuscript2020(animalIDs{1,b},AnalysisResults);
%     multiWaitbar_Manuscript2020('Analyzing cross correlation','Value',b/length(animalIDs));
% end

% %% BLOCK PURPOSE: [7] 
% for b = 1:length(animalIDs)
%     [AnalysisResults] = AnalyzeXCorr_Manuscript2020(animalIDs{1,b},AnalysisResults);
%     multiWaitbar_Manuscript2020('Analyzing cross correlation','Value',b/length(animalIDs));
% end

% %% BLOCK PURPOSE: [8] 
% for b = 1:length(animalIDs)
%     [AnalysisResults] = AnalyzeXCorr_Manuscript2020(animalIDs{1,b},AnalysisResults);
%     multiWaitbar_Manuscript2020('Analyzing cross correlation','Value',b/length(animalIDs));
% end

answer = questdlg('Would you like to save the analysis results structure?','','yes','no','yes');
if strcmp(answer,'yes')
    save('ComparisonData.mat','ComparisonData')
end
disp('Loading analysis results and generating figures...'); disp(' ')

end
