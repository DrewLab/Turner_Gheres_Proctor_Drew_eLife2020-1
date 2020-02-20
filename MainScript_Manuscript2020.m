function [] = MainScript_Manuscript2020()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generates KLT's main and supplemental figs for the 2020 sleep paper.
%
%            Scripts used to pre-process the original data are located in the folder "Pre-processing-scripts".
%            Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________

clear; clc;
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
    multiWaitbar_Manuscript2020('Analyzing coherence',0,'Color','K'); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing power spectra',0,'Color','Y'); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing cross correlation',0,'Color','G'); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing evoked responses',0,'Color','P'); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing Pearson''s correlation coefficients',0,'Color','O'); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing behavioral hemodynamics',0,'Color','C'); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing behavioral vessel diameter',0,'Color','R'); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing behavioral heart rate' ,0,'Color','W'); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing behavioral transitions',0,'Color','M'); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing hemodynamic response functions',0,'Color','B'); pause(0.25);
    multiWaitbar_Manuscript2020('Analyzing laser doppler flow',0,'Color','A'); pause(0.25);
    % Run analysis and output a structure containing all the analyzed data.
    [AnalysisResults] = AnalyzeData_Manuscript2020(rootFolder);
    multiWaitbar_Manuscript2020('CloseAll');
end

%% Individual figures can be re-run after the analysis has completed.
% AvgCoherence_Manuscript2020(AnalysisResults)
% AvgPowerSpectra_Manuscript2020(AnalysisResults)
% AvgXCorr_Manuscript2020(AnalysisResults)
% AvgStim_Manuscript2020(AnalysisResults)
% AvgWhisk_Manuscript2020(AnalysisResults)
% AvgCorrCoeff_Manuscript2020(AnalysisResults)
% AvgCBVandHeartRate_Manuscript2020(AnalysisResults)
% AvgResponseFunctionPredictions_Manuscript2020(AnalysisResults)
disp('MainScript Analysis - Complete'); disp(' ')

%% Informational figures with function dependencies for the various analysis and the time per vessel.
% To view individual summary figures, change the value of line 72 to false. You will then be prompted to manually select
% any number of figures (CTL-A for all) inside any of the five folders. You can only do one animal at a time.
% functionNames = {'MainScript_Manuscript2020','StageOneProcessing_Manuscript2020','StageTwoProcessing_Manuscript2020','StageThreeProcessing_Manuscript2020'};
% functionList = {};
% for a = 1:length(functionNames)
%     [functionList] = GetFuncDependencies_Manuscript2020(a,functionNames{1,a},functionNames,functionList);
% end

end

function [AnalysisResults] = AnalyzeData_Manuscript2020(rootFolder)
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120'};
saveFigs = 'y';
if exist('AnalysisResults.mat') == 2
    load('AnalysisResults.mat')
else
    AnalysisResults = [];
end

%% Block [1] Analyze the coherence between bilateral hemispheres (IOS)
runFromStart = 'n';
for a = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,a})) == false || isfield(AnalysisResults.(animalIDs{1,a}),'Coherence') == false || strcmp(runFromStart,'y') == true 
        [AnalysisResults] = AnalyzeCoherence_Manuscript2020(animalIDs{1,a},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing coherence','Value',a/length(animalIDs));
end

%% Block [2] Analyze the power spectra of each single hemisphere (IOS)
runFromStart = 'n';
for b = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,b})) == false || isfield(AnalysisResults.(animalIDs{1,b}),'PowerSpectra') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzePowerSpectrum_Manuscript2020(animalIDs{1,b},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing power spectra','Value',b/length(animalIDs));
end

%% Block [3] Analyze the cross-correlation between local neural activity and hemodynamics (IOS)
runFromStart = 'n';
for c = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,c})) == false || isfield(AnalysisResults.(animalIDs{1,c}),'XCorr') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeXCorr_Manuscript2020(animalIDs{1,c},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing cross correlation','Value',c/length(animalIDs));
end

%% Block [4] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
runFromStart = 'n';
for d = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,d})) == false || isfield(AnalysisResults.(animalIDs{1,d}),'EvokedAvgs') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeEvokedResponses_Manuscript2020(animalIDs{1,d},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing evoked responses','Value',d/length(animalIDs));
end

%% Block [5] Analyze the Pearson's correlation coefficient between neural/hemodynamic signals (IOS)
runFromStart = 'n';
for e = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,e})) == false || isfield(AnalysisResults.(animalIDs{1,e}),'CorrCoeff') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeCorrCoeffs_Manuscript2020(animalIDs{1,e},rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing Pearson''s correlation coefficients','Value',e/length(animalIDs));
end

%% Block [6] Analyze the mean HbT during different behaviors
runFromStart = 'n';
for f = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,f})) == false || isfield(AnalysisResults.(animalIDs{1,f}),'MeanCBV') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeMeanCBV_Manuscript2020(animalIDs{1,f},rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing behavioral hemodynamics','Value',f/length(animalIDs));
end

%% Block [7] Analyze the mean vessel diameter during different behaviors
% runFromStart = 'n';
% for g = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,g})) == false || isfield(AnalysisResults.(animalIDs{1,g}),'MeanCBV') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeMeanCBV_Manuscript2020(animalIDs{1,f},rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing behavioral vessel diameter','Value',f/length(animalIDs));
% end

%% Block [8] Analyze the mean heart rate during different behaviors
runFromStart = 'n';
for h = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,h})) == false || isfield(AnalysisResults.(animalIDs{1,h}),'MeanHR') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeMeanHeartRate_Manuscript2020(animalIDs{1,h},rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing behavioral heart rate','Value',h/length(animalIDs));
end

%% Block [9] Analyzing behavioral transitions
% runFromStart = 'n';
% for i = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,i})) == false || isfield(AnalysisResults.(animalIDs{1,i}),'Transitions') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeTransitionalAverages_IOS_Manuscript2020(animalIDs{1,i},saveFigs,rootFolder,AnalysisResults);
%     multiWaitbar_Manuscript2020('Analyzing behavioral transitions','Value',i/length(animalIDs));
% end

%% Block [10] Analyze the impulse/gamma response functions and calculate prediction accuracy
runFromStart = 'n';
for j = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,j})) == false || isfield(AnalysisResults.(animalIDs{1,j}),'HRFs') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeHRF_Manuscript2020(animalIDs{1,j},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing hemodynamic response functions','Value',j/length(animalIDs));
end

%% Block [11] Analyze mean laser doppler flow during different behaviors
runFromStart = 'n';
for k = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,k})) == false || isfield(AnalysisResults.(animalIDs{1,k}),'LDFlow') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeLaserDoppler_Manuscript2020(animalIDs{1,k},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing laser doppler flow','Value',k/length(animalIDs));
end

%% Fin.
disp('Loading analysis results and generating figures...'); disp(' ')

end

