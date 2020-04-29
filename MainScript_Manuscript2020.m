function [] = MainScript_Manuscript2020()
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
% Purpose: Generates KLT's main and supplemental figs for the 2020 sleep paper.
%
% Scripts used to pre-process the original data are located in the folder "Pre-processing-scripts".
% Functions that are used in both the analysis and pre-processing are located in the analysis folder.
%________________________________________________________________________________________________________________________

% T108
% 190822_11_52_51
% 
% T122
% 200218_10_40_55
% 
% T108
% 190822_13_59_30
% 
% T123
% 200224_16_27_59
%
% T126
% RH_200310_12_46_04_011_A3
%
% T115
% RH_191119_15_47_21_024_A2
%
% T115
% RH_191122_16_39_22_011_P1
%
% T116
% RH_191120_11_18_02_009_A2


clear; clc; close all;
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
multiWaitbar_Manuscript2020('Analyzing model cross validation distribution',0,'Color','B'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing sleep probability',0,'Color','W'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing behavioral transitions',0,'Color','B'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing vessel behavioral transitions',0,'Color','B'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing behavioral hemodynamics',0,'Color','W'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing behavioral vessel diameter',0,'Color','B'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing behavioral heart rate' ,0,'Color','W'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing laser doppler flow',0,'Color','B'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing coherence',0,'Color','B'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing power spectra',0,'Color','W'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing vessel power spectra',0,'Color','W'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing Pearson''s correlation coefficients',0,'Color','B'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing cross correlation',0,'Color','B'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing evoked responses',0,'Color','W'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing vessel evoked responses',0,'Color','W'); pause(0.25);
multiWaitbar_Manuscript2020('Analyzing hemodynamic response functions',0,'Color','W'); pause(0.25);
% Run analysis and output a structure containing all the analyzed data.
[AnalysisResults] = AnalyzeData_Manuscript2020(rootFolder);
multiWaitbar_Manuscript2020('CloseAll');

%% main figure panels
% FigurePanelOne_Manuscript2020(rootFolder)
% FigurePanelTwo_Manuscript2020(rootFolder,AnalysisResults)
% FigurePanelThree_Manuscript2020(rootFolder)
% FigurePanelFour_Manuscript2020(rootFolder,AnalysisResults)
% FigurePanelFive_Manuscript2020(rootFolder,AnalysisResults)
% FigurePanelSix_Manuscript2020(rootFolder,AnalysisResults)
% FigurePanelEight_Manuscript2020(rootFolder,AnalysisResults)
% FigurePanelNine_Manuscript2020(rootFolder,AnalysisResults)
%% supplemental figure panels
% SupplementalFigurePanelOne_Manuscript2020(rootFolder,AnalysisResults)
% SupplementalFigurePanelTwo_Manuscript2020(rootFolder)
% SupplementalFigurePanelThree_Manuscript2020(rootFolder,AnalysisResults)
% SupplementalFigurePanelFour_Manuscript2020(rootFolder,AnalysisResults)
% SupplementalFigurePanelFive_Manuscript2020(rootFolder,AnalysisResults)
% SupplementalFigurePanelSix_Manuscript2020(rootFolder,AnalysisResults)
% SupplementalFigurePanelSeven_Manuscript2020(rootFolder,AnalysisResults)

%%  
% AvgConfMatrixAndCrossValidations_IOS_Manuscript2020(rootFolder,AnalysisResults)
% AvgSleepProbability_Manuscript2020(rootFolder,AnalysisResults)
% AvgBehaviorTransitions_Manuscript2020(rootFolder,AnalysisResults)
% AvgStim_Manuscript2020(rootFolder,AnalysisResults)
% AvgWhisk_Manuscript2020(rootFolder,AnalysisResults)
% AvgCorrCoeff_Manuscript2020(rootFolder,AnalysisResults)
% AvgCBVandHeartRate_Manuscript2020(rootFolder,AnalysisResults)
% [SimulationData] = AvgVesselDiameter_Manuscript2020(rootFolder,AnalysisResults);
% [SimulationData] = AvgVesselPowerSpectra_Manuscript2020(rootFolder,AnalysisResults,SimulationData);
% AvgVesselEvokedAvgs_Manuscript2020(rootFolder,AnalysisResults,SimulationData)
% AvgLaserDopplerFlow_Manuscript2020(rootFolder,AnalysisResults)
% AvgCoherence_Manuscript2020(rootFolder,AnalysisResults)
% AvgPowerSpectra_Manuscript2020(rootFolder,AnalysisResults)
% AvgXCorr_Manuscript2020(rootFolder,AnalysisResults)
% AvgResponseFunctionPredictions_Manuscript2020(rootFolder,AnalysisResults)
% PixelDriftExample_Manuscript2020(rootFolder,AnalysisResults)
% CrossCorrelationROIExample_Manuscript2020(rootFolder,AnalysisResults)
% disp('MainScript Analysis - Complete'); disp(' ')
% sendmail('kevinlturnerjr@gmail.com','Manuscript2020 Analysis Complete');

% %% Informational figures with function dependencies for the various analysis and the time per vessel.
% functionNames = {'MainScript_Manuscript2020','StageOneProcessing_IOS_Manuscript2020','StageTwoProcessing_IOS_Manuscript2020','StageThreeProcessing_IOS_Manuscript2020','SleepScoreMainScript_IOS_Manuscript2020',...
%     'StageOneProcessing_2P_Manuscript2020','StageTwoProcessing_2P_Manuscript2020','StageThreeProcessing_2P_Manuscript2020','SleepScoreMainScript_2P_Manuscript2020'};
% functionList = {};
% for a = 1:length(functionNames)
%     [functionList] = GetFuncDependencies_Manuscript2020(a,functionNames{1,a},functionNames,functionList);
% end

end

function [AnalysisResults] = AnalyzeData_Manuscript2020(rootFolder)
% IOS animal IDs
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
% two photon animal IDs
animalIDs2 = {'T115','T116','T117','T118','T125','T126'};
saveFigs = 'y';
if exist('AnalysisResults.mat') == 2
    load('AnalysisResults.mat')
else
    AnalysisResults = [];
end

% %% Block [0] Analyze the cross validation distribution of 100 iterations of real data and shuffled data
% runFromStart = 'n';
% for aa = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,aa})) == false || isfield(AnalysisResults.(animalIDs{1,aa}),'ModelCrossValidation') == false || strcmp(runFromStart,'y') == true 
%         [AnalysisResults] = AnalyzeModelAccuracy_Manuscript2020(animalIDs{1,aa},saveFigs,rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing model cross validation distribution','Value',aa/length(animalIDs));
% end
% 
% %% Block [1] Analyze the probability of an animal being awake or asleep based on duration of the trial
% runFromStart = 'n';
% for bb = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,bb})) == false || isfield(AnalysisResults.(animalIDs{1,bb}),'SleepProbability') == false || strcmp(runFromStart,'y') == true 
%         [AnalysisResults] = AnalyzeAwakeProbability_Manuscript2020(animalIDs{1,bb},saveFigs,rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing sleep probability','Value',bb/length(animalIDs));
% end

%% Block [2] Analyze the coherence between bilateral hemispheres (IOS)
runFromStart = 'y';
for cc = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,cc})) == false || isfield(AnalysisResults.(animalIDs{1,cc}),'Coherence') == false || strcmp(runFromStart,'y') == true 
        [AnalysisResults] = AnalyzeCoherence_Manuscript2020(animalIDs{1,cc},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing coherence','Value',cc/length(animalIDs));
end

%% Block [3] Analyze the power spectra of each single hemisphere (IOS)
runFromStart = 'y';
for dd = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,dd})) == false || isfield(AnalysisResults.(animalIDs{1,dd}),'PowerSpectra') == false || strcmp(runFromStart,'y') == true
            [AnalysisResults] = AnalyzePowerSpectrum_Manuscript2020(animalIDs{1,dd},saveFigs,rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing power spectra','Value',dd/length(animalIDs));
end

% %% Block [4] Analyze the cross-correlation between local neural activity and hemodynamics (IOS)
% runFromStart = 'n';
% for ee = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,ee})) == false || isfield(AnalysisResults.(animalIDs{1,ee}),'XCorr') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeXCorr_Manuscript2020(animalIDs{1,ee},saveFigs,rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing cross correlation','Value',ee/length(animalIDs));
% end
% 
% %% Block [5] Analyze the stimulus-evoked and whisking-evoked neural/hemodynamic responses (IOS)
% runFromStart = 'n';
% for ff = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,ff})) == false || isfield(AnalysisResults.(animalIDs{1,ff}),'EvokedAvgs') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeEvokedResponses_Manuscript2020(animalIDs{1,ff},saveFigs,rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing evoked responses','Value',ff/length(animalIDs));
% end
% 
% %% Block [6] Analyze the Pearson's correlation coefficient between neural/hemodynamic signals (IOS)
% runFromStart = 'n';
% for gg = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,gg})) == false || isfield(AnalysisResults.(animalIDs{1,gg}),'CorrCoeff') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeCorrCoeffs_Manuscript2020(animalIDs{1,gg},rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing Pearson''s correlation coefficients','Value',gg/length(animalIDs));
% end
% 
% %% Block [7] Analyze the mean HbT during different behaviors
% runFromStart = 'n';
% for hh = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,hh})) == false || isfield(AnalysisResults.(animalIDs{1,hh}),'MeanCBV') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeMeanCBV_Manuscript2020(animalIDs{1,hh},rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing behavioral hemodynamics','Value',hh/length(animalIDs));
% end
% 
%% Block [8] Analyze the mean vessel diameter during different behaviors
runFromStart = 'n';
for ii = 1:length(animalIDs2)
    if isfield(AnalysisResults,(animalIDs2{1,ii})) == false || isfield(AnalysisResults.(animalIDs2{1,ii}),'MeanVesselDiameter') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeMeanVesselDiameter_Manuscript2020(animalIDs2{1,ii},rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing behavioral vessel diameter','Value',ii/length(animalIDs2));
end
% 
% %% Block [9] Analyze the mean heart rate during different behaviors
% runFromStart = 'n';
% for jj = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,jj})) == false || isfield(AnalysisResults.(animalIDs{1,jj}),'MeanHR') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeMeanHeartRate_Manuscript2020(animalIDs{1,jj},rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing behavioral heart rate','Value',jj/length(animalIDs));
% end
% 
% %% Block [10] Analyze behavioral transitions
% runFromStart = 'n';
% for kk = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,kk})) == false || isfield(AnalysisResults.(animalIDs{1,kk}),'Transitions') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeTransitionalAverages_Manuscript2020(animalIDs{1,kk},saveFigs,rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing behavioral transitions','Value',kk/length(animalIDs));
% end
% 
% %% Block [11] Analyze the impulse/gamma response functions and calculate prediction accuracy
% runFromStart = 'n';
% for ll = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,ll})) == false || isfield(AnalysisResults.(animalIDs{1,ll}),'HRFs') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeHRF_Manuscript2020(animalIDs{1,ll},saveFigs,rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing hemodynamic response functions','Value',ll/length(animalIDs));
% end
% 
% %% Block [12] Analyze mean laser doppler flow during different behaviors
% runFromStart = 'n';
% for mm = 1:length(animalIDs)
%     if isfield(AnalysisResults,(animalIDs{1,mm})) == false || isfield(AnalysisResults.(animalIDs{1,mm}),'LDFlow') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeLaserDoppler_Manuscript2020(animalIDs{1,mm},rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing laser doppler flow','Value',mm/length(animalIDs));
% end
% 
% %% Block [13] Analyze vessel power and evoked responses during different behaviors
% runFromStart = 'n';
% for nn = 1:length(animalIDs2)
%     if isfield(AnalysisResults,(animalIDs2{1,nn})) == false || isfield(AnalysisResults.(animalIDs2{1,nn}),'PowerSpectra') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeVesselPowerSpectrum_Manuscript2020(animalIDs2{1,nn},saveFigs,rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing vessel power spectra','Value',nn/length(animalIDs2));
% end
%% Block [14] Anayze behavioral transitions
% runFromStart = 'n';
% for nn = 1:length(animalIDs2)
%     if isfield(AnalysisResults,(animalIDs2{1,nn})) == false || isfield(AnalysisResults.(animalIDs2{1,nn}),'Transitions') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeVesselTransitionalAverages_Manuscript2020(animalIDs2{1,nn},rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing vessel behavioral transitions','Value',nn/length(animalIDs2));
% end
% %% Block [15] Analyze the whisking-evoked arteriole response (2P)
% runFromStart = 'n';
% for nn = 1:length(animalIDs2)
%     if isfield(AnalysisResults,(animalIDs2{1,nn})) == false || isfield(AnalysisResults.(animalIDs2{1,nn}),'EvokedAvgs') == false || strcmp(runFromStart,'y') == true
%         [AnalysisResults] = AnalyzeVesselEvokedResponses_Manuscript2020(animalIDs2{1,nn},saveFigs,rootFolder,AnalysisResults);
%     end
%     multiWaitbar_Manuscript2020('Analyzing vessel evoked responses','Value',nn/length(animalIDs2));
% end
%% Block [15] Analyze the whisking-evoked arteriole response (2P)
runFromStart = 'n';
for nn = 1:length(animalIDs)
    if isfield(AnalysisResults,(animalIDs{1,nn})) == false || isfield(AnalysisResults.(animalIDs{1,nn}),'BehaviorDistributions') == false || strcmp(runFromStart,'y') == true
        [AnalysisResults] = AnalyzeBehavioralDistributions_Manuscript2020(animalIDs{1,nn},rootFolder,AnalysisResults);
    end
    multiWaitbar_Manuscript2020('Analyzing behavioral distributions','Value',nn/length(animalIDs));
end
%% Fin.
disp('Loading analysis results and generating figures...'); disp(' ')

end
