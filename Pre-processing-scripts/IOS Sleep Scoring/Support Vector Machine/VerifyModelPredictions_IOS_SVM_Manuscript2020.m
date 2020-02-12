function [modelAccuracy] = VerifyModelPredictions_IOS_SVM_Manuscript2020(animalID,RestingBaselines,startingDirectory,validationDirectory,saveFigs,baselineType,modelAccuracy)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https:\\github.com\KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose:
%________________________________________________________________________________________________________________________

% cd to model location and load it
modelDirectory = [startingDirectory '\Support Files\'];
cd(modelDirectory)
load('IOS_SVM_SleepScoringModel.mat')
% cd back to validation data set
cd(startingDirectory)
cd(validationDirectory)
% get model data file names
modelDataFileStruct = dir('*_ModelData.mat');
modelDataFiles = {modelDataFileStruct.name}';
modelDataFileIDs = char(modelDataFiles);
% predict behavior events of the model files
PredictBehaviorEvents_IOS_SVM_Manuscript2020(modelDataFileIDs,SVMModel)
% get scoring results file name and load it into workspace
scoringResultsDataFileStruct = dir('*_SleepScoringResults.mat');
scoringResultsDataFiles = {scoringResultsDataFileStruct.name}';
scoringResultsDataFileID = char(scoringResultsDataFiles);
load(scoringResultsDataFileID)
% get proc data file names
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% get the filenames of all the training data files to compare the manual scores to the model results
trainingDataFileStruct = dir('*_TrainingData.mat');
trainingDataFiles = {trainingDataFileStruct.name}';
trainingDataFileIDs = char(trainingDataFiles);
for b = 1:size(trainingDataFileIDs,1)
    modelDataFileID = modelDataFileIDs(b,:);
    procDataFileID = procDataFileIDs(b,:);
    trainingDataFileID = trainingDataFileIDs(b,:);
    load(trainingDataFileID)
    c = 1;
    clear fileSVMScores
    % pull out each file's manual scores
    for d = 1:length(SVMResults.fileIDs)
        if strcmp(modelDataFileID,SVMResults.fileIDs{d,1}) == true
            fileSVMScores{c,1} = SVMResults.labels{d,1}; %#ok<*AGROW>
            c = c + 1;
        end
    end
    fileManualScores = trainingTable.behavState;
    if b == 1
        allSVMScores = fileSVMScores;
        allManualScores = fileManualScores;
    else
        allSVMScores = vertcat(allSVMScores,fileSVMScores);
        allManualScores = vertcat(allManualScores,fileManualScores);
    end    
    % show time indece comparison of each score if desired
    if strcmp(saveFigs,'y') == true
        % Not sleep SVM index
        svmNotSleepInds = double(strcmp(fileSVMScores,'Not Sleep'));
        svmNotSleepInds(svmNotSleepInds~=1) = NaN;
        % NREM SVM index
        svmNREMSleepInds = double(strcmp(fileSVMScores,'NREM Sleep'));
        svmNREMSleepInds(svmNREMSleepInds~=1) = NaN;
        % REM SVM index
        svmREMSleepInds = double(strcmp(fileSVMScores,'REM Sleep'));
        svmREMSleepInds(svmREMSleepInds~=1) = NaN;
        % Not sleep manual index
        manualNotSleepInds = double(strcmp(fileManualScores,'Not Sleep'));
        manualNotSleepInds(manualNotSleepInds~=1) = NaN;
        % NREM manual index
        manualNREMSleepInds = double(strcmp(fileManualScores,'NREM Sleep'));
        manualNREMSleepInds(manualNREMSleepInds~=1) = NaN;
        % REM manual index
        manualREMSleepInds = double(strcmp(fileManualScores,'REM Sleep'));
        manualREMSleepInds(manualREMSleepInds~=1) = NaN;
        xInds = 1:5:900;
        [figHandle] = GenerateSingleFigures_IOS_SVM_Manuscript2020(procDataFileID,RestingBaselines,baselineType);
        % show results on figure
        subplot(6,1,1)
        yyaxis left
        ylimits1 = ylim;
        yMax1 = ylimits1(2);
        yInds_svmNotSleep = svmNotSleepInds*yMax1*1.2;
        yInds_svmNREMSleep = svmNREMSleepInds*yMax1*1.2;
        yInds_svmREMSleep = svmREMSleepInds*yMax1*1.2;
        yInds_manualNotSleep = manualNotSleepInds*yMax1*1.4;
        yInds_manualNREMSleep = manualNREMSleepInds*yMax1*1.4;
        yInds_manualREMSleep = manualREMSleepInds*yMax1*1.4;
        hold on
        scatter(xInds,yInds_svmNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_svmNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_svmREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');
        scatter(xInds,yInds_manualNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_manualNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_manualREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');
        subplot(6,1,2)
        yyaxis left
        ylimits2 = ylim;
        yMax2 = ylimits2(2);
        yInds_svmNotSleep = svmNotSleepInds*yMax2*1.1;
        yInds_svmNREMSleep = svmNREMSleepInds*yMax2*1.1;
        yInds_svmREMSleep = svmREMSleepInds*yMax2*1.1;
        yInds_manualNotSleep = manualNotSleepInds*yMax2*1.2;
        yInds_manualNREMSleep = manualNREMSleepInds*yMax2*1.2;
        yInds_manualREMSleep = manualREMSleepInds*yMax2*1.2;
        hold on
        scatter(xInds,yInds_svmNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_svmNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_svmREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');
        scatter(xInds,yInds_manualNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_manualNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_manualREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');
        ylim([-20,80])
        yyaxis right
        ylim([5,15])
        subplot(6,1,3)
        yyaxis left
        ylimits3 = ylim;
        yMax3 = ylimits3(2);
        yInds_svmNotSleep = svmNotSleepInds*yMax3*1.2;
        yInds_svmNREMSleep = svmNREMSleepInds*yMax3*1.2;
        yInds_svmREMSleep = svmREMSleepInds*yMax3*1.2;
        yInds_manualNotSleep = manualNotSleepInds*yMax3*1.4;
        yInds_manualNREMSleep = manualNREMSleepInds*yMax3*1.4;
        yInds_manualREMSleep = manualREMSleepInds*yMax3*1.4;
        hold on
        scatter(xInds,yInds_svmNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_svmNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_svmREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');
        scatter(xInds,yInds_manualNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_manualNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_manualREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');       
        subplot(6,1,4)
        yyaxis left
        yMax4 = 100;
        yInds_svmNotSleep = svmNotSleepInds*yMax4*1.2;
        yInds_svmNREMSleep = svmNREMSleepInds*yMax4*1.2;
        yInds_svmREMSleep = svmREMSleepInds*yMax4*1.2;
        yInds_manualNotSleep = manualNotSleepInds*yMax4*1.6;
        yInds_manualNREMSleep = manualNREMSleepInds*yMax4*1.6;
        yInds_manualREMSleep = manualREMSleepInds*yMax4*1.6;
        hold on
        scatter(xInds,yInds_svmNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_svmNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_svmREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');
        scatter(xInds,yInds_manualNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_manualNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_manualREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');
        subplot(6,1,5)
        yyaxis left
        yMax5 = 100;
        yInds_svmNotSleep = svmNotSleepInds*yMax5*1.2;
        yInds_svmNREMSleep = svmNREMSleepInds*yMax5*1.2;
        yInds_svmREMSleep = svmREMSleepInds*yMax5*1.2;
        yInds_manualNotSleep = manualNotSleepInds*yMax5*1.6;
        yInds_manualNREMSleep = manualNREMSleepInds*yMax5*1.6;
        yInds_manualREMSleep = manualREMSleepInds*yMax5*1.6;
        hold on
        scatter(xInds,yInds_svmNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_svmNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_svmREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');
        scatter(xInds,yInds_manualNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_manualNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_manualREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');       
        subplot(6,1,6)
        yyaxis left
        yMax6 = 100;
        yInds_svmNotSleep = svmNotSleepInds*yMax6*1.2;
        yInds_svmNREMSleep = svmNREMSleepInds*yMax6*1.2;
        yInds_svmREMSleep = svmREMSleepInds*yMax6*1.2;
        yInds_manualNotSleep = manualNotSleepInds*yMax6*1.6;
        yInds_manualNREMSleep = manualNREMSleepInds*yMax6*1.6;
        yInds_manualREMSleep = manualREMSleepInds*yMax6*1.6;
        hold on
        scatter(xInds,yInds_svmNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_svmNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_svmREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');
        scatter(xInds,yInds_manualNotSleep,'MarkerEdgeColor','k','MarkerFaceColor','k');
        scatter(xInds,yInds_manualNREMSleep,'MarkerEdgeColor','b','MarkerFaceColor','b');
        scatter(xInds,yInds_manualREMSleep,'MarkerEdgeColor','r','MarkerFaceColor','r');
        % save figures to proper directory
        dirpath = [startingDirectory '\' animalID '\Figures\SVM Validation\'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(figHandle,[dirpath procDataFileID(1:end-12) 'SVM_ModelValidation']);
        close(figHandle)
    end
end
% confusion matrix
confMat = figure;
cm = confusionchart(allManualScores,allSVMScores);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
cm.Title = 'T109 SVM Confusion Matrix';
dirpath = [startingDirectory '\' animalID '\Figures\SVM Validation\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(confMat,[dirpath animalID '_IOS_SVM_ConfusionMatrix']);
close(confMat)
% "Kevin accuracy" for rough estimation of success 
totalScores = length(allSVMScores);
positiveScores = sum(strcmp(allSVMScores,allManualScores));
modelAccuracy.(animalID) = (positiveScores/totalScores)*100;
disp([animalID ' SVM model prediction success: ' num2str((positiveScores/totalScores)*100) '%']); disp(' ')

end
