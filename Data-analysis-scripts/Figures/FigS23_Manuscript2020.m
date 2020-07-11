function [AnalysisResults] = FigS23_Manuscript2020(rootFolder,saveFigs,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel S23 for Turner_Kederasetti_Gheres_Proctorostanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% set-up and process data
% load in all the ConfusionData.mat structure
startingDirectory = cd;
confusionDataDirectory = [startingDirectory '\Summary Figures and Structures\Confusion Matricies\'];
cd(confusionDataDirectory)
load('ConfusionData.mat','-mat')
% pull out confusion matrix values
modelName = 'RF';
holdYlabels = [];
holdXlabels = [];
for bb = 1:length(ConfusionData.(modelName).testYlabels)
    holdYlabels = vertcat(holdYlabels,ConfusionData.(modelName).testYlabels{bb,1}); %#ok<*AGROW>
    holdXlabels = vertcat(holdXlabels,ConfusionData.(modelName).testXlabels{bb,1});
end
%% Figure panel S23
confMat = figure('Name','FigS23 (a)'); %#ok<*NASGU>
%% confusion matrix
cm = confusionchart(holdYlabels,holdXlabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Supplemental Figure S23 Turner Manuscript 2020','','[S23a] Random forest unseen data confusion matrix',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
%% save location
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder '\Summary Figures and Structures\'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(confMat,[dirpath 'FigS23']);
    set(confMat,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'FigS23'])
end

end
