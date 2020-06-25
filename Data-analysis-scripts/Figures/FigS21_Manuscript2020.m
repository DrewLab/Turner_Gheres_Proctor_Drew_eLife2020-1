function [AnalysisResults] = FigS21_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel S21 for Turner_Kederasetti_Gheres_Proctorostanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% Set-up and process data for Fig S21 (a)
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
%% Figure panel S15
confMat = figure('Name','FigS21 (a)');
%% confusion matrix
cm = confusionchart(holdYlabels,holdXlabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Supplemental Figure S21 Turner Manuscript 2020','','[S21a] Random forest unseen data confusion matrix',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
%% save location
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
savefig(confMat,[dirpath 'FigS21']);
set(confMat,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'FigS21'])

end
