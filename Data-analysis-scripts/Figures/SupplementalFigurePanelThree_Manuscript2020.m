function [] = SupplementalFigurePanelThree_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: 
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
%% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    data.oobError(aa,1) = round(AnalysisResults.(animalID).ModelCrossValidation.oobErr,1);
    data.shuffMean(aa,1) = round(mean(AnalysisResults.(animalID).ModelCrossValidation.shuff_oobErr),1);
    data.shuffStD(aa,1) = round(std(AnalysisResults.(animalID).ModelCrossValidation.shuff_oobErr),1);
    data.shuffMax(aa,1) = round(max(AnalysisResults.(animalID).ModelCrossValidation.shuff_oobErr),1);
    data.shuffMin(aa,1) = round(min(AnalysisResults.(animalID).ModelCrossValidation.shuff_oobErr),1);
end
%% load in all the ConfusionData.mat structure
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
%% Supplemental figure panel three
% confusion matrix
confMat = figure;
cm = confusionchart(holdYlabels,holdXlabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
modelAccuracy = round((sum(confVals([1,5,9])/totalScores))*100,1);
cm.Title = {'Supplemental Figure Panel 3 - Turner Manuscript 2020','','[A] Random forest unseen data confusion matrix',['total accuracy: ' num2str(modelAccuracy) ' (%)']};
% save location
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
savefig(confMat,[dirpath 'Supplemental Figure Panel 3 (A)']);
set(confMat,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Supplemental Figure Panel 3 (A)'])
%% table of shuffled data and individual animal accuracy
oobTable = figure;
variableNames = {'oobErr','shuff_oobErr_Mean','shuff_oobErr_StD','shuff_oobErr_Max','shuff_oobErr_Min'};
T = table(data.oobError,data.shuffMean,data.shuffStD,data.shuffMax,data.shuffMin,'RowNames',animalIDs,'VariableNames',variableNames);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
uicontrol('Style','text','Position',[20,100,200,20],'String',{'Supplemental Figure Panel Three (Continued)','to be made in word at a later point'});
savefig(oobTable,[dirpath 'Supplemental Figure Panel 3 (B)']);
set(oobTable,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Supplemental Figure Panel 3 (B)'])
cd(startingDirectory)

end
