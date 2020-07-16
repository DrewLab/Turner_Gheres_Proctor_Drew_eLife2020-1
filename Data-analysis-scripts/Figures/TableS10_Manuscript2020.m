function [AnalysisResults] = TableS10_Manuscript2020(rootFolder,saveFigs,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table S10 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% Set-up and process data
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
% cd through each animal's directory and extract the appropriate analysis results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    data.oobError(aa,1) = round(AnalysisResults.(animalID).ModelAccuracy.oobErr,1);
    data.shuffMean(aa,1) = round(mean(AnalysisResults.(animalID).ModelAccuracy.shuff_oobErr),1);
    data.shuffStD(aa,1) = round(std(AnalysisResults.(animalID).ModelAccuracy.shuff_oobErr),1);
end
meanOOB = mean(data.oobError,1);
stdOOB = std(data.oobError,0,1);
shuffMeanOOB = mean(data.shuffMean,1);
shuffStDOOB = std(data.shuffMean,0,1);
%% Table S10
summaryTable = figure('Name','TableS10'); %#ok<*NASGU>
variableNames = {'oobErr','shuff_oobErr_Mean'};
T = table(data.oobError,data.shuffMean,'RowNames',animalIDs,'VariableNames',variableNames);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
uicontrol('Style','text','Position',[700,600,100,150],'String',{'Mean OOBerror (%): ' num2str(meanOOB) ' +/- ' num2str(stdOOB),'Mean Shuffled OOBerror (%): ' num2str(shuffMeanOOB) ' +/- ' num2str(shuffStDOOB)});
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder '\Summary Figures and Structures\MATLAB Analysis Figures\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryTable,[dirpath 'TableS10']);
end

end
