function [AnalysisResults] = Table1_Manuscript2020(rootFolder,AnalysisResults) %#ok<INUSL>
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table 1 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% set-up and process data
columnNames = AnalysisResults.PSD.columnNames;
rowNames = {'Gamma_S01_meanStD','Gamma_S01_pVal','HbT_S01_meanStD','HbT_S01_pVal','TwoP_S01_meanStD','TwoP_S01_pVal'...
    'Gamma_S001_meanStD','Gamma_S001_pVal','HbT_S001_meanStD','HbT_S001_pVal','TwoP_S001_meanStD','TwoP_S001_pVal'};
T(1,:) = cell2table(AnalysisResults.PSD.gammaBandPower.meanStD01);
T(2,:) = cell2table(AnalysisResults.PSD.gammaBandPower.p01);
T(3,:) = cell2table(AnalysisResults.PSD.CBV_HbT.meanStD01);
T(4,:) = cell2table(AnalysisResults.PSD.CBV_HbT.p01);
T(5,:) = cell2table(AnalysisResults.PSD.TwoP.meanStD01);
T(6,:) = cell2table(AnalysisResults.PSD.TwoP.p01);
T(7,:) = cell2table(AnalysisResults.PSD.gammaBandPower.meanStD001);
T(8,:) = cell2table(AnalysisResults.PSD.gammaBandPower.p001);
T(9,:) = cell2table(AnalysisResults.PSD.CBV_HbT.meanStD001);
T(10,:) = cell2table(AnalysisResults.PSD.CBV_HbT.p001);
T(11,:) = cell2table(AnalysisResults.PSD.TwoP.meanStD001);
T(12,:) = cell2table(AnalysisResults.PSD.TwoP.p001);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% Table 1
summaryTable = figure('Name','Table1'); %#ok<*NASGU>
sgtitle('Table 1 Turner Manuscript 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
%% save figure(s)
% dirpath = [rootFolder '\Summary Figures and Structures\'];
% if ~exist(dirpath,'dir')
%     mkdir(dirpath);
% end
% savefig(summaryTable,[dirpath 'Table1']);

end
