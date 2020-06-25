function [AnalysisResults] = Table2_Manuscript2020(rootFolder,AnalysisResults) %#ok<INUSL>
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table 2 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% set-up and process data
columnNames = AnalysisResults.Coherr.columnNames;
rowNames = {'Gamma_C01_meanStD','Gamma_C01_pVal','HbT_C01_meanStD','HbT_C01_pVal'...
    'Gamma_C001_meanStD','Gamma_C001_pVal','HbT_C001_meanStD','HbT_C001_pVal'};
T(1,:) = cell2table(AnalysisResults.Coherr.gammaBandPower.meanStD01);
T(2,:) = cell2table(AnalysisResults.Coherr.gammaBandPower.p01);
T(3,:) = cell2table(AnalysisResults.Coherr.CBV_HbT.meanStD01);
T(4,:) = cell2table(AnalysisResults.Coherr.CBV_HbT.p01);
T(5,:) = cell2table(AnalysisResults.Coherr.gammaBandPower.meanStD001);
T(6,:) = cell2table(AnalysisResults.Coherr.gammaBandPower.p001);
T(7,:) = cell2table(AnalysisResults.Coherr.CBV_HbT.meanStD001);
T(8,:) = cell2table(AnalysisResults.Coherr.CBV_HbT.p001);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% Table 2
summaryTable = figure('Name','Table2'); %#ok<*NASGU>
sgtitle('Table 2 Turner Manuscript 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
%% save figure(s)
% dirpath = [rootFolder '\Summary Figures and Structures\'];
% if ~exist(dirpath,'dir')
%     mkdir(dirpath);
% end
% savefig(summaryTable,[dirpath 'Table2']);

end
