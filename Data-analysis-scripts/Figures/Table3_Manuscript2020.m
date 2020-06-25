function [AnalysisResults] = Table3_Manuscript2020(rootFolder,AnalysisResults) %#ok<INUSL>
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table 3 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% set-up and process data
columnNames = AnalysisResults.CorrCoef.columnNames;
rowNames = {'Gamma_R_meanStD','Gamma_R_pVal','HbT_R_meanStD','HbT_R_pVal'};
T(1,:) = cell2table(AnalysisResults.CorrCoef.gammaBandPower.meanStD);
T(2,:) = cell2table(AnalysisResults.CorrCoef.gammaBandPower.p);
T(3,:) = cell2table(AnalysisResults.CorrCoef.CBV_HbT.meanStD);
T(4,:) = cell2table(AnalysisResults.CorrCoef.CBV_HbT.p);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% Table 3
summaryTable = figure('Name','Table3'); %#ok<*NASGU>
sgtitle('Table 3 Turner Manuscript 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
%% save figure(s)
% dirpath = [rootFolder '\Summary Figures and Structures\'];
% if ~exist(dirpath,'dir')
%     mkdir(dirpath);
% end
% savefig(summaryTable,[dirpath 'Table3']);

end
