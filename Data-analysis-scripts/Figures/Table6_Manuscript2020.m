function [] = Table6_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table S6 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% set-up and process data for Table S6
columnNames = AnalysisResults.CorrCoef.columnNames;
rowNames = {'Delta_R_meanStD','Delta_R_pVal','Theta_R_meanStD','Theta_R_pVal'...
    'Alpha_R_meanStD','Alpha_R_pVal','Beta_R_meanStD','Beta_R_pVal'};
T(1,:) = cell2table(AnalysisResults.CorrCoef.deltaBandPower.meanStD);
T(2,:) = cell2table(AnalysisResults.CorrCoef.deltaBandPower.p);
T(3,:) = cell2table(AnalysisResults.CorrCoef.thetaBandPower.meanStD);
T(4,:) = cell2table(AnalysisResults.CorrCoef.thetaBandPower.p);
T(5,:) = cell2table(AnalysisResults.CorrCoef.alphaBandPower.meanStD);
T(6,:) = cell2table(AnalysisResults.CorrCoef.alphaBandPower.p);
T(7,:) = cell2table(AnalysisResults.CorrCoef.betaBandPower.meanStD);
T(8,:) = cell2table(AnalysisResults.CorrCoef.betaBandPower.p);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
%% Table S6
summaryTable = figure('Name','TableS6');
sgtitle('Table S6 Turner Manuscript 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
savefig(summaryTable,[dirpath 'TableS6']);

end
