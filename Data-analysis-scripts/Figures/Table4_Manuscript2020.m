function [AnalysisResults] = Table4_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table S4 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% set-up and process data for Table S4
columnNames = AnalysisResults.PSD.columnNames;
rowNames = {'Delta_S01_meanStD','Delta_S01_pVal','Theta_S01_meanStD','Theta_S01_pVal'...
    'Alpha_S01_meanStD','Alpha_S01_pVal','Beta_S01_meanStD','Beta_S01_pVal'...
    'Delta_S001_meanStD','Delta_S001_pVal','Theta_S001_meanStD','Theta_S001_pVal'...
    'Alpha_S001_meanStD','Alpha_S001_pVal','Beta_S001_meanStD','Beta_S001_pVal'};
T(1,:) = cell2table(AnalysisResults.PSD.deltaBandPower.meanStD01);
T(2,:) = cell2table(AnalysisResults.PSD.deltaBandPower.p01);
T(3,:) = cell2table(AnalysisResults.PSD.thetaBandPower.meanStD01);
T(4,:) = cell2table(AnalysisResults.PSD.thetaBandPower.p01);
T(5,:) = cell2table(AnalysisResults.PSD.alphaBandPower.meanStD01);
T(6,:) = cell2table(AnalysisResults.PSD.alphaBandPower.p01);
T(7,:) = cell2table(AnalysisResults.PSD.betaBandPower.meanStD01);
T(8,:) = cell2table(AnalysisResults.PSD.betaBandPower.p01);
T(9,:) = cell2table(AnalysisResults.PSD.deltaBandPower.meanStD001);
T(10,:) = cell2table(AnalysisResults.PSD.deltaBandPower.p001);
T(11,:) = cell2table(AnalysisResults.PSD.thetaBandPower.meanStD001);
T(12,:) = cell2table(AnalysisResults.PSD.thetaBandPower.p001);
T(13,:) = cell2table(AnalysisResults.PSD.alphaBandPower.meanStD001);
T(14,:) = cell2table(AnalysisResults.PSD.alphaBandPower.p001);
T(15,:) = cell2table(AnalysisResults.PSD.betaBandPower.meanStD001);
T(16,:) = cell2table(AnalysisResults.PSD.betaBandPower.p001);
T.Properties.RowNames = rowNames;
T.Properties.VariableNames = columnNames;
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
%% Table S4
summaryTable = figure('Name','TableS4');
sgtitle('Table S4 Turner Manuscript 2020')
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
savefig(summaryTable,[dirpath 'TableS4']);

end
