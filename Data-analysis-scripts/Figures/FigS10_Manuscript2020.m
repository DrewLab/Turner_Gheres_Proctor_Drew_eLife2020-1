function [] = FigS10_Manuscript2020(rootFolder)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
allCatLabels = [];
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
%% extract data from each animal's sleep scoring results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    dataLoc = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLoc)
    scoringResults = 'Forest_ScoringResults.mat';
    load(scoringResults,'-mat')
    numberOfScores(a,1) = length(ScoringResults.alllabels); %#ok<*AGROW,*NASGU>
    indAwakePerc(a,1) = round((sum(strcmp(ScoringResults.alllabels,'Not Sleep'))/length(ScoringResults.alllabels))*100,1);
    indNremPerc(a,1) = round((sum(strcmp(ScoringResults.alllabels,'NREM Sleep'))/length(ScoringResults.alllabels))*100,1);
    indRemPerc(a,1) = round((sum(strcmp(ScoringResults.alllabels,'REM Sleep'))/length(ScoringResults.alllabels))*100,1);
    allCatLabels = vertcat(allCatLabels,ScoringResults.alllabels); 
end
labels = {'Awake','NREM','REM'};
% mean percentage of each state between animals
meanAwakePerc = mean(indAwakePerc,1);
stdAwakePerc = std(indAwakePerc,0,1);
meanNremPerc = mean(indNremPerc,1);
stdNremPerc = std(indNremPerc,0,1);
meanRemPerc = mean(indRemPerc,1);
stdRemPerc = std(indRemPerc,0,1);
meanPercs = horzcat(meanAwakePerc,meanNremPerc,meanRemPerc);
% percentage of each state for all labels together
allAwakePerc = round((sum(strcmp(allCatLabels,'Not Sleep'))/length(allCatLabels))*100,1);
allNremPerc = round((sum(strcmp(allCatLabels,'NREM Sleep'))/length(allCatLabels))*100,1);
allRemPerc = round((sum(strcmp(allCatLabels,'REM Sleep'))/length(allCatLabels))*100,1);
meanAllPercs = horzcat(allAwakePerc,allNremPerc,allRemPerc);
% total time per animal behavioral states
labelTime = 5;   % seconds
IOS_indTotalTimeHours = ((numberOfScores*labelTime)/60)/60;
IOS_allTimeHours = sum(IOS_indTotalTimeHours);
IOS_meanTimeHours = mean(IOS_indTotalTimeHours,1);
IOS_stdTimeHours = std(IOS_indTotalTimeHours,0,1);
allTimeDays = sum(IOS_indTotalTimeHours)/24;
totalTimeAwake = IOS_indTotalTimeHours.*(indAwakePerc/100);
meanAwakeHours = mean(totalTimeAwake,1);
stdAwakeHours = std(totalTimeAwake,0,1);
totalTimeNREM = IOS_indTotalTimeHours.*(indNremPerc/100);
meanNREMHours = mean(totalTimeNREM);
stdNREMHours = std(totalTimeNREM,0,1);
totalTimeREM = IOS_indTotalTimeHours.*(indRemPerc/100);
meanREMHours = mean(totalTimeREM);
stdREMHours = std(totalTimeREM,0,1);
%% 2p data
animalIDs2 = {'T115','T116','T117','T118','T125','T126'};
allFileIDs = [];
% extract data from each animal's sleep scoring results
for a = 1:length(animalIDs2)
    animalID = animalIDs2{1,a};
    dataLoc = [rootFolder '/' animalID '/2P Data/'];
    cd(dataLoc)
     % Character list of all MergedData files
    mergedDirectory = dir('*_MergedData.mat');
    mergedDataFiles = {mergedDirectory.name}';
    mergedDataFileIDs = char(mergedDataFiles);
    allFileIDs(a,1) = size(mergedDataFileIDs,1);
end
PLSM_indTotalTimeHours = (allFileIDs.*15)./60;
PLSM_allTimeHours = sum(PLSM_indTotalTimeHours);
PLSM_meanTimeHours = mean(PLSM_indTotalTimeHours,1);
PLSM_stdTimeHours = std(PLSM_indTotalTimeHours,0,1);
allHours = IOS_allTimeHours + PLSM_allTimeHours;
%% supplemental figure panel
summaryFigure = figure;
sgtitle('Temporary Supplemental Figure Panel - Turner Manuscript 2020')
%% [A] Perc of behav states scores
ax1 = subplot(2,2,1);
xInds = ones(1,length(animalIDs));
s1 = scatter(xInds*1,indAwakePerc,75,'MarkerEdgeColor','k','MarkerFaceColor',colors_Manuscript2020('rich black'),'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,meanAwakePerc,stdAwakePerc,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(xInds*2,indNremPerc,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,meanNremPerc,stdNremPerc,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(xInds*3,indRemPerc,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,meanRemPerc,stdRemPerc,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[A] Arousal state probability','sleep scoring labels',''})
ylabel('Probability (%)')
legend([s1,s2,s3],'Awake','NREM','REM')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(labels) + 1])
ylim([0,100])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [B] Ternary
ax2 = subplot(2,2,2);
terplot_Manuscript2020();
[hd] = ternaryc_Manuscript2020(indAwakePerc/100,indNremPerc/100,indRemPerc/100);
hlabels = terlabel_Manuscript2020('Not asleep','NREM sleep','REM sleep');
%% [C] Perc of behav states scores
ax3 = subplot(2,2,3);
p1 = pie(meanPercs);
pText = findobj(p1,'Type','text');
percentValues = get(pText,'String'); 
txt = {'Awake: ';'NREM: ';'REM: '}; 
combinedtxt = strcat(txt,percentValues); 
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
title({'[C] Sleep scoring label probability','Mean animal sleep scoring labels',''})
%% [D] Perc of behav states scores (All labels together)
ax4 = subplot(2,2,4);
p2 = pie(meanAllPercs);
pText = findobj(p2,'Type','text');
percentValues = get(pText,'String'); 
txt = {'Awake: ';'NREM: ';'REM: '}; 
combinedtxt = strcat(txt,percentValues); 
pText(1).String = combinedtxt(1);
pText(2).String = combinedtxt(2);
pText(3).String = combinedtxt(3);
title({'[D] Sleep scoring label probability','All sleep scoring labels',''})
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Test Supplemental Figure Panel']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Test Supplemental Figure Pane'])
%% table of IOS random forest scores
arousalStateTime = figure;
variableNames = {'TotalTimeHrs','AwakePerc','AwakeTimeHrs','NREMPerc','NREMTimeHrs','REMPerc','REMTimeHrs'};
T = table(IOS_indTotalTimeHours,indAwakePerc,totalTimeAwake,indNremPerc,totalTimeNREM,indRemPerc,totalTimeREM,'RowNames',animalIDs,'VariableNames',variableNames);
uitable('Data',T{:,:},'ColumnName',T.Properties.VariableNames,'RowName',T.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
uicontrol('Style','text','Position',[200,300,200,200],'String',{'Total Time (days): ' num2str(allTimeDays),'Mean Awake time per animal (Hrs): ' num2str(meanAwakeHours) ' +/- ' num2str(stdAwakeHours),'Mean NREM time per animal (Hrs): ' num2str(meanNREMHours) ' +/- ' num2str(stdNREMHours),'Mean REM time per animal (Hrs): ' num2str(meanREMHours) ' +/- ' num2str(stdREMHours)});
savefig(arousalStateTime,[dirpath 'Test Supplemental Figure Panel Table']);
% set(arousalStateTime,'PaperPositionMode','auto');
% print('-painters','-dpdf','-fillpage',[dirpath 'Test Supplemental Figure Panel Table'])
end