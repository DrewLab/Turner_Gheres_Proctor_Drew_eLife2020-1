function [] = FigurePanelFive_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T115','T117','T118','T125','T126'};   % T116 has no REM events
%% REM dilation and trnasition REM to Awake
% cd through each animal's directory and extract the appropriate analysis results
data.NREMtoREM.data = []; data.REMtoAwake.data = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    evokedBehavFields = fieldnames(AnalysisResults.(animalID).Transitions);
    for bb = 1:length(evokedBehavFields)
        behavField = evokedBehavFields{bb,1};
        vesselIDs = fieldnames(AnalysisResults.(animalID).Transitions.(behavField));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            data.(behavField).data = vertcat(data.(behavField).data,AnalysisResults.(animalID).Transitions.(behavField).(vesselID).mean);
            data.(behavField).timeVector = AnalysisResults.(animalID).Transitions.(behavField).(vesselID).timeVector;
        end
    end
end
% take the average of the vessels for each behavior
evokedBehavFields = {'NREMtoREM','REMtoAwake'};
for dd = 1:length(evokedBehavFields)
    behavField = evokedBehavFields{1,dd};
    data.(behavField).mean = mean(data.(behavField).data,1);
    data.(behavField).StD = std(data.(behavField).data,0,1);
end
%% Figure panel fiveB
summaryFigure = figure;
sgtitle('Figure panel 5 - Turner Manuscript 2020')
%% [A] NREM to REM transition
ax1 = subplot(1,2,1);
plot(data.NREMtoREM.timeVector,data.NREMtoREM.mean,'color',colors_Manuscript2020('rich black'),'LineWidth',2)
hold on;
plot(data.NREMtoREM.timeVector,data.NREMtoREM.mean + data.NREMtoREM.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.NREMtoREM.timeVector,data.NREMtoREM.mean - data.NREMtoREM.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[A] NREM to REM transition')
xlabel('Time (s)')
ylabel('\DeltaD/D (%)')
axis tight
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [B] REM to Awake transition
ax2 = subplot(1,2,2);
plot(data.REMtoAwake.timeVector,data.REMtoAwake.mean,'color',colors_Manuscript2020('rich black'),'LineWidth',2)
hold on;
plot(data.REMtoAwake.timeVector,data.REMtoAwake.mean + data.REMtoAwake.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.REMtoAwake.timeVector,data.REMtoAwake.mean - data.REMtoAwake.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[B] REM to Awake transition')
xlabel('Time (s)')
ylabel('\DeltaD/D (%)')
axis tight
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Figure Panel 5']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Figure Panel 5'])

end
