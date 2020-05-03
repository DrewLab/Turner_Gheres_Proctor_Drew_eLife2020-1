function [] = SupplementalFigurePanelSevenB_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T115','T116','T117','T118','T125','T126'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
%% cd through each animal's directory and extract the appropriate analysis results
data.EvokedAvgs.ShortWhisks.means = [];
data.EvokedAvgs.IntermediateWhisks.means = [];
data.EvokedAvgs.LongWhisks.means = [];
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for c = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,c};
        vesselIDs = fieldnames(AnalysisResults.(animalID).EvokedAvgs.(whiskDataType));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            if strcmp(vesselID(1),'V') == false
                % LH cortical
                data.EvokedAvgs.(whiskDataType).means = vertcat(data.EvokedAvgs.(whiskDataType).means,AnalysisResults.(animalID).EvokedAvgs.(whiskDataType).(vesselID).mean);
                data.EvokedAvgs.(whiskDataType).timeVector = AnalysisResults.(animalID).EvokedAvgs.(whiskDataType).(vesselID).timeVector;
            end
        end
    end
end
% mean/std of the data
for e = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,e};
    data.EvokedAvgs.(whiskDataType).mean = mean(data.EvokedAvgs.(whiskDataType).means,1);
    data.EvokedAvgs.(whiskDataType).StD = std(data.EvokedAvgs.(whiskDataType).means,0,1);
end
%% summary figure(s)
summaryFigure = figure;
sgtitle('Supplemental Figure Panel 8 - Turner Manuscript 2020')
%% [A] Short whisks
ax1 = subplot(1,3,1);
plot(data.EvokedAvgs.ShortWhisks.timeVector,data.EvokedAvgs.ShortWhisks.mean,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.EvokedAvgs.ShortWhisks.timeVector,data.EvokedAvgs.ShortWhisks.mean + data.EvokedAvgs.ShortWhisks.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.EvokedAvgs.ShortWhisks.timeVector,data.EvokedAvgs.ShortWhisks.mean - data.EvokedAvgs.ShortWhisks.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[A] Short whisk response')
ylabel('\DeltaD/D (%)')
xlabel('Peri-whisk time (s)')   
axis square
xlim([-2,10])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [B] Intermediate whisks
ax2 = subplot(1,3,2);
plot(data.EvokedAvgs.IntermediateWhisks.timeVector,data.EvokedAvgs.IntermediateWhisks.mean,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.EvokedAvgs.IntermediateWhisks.timeVector,data.EvokedAvgs.IntermediateWhisks.mean + data.EvokedAvgs.IntermediateWhisks.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.EvokedAvgs.IntermediateWhisks.timeVector,data.EvokedAvgs.IntermediateWhisks.mean - data.EvokedAvgs.IntermediateWhisks.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[B] Intermed whisk response')
ylabel('\DeltaD/D (%)')
xlabel('Peri-whisk time (s)')   
axis square
xlim([-2,10])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [C] Long whisks
ax3 = subplot(1,3,3);
plot(data.EvokedAvgs.LongWhisks.timeVector,data.EvokedAvgs.LongWhisks.mean,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.EvokedAvgs.LongWhisks.timeVector,data.EvokedAvgs.LongWhisks.mean + data.EvokedAvgs.LongWhisks.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.EvokedAvgs.LongWhisks.timeVector,data.EvokedAvgs.LongWhisks.mean - data.EvokedAvgs.LongWhisks.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[C] Long whisk response')
ylabel('\DeltaD/D (%)')
xlabel('Peri-whisk time (s)')   
axis square
xlim([-2,10])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
linkaxes([ax1,ax2,ax3],'xy')
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Supplemental Figure Panel 7b']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Supplemental Figure Panel 7b'])

end
