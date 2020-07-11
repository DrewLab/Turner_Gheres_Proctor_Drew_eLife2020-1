function [AnalysisResults] = FigS20_Manuscript2020(rootFolder,saveFigs,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel S20 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','Whisk','NREM','REM'};
dataTypes = {'gammaBandPower'};
colorRest = [(51/256),(160/256),(44/256)];
colorNREM = [(192/256),(0/256),(256/256)];
colorREM = [(255/256),(140/256),(0/256)];
% colorAwake = [(256/256),(192/256),(0/256)];
% colorSleep = [(0/256),(128/256),(256/256)];
% colorAll = [(184/256),(115/256),(51/256)];
colorWhisk = [(31/256),(120/256),(180/256)];
% colorStim = [(256/256),(28/256),(207/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% extract data from each animal's sleep scoring results
TwoPallCatMeans = AnalysisResults.TwoPSleepProbability.TwoPCatMeans;
awakeProbPerc = AnalysisResults.TwoPSleepProbability.awakeProbPerc;
nremProbPerc = AnalysisResults.TwoPSleepProbability.nremProbPerc;
remProbPerc = AnalysisResults.TwoPSleepProbability.remProbPerc;
%% take data from each animal corresponding to the CBV-gamma relationship
catHbT = []; catGam = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        if isfield(catHbT,behavField) == false
            catHbT.(behavField) = []; 
            catGam.(behavField) = [];
        end
        catHbT.(behavField) = cat(1,catHbT.(behavField),AnalysisResults.(animalID).HbTvsGamma.(behavField).HbT.MeanAdjLH,AnalysisResults.(animalID).HbTvsGamma.(behavField).HbT.MeanAdjRH);
        catGam.(behavField) = cat(1,catGam.(behavField),AnalysisResults.(animalID).HbTvsGamma.(behavField).Gamma.MeanAdjLH,AnalysisResults.(animalID).HbTvsGamma.(behavField).Gamma.MeanAdjRH);
    end
end
%% Fig. S20
summaryFigure = figure('Name','FigS20 (a,b)'); %#ok<*NASGU>
sgtitle('Figure Panel S20 (a,b) Turner Manuscript 2020')
%% [S20a] HbT vs. arousal state probability
ax1 = subplot(1,2,1);
edges = -20:1:50;
yyaxis right
h1 = histogram(TwoPallCatMeans,edges,'Normalization','probability','EdgeColor','k','FaceColor',colors_Manuscript2020('dark candy apple red'));
ylabel({'5-sec Mean \DeltaD/D (%)','Probability distribution'},'rotation',-90,'VerticalAlignment','bottom')
yyaxis left
p1 = plot(edges,sgolayfilt(medfilt1(awakeProbPerc,10,'truncate'),3,17),'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
p2 = plot(edges,sgolayfilt(medfilt1(nremProbPerc,10,'truncate'),3,17),'-','color',colorNREM,'LineWidth',2);
p3 = plot(edges,sgolayfilt(medfilt1(remProbPerc,10,'truncate'),3,17),'-','color',colorREM,'LineWidth',2);
ylabel({'Arousal-state probability (%)'})
xlim([-20,50])
ylim([0,90])
legend([p1,p2,p3,h1],'Awake','NREM','REM','\DeltaD/D','Location','NorthEast')
title({'[S20a] 5-sec mean \DeltaD/D (%)','vs. arousal state probability',''})
xlabel({'\DeltaD/D (%)','1 (%) bins'})
axis square
set(gca,'box','off')
set(gca,'TickLength',[0.03,0.03]);
ylim([0,90])
xlim([-20,50])
set(h1,'facealpha',0.2);
ax1.TickLength = [0.03,0.03];
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = colors_Manuscript2020('dark candy apple red');
%% [S20b]
ax2 = subplot(1,2,2);
s1 = scatter(catGam.NREM,catHbT.NREM,'MarkerFaceColor',colorNREM,'MarkerEdgeColor','k');
hold on
s2 = scatter(catGam.REM,catHbT.REM,'MarkerFaceColor',colorREM,'MarkerEdgeColor','k');
s3 = scatter(catGam.Whisk,catHbT.Whisk,'MarkerFaceColor',colorWhisk,'MarkerEdgeColor','k');
s4 = scatter(catGam.Rest,catHbT.Rest,'MarkerFaceColor',colorRest,'MarkerEdgeColor','k');
xlabel('Gamma-band DeltaP/P (%)')
ylabel('\Delta[HbT] (\muM)')
title({'[S20b] Gamma-band \DeltaP/P (%)','vs. \Delta[HbT] (\muM)',''})
legend([s1,s2,s3,s4],'Rest','Whisk','NREM','REM','Location','NorthEast')
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder '\Summary Figures and Structures\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'FigS20']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'FigS20'])
end

end
