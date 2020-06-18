function [] = test_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
colorA = [(0/256),(166/256),(81/256)];   % rest color
colorB = [(191/256),(0/256),(255/256)];   % NREM color
colorC = [(254/256),(139/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(179/256)];  % whisk color
%% Mean HbT comparison between behaviors
% pre-allocate the date for each day
behavFields = {'Rest','Whisk','NREM','REM'};
for aa = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        data.(animalID).(behavField).HbTindLH = AnalysisResults.(animalID).HbTvsGamma.(behavField).HbT.IndAdjLH;
        data.(animalID).(behavField).HbTindRH = AnalysisResults.(animalID).HbTvsGamma.(behavField).HbT.IndAdjRH;
        data.(animalID).(behavField).GamIndLH = AnalysisResults.(animalID).HbTvsGamma.(behavField).Gamma.IndAdjLH;
        data.(animalID).(behavField).GamIndRH = AnalysisResults.(animalID).HbTvsGamma.(behavField).Gamma.IndAdjRH;
        data.(animalID).(behavField).HbTmeanLH = AnalysisResults.(animalID).HbTvsGamma.(behavField).HbT.MeanAdjLH;
        data.(animalID).(behavField).HbTmeanRH = AnalysisResults.(animalID).HbTvsGamma.(behavField).HbT.MeanAdjRH;
        data.(animalID).(behavField).GamMeanLH = AnalysisResults.(animalID).HbTvsGamma.(behavField).Gamma.MeanAdjLH;
        data.(animalID).(behavField).GamMeanRH = AnalysisResults.(animalID).HbTvsGamma.(behavField).Gamma.MeanAdjRH;
    end
end
% concatenate into arrays
data.Rest.IndHbT = []; data.Rest.IndGam = []; data.Rest.MeanHbT = []; data.Rest.MeanGam = [];
data.Whisk.IndHbT = []; data.Whisk.IndGam = []; data.Whisk.MeanHbT = []; data.Whisk.MeanGam = [];
data.NREM.IndHbT = []; data.NREM.IndGam = []; data.NREM.MeanHbT = []; data.NREM.MeanGam = [];
data.REM.IndHbT = []; data.REM.IndGam = []; data.REM.MeanHbT = []; data.REM.MeanGam = [];
for cc = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,cc};
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        for ee = 1:length(data.(animalID).(behavField).HbTindLH)
            data.(behavField).IndHbT = cat(2,data.(behavField).IndHbT,data.(animalID).(behavField).HbTindLH{ee,1},data.(animalID).(behavField).HbTindRH{ee,1});
            data.(behavField).IndGam = cat(2,data.(behavField).IndGam,data.(animalID).(behavField).GamIndLH{ee,1},data.(animalID).(behavField).GamIndRH{ee,1});
            data.(behavField).MeanHbT = cat(2,data.(behavField).MeanHbT,data.(animalID).(behavField).HbTmeanLH(ee,1),data.(animalID).(behavField).HbTmeanRH(ee,1));
            data.(behavField).MeanGam = cat(2,data.(behavField).MeanGam,data.(animalID).(behavField).GamMeanLH(ee,1),data.(animalID).(behavField).GamMeanRH(ee,1));
        end
    end
end
%%
figure;
% scatter(data.NREM.MeanGam*100,data.NREM.MeanHbT,'MarkerFaceColor',colorB)
hold on
% scatter(data.REM.MeanGam*100,data.REM.MeanHbT,'MarkerFaceColor',colorC)
% scatter(data.Whisk.MeanGam*100,data.Whisk.MeanHbT,'MarkerFaceColor',colorD)
scatter(data.Rest.MeanGam*100,data.Rest.MeanHbT,'MarkerFaceColor',colorA)
% ylim([-20,120])
% xlim([-50,200])
xlabel('Gamma-band power (%)')
ylabel('\Delta[HbT] (\muM)')

%%
figure; 
histogram2(data.Rest.MeanHbT,data.Rest.MeanGam,'Normalization','pdf','XBinLimits',[-20,120],'YBinLimits',[-0.5,10]);
hold on
histogram2(data.NREM.MeanHbT,data.NREM.MeanGam,'Normalization','pdf','XBinLimits',[-20,120],'YBinLimits',[-0.5,10]);
histogram2(data.REM.MeanHbT,data.REM.MeanGam,'Normalization','pdf','XBinLimits',[-20,120],'YBinLimits',[-0.5,10]);
histogram2(data.Whisk.MeanHbT,data.Whisk.MeanGam,'Normalization','pdf','XBinLimits',[-20,120],'YBinLimits',[-0.5,10]);

%%
figure; 
histogram2(data.Rest.MeanHbT,data.Rest.MeanGam,'Normalization','countdensity')%,'Normalization','pdf','BinWidth',[10,0.25],'XBinLimits',[-50,150],'YBinLimits',[-.5,2]);
hold on
histogram2(data.NREM.MeanHbT,data.NREM.MeanGam,'Normalization','countdensity')%,'Normalization','pdf','BinWidth',[10,0.25],'XBinLimits',[-50,150],'YBinLimits',[-.5,2]);
histogram2(data.REM.MeanHbT,data.REM.MeanGam,'Normalization','pdf')%,'Normalization','pdf','BinWidth',[10,0.25],'XBinLimits',[-50,150],'YBinLimits',[-.5,2]);
histogram2(data.Whisk.MeanHbT,data.Whisk.MeanGam,'Normalization','pdf')%,'Normalization','pdf','BinWidth',[10,0.25],'XBinLimits',[-50,150],'YBinLimits',[-.5,2]);

%%
figure; 
histogram2(data.Rest.IndHbT,data.Rest.IndGam)%,'Normalization','pdf','BinWidth',[10,0.25],'XBinLimits',[-50,150],'YBinLimits',[-.5,2]);
hold on
histogram2(data.NREM.IndHbT,data.NREM.IndGam)%,'Normalization','pdf','BinWidth',[10,0.25],'XBinLimits',[-50,150],'YBinLimits',[-.5,2]);
histogram2(data.REM.IndHbT,data.REM.IndGam)%,'Normalization','pdf','BinWidth',[10,0.25],'XBinLimits',[-50,150],'YBinLimits',[-.5,2]);
histogram2(data.Whisk.IndHbT,data.Whisk.IndGam)%,'Normalization','countdensity')%,'Normalization','pdf','BinWidth',[10,0.25],'XBinLimits',[-50,150],'YBinLimits',[-.5,2]);


%% Pixel panel 5
summaryFigure = figure('Name','Fig5 (a-f)');
sgtitle('Figure panel 5 (a-f) Turner Manuscript 2020')
%% [5a] Mean HbT during different behaviors
ax1 = subplot(2,3,1);
HbT_xInds = ones(1,length(IOS_animalIDs)*2);
s1 = scatter(HbT_xInds*1,procData.HbT.Rest.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,procData.HbT.Rest.MeanCBV,procData.HbT.Rest.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(HbT_xInds*2,procData.HbT.Whisk.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,procData.HbT.Whisk.MeanCBV,procData.HbT.Whisk.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(HbT_xInds*3,procData.HbT.Stim.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorE,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,procData.HbT.Stim.MeanCBV,procData.HbT.Stim.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(HbT_xInds*4,procData.HbT.NREM.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,procData.HbT.NREM.MeanCBV,procData.HbT.NREM.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(HbT_xInds*5,procData.HbT.REM.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,procData.HbT.REM.MeanCBV,procData.HbT.REM.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title({'[5a] Mean \DeltaHbT (\muM)','during arousal-states',''})
ylabel('\DeltaHbT (\muM)')
legend([s1,s2,s3,s4,s5],'Awake Rest','Whisk','Stim','NREM','REM','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
ylim([-10,100])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [5b] Mean vessel diameter during different behaviors
ax2 = subplot(2,3,2);
TwoP_xIndsRest = ones(1,length(procData.TwoP.Rest.IndMeanDiam));
TwoP_xIndsWhisk = ones(1,length(procData.TwoP.Whisk.IndMeanDiam));
TwoP_xIndsNREM = ones(1,length(procData.TwoP.NREM.IndMeanDiam));
TwoP_xIndsREM = ones(1,length(procData.TwoP.REM.IndMeanDiam));
scatter(TwoP_xIndsRest*1,procData.TwoP.Rest.IndMeanDiam,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,procData.TwoP.Rest.MeanDiam,procData.TwoP.Rest.StdMeanDiam,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(TwoP_xIndsWhisk*2,procData.TwoP.Whisk.IndMeanDiam,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,procData.TwoP.Whisk.MeanDiam,procData.TwoP.Whisk.StdMeanDiam,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(TwoP_xIndsNREM*3,procData.TwoP.NREM.IndMeanDiam,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,procData.TwoP.NREM.MeanDiam,procData.TwoP.NREM.StdMeanDiam,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(TwoP_xIndsREM*4,procData.TwoP.REM.IndMeanDiam,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,procData.TwoP.REM.MeanDiam,procData.TwoP.REM.StdMeanDiam,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[5b] Mean \DeltaD/D (%)','during arousal-states',''})
ylabel('\DeltaD/D (%)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(TwoP_behavFields) + 1])
ylim([-10,60])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [5c] Mean vessel diameter during different behaviors
ax3 = subplot(2,3,3);
LDF_xInds = ones(1,length(LDF_animalIDs));
scatter(LDF_xInds*1,procData.LDF.Rest.IndMeanLDF,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,procData.LDF.Rest.MeanLDF,procData.LDF.Rest.StdLDF,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(LDF_xInds*2,procData.LDF.Whisk.IndMeanLDF,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,procData.LDF.Whisk.MeanLDF,procData.LDF.Whisk.StdLDF,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(LDF_xInds*3,procData.LDF.NREM.IndMeanLDF,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,procData.LDF.NREM.MeanLDF,procData.LDF.NREM.StdLDF,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(LDF_xInds*4,procData.LDF.REM.IndMeanLDF,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,procData.LDF.REM.MeanLDF,procData.LDF.REM.StdLDF,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[5c] Mean \DeltaQ/Q (%)','during arousal-states',''})
ylabel('\DeltaQ/Q (%)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(LDF_behavFields) + 1])
ylim([-10,80])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [5a bottom] Mean HbT distribution during different behaviors
ax4 = subplot(2,3,4);
edges = -35:15:150;
[curve1] = SmoothHistogramBins_Manuscript2020(procData.HbT.Rest.CatCBV,edges);
[curve2] = SmoothHistogramBins_Manuscript2020(procData.HbT.Whisk.CatCBV,edges);
[curve3] = SmoothHistogramBins_Manuscript2020(procData.HbT.Stim.CatCBV,edges);
[curve4] = SmoothHistogramBins_Manuscript2020(procData.HbT.NREM.CatCBV,edges);
[curve5] = SmoothHistogramBins_Manuscript2020(procData.HbT.REM.CatCBV,edges);
before = findall(gca);
fnplt(curve1);
added = setdiff(findall(gca),before);
set(added,'Color',colorA)
hold on
before = findall(gca);
fnplt(curve2);
added = setdiff(findall(gca),before);
set(added,'Color',colorD)
before = findall(gca);
fnplt(curve3);
added = setdiff(findall(gca),before);
set(added,'Color',colorE)
before = findall(gca);
fnplt(curve4);
added = setdiff(findall(gca),before);
set(added,'Color',colorB)
before = findall(gca);
fnplt(curve5);
added = setdiff(findall(gca),before);
set(added,'Color',colorC)
title({'\DeltaHbT (\muM)','arousal-state distribution',''})
xlabel('\DeltaHbT (\muM)')
ylabel('Probability')
axis square
set(gca,'box','off')
ylim([0,1])
ax4.TickLength = [0.03,0.03];
%% [5b bottom] vessel diameter distribution during different behaviors
ax5 = subplot(2,3,5);
edges = -20:10:70;
[curve1] = SmoothHistogramBins_Manuscript2020(procData.TwoP.Rest.CatIndDiam,edges);
[curve2] = SmoothHistogramBins_Manuscript2020(procData.TwoP.Whisk.CatIndDiam,edges);
[curve3] = SmoothHistogramBins_Manuscript2020(procData.TwoP.NREM.CatIndDiam,edges);
[curve4] = SmoothHistogramBins_Manuscript2020(procData.TwoP.REM.CatIndDiam,edges);
before = findall(gca);
fnplt(curve1);
added = setdiff(findall(gca),before);
set(added,'Color',colorA)
hold on
before = findall(gca);
fnplt(curve2);
added = setdiff(findall(gca),before);
set(added,'Color',colorD)
before = findall(gca);
fnplt(curve3);
added = setdiff(findall(gca),before);
set(added,'Color',colorB)
before = findall(gca);
fnplt(curve4);
added = setdiff(findall(gca),before);
set(added,'Color',colorC)
title({'\DeltaD/D (%)','arousal-state distribution',''})
xlabel('\DeltaD/D (%)')
ylabel('Probability')
axis square
set(gca,'box','off')
xlim([-20,70])
ylim([0,1])
ax5.TickLength = [0.03,0.03];
%% [5c bottom] LDF arousal-state vessel distribution
ax6 = subplot(2,3,6);
edgesA = -30:3:80;
edgesB = -30:20:80;
[curve1] = SmoothHistogramBins_Manuscript2020(procData.LDF.Rest.CatLDF,edgesA);
[curve2] = SmoothHistogramBins_Manuscript2020(procData.LDF.Whisk.CatLDF,edgesB);
[curve3] = SmoothHistogramBins_Manuscript2020(procData.LDF.NREM.CatLDF,edgesB);
[curve4] = SmoothHistogramBins_Manuscript2020(procData.LDF.REM.CatLDF,edgesB);
before = findall(gca);
fnplt(curve1);
added = setdiff(findall(gca),before);
set(added,'Color',colorA)
hold on
before = findall(gca);
fnplt(curve2);
added = setdiff(findall(gca),before);
set(added,'Color',colorD)
before = findall(gca);
fnplt(curve3);
added = setdiff(findall(gca),before);
set(added,'Color',colorB)
before = findall(gca);
fnplt(curve4);
added = setdiff(findall(gca),before);
set(added,'Color',colorC)
title({'\DeltaQ/Q (%)','arousal-state distribution',''})
xlabel('\DeltaQ/Q (%)')
ylabel('Probability')
axis square
axis tight
ylim([0,1])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
set(summaryFigure,'PaperPositionMode','auto');
savefig(summaryFigure,[dirpath 'Fig5']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Fig5'])
%% statistical diary
diaryFile = [dirpath 'Fig5_Statistics.txt'];
if exist(diaryFile,'file') == 2
    delete(diaryFile)
end
diary(diaryFile)
diary on
% HbT statistical diary
disp('======================================================================================================================')
disp('[5a] Generalized linear mixed-effects model statistics for mean HbT during Rest, Whisk, Stim, NREM, and REM')
disp('======================================================================================================================')
disp(HbTStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.05 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HbTCI{1,1}(1,:))])
disp(['Whisk: ' num2str(HbTCI{1,1}(2,:))])
disp(['Stim: ' num2str(HbTCI{1,1}(3,:))])
disp(['NREM: ' num2str(HbTCI{1,1}(4,:))])
disp(['REM: ' num2str(HbTCI{1,1}(5,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.01 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HbTCI{2,1}(1,:))])
disp(['Whisk: ' num2str(HbTCI{2,1}(2,:))])
disp(['Stim: ' num2str(HbTCI{2,1}(3,:))])
disp(['NREM: ' num2str(HbTCI{2,1}(4,:))])
disp(['REM: ' num2str(HbTCI{2,1}(5,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.001 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HbTCI{3,1}(1,:))])
disp(['Whisk: ' num2str(HbTCI{3,1}(2,:))])
disp(['Stim: ' num2str(HbTCI{3,1}(3,:))])
disp(['NREM: ' num2str(HbTCI{3,1}(4,:))])
disp(['REM: ' num2str(HbTCI{3,1}(5,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
% Peak vessel diameter statistical diary
disp('======================================================================================================================')
disp('[5b] Generalized linear mixed-effects model statistics for mean vessel diameter during Rest, Whisk, NREM, and REM')
disp('======================================================================================================================')
disp(vesselStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.05 confidence intervals with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(vesselCI{1,1}(1,:))])
disp(['Whisk: ' num2str(vesselCI{1,1}(2,:))])
disp(['NREM: ' num2str(vesselCI{1,1}(3,:))])
disp(['REM: ' num2str(vesselCI{1,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.01 confidence intervals with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(vesselCI{2,1}(1,:))])
disp(['Whisk: ' num2str(vesselCI{2,1}(2,:))])
disp(['NREM: ' num2str(vesselCI{2,1}(3,:))])
disp(['REM: ' num2str(vesselCI{2,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.001 confidence intervals with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(vesselCI{3,1}(1,:))])
disp(['Whisk: ' num2str(vesselCI{3,1}(2,:))])
disp(['NREM: ' num2str(vesselCI{3,1}(3,:))])
disp(['REM: ' num2str(vesselCI{3,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
% LDF flow statistical diary
disp('======================================================================================================================')
disp('[5c] Generalized linear mixed-effects model statistics for mean doppler flow during Rest, Whisk, NREM, and REM')
disp('======================================================================================================================')
disp(flowStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.05 confidence intervals with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(flowCI{1,1}(1,:))])
disp(['Whisk: ' num2str(flowCI{1,1}(2,:))])
disp(['NREM: ' num2str(flowCI{1,1}(3,:))])
disp(['REM: ' num2str(flowCI{1,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.01 confidence intervals with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(flowCI{2,1}(1,:))])
disp(['Whisk: ' num2str(flowCI{2,1}(2,:))])
disp(['NREM: ' num2str(flowCI{2,1}(3,:))])
disp(['REM: ' num2str(flowCI{2,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.001 confidence intervals with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(flowCI{3,1}(1,:))])
disp(['Whisk: ' num2str(flowCI{3,1}(2,:))])
disp(['NREM: ' num2str(flowCI{3,1}(3,:))])
disp(['REM: ' num2str(flowCI{3,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
diary off

end
