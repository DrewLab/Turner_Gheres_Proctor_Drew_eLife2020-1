function [] = FigurePanelFive_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

IOSanimalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
peakVesselAnimalIDs = {'T115','T116','T117','T118','T125','T126'};
dopplerAnimalIDs = {'T109','T110','T111','T119','T120','T121'};
modelType = 'Forest';
numComparisons = 3;
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(180/256)];  % whisk color
%% Mean HbT and heart rate comparison between behaviors
% cd through each animal's directory and extract the appropriate analysis results
IOS_behavFields = {'Rest','Whisk','NREM','REM'};
for a = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,a};
    for b = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,b};
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
            data.(behavField).CBV_HbT.meanLH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjLH);
            data.(behavField).CBV_HbT.meanRH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjRH);
            data.(behavField).CBV_HbT.allLH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjLH;
            data.(behavField).CBV_HbT.allRH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjRH;
        elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            data.(behavField).CBV_HbT.meanLH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjLH);
            data.(behavField).CBV_HbT.meanRH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjRH);
            data.(behavField).CBV_HbT.allLH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjLH;
            data.(behavField).CBV_HbT.allRH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjRH;
        end
        data.(behavField).CBV_HbT.animalID{a,1} = animalID;
        data.(behavField).CBV_HbT.behavior{a,1} = behavField;
        data.(behavField).CBV_HbT.LH{a,1} = 'LH';
        data.(behavField).CBV_HbT.RH{a,1} = 'RH';
    end
end
% take the mean and standard deviation of each set of signals
for e = 1:length(IOS_behavFields)
    behavField = IOS_behavFields{1,e};
    data.(behavField).CBV_HbT.Comb = cat(1,data.(behavField).CBV_HbT.meanLH,data.(behavField).CBV_HbT.meanRH);
    data.(behavField).CBV_HbT.catAllLH = [];
    data.(behavField).CBV_HbT.catAllRH = [];
    for h = 1:length(data.(behavField).CBV_HbT.allLH)
        data.(behavField).CBV_HbT.catAllLH = cat(1,data.(behavField).CBV_HbT.catAllLH,data.(behavField).CBV_HbT.allLH{h,1});
        data.(behavField).CBV_HbT.catAllRH = cat(1,data.(behavField).CBV_HbT.catAllRH,data.(behavField).CBV_HbT.allRH{h,1});
    end
    data.(behavField).CBV_HbT.allComb = cat(1,data.(behavField).CBV_HbT.catAllLH,data.(behavField).CBV_HbT.catAllRH);
    data.(behavField).CBV_HbT.meanCBV = mean(data.(behavField).CBV_HbT.Comb);
    data.(behavField).CBV_HbT.stdCBV = std(data.(behavField).CBV_HbT.Comb,0,1);
end
% statistics - linear mixed effects model
% HbT
HbT_alphaConf = 0.005;
HbTtableSize = cat(1,data.Rest.CBV_HbT.meanLH,data.Rest.CBV_HbT.meanRH,data.Whisk.CBV_HbT.meanLH,data.Whisk.CBV_HbT.meanRH,...
    data.NREM.CBV_HbT.meanLH,data.NREM.CBV_HbT.meanRH,data.REM.CBV_HbT.meanLH,data.REM.CBV_HbT.meanRH);
HbTTable = table('Size',[size(HbTtableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','HbT','Behavior','Hemisphere'});
HbTTable.Mouse = cat(1,data.Rest.CBV_HbT.animalID,data.Rest.CBV_HbT.animalID,data.Whisk.CBV_HbT.animalID,data.Whisk.CBV_HbT.animalID,...
    data.NREM.CBV_HbT.animalID,data.NREM.CBV_HbT.animalID,data.REM.CBV_HbT.animalID,data.REM.CBV_HbT.animalID);
HbTTable.HbT = cat(1,data.Rest.CBV_HbT.meanLH,data.Rest.CBV_HbT.meanRH,data.Whisk.CBV_HbT.meanLH,data.Whisk.CBV_HbT.meanRH,...
    data.NREM.CBV_HbT.meanLH,data.NREM.CBV_HbT.meanRH,data.REM.CBV_HbT.meanLH,data.REM.CBV_HbT.meanRH);
HbTTable.Behavior = cat(1,data.Rest.CBV_HbT.behavior,data.Rest.CBV_HbT.behavior,data.Whisk.CBV_HbT.behavior,data.Whisk.CBV_HbT.behavior,...
    data.NREM.CBV_HbT.behavior,data.NREM.CBV_HbT.behavior,data.REM.CBV_HbT.behavior,data.REM.CBV_HbT.behavior);
HbTTable.Hemisphere = cat(1,data.Rest.CBV_HbT.LH,data.Rest.CBV_HbT.RH,data.Whisk.CBV_HbT.LH,data.Whisk.CBV_HbT.RH,...
    data.NREM.CBV_HbT.LH,data.NREM.CBV_HbT.RH,data.REM.CBV_HbT.LH,data.REM.CBV_HbT.RH);
HbTFitFormula = 'HbT ~ 1 + Behavior + (1|Mouse) + (1|Hemisphere)';
HbTStats = fitglme(HbTTable,HbTFitFormula);
HbTCI = coefCI(HbTStats,'Alpha',(HbT_alphaConf/numComparisons));
%% Peak vessel diameter comparison between behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Rest.PeakVD.data = []; data.Rest.PeakVD.indData = []; data.Rest.PeakVD.animalID = {}; data.Rest.PeakVD.behavior = {}; data.Rest.PeakVD.vID = {};
data.Whisk.PeakVD.data = []; data.Whisk.PeakVD.indData = []; data.Whisk.PeakVD.animalID = {}; data.Whisk.PeakVD.behavior = {}; data.Whisk.PeakVD.vID = {};
data.NREM.PeakVD.data = []; data.NREM.PeakVD.indData = []; data.NREM.PeakVD.animalID = {}; data.NREM.PeakVD.behavior = {}; data.NREM.PeakVD.vID = {};
data.REM.PeakVD.data = []; data.REM.PeakVD.indData = []; data.REM.PeakVD.animalID = {}; data.REM.PeakVD.behavior = {}; data.REM.PeakVD.vID = {};
for aa = 1:length(peakVesselAnimalIDs)
    animalID = peakVesselAnimalIDs{1,aa};
    PeakVD_behavFields = fieldnames(AnalysisResults.(animalID).MeanVesselDiameter);
    for bb = 1:length(PeakVD_behavFields)
        behavField = PeakVD_behavFields{bb,1};
        vesselIDs = fieldnames(AnalysisResults.(animalID).MeanVesselDiameter.(behavField));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            if strcmp(vesselID(1),'V') == false
                data.(behavField).PeakVD.data = vertcat(data.(behavField).PeakVD.data,AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vesselID).mean);
                data.(behavField).PeakVD.indData = vertcat(data.(behavField).PeakVD.indData,AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vesselID).indEvents);
                data.(behavField).PeakVD.animalID = vertcat(data.(behavField).PeakVD.animalID,animalID);
                data.(behavField).PeakVD.behavior = vertcat(data.(behavField).PeakVD.behavior,behavField);
                data.(behavField).PeakVD.vID = vertcat(data.(behavField).PeakVD.vID,vesselID);
            end
        end
    end
end
% take the average of the vessels for each behavior
for dd = 1:length(PeakVD_behavFields)
    behavField = PeakVD_behavFields{dd,1};
    data.(behavField).PeakVD.mean = mean(data.(behavField).PeakVD.data);
    data.(behavField).PeakVD.StD = std(data.(behavField).PeakVD.data,0,1);
end
% statistics - linear mixed effects model
peakVD_alphaConf = 0.001;
tableSize = cat(1,data.Rest.PeakVD.animalID,data.Whisk.PeakVD.animalID,data.NREM.PeakVD.animalID,data.REM.PeakVD.animalID);
vesselDiameterTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Diameter','Behavior','Vessel'});
vesselDiameterTable.Mouse = cat(1,data.Rest.PeakVD.animalID,data.Whisk.PeakVD.animalID,data.NREM.PeakVD.animalID,data.REM.PeakVD.animalID);
vesselDiameterTable.Diameter = cat(1,data.Rest.PeakVD.data,data.Whisk.PeakVD.data,data.NREM.PeakVD.data,data.REM.PeakVD.data);
vesselDiameterTable.Behavior = cat(1,data.Rest.PeakVD.behavior,data.Whisk.PeakVD.behavior,data.NREM.PeakVD.behavior,data.REM.PeakVD.behavior);
vesselDiameterTable.Vessel = cat(1,data.Rest.PeakVD.vID,data.Whisk.PeakVD.vID,data.NREM.PeakVD.vID,data.REM.PeakVD.vID);
vesselFitFormula = 'Diameter ~ 1 + Behavior + (1|Mouse) + (1|Vessel)';
vesselStats = fitglme(vesselDiameterTable,vesselFitFormula);
vesselCI = coefCI(vesselStats,'Alpha',(peakVD_alphaConf/numComparisons));
%% LDf comparison between behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Rest.LDFlow.indFlowMeans = []; data.Whisk.LDFlow.indFlowMeans = [];
data.NREM.LDFlow.indFlowMeans = []; data.REM.LDFlow.indFlowMeans = [];
for a = 1:length(dopplerAnimalIDs)
    animalID = dopplerAnimalIDs{1,a};
    for b = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,b};
        data.(behavField).LDFlow.indFlowMeans = vertcat(data.(behavField).LDFlow.indFlowMeans,AnalysisResults.(animalID).LDFlow.(behavField));
        data.(behavField).LDFlow.flowMeans(a,1) = mean(AnalysisResults.(animalID).LDFlow.(behavField));
        data.(behavField).LDFlow.animalID{a,1} = animalID;
        data.(behavField).LDFlow.behavior{a,1} = behavField;
    end
end
% take average of the flow for each behavior
for c = 1:length(IOS_behavFields)
    behavField = IOS_behavFields{1,c};
    data.(behavField).LDFlow.behavMean = mean(data.(behavField).LDFlow.flowMeans);
    data.(behavField).LDFlow.behavStD = std(data.(behavField).LDFlow.flowMeans,0,1);
end
% statistics - linear mixed effects model
LDflow_alphaConf = 0.001;
tableSize = cat(1,data.Rest.LDFlow.flowMeans,data.Whisk.LDFlow.flowMeans,data.NREM.LDFlow.flowMeans,data.REM.LDFlow.flowMeans);
flowTable = table('Size',[size(tableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','Flow','Behavior'});
flowTable.Mouse = cat(1,data.Rest.LDFlow.animalID,data.Whisk.LDFlow.animalID,data.NREM.LDFlow.animalID,data.REM.LDFlow.animalID);
flowTable.Flow = cat(1,data.Rest.LDFlow.flowMeans,data.Whisk.LDFlow.flowMeans,data.NREM.LDFlow.flowMeans,data.REM.LDFlow.flowMeans);
flowTable.Behavior = cat(1,data.Rest.LDFlow.behavior,data.Whisk.LDFlow.behavior,data.NREM.LDFlow.behavior,data.REM.LDFlow.behavior);
flowFitFormula = 'Flow ~ 1 + Behavior + (1|Mouse)';
flowStats = fitglme(flowTable,flowFitFormula);
flowCI = coefCI(flowStats,'Alpha',(LDflow_alphaConf/numComparisons));
%% Pixel panel 5
summaryFigure = figure;
sgtitle('Figure panel 5 - Turner Manuscript 2020')
%% [A] Mean HbT during different behaviors
ax1 = subplot(2,3,1);
HbT_xInds = ones(1,length(IOSanimalIDs)*2);
s1 = scatter(HbT_xInds*1,data.Rest.CBV_HbT.Comb,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.CBV_HbT.meanCBV,data.Rest.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(HbT_xInds*2,data.Whisk.CBV_HbT.Comb,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Whisk.CBV_HbT.meanCBV,data.Whisk.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(HbT_xInds*3,data.NREM.CBV_HbT.Comb,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.NREM.CBV_HbT.meanCBV,data.NREM.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(HbT_xInds*4,data.REM.CBV_HbT.Comb,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.REM.CBV_HbT.meanCBV,data.REM.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[A] Mean \DeltaHbT (\muM)','during arousal-states',''})
ylabel('\DeltaHbT (\muM)')
legend([s1,s2,s3,s4],'Awake Rest','Whisking','NREM','REM','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(IOS_behavFields) + 1])
ylim([-10,100])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [B] Peak vessel diameter during different behaviors
ax4 = subplot(2,3,2);
peakVD_xIndsRest = ones(1,length(data.Rest.PeakVD.data));
peakVD_xIndsWhisk = ones(1,length(data.Whisk.PeakVD.data));
peakVD_xIndsNREM = ones(1,length(data.NREM.PeakVD.data));
peakVD_xIndsREM = ones(1,length(data.REM.PeakVD.data));
scatter(peakVD_xIndsRest*1,data.Rest.PeakVD.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.PeakVD.mean,data.Rest.PeakVD.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(peakVD_xIndsWhisk*2,data.Whisk.PeakVD.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Whisk.PeakVD.mean,data.Whisk.PeakVD.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(peakVD_xIndsNREM*3,data.NREM.PeakVD.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.NREM.PeakVD.mean,data.NREM.PeakVD.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(peakVD_xIndsREM*4,data.REM.PeakVD.data,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.REM.PeakVD.mean,data.REM.PeakVD.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[B] Peak \DeltaD/D (%)','during arousal-states',''})
ylabel('\DeltaD/D (%)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(IOS_behavFields) + 1])
ylim([-5,70])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [C] Mean arousal-state LDF
ax3 = subplot(2,3,3);
LDF_xInds = ones(1,length(dopplerAnimalIDs));
scatter(LDF_xInds*1,data.Rest.LDFlow.flowMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.LDFlow.behavMean,data.Rest.LDFlow.behavStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(LDF_xInds*2,data.Whisk.LDFlow.flowMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Whisk.LDFlow.behavMean,data.Whisk.LDFlow.behavStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(LDF_xInds*3,data.NREM.LDFlow.flowMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.NREM.LDFlow.behavMean,data.NREM.LDFlow.behavStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(LDF_xInds*4,data.REM.LDFlow.flowMeans,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.REM.LDFlow.behavMean,data.REM.LDFlow.behavStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[F] Mean LDF','during arousal-states',''})
ylabel('/DeltaQ/Q (%)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(IOS_behavFields) + 1])
ylim([-10,70])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [D] Mean HbT distribution during different behaviors
ax4 = subplot(2,3,4);
edges = -25:15:130;
[curve1] = SmoothHistogramBins_Manuscript2020(data.Rest.CBV_HbT.allComb,edges);
[curve2] = SmoothHistogramBins_Manuscript2020(data.Whisk.CBV_HbT.allComb,edges);
[curve3] = SmoothHistogramBins_Manuscript2020(data.NREM.CBV_HbT.allComb,edges);
[curve4] = SmoothHistogramBins_Manuscript2020(data.REM.CBV_HbT.allComb,edges);
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
title({'[D] Mean \DeltaHbT (\muM)','arousal-state distribution',''})
xlabel('\DeltaHbT (\muM)')
ylabel('Probability')
axis square
set(gca,'box','off')
axis tight
ax4.TickLength = [0.03,0.03];
%% [E] Peak vessel diameter arousal-state vessel distribution
ax5 = subplot(2,3,5);
edges = -5:10:80;
[curve1] = SmoothHistogramBins_Manuscript2020(data.Rest.PeakVD.indData,edges);
[curve2] = SmoothHistogramBins_Manuscript2020(data.Whisk.PeakVD.indData,edges);
[curve3] = SmoothHistogramBins_Manuscript2020(data.NREM.PeakVD.indData,edges);
[curve4] = SmoothHistogramBins_Manuscript2020(data.REM.PeakVD.indData,edges);
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
title({'[E] Peak \DeltaD/D (%)','arousal-state distribution',''})
xlabel('Peak \DeltaD/D (%)')
ylabel('Probability')
axis square
axis tight
y1 = ylim;
ylim([0,y1(2)])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [F] Peak vessel diameter arousal-state vessel distribution
ax6 = subplot(2,3,6);
edges = -25:10:75;
[curve1] = SmoothHistogramBins_Manuscript2020(data.Rest.LDFlow.indFlowMeans,edges);
[curve2] = SmoothHistogramBins_Manuscript2020(data.Whisk.LDFlow.indFlowMeans,edges);
[curve3] = SmoothHistogramBins_Manuscript2020(data.NREM.LDFlow.indFlowMeans,edges);
[curve4] = SmoothHistogramBins_Manuscript2020(data.REM.LDFlow.indFlowMeans,edges);
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
title({'[F] Mean LDF','arousal-state distribution',''})
xlabel('/DeltaQ/Q (%)')
ylabel('Probability')
axis square
axis tight
y1 = ylim;
ylim([0,y1(2)])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
set(summaryFigure,'PaperPositionMode','auto');
savefig(summaryFigure,[dirpath 'Figure Panel 5']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Figure Panel 5'])
%% statistical diary
diaryFile = [dirpath 'FigurePanel5_Statistics.txt'];
if exist(diaryFile,'file') == true
    delete(diaryFile)
end
diary(diaryFile)
diary on
% HbT statistical diary
disp('======================================================================================================================')
disp('[A] Generalized linear mixed-effects model statistics for mean HbT during Rest, Whisking, NREM, and REM')
disp('======================================================================================================================')
disp(HbTStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.005 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HbTCI(1,:))])
disp(['Whisk: ' num2str(HbTCI(2,:))])
disp(['NREM: ' num2str(HbTCI(3,:))])
disp(['REM: ' num2str(HbTCI(4,:))])
% Peak vessel diameter statistical diary
disp('======================================================================================================================')
disp('[B] Generalized linear mixed-effects model statistics for mean peak vessel diameter during Rest, Whisking, NREM, and REM')
disp('======================================================================================================================')
disp(vesselStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.001 confidence intervals with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(vesselCI(1,:))])
disp(['Whisk: ' num2str(vesselCI(2,:))])
disp(['NREM: ' num2str(vesselCI(3,:))])
disp(['REM: ' num2str(vesselCI(4,:))])
% LDF flow statistical diary
disp('======================================================================================================================')
disp('[C] Generalized linear mixed-effects model statistics for mean doppler flow during Rest, Whisking, NREM, and REM')
disp('======================================================================================================================')
disp(flowStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.05 confidence intervals with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(flowCI(1,:))])
disp(['Whisk: ' num2str(flowCI(2,:))])
disp(['NREM: ' num2str(flowCI(3,:))])
disp(['REM: ' num2str(flowCI(4,:))])
disp('======================================================================================================================')
diary off

end
