function [] = AvgCBVandHeartRate_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Calculate the average hemodynamics and heart rate during different behavioral states
%________________________________________________________________________________________________________________________

IOSanimalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
isoAnimalIDs = {'T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Whisk','Rest','NREM','REM','Iso'};
modelType = 'Forest';
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(180/256)];  % whisk color
colorE = [(0/256),(256/256),(256/256)];   % Isoflurane color

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
            data.(behavField).CBV_HbT.meanLH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjLH);
            data.(behavField).CBV_HbT.meanRH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjRH);
            data.(behavField).CBV_HbT.allLH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjLH;
            data.(behavField).CBV_HbT.allRH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjRH;
            data.(behavField).HR(a,1) = mean(AnalysisResults.(animalID).MeanHR.(behavField));
        elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            data.(behavField).CBV_HbT.meanLH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjLH);
            data.(behavField).CBV_HbT.meanRH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjRH);
            data.(behavField).CBV_HbT.allLH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjLH;
            data.(behavField).CBV_HbT.allRH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjRH;
            data.(behavField).HR(a,1) = mean(AnalysisResults.(animalID).MeanHR.(behavField));
        end
        data.(behavField).animalID{a,1} = animalID;
        data.(behavField).behavior{a,1} = behavField;
        data.(behavField).LH{a,1} = 'LH';
        data.(behavField).RH{a,1} = 'RH';
    end
end
% pull out isoflurane data
for a = 1:length(isoAnimalIDs)
    isoAnimalID = isoAnimalIDs{1,a};
    data.(behavField).CBV_HbT.meanLH(a,1) = AnalysisResults.(isoAnimalID).MeanCBV.(behavField).CBV_HbT.adjLH;
    data.(behavField).CBV_HbT.meanRH(a,1) = AnalysisResults.(isoAnimalID).MeanCBV.(behavField).CBV_HbT.adjRH;
end
% take the mean and standard deviation of each set of signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    data.(behavField).CBV_HbT.Comb = cat(1,data.(behavField).CBV_HbT.meanLH,data.(behavField).CBV_HbT.meanRH);
    data.(behavField).CBV_HbT.catAllLH = [];
    data.(behavField).CBV_HbT.catAllRH = [];
    if strcmp(behavField,'Iso') == false
        for h = 1:length(data.(behavField).CBV_HbT.allLH)
            data.(behavField).CBV_HbT.catAllLH = cat(1,data.(behavField).CBV_HbT.catAllLH,data.(behavField).CBV_HbT.allLH{h,1});
            data.(behavField).CBV_HbT.catAllRH = cat(1,data.(behavField).CBV_HbT.catAllRH,data.(behavField).CBV_HbT.allRH{h,1});
        end
        data.(behavField).CBV_HbT.allComb = cat(1,data.(behavField).CBV_HbT.catAllLH,data.(behavField).CBV_HbT.catAllRH);
        data.(behavField).meanHR = mean(data.(behavField).HR);
        data.(behavField).stdHR = std(data.(behavField).HR,0,1);
    end
    data.(behavField).CBV_HbT.meanCBV = mean(data.(behavField).CBV_HbT.Comb);
    data.(behavField).CBV_HbT.stdCBV = std(data.(behavField).CBV_HbT.Comb,0,1);
end

%% statistics - linear mixed effects model
alphaConf1 = 0.05;
alphaConf2 = 0.005;
numComparisons = 3;
% heart rate
HRtableSize = cat(1,data.Rest.HR,data.Whisk.HR,data.NREM.HR,data.REM.HR);
HRTable = table('Size',[size(HRtableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','HR','Behavior'});
HRTable.Mouse = cat(1,data.Rest.animalID,data.Whisk.animalID,data.NREM.animalID,data.REM.animalID);
HRTable.HR = cat(1,data.Rest.HR,data.Whisk.HR,data.NREM.HR,data.REM.HR);
HRTable.Behavior = cat(1,data.Rest.behavior,data.Whisk.behavior,data.NREM.behavior,data.REM.behavior);
HRFitFormula = 'HR ~ 1 + Behavior + (1|Mouse)';
HRStats = fitglme(HRTable,HRFitFormula);
HRCI = coefCI(HRStats,'Alpha',(alphaConf1/numComparisons));
% HbT
HbTtableSize = cat(1,data.Rest.CBV_HbT.meanLH,data.Rest.CBV_HbT.meanRH,data.Whisk.CBV_HbT.meanLH,data.Whisk.CBV_HbT.meanRH,...
    data.NREM.CBV_HbT.meanLH,data.NREM.CBV_HbT.meanRH,data.REM.CBV_HbT.meanLH,data.REM.CBV_HbT.meanRH);
HbTTable = table('Size',[size(HbTtableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','HbT','Behavior','Hemisphere'});
HbTTable.Mouse = cat(1,data.Rest.animalID,data.Rest.animalID,data.Whisk.animalID,data.Whisk.animalID,...
    data.NREM.animalID,data.NREM.animalID,data.REM.animalID,data.REM.animalID);
HbTTable.HbT = cat(1,data.Rest.CBV_HbT.meanLH,data.Rest.CBV_HbT.meanRH,data.Whisk.CBV_HbT.meanLH,data.Whisk.CBV_HbT.meanRH,...
    data.NREM.CBV_HbT.meanLH,data.NREM.CBV_HbT.meanRH,data.REM.CBV_HbT.meanLH,data.REM.CBV_HbT.meanRH);
HbTTable.Behavior = cat(1,data.Rest.behavior,data.Rest.behavior,data.Whisk.behavior,data.Whisk.behavior,...
    data.NREM.behavior,data.NREM.behavior,data.REM.behavior,data.REM.behavior);
HbTTable.Hemisphere = cat(1,data.Rest.LH,data.Rest.RH,data.Whisk.LH,data.Whisk.RH,...
    data.NREM.LH,data.NREM.RH,data.REM.LH,data.REM.RH);
HbTFitFormula = 'HbT ~ 1 + Behavior + (1|Mouse) + (1|Hemisphere)';
HbTStats = fitglme(HbTTable,HbTFitFormula);
HbTCI = coefCI(HbTStats,'Alpha',(alphaConf2/numComparisons));

%% summary figure(s)
summaryFigure = figure;
sgtitle('Mean Hemodynamics and Heart Rate')
HbT_xInds = ones(1,length(IOSanimalIDs)*2);
xIndsB = ones(1,length(isoAnimalIDs)*2);
%% CBV HbT
% scatter plot of mean HbT per behavior
subplot(1,3,1);
s1 = scatter(HbT_xInds*1,data.Whisk.CBV_HbT.Comb,100,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.CBV_HbT.meanCBV,data.Whisk.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 15;
e1.CapSize = 15;
s2 = scatter(HbT_xInds*2,data.Rest.CBV_HbT.Comb,100,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Rest.CBV_HbT.meanCBV,data.Rest.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 15;
e2.CapSize = 15;
s3 = scatter(HbT_xInds*3,data.NREM.CBV_HbT.Comb,100,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.NREM.CBV_HbT.meanCBV,data.NREM.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 15;
e3.CapSize = 15;
s4 = scatter(HbT_xInds*4,data.REM.CBV_HbT.Comb,100,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.REM.CBV_HbT.meanCBV,data.REM.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 15;
e4.CapSize = 15;
s5 = scatter(xIndsB*5,data.Iso.CBV_HbT.Comb,100,'MarkerEdgeColor','k','MarkerFaceColor',colorE,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.Iso.CBV_HbT.meanCBV,data.Iso.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 15;
e5.CapSize = 15;
title('Mean \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
legend([s1,s2,s3,s4,s5],'Whisking','Awake Rest','NREM','REM','Isoflurane','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')
% histogram of all events
subplot(1,3,2);
edges = -25:15:130;
[curve1] = SmoothHistogramBins_Manuscript2020(data.Whisk.CBV_HbT.allComb,edges);
[curve2] = SmoothHistogramBins_Manuscript2020(data.Rest.CBV_HbT.allComb,edges);
[curve3] = SmoothHistogramBins_Manuscript2020(data.NREM.CBV_HbT.allComb,edges);
[curve4] = SmoothHistogramBins_Manuscript2020(data.REM.CBV_HbT.allComb,edges);
before = findall(gca);
fnplt(curve1);
added = setdiff(findall(gca),before);
set(added,'Color',colorD)
hold on
before = findall(gca);
fnplt(curve2);
added = setdiff(findall(gca),before);
set(added,'Color',colorA)
before = findall(gca);
fnplt(curve3);
added = setdiff(findall(gca),before);
set(added,'Color',colorB)
before = findall(gca);
fnplt(curve4);
added = setdiff(findall(gca),before);
set(added,'Color',colorC)
title('\DeltaHbT (\muM) Distribution')
xlabel('\DeltaHbT (\muM)')
ylabel('Probability')
axis square
set(gca,'box','off')
axis tight
% scatter plot of mean heart rate per behavior
subplot(1,3,3)
HR_xInds = ones(1,length(IOSanimalIDs));
scatter(HR_xInds*1,data.Whisk.HR,100,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
hold on
e6 = errorbar(1,data.Whisk.meanHR,data.Whisk.stdHR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 15;
e6.CapSize = 15;
scatter(HR_xInds*2,data.Rest.HR,100,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
e7 = errorbar(2,data.Rest.meanHR,data.Rest.stdHR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 15;
e7.CapSize = 15;
scatter(HR_xInds*3,data.NREM.HR,100,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e8 = errorbar(3,data.NREM.meanHR,data.NREM.stdHR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 15;
e8.CapSize = 15;
scatter(HR_xInds*4,data.REM.HR,100,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e9 = errorbar(4,data.REM.meanHR,data.REM.stdHR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 15;
e9.CapSize = 15;
title('Mean Heart Rate')
ylabel('Freq (Hz)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields)])
set(gca,'box','off')
% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Hemodynamics and Heart Rate']);
% heart rate statistical diary
diary([dirpath 'Behavior_MeanHeartRate_Stats.txt'])
diary on
disp('Generalized linear mixed-effects model statistics for mean heart rate during Rest, Whisking, NREM, and REM')
disp('======================================================================================================================')
disp(HRStats)
disp('======================================================================================================================')
disp('Alpha = 0.05 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HRCI(1,:))])
disp(['Whisk: ' num2str(HRCI(2,:))])
disp(['NREM: ' num2str(HRCI(3,:))])
disp(['REM: ' num2str(HRCI(4,:))])
diary off
% HbT statistical diary
diary([dirpath 'Behavior_MeanHbT_Stats.txt'])
diary on
disp('Generalized linear mixed-effects model statistics for mean HbT during Rest, Whisking, NREM, and REM')
disp('======================================================================================================================')
disp(HbTStats)
disp('======================================================================================================================')
disp('Alpha = 0.005 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HbTCI(1,:))])
disp(['Whisk: ' num2str(HbTCI(2,:))])
disp(['NREM: ' num2str(HbTCI(3,:))])
disp(['REM: ' num2str(HbTCI(4,:))])
diary off

end

