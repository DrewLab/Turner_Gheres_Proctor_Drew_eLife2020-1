function [] = Fig5_testManuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
colorA = [(0/256),(166/256),(81/256)];   % rest color
colorB = [(191/256),(0/256),(255/256)];   % NREM color
colorC = [(254/256),(139/256),(0/256)];   % REM color
%% Mean HbT comparison between behaviors
% pre-allocate the date for each day
behavFields = {'Rest','NREM','REM'};
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    for bb = 1:length(behavFields)
        behavField = behavFields{1,bb};
        for cc = 1:length(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH)
            data.HbT.(animalID).(behavField).maxLH(cc,1) = max(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{cc,1});
            data.HbT.(animalID).(behavField).maxRH(cc,1) = max(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{cc,1});
            data.HbT.(animalID).(behavField).p2pLH(cc,1) = abs(max(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{cc,1}));
            data.HbT.(animalID).(behavField).p2pRH(cc,1) = abs(max(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{cc,1})) + abs(min(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{cc,1}));
        end
    end
end
% 
for dd = 1:length(animalIDs)
    animalID = animalIDs{1,dd};
    for ee = 1:length(behavFields)
        behavField = behavFields{1,ee};
        if isfield(data.HbT,behavField) == false
            data.HbT.(behavField).catMax = [];
            data.HbT.(behavField).catP2P = [];
        end
        data.HbT.(behavField).catMax = cat(1,data.HbT.(behavField).catMax,mean(data.HbT.(animalID).(behavField).maxLH),mean(data.HbT.(animalID).(behavField).maxRH));
        data.HbT.(behavField).catP2P = cat(1,data.HbT.(behavField).catP2P,mean(data.HbT.(animalID).(behavField).p2pLH),mean(data.HbT.(animalID).(behavField).p2pRH));
    end
end
%
for ff = 1:length(behavFields)
    behavField = behavFields{1,ff};
    data.HbT.(behavField).meanMax = mean(data.HbT.(behavField).catMax,1);
    data.HbT.(behavField).stdMax = std(data.HbT.(behavField).catMax,0,1);
    data.HbT.(behavField).meanP2P = mean(data.HbT.(behavField).catP2P,1);
    data.HbT.(behavField).stdP2P = std(data.HbT.(behavField).catP2P,0,1);
end
%% text values
round(data.HbT.NREM.meanP2P,1)
round(data.HbT.NREM.stdP2P,1)
round(data.HbT.Rest.meanMax,1)
round(data.HbT.Rest.stdMax,1)
round(data.HbT.REM.meanMax,1)
round(data.HbT.REM.stdMax,1)
%% [5c] Mean vessel diameter during different behaviors
figure
ax1 = subplot(1,2,1);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.HbT.Rest.catP2P,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.HbT.Rest.meanP2P,data.HbT.Rest.stdP2P,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.HbT.NREM.catP2P,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(2,data.HbT.NREM.meanP2P,data.HbT.NREM.stdP2P,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*3,data.HbT.REM.catP2P,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(3,data.HbT.REM.meanP2P,data.HbT.REM.stdP2P,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[] Mean Peak-to-Peak \DeltaHbT (\muM)','during arousal-states',''})
ylabel('\DeltaHbT (\muM)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
% ylim([-10,80])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [5a bottom] Mean HbT distribution during different behaviors
ax2 = subplot(1,2,2);
xInds = ones(1,length(animalIDs)*2);
scatter(xInds*1,data.HbT.Rest.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.HbT.Rest.meanMax,data.HbT.Rest.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.HbT.NREM.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(2,data.HbT.NREM.meanMax,data.HbT.NREM.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*3,data.HbT.REM.catMax,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(3,data.HbT.REM.meanMax,data.HbT.REM.stdMax,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
title({'[] Mean Peak \DeltaHbT (\muM)','during arousal-states',''})
ylabel('\DeltaHbT (\muM)')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
% ylim([-10,80])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
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
