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
% shift the HbT paired samples
for cc = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,cc};
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        for ee = 1:length(data.(animalID).(behavField).HbTindLH)
            test.(animalID).(behavField).HbTindLH{ee,1} = data.(animalID).(behavField).HbTindLH{ee,1}(31:end);
            test.(animalID).(behavField).HbTindRH{ee,1} = data.(animalID).(behavField).HbTindRH{ee,1}(31:end);
            test.(animalID).(behavField).GamIndLH{ee,1} = data.(animalID).(behavField).GamIndLH{ee,1}(1:end-30);
            test.(animalID).(behavField).GamIndRH{ee,1} = data.(animalID).(behavField).GamIndRH{ee,1}(1:end-30);
        end
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
            data.(behavField).IndHbT = cat(2,data.(behavField).IndHbT,test.(animalID).(behavField).HbTindLH{ee,1},test.(animalID).(behavField).HbTindRH{ee,1});
            data.(behavField).IndGam = cat(2,data.(behavField).IndGam,test.(animalID).(behavField).GamIndLH{ee,1},test.(animalID).(behavField).GamIndRH{ee,1});
            data.(behavField).MeanHbT = cat(2,data.(behavField).MeanHbT,data.(animalID).(behavField).HbTmeanLH(ee,1),data.(animalID).(behavField).HbTmeanRH(ee,1));
            data.(behavField).MeanGam = cat(2,data.(behavField).MeanGam,data.(animalID).(behavField).GamMeanLH(ee,1),data.(animalID).(behavField).GamMeanRH(ee,1));
        end
    end
end
%% paired event means
figure;
scatter(data.NREM.MeanGam*100,data.NREM.MeanHbT,'MarkerFaceColor',colorB,'MarkerEdgeColor','k')
hold on
scatter(data.REM.MeanGam*100,data.REM.MeanHbT,'MarkerFaceColor',colorC,'MarkerEdgeColor','k')
scatter(data.Whisk.MeanGam*100,data.Whisk.MeanHbT,'MarkerFaceColor',colorD,'MarkerEdgeColor','k')
scatter(data.Rest.MeanGam*100,data.Rest.MeanHbT,'MarkerFaceColor',colorA,'MarkerEdgeColor','k')
ylim([-20,120])
xlim([-50,200])
title('Mean event paired samples')
xlabel('Gamma-band power (%)')
ylabel('\Delta[HbT] (\muM)')
%% paired individual samples
figure;
histogram2(data.NREM.IndGam*100,data.NREM.IndHbT,'Normalization','probability','FaceColor',colorB,'EdgeColor','k')
hold on
histogram2(data.REM.IndGam*100,data.REM.IndHbT,'Normalization','probability','FaceColor',colorC,'EdgeColor','k')
histogram2(data.Whisk.IndGam*100,data.Whisk.IndHbT,'Normalization','probability','FaceColor',colorD,'EdgeColor','k')
histogram2(data.Rest.IndGam*100,data.Rest.IndHbT,'Normalization','probability','FaceColor',colorA,'EdgeColor','k')
ylim([-20,120])
xlim([-50,200])
title('1-sec shifted (HbT) individual paired samples')
xlabel('Gamma-band power (%)')
ylabel('\Delta[HbT] (\muM)')

end
