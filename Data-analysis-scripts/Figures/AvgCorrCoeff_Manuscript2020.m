function [] = AvgCorrCoeff_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Average behavior-dependent Pearson's correlation coefficients
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','NREM','REM','Whisk'};
corrCoeff_dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(180/256)];  % whisk color

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(corrCoeff_dataTypes)
            corrCoeff_dataType = corrCoeff_dataTypes{1,c};
            data.CorrCoef.(behavField).(corrCoeff_dataType).R{a,1} = AnalysisResults.(animalID).CorrCoeff.(behavField).(corrCoeff_dataType).R;
            data.CorrCoef.(behavField).(corrCoeff_dataType).meanRs(a,1) = AnalysisResults.(animalID).CorrCoeff.(behavField).(corrCoeff_dataType).meanR;
            data.CorrCoef.(behavField).animalID{a,1} = animalID;
            data.CorrCoef.(behavField).behavior{a,1} = behavField;
        end
    end
end
% concatenate the data and take mean/STD
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(corrCoeff_dataTypes)
        corrCoeff_dataType = corrCoeff_dataTypes{1,f};
        data.CorrCoef.(behavField).(corrCoeff_dataType).catR = [];
        for z = 1:length(data.CorrCoef.(behavField).(corrCoeff_dataType).R) 
            data.CorrCoef.(behavField).(corrCoeff_dataType).catR = cat(1,data.CorrCoef.(behavField).(corrCoeff_dataType).catR,data.CorrCoef.(behavField).(corrCoeff_dataType).R{z,1});
        end
        data.CorrCoef.(behavField).(corrCoeff_dataType).meanR = mean(data.CorrCoef.(behavField).(corrCoeff_dataType).meanRs,1);
        data.CorrCoef.(behavField).(corrCoeff_dataType).stdR = std(data.CorrCoef.(behavField).(corrCoeff_dataType).meanRs,0,1);
    end
end

%% statistics - linear mixed effects model
alphaConf1 = 0.001;
alphaConf2 = 0.05;
numComparisons = 3;
% HbT
HbTtableSize = cat(1,data.Rest.CBV_HbT.meanRs,data.Whisk.CBV_HbT.meanRs,data.NREM.CBV_HbT.meanRs,data.REM.CBV_HbT.meanRs);
HbTTable = table('Size',[size(HbTtableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
HbTTable.Mouse = cat(1,data.Rest.animalID,data.Whisk.animalID,data.NREM.animalID,data.REM.animalID);
HbTTable.CorrCoef = cat(1,data.Rest.CBV_HbT.meanRs,data.Whisk.CBV_HbT.meanRs,data.NREM.CBV_HbT.meanRs,data.REM.CBV_HbT.meanRs);
HbTTable.Behavior = cat(1,data.Rest.behavior,data.Whisk.behavior,data.NREM.behavior,data.REM.behavior);
HbTFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
HbTStats = fitglme(HbTTable,HbTFitFormula);
HbTCI = coefCI(HbTStats,'Alpha',(alphaConf1/numComparisons));
% gamma-band power
gammaTableSize = cat(1,data.Rest.gammaBandPower.meanRs,data.Whisk.gammaBandPower.meanRs,data.NREM.gammaBandPower.meanRs,data.REM.gammaBandPower.meanRs);
gammaTable = table('Size',[size(gammaTableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
gammaTable.Mouse = cat(1,data.Rest.animalID,data.Whisk.animalID,data.NREM.animalID,data.REM.animalID);
gammaTable.CorrCoef = cat(1,data.Rest.gammaBandPower.meanRs,data.Whisk.gammaBandPower.meanRs,data.NREM.gammaBandPower.meanRs,data.REM.gammaBandPower.meanRs);
gammaTable.Behavior = cat(1,data.Rest.behavior,data.Whisk.behavior,data.NREM.behavior,data.REM.behavior);
gammaFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
gammaStats = fitglme(gammaTable,gammaFitFormula);
gammaCI = coefCI(gammaStats,'Alpha',(alphaConf2/numComparisons));

%% summary figure(s)
summaryFigure = figure;
sgtitle('Pearson''s Correlation Coefficients')
CC_xInds = ones(1,length(animalIDs));
%% CBV HbT
p1 = subplot(2,3,1);
s1 = scatter(CC_xInds*1,data.Whisk.CBV_HbT.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.CBV_HbT.meanR,data.Whisk.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 15;
e1.CapSize = 15;
s2 = scatter(CC_xInds*2,data.Rest.CBV_HbT.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.CBV_HbT.meanR,data.Rest.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 15;
e2.CapSize = 15;
s3 = scatter(CC_xInds*3,data.NREM.CBV_HbT.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.CBV_HbT.meanR,data.NREM.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 15;
e3.CapSize = 15;
s4 = scatter(CC_xInds*4,data.REM.CBV_HbT.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.CBV_HbT.meanR,data.REM.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 15;
e4.CapSize = 15;
title('\DeltaHbT (\muM)')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
legend([s1,s2,s3,s4],'Whisking','Awake Rest','NREM','REM','Location','SouthEast')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

%% Delta-band power
p2 = subplot(2,3,2);
scatter(CC_xInds*1,data.Whisk.deltaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.deltaBandPower.meanR,data.Whisk.deltaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 15;
e1.CapSize = 15;
scatter(CC_xInds*2,data.Rest.deltaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.deltaBandPower.meanR,data.Rest.deltaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 15;
e2.CapSize = 15;
scatter(CC_xInds*3,data.NREM.deltaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.deltaBandPower.meanR,data.NREM.deltaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 15;
e3.CapSize = 15;
scatter(CC_xInds*4,data.REM.deltaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.deltaBandPower.meanR,data.REM.deltaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 15;
e4.CapSize = 15;
title('Delta-band [1-4 Hz]')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

%% Theta-band power
p3 = subplot(2,3,3);
scatter(CC_xInds*1,data.Whisk.thetaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.thetaBandPower.meanR,data.Whisk.thetaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 15;
e1.CapSize = 15;
scatter(CC_xInds*2,data.Rest.thetaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.thetaBandPower.meanR,data.Rest.thetaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 15;
e2.CapSize = 15;
scatter(CC_xInds*3,data.NREM.thetaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.thetaBandPower.meanR,data.NREM.thetaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 15;
e3.CapSize = 15;
scatter(CC_xInds*4,data.REM.thetaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.thetaBandPower.meanR,data.REM.thetaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 15;
e4.CapSize = 15;
title('Theta-band [4-10 Hz]')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

%% Alpha-band power
p4 = subplot(2,3,4);
scatter(CC_xInds*1,data.Whisk.alphaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.alphaBandPower.meanR,data.Whisk.alphaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 15;
e1.CapSize = 15;
scatter(CC_xInds*2,data.Rest.alphaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.alphaBandPower.meanR,data.Rest.alphaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 15;
e2.CapSize = 15;
scatter(CC_xInds*3,data.NREM.alphaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.alphaBandPower.meanR,data.NREM.alphaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 15;
e3.CapSize = 15;
scatter(CC_xInds*4,data.REM.alphaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.alphaBandPower.meanR,data.REM.alphaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 15;
e4.CapSize = 15;
title('Alpha-band [10-13 Hz]')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

%% Beta-band power
p5 = subplot(2,3,5);
scatter(CC_xInds*1,data.Whisk.betaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.betaBandPower.meanR,data.Whisk.betaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 15;
e1.CapSize = 15;
scatter(CC_xInds*2,data.Rest.betaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.betaBandPower.meanR,data.Rest.betaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 15;
e2.CapSize = 15;
scatter(CC_xInds*3,data.NREM.betaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.betaBandPower.meanR,data.NREM.betaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 15;
e3.CapSize = 15;
scatter(CC_xInds*4,data.REM.betaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.betaBandPower.meanR,data.REM.betaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 15;
e4.CapSize = 15;
title('Beta-band [13-30 Hz]')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')

%% Gamma-band power
p6 = subplot(2,3,6);
scatter(CC_xInds*1,data.Whisk.gammaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.gammaBandPower.meanR,data.Whisk.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 15;
e1.CapSize = 15;
scatter(CC_xInds*2,data.Rest.gammaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Rest.gammaBandPower.meanR,data.Rest.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 15;
e2.CapSize = 15;
scatter(CC_xInds*3,data.NREM.gammaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.NREM.gammaBandPower.meanR,data.NREM.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 15;
e3.CapSize = 15;
scatter(CC_xInds*4,data.REM.gammaBandPower.meanRs,100,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.REM.gammaBandPower.meanR,data.REM.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 15;
e4.CapSize = 15;
title('Gamma-band [30-100 Hz]')
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0 length(behavFields)+1])
set(gca,'box','off')
linkaxes([p1,p2,p3,p4,p5,p6],'xy')

%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Peason''s Correlation Coefficients']);
% HbT statistical diary
diary([dirpath 'Behavior_MeanHbTCorrCoef_Stats.txt'])
diary on
disp('Generalized linear mixed-effects model statistics for mean HbT corr. coef during Rest, Whisking, NREM, and REM')
disp('======================================================================================================================')
disp(HbTStats)
disp('======================================================================================================================')
disp('Alpha = 0.001 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HbTCI(1,:))])
disp(['Whisk: ' num2str(HbTCI(2,:))])
disp(['NREM: ' num2str(HbTCI(3,:))])
disp(['REM: ' num2str(HbTCI(4,:))])
diary off
% gamma statistical diary
diary([dirpath 'Behavior_MeanGammaBandCorrCoef_Stats.txt'])
diary on
disp('Generalized linear mixed-effects model statistics for mean gamma-band corr. coef during Rest, Whisking, NREM, and REM')
disp('======================================================================================================================')
disp(gammaStats)
disp('======================================================================================================================')
disp('Alpha = 0.05 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(gammaCI(1,:))])
disp(['Whisk: ' num2str(gammaCI(2,:))])
disp(['NREM: ' num2str(gammaCI(3,:))])
disp(['REM: ' num2str(gammaCI(4,:))])
diary off

end
