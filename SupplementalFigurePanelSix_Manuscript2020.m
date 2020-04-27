function [] = SupplementalFigurePanelSix_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','NREM','REM'};
modelType = 'Forest';
dataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower'};
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
%% Average coherence during different behaviors 
% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        if strcmp(behavField,'Rest') == true
            for c = 1:length(dataTypes)
                dataType = dataTypes{1,c};
                data.Coherr.(behavField).(dataType).C(:,a) = (AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                data.Coherr.(behavField).(dataType).f(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(dataType).f;
                data.Coherr.(behavField).(dataType).confC(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC;
            end
        elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            for d = 1:length(dataTypes)
                dataType = dataTypes{1,d};
                data.Coherr.(behavField).(modelType).(dataType).C(:,a) = (AnalysisResults.(animalID).Coherence.(behavField).(modelType).(dataType).C);
                data.Coherr.(behavField).(modelType).(dataType).f(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(modelType).(dataType).f;
                data.Coherr.(behavField).(modelType).(dataType).confC(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(modelType).(dataType).confC;
            end
        end
    end
end
% take the mean and standard deviation of each set of signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    if strcmp(behavField,'Rest') == true
        for f = 1:length(dataTypes)
            dataType = dataTypes{1,f};
            data.Coherr.(behavField).(dataType).meanC = mean(data.Coherr.(behavField).(dataType).C,2);
            data.Coherr.(behavField).(dataType).stdC = std(data.Coherr.(behavField).(dataType).C,0,2);
            data.Coherr.(behavField).(dataType).meanf = mean(data.Coherr.(behavField).(dataType).f,2);
            data.Coherr.(behavField).(dataType).maxConfC = max(data.Coherr.(behavField).(dataType).confC);
            data.Coherr.(behavField).(dataType).maxConfC_Y = ones(length(data.Coherr.(behavField).(dataType).meanf),1)*data.Coherr.(behavField).(dataType).maxConfC;
        end
    elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
        for h = 1:length(dataTypes)
            dataType = dataTypes{1,h};
            data.Coherr.(behavField).(modelType).(dataType).meanC = mean(data.Coherr.(behavField).(modelType).(dataType).C,2);
            data.Coherr.(behavField).(modelType).(dataType).stdC = std(data.Coherr.(behavField).(modelType).(dataType).C,0,2);
            data.Coherr.(behavField).(modelType).(dataType).meanf = mean(data.Coherr.(behavField).(modelType).(dataType).f,2);
            data.Coherr.(behavField).(modelType).(dataType).maxConfC = max(data.Coherr.(behavField).(modelType).(dataType).confC);
            data.Coherr.(behavField).(modelType).(dataType).maxConfC_Y = ones(length(data.Coherr.(behavField).(modelType).(dataType).meanf),1)*data.Coherr.(behavField).(modelType).(dataType).maxConfC;
        end
    end
end
%% Power spectra during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            data.PowerSpec.(behavField).(dataType).adjLH.S(:,a) = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.S;
            data.PowerSpec.(behavField).(dataType).adjLH.f(:,a) = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.f;
            data.PowerSpec.(behavField).(dataType).adjRH.S(:,a) = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.S;
            data.PowerSpec.(behavField).(dataType).adjRH.f(:,a) = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.f;
        end
    end
end
% concatenate the data from the left and right hemispheres
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        data.PowerSpec.(behavField).(dataType).cat_S = cat(2,data.PowerSpec.(behavField).(dataType).adjLH.S,data.PowerSpec.(behavField).(dataType).adjRH.S);
        data.PowerSpec.(behavField).(dataType).cat_f = cat(2,data.PowerSpec.(behavField).(dataType).adjLH.f,data.PowerSpec.(behavField).(dataType).adjRH.f);
    end
end
% take the mean and standard deviation of each set of signals before normalizing
for h = 1:length(behavFields)
    behavField = behavFields{1,h};
    for j = 1:length(dataTypes)
        dataType = dataTypes{1,j};
        data.PowerSpec.(behavField).(dataType).preMeanCortS = mean(data.PowerSpec.(behavField).(dataType).cat_S,2);
    end
end
% normalize by the peak average power during the restng trace
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    for j = 1:length(dataTypes)
        dataType = dataTypes{1,j};
        for ee = 1:size(data.PowerSpec.(behavField).(dataType).cat_S,2)
            data.PowerSpec.(behavField).(dataType).normCat_S(:,ee) = (data.PowerSpec.(behavField).(dataType).cat_S(:,ee))*(1/max(data.PowerSpec.Rest.(dataType).preMeanCortS));
        end
    end
end
% take the mean and standard deviation of each set of signals
for h = 1:length(behavFields)
    behavField = behavFields{1,h};
    for j = 1:length(dataTypes)
        dataType = dataTypes{1,j};
        data.PowerSpec.(behavField).(dataType).meanCortS = mean(data.PowerSpec.(behavField).(dataType).normCat_S,2);
        data.PowerSpec.(behavField).(dataType).stdCortS = std(data.PowerSpec.(behavField).(dataType).normCat_S,0,2);
        data.PowerSpec.(behavField).(dataType).meanCortf = mean(data.PowerSpec.(behavField).(dataType).cat_f,2);
    end
end
%% Pearson's correlations during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            data.CorrCoef.(behavField).(dataType).R{a,1} = AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).R;
            data.CorrCoef.(behavField).(dataType).meanRs(a,1) = AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR;
            data.CorrCoef.(behavField).animalID{a,1} = animalID;
            data.CorrCoef.(behavField).behavior{a,1} = behavField;
        end
    end
end
% concatenate the data and take mean/STD
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        data.CorrCoef.(behavField).(dataType).catR = [];
        for z = 1:length(data.CorrCoef.(behavField).(dataType).R) 
            data.CorrCoef.(behavField).(dataType).catR = cat(1,data.CorrCoef.(behavField).(dataType).catR,data.CorrCoef.(behavField).(dataType).R{z,1});
        end
        data.CorrCoef.(behavField).(dataType).meanR = mean(data.CorrCoef.(behavField).(dataType).meanRs,1);
        data.CorrCoef.(behavField).(dataType).stdR = std(data.CorrCoef.(behavField).(dataType).meanRs,0,1);
    end
end
%% statistics - linear mixed effects model
numComparisons = 3;
% delta-band power
CCDelta_alphaConf = 0.05;
deltaTableSize = cat(1,data.CorrCoef.Rest.deltaBandPower.meanRs,data.CorrCoef.NREM.deltaBandPower.meanRs,data.CorrCoef.REM.deltaBandPower.meanRs);
deltaTable = table('Size',[size(deltaTableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
deltaTable.Mouse = cat(1,data.CorrCoef.Rest.animalID,data.CorrCoef.NREM.animalID,data.CorrCoef.REM.animalID);
deltaTable.CorrCoef = cat(1,data.CorrCoef.Rest.deltaBandPower.meanRs,data.CorrCoef.NREM.deltaBandPower.meanRs,data.CorrCoef.REM.deltaBandPower.meanRs);
deltaTable.Behavior = cat(1,data.CorrCoef.Rest.behavior,data.CorrCoef.NREM.behavior,data.CorrCoef.REM.behavior);
deltaFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
deltaStats = fitglme(deltaTable,deltaFitFormula);
deltaCI = coefCI(deltaStats,'Alpha',(CCDelta_alphaConf/numComparisons));
% theta-band power
CCTheta_alphaConf = 0.05;
thetaTableSize = cat(1,data.CorrCoef.Rest.thetaBandPower.meanRs,data.CorrCoef.NREM.thetaBandPower.meanRs,data.CorrCoef.REM.thetaBandPower.meanRs);
thetaTable = table('Size',[size(thetaTableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
thetaTable.Mouse = cat(1,data.CorrCoef.Rest.animalID,data.CorrCoef.NREM.animalID,data.CorrCoef.REM.animalID);
thetaTable.CorrCoef = cat(1,data.CorrCoef.Rest.thetaBandPower.meanRs,data.CorrCoef.NREM.thetaBandPower.meanRs,data.CorrCoef.REM.thetaBandPower.meanRs);
thetaTable.Behavior = cat(1,data.CorrCoef.Rest.behavior,data.CorrCoef.NREM.behavior,data.CorrCoef.REM.behavior);
thetaFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
thetaStats = fitglme(thetaTable,thetaFitFormula);
thetaCI = coefCI(thetaStats,'Alpha',(CCTheta_alphaConf/numComparisons));
% alpha-band power
CCAlpha_alphaConf = 0.05;
alphaTableSize = cat(1,data.CorrCoef.Rest.alphaBandPower.meanRs,data.CorrCoef.NREM.alphaBandPower.meanRs,data.CorrCoef.REM.alphaBandPower.meanRs);
alphaTable = table('Size',[size(alphaTableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
alphaTable.Mouse = cat(1,data.CorrCoef.Rest.animalID,data.CorrCoef.NREM.animalID,data.CorrCoef.REM.animalID);
alphaTable.CorrCoef = cat(1,data.CorrCoef.Rest.alphaBandPower.meanRs,data.CorrCoef.NREM.alphaBandPower.meanRs,data.CorrCoef.REM.alphaBandPower.meanRs);
alphaTable.Behavior = cat(1,data.CorrCoef.Rest.behavior,data.CorrCoef.NREM.behavior,data.CorrCoef.REM.behavior);
alphaFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
alphaStats = fitglme(alphaTable,alphaFitFormula);
alphaCI = coefCI(alphaStats,'Alpha',(CCAlpha_alphaConf/numComparisons));
% beta-band power
CCBeta_alphaConf = 0.05;
betaTableSize = cat(1,data.CorrCoef.Rest.betaBandPower.meanRs,data.CorrCoef.NREM.betaBandPower.meanRs,data.CorrCoef.REM.betaBandPower.meanRs);
betaTable = table('Size',[size(betaTableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
betaTable.Mouse = cat(1,data.CorrCoef.Rest.animalID,data.CorrCoef.NREM.animalID,data.CorrCoef.REM.animalID);
betaTable.CorrCoef = cat(1,data.CorrCoef.Rest.betaBandPower.meanRs,data.CorrCoef.NREM.betaBandPower.meanRs,data.CorrCoef.REM.betaBandPower.meanRs);
betaTable.Behavior = cat(1,data.CorrCoef.Rest.behavior,data.CorrCoef.NREM.behavior,data.CorrCoef.REM.behavior);
betaFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
betaStats = fitglme(betaTable,betaFitFormula);
betaCI = coefCI(betaStats,'Alpha',(CCBeta_alphaConf/numComparisons));
%% summary figure(s)
summaryFigure = figure;
sgtitle('Turner Manuscript 2020 - Supplemental Figure Panel Six')
CC_xInds = ones(1,length(animalIDs));
%% [A] Coherence between bilateral delta-band power during different arousal-states
ax1 = subplot(4,3,1);
s1 = semilogx(data.Coherr.Rest.deltaBandPower.meanf,data.Coherr.Rest.deltaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Coherr.Rest.deltaBandPower.meanf,data.Coherr.Rest.deltaBandPower.meanC + data.Coherr.Rest.deltaBandPower.stdC,'color',colorA,'LineWidth',0.5)
semilogx(data.Coherr.Rest.deltaBandPower.meanf,data.Coherr.Rest.deltaBandPower.meanC - data.Coherr.Rest.deltaBandPower.stdC,'color',colorA,'LineWidth',0.5)
s2 = semilogx(data.Coherr.NREM.(modelType).deltaBandPower.meanf,data.Coherr.NREM.(modelType).deltaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.Coherr.NREM.(modelType).deltaBandPower.meanf,data.Coherr.NREM.(modelType).deltaBandPower.meanC + data.Coherr.NREM.(modelType).deltaBandPower.stdC,'color',colorB,'LineWidth',0.5)
semilogx(data.Coherr.NREM.(modelType).deltaBandPower.meanf,data.Coherr.NREM.(modelType).deltaBandPower.meanC - data.Coherr.NREM.(modelType).deltaBandPower.stdC,'color',colorB,'LineWidth',0.5)
s3 = semilogx(data.Coherr.REM.(modelType).deltaBandPower.meanf,data.Coherr.REM.(modelType).deltaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.Coherr.REM.(modelType).deltaBandPower.meanf,data.Coherr.REM.(modelType).deltaBandPower.meanC + data.Coherr.REM.(modelType).deltaBandPower.stdC,'color',colorC,'LineWidth',0.5)
semilogx(data.Coherr.REM.(modelType).deltaBandPower.meanf,data.Coherr.REM.(modelType).deltaBandPower.meanC - data.Coherr.REM.(modelType).deltaBandPower.stdC,'color',colorC,'LineWidth',0.5)
semilogx(data.Coherr.Rest.deltaBandPower.meanf,data.Coherr.Rest.deltaBandPower.maxConfC_Y,'-','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.(modelType).deltaBandPower.meanf,data.Coherr.NREM.(modelType).deltaBandPower.maxConfC_Y,'-','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.(modelType).deltaBandPower.meanf,data.Coherr.REM.(modelType).deltaBandPower.maxConfC_Y,'-','color',colorC,'LineWidth',1);
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[A] Bilateral coherence','Delta-band [1-4 Hz]'})
legend([s1,s2,s3],'Rest','NREM','REM','Location','SouthEast')
axis square
xlim([0.1,0.5])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [B] Power spectra of delta-band power during different arousal-states
ax2 = subplot(4,3,2);
loglog(data.PowerSpec.Rest.deltaBandPower.meanCortf,data.PowerSpec.Rest.deltaBandPower.meanCortS,'color',colorA,'LineWidth',2);
hold on
loglog(data.PowerSpec.NREM.deltaBandPower.meanCortf,data.PowerSpec.NREM.deltaBandPower.meanCortS,'color',colorB,'LineWidth',2);
loglog(data.PowerSpec.REM.deltaBandPower.meanCortf,data.PowerSpec.REM.deltaBandPower.meanCortS,'color',colorC,'LineWidth',2);
title({'[B] Cortical power','Delta-band [1-4 Hz]'})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
y2 = ylim(ax2);
ylim([y2(1)/2,y2(2)*2])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [C] Pearson's correlations between bilateral delta-band power during different arousal-states
ax3 = subplot(4,3,3);
scatter(CC_xInds*1,data.CorrCoef.Rest.deltaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.CorrCoef.Rest.deltaBandPower.meanR,data.CorrCoef.Rest.deltaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(CC_xInds*2,data.CorrCoef.NREM.deltaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.CorrCoef.NREM.deltaBandPower.meanR,data.CorrCoef.NREM.deltaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(CC_xInds*3,data.CorrCoef.REM.deltaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.CorrCoef.REM.deltaBandPower.meanR,data.CorrCoef.REM.deltaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[C] Cortical Pearson''s corr. coef','Delta-band [1-4 Hz]'})
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [D] Coherence between bilateral theta-band power during different arousal-states
ax4 = subplot(4,3,4);
semilogx(data.Coherr.Rest.thetaBandPower.meanf,data.Coherr.Rest.thetaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Coherr.Rest.thetaBandPower.meanf,data.Coherr.Rest.thetaBandPower.meanC + data.Coherr.Rest.thetaBandPower.stdC,'color',colorA,'LineWidth',0.5)
semilogx(data.Coherr.Rest.thetaBandPower.meanf,data.Coherr.Rest.thetaBandPower.meanC - data.Coherr.Rest.thetaBandPower.stdC,'color',colorA,'LineWidth',0.5)
semilogx(data.Coherr.NREM.(modelType).thetaBandPower.meanf,data.Coherr.NREM.(modelType).thetaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.Coherr.NREM.(modelType).thetaBandPower.meanf,data.Coherr.NREM.(modelType).thetaBandPower.meanC + data.Coherr.NREM.(modelType).thetaBandPower.stdC,'color',colorB,'LineWidth',0.5)
semilogx(data.Coherr.NREM.(modelType).thetaBandPower.meanf,data.Coherr.NREM.(modelType).thetaBandPower.meanC - data.Coherr.NREM.(modelType).thetaBandPower.stdC,'color',colorB,'LineWidth',0.5)
semilogx(data.Coherr.REM.(modelType).thetaBandPower.meanf,data.Coherr.REM.(modelType).thetaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.Coherr.REM.(modelType).thetaBandPower.meanf,data.Coherr.REM.(modelType).thetaBandPower.meanC + data.Coherr.REM.(modelType).thetaBandPower.stdC,'color',colorC,'LineWidth',0.5)
semilogx(data.Coherr.REM.(modelType).thetaBandPower.meanf,data.Coherr.REM.(modelType).thetaBandPower.meanC - data.Coherr.REM.(modelType).thetaBandPower.stdC,'color',colorC,'LineWidth',0.5)
semilogx(data.Coherr.Rest.thetaBandPower.meanf,data.Coherr.Rest.thetaBandPower.maxConfC_Y,'-','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.(modelType).thetaBandPower.meanf,data.Coherr.NREM.(modelType).thetaBandPower.maxConfC_Y,'-','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.(modelType).thetaBandPower.meanf,data.Coherr.REM.(modelType).thetaBandPower.maxConfC_Y,'-','color',colorC,'LineWidth',1);
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[D] Bilateral coherence','Theta-band [4-10 Hz]'})
axis square
xlim([0.1,0.5])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [E] Power spectra of theta-band power during different arousal-states
ax5 = subplot(4,3,5);
loglog(data.PowerSpec.Rest.thetaBandPower.meanCortf,data.PowerSpec.Rest.thetaBandPower.meanCortS,'color',colorA,'LineWidth',2)
hold on
loglog(data.PowerSpec.NREM.thetaBandPower.meanCortf,data.PowerSpec.NREM.thetaBandPower.meanCortS,'color',colorB,'LineWidth',2);
loglog(data.PowerSpec.REM.thetaBandPower.meanCortf,data.PowerSpec.REM.thetaBandPower.meanCortS,'color',colorC,'LineWidth',2);
title({'[E] Cortical power','Theta-band [4-10 Hz]'})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
y5 = ylim(ax5);
ylim([y5(1)/2,y5(2)*2])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [F] Pearson's correlations between bilateral theta-band power during different arousal-states
ax6 = subplot(4,3,6);
scatter(CC_xInds*1,data.CorrCoef.Rest.thetaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.CorrCoef.Rest.thetaBandPower.meanR,data.CorrCoef.Rest.thetaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(CC_xInds*2,data.CorrCoef.NREM.thetaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.CorrCoef.NREM.thetaBandPower.meanR,data.CorrCoef.NREM.thetaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(CC_xInds*3,data.CorrCoef.REM.thetaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.CorrCoef.REM.thetaBandPower.meanR,data.CorrCoef.REM.thetaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[F] Cortical Pearson''s corr. coef','Theta-band [4-10 Hz]'})
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [G] Coherence between bilateral alpha-band power during different arousal-states
ax7 = subplot(4,3,7);
semilogx(data.Coherr.Rest.alphaBandPower.meanf,data.Coherr.Rest.alphaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Coherr.Rest.alphaBandPower.meanf,data.Coherr.Rest.alphaBandPower.meanC + data.Coherr.Rest.alphaBandPower.stdC,'color',colorA,'LineWidth',0.5)
semilogx(data.Coherr.Rest.alphaBandPower.meanf,data.Coherr.Rest.alphaBandPower.meanC - data.Coherr.Rest.alphaBandPower.stdC,'color',colorA,'LineWidth',0.5)
semilogx(data.Coherr.NREM.(modelType).alphaBandPower.meanf,data.Coherr.NREM.(modelType).alphaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.Coherr.NREM.(modelType).alphaBandPower.meanf,data.Coherr.NREM.(modelType).alphaBandPower.meanC + data.Coherr.NREM.(modelType).alphaBandPower.stdC,'color',colorB,'LineWidth',0.5)
semilogx(data.Coherr.NREM.(modelType).alphaBandPower.meanf,data.Coherr.NREM.(modelType).alphaBandPower.meanC - data.Coherr.NREM.(modelType).alphaBandPower.stdC,'color',colorB,'LineWidth',0.5)
semilogx(data.Coherr.REM.(modelType).alphaBandPower.meanf,data.Coherr.REM.(modelType).alphaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.Coherr.REM.(modelType).alphaBandPower.meanf,data.Coherr.REM.(modelType).alphaBandPower.meanC + data.Coherr.REM.(modelType).alphaBandPower.stdC,'color',colorC,'LineWidth',0.5)
semilogx(data.Coherr.REM.(modelType).alphaBandPower.meanf,data.Coherr.REM.(modelType).alphaBandPower.meanC - data.Coherr.REM.(modelType).alphaBandPower.stdC,'color',colorC,'LineWidth',0.5)
semilogx(data.Coherr.Rest.alphaBandPower.meanf,data.Coherr.Rest.alphaBandPower.maxConfC_Y,'-','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.(modelType).alphaBandPower.meanf,data.Coherr.NREM.(modelType).alphaBandPower.maxConfC_Y,'-','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.(modelType).alphaBandPower.meanf,data.Coherr.REM.(modelType).alphaBandPower.maxConfC_Y,'-','color',colorC,'LineWidth',1);
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[G] Bilateral coherence','Alpha-band [10-13 Hz]'})
axis square
xlim([0.1,0.5])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [H] Power spectra of alpha-band power during different arousal-states
ax8 = subplot(4,3,8);
loglog(data.PowerSpec.Rest.alphaBandPower.meanCortf,data.PowerSpec.Rest.alphaBandPower.meanCortS,'color',colorA,'LineWidth',2)
hold on
loglog(data.PowerSpec.NREM.alphaBandPower.meanCortf,data.PowerSpec.NREM.alphaBandPower.meanCortS,'color',colorB,'LineWidth',2);
loglog(data.PowerSpec.REM.alphaBandPower.meanCortf,data.PowerSpec.REM.alphaBandPower.meanCortS,'color',colorC,'LineWidth',2);
title({'[H] Cortical power','Alpha-band [10-13 Hz]'})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
y8 = ylim(ax5);
ylim([y8(1)/2,y8(2)*2])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% [I] Pearson's correlations between bilateral alpha-band power during different arousal-states
ax9 = subplot(4,3,9);
scatter(CC_xInds*1,data.CorrCoef.Rest.alphaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.CorrCoef.Rest.alphaBandPower.meanR,data.CorrCoef.Rest.alphaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(CC_xInds*2,data.CorrCoef.NREM.alphaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.CorrCoef.NREM.alphaBandPower.meanR,data.CorrCoef.NREM.alphaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(CC_xInds*3,data.CorrCoef.REM.alphaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.CorrCoef.REM.alphaBandPower.meanR,data.CorrCoef.REM.alphaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[I] Cortical Pearson''s corr. coef','Alpha-band [10-13 Hz]'})
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [J] Coherence between bilateral beta-band power during different arousal-states
ax10 = subplot(4,3,10);
semilogx(data.Coherr.Rest.betaBandPower.meanf,data.Coherr.Rest.betaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Coherr.Rest.betaBandPower.meanf,data.Coherr.Rest.betaBandPower.meanC + data.Coherr.Rest.betaBandPower.stdC,'color',colorA,'LineWidth',0.5)
semilogx(data.Coherr.Rest.betaBandPower.meanf,data.Coherr.Rest.betaBandPower.meanC - data.Coherr.Rest.betaBandPower.stdC,'color',colorA,'LineWidth',0.5)
semilogx(data.Coherr.NREM.(modelType).betaBandPower.meanf,data.Coherr.NREM.(modelType).betaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.Coherr.NREM.(modelType).betaBandPower.meanf,data.Coherr.NREM.(modelType).betaBandPower.meanC + data.Coherr.NREM.(modelType).betaBandPower.stdC,'color',colorB,'LineWidth',0.5)
semilogx(data.Coherr.NREM.(modelType).betaBandPower.meanf,data.Coherr.NREM.(modelType).betaBandPower.meanC - data.Coherr.NREM.(modelType).betaBandPower.stdC,'color',colorB,'LineWidth',0.5)
semilogx(data.Coherr.REM.(modelType).betaBandPower.meanf,data.Coherr.REM.(modelType).betaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.Coherr.REM.(modelType).betaBandPower.meanf,data.Coherr.REM.(modelType).betaBandPower.meanC + data.Coherr.REM.(modelType).betaBandPower.stdC,'color',colorC,'LineWidth',0.5)
semilogx(data.Coherr.REM.(modelType).betaBandPower.meanf,data.Coherr.REM.(modelType).betaBandPower.meanC - data.Coherr.REM.(modelType).betaBandPower.stdC,'color',colorC,'LineWidth',0.5)
semilogx(data.Coherr.Rest.betaBandPower.meanf,data.Coherr.Rest.betaBandPower.maxConfC_Y,'-','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.(modelType).betaBandPower.meanf,data.Coherr.NREM.(modelType).betaBandPower.maxConfC_Y,'-','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.(modelType).betaBandPower.meanf,data.Coherr.REM.(modelType).betaBandPower.maxConfC_Y,'-','color',colorC,'LineWidth',1);
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[J] Bilateral coherence','Beta-band [13-30 Hz]'})
axis square
xlim([0.1,0.5])
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];
%% [K] Power spectra of beta-band power during different arousal-states
ax11 = subplot(4,3,11);
loglog(data.PowerSpec.Rest.betaBandPower.meanCortf,data.PowerSpec.Rest.betaBandPower.meanCortS,'color',colorA,'LineWidth',2)
hold on
loglog(data.PowerSpec.NREM.betaBandPower.meanCortf,data.PowerSpec.NREM.betaBandPower.meanCortS,'color',colorB,'LineWidth',2);
loglog(data.PowerSpec.REM.betaBandPower.meanCortf,data.PowerSpec.REM.betaBandPower.meanCortS,'color',colorC,'LineWidth',2);
title({'[K] Cortical power','Beta-band [13-30 Hz]'})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
y11 = ylim(ax5);
ylim([y11(1)/2,y11(2)*2])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
%% [L] Pearson's correlations between bilateral beta-band power during different arousal-states
ax12 = subplot(4,3,12);
scatter(CC_xInds*1,data.CorrCoef.Rest.betaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.CorrCoef.Rest.betaBandPower.meanR,data.CorrCoef.Rest.betaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(CC_xInds*2,data.CorrCoef.NREM.betaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.CorrCoef.NREM.betaBandPower.meanR,data.CorrCoef.NREM.betaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(CC_xInds*3,data.CorrCoef.REM.betaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.CorrCoef.REM.betaBandPower.meanR,data.CorrCoef.REM.betaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
title({'[L] Cortical Pearson''s corr. coef','Beta-band [13-30 Hz]'})
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Supplemental Figure Panel 6']);
print('-painters','-dpdf','-fillpage',[dirpath 'Supplemental Figure Panel 6'])
%% statistical diary
diaryFile = [dirpath 'SupplementalFigurePanel6_Statistics.txt'];
if exist(diaryFile,'file') == true
    delete(diaryFile)
end
diary(diaryFile)
diary on
% delta statistical diary
disp('======================================================================================================================')
disp('[C] Generalized linear mixed-effects model statistics for mean delta-band corr. coef during Rest, NREM, and REM')
disp('======================================================================================================================')
disp(deltaStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.05 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(deltaCI(1,:))])
disp(['NREM: ' num2str(deltaCI(2,:))])
disp(['REM: ' num2str(deltaCI(3,:))])
% theta statistical diary
disp('======================================================================================================================')
disp('[F] Generalized linear mixed-effects model statistics for mean theta-band corr. coef during Rest, NREM, and REM')
disp('======================================================================================================================')
disp(thetaStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.05 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(thetaCI(1,:))])
disp(['NREM: ' num2str(thetaCI(2,:))])
disp(['REM: ' num2str(thetaCI(3,:))])
disp('======================================================================================================================')
% alpha statistical diary
disp('======================================================================================================================')
disp('[I] Generalized linear mixed-effects model statistics for mean alpha-band corr. coef during Rest, NREM, and REM')
disp('======================================================================================================================')
disp(alphaStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.05 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(alphaCI(1,:))])
disp(['NREM: ' num2str(alphaCI(2,:))])
disp(['REM: ' num2str(alphaCI(3,:))])
disp('======================================================================================================================')
% beta statistical diary
disp('======================================================================================================================')
disp('[L] Generalized linear mixed-effects model statistics for mean beta-band corr. coef during Rest, NREM, and REM')
disp('======================================================================================================================')
disp(betaStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.05 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(betaCI(1,:))])
disp(['NREM: ' num2str(betaCI(2,:))])
disp(['REM: ' num2str(betaCI(3,:))])
disp('======================================================================================================================')
diary off

end

