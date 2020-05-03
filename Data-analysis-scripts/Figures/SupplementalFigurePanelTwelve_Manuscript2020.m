function [] = SupplementalFigurePanelTwelve_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose:
%________________________________________________________________________________________________________________________

%% valid animals, behaviors, data types, and associated colors
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Contra','Whisk','Rest','NREM','REM'};
neuralBands = {'gammaBandPower','muaPower'};
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(180/256)];  % whisk color
colorE = [(256/256),(28/256),(207/256)];  % contra color
%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for c = 1:length(neuralBands)
        neuralBand = neuralBands{1,c};
        for b = 1:length(behavFields)
            behavior = behavFields{1,b};
                % Gamma function kernels
                data.(neuralBand).(behavior).adjLH.gammaTimeVec{a,1} = AnalysisResults.(animalID).HRFs.(neuralBand).adjLH.(behavior).gammaTimeVec;
                data.(neuralBand).(behavior).adjRH.gammaTimeVec{a,1}  = AnalysisResults.(animalID).HRFs.(neuralBand).adjLH.(behavior).gammaTimeVec;
                data.(neuralBand).(behavior).adjLH.gammaFunc{a,1} = AnalysisResults.(animalID).HRFs.(neuralBand).adjRH.(behavior).gammaFunc;
                data.(neuralBand).(behavior).adjRH.gammaFunc{a,1} = AnalysisResults.(animalID).HRFs.(neuralBand).adjRH.(behavior).gammaFunc;
                % Impulse response function kernels
                data.(neuralBand).(behavior).adjLH.IRtimeVec{a,1} = AnalysisResults.(animalID).HRFs.(neuralBand).adjLH.(behavior).IRtimeVec;
                data.(neuralBand).(behavior).adjRH.IRtimeVec{a,1}  = AnalysisResults.(animalID).HRFs.(neuralBand).adjLH.(behavior).IRtimeVec;
                data.(neuralBand).(behavior).adjLH.IR{a,1} = AnalysisResults.(animalID).HRFs.(neuralBand).adjRH.(behavior).IR;
                data.(neuralBand).(behavior).adjRH.IR{a,1} = AnalysisResults.(animalID).HRFs.(neuralBand).adjRH.(behavior).IR;
                % Behavior-derived R2 predictions
                data.(neuralBand).(behavior).adjLH.medR2(a,1) = AnalysisResults.(animalID).HRFs.(neuralBand).adjLH.(behavior).medIndR2;
                data.(neuralBand).(behavior).adjRH.medR2(a,1) = AnalysisResults.(animalID).HRFs.(neuralBand).adjRH.(behavior).medIndR2;
                data.(neuralBand).(behavior).adjLH.meanR2(a,1) = AnalysisResults.(animalID).HRFs.(neuralBand).adjLH.(behavior).meanR2;
                data.(neuralBand).(behavior).adjRH.meanR2(a,1) = AnalysisResults.(animalID).HRFs.(neuralBand).adjRH.(behavior).meanR2;
        end
    end
end
%% take the mean and standard deviation of each set of signals
for c = 1:length(neuralBands)
    neuralBand = neuralBands{1,c};
    for b = 1:length(behavFields)
        behavior = behavFields{1,b};
        % average IR functions
        data.(neuralBand).(behavior).Comb.IR = cat(1,cell2mat(data.(neuralBand).(behavior).adjLH.IR),cell2mat(data.(neuralBand).(behavior).adjRH.IR));
        data.(neuralBand).(behavior).Comb.IRtimeVec = cat(1,cell2mat(data.(neuralBand).(behavior).adjLH.IRtimeVec),cell2mat(data.(neuralBand).(behavior).adjRH.IRtimeVec));
        data.(neuralBand).(behavior).meanIR = mean(data.(neuralBand).(behavior).Comb.IR,1);
        data.(neuralBand).(behavior).stdIR = std(data.(neuralBand).(behavior).Comb.IR,0,1);
        data.(neuralBand).(behavior).meanIRtimeVec = mean(data.(neuralBand).(behavior).Comb.IRtimeVec);
        % average gamme functions
        data.(neuralBand).(behavior).Comb.Gamma = cat(1,cell2mat(data.(neuralBand).(behavior).adjLH.gammaFunc),cell2mat(data.(neuralBand).(behavior).adjRH.gammaFunc));
        data.(neuralBand).(behavior).Comb.gammaTimeVec = cat(1,cell2mat(data.(neuralBand).(behavior).adjLH.gammaTimeVec),cell2mat(data.(neuralBand).(behavior).adjRH.gammaTimeVec));
        data.(neuralBand).(behavior).meanGamma = mean(data.(neuralBand).(behavior).Comb.Gamma,1);
        data.(neuralBand).(behavior).stdGamma = std(data.(neuralBand).(behavior).Comb.Gamma,0,1);
        data.(neuralBand).(behavior).meanGammaTimeVec = mean(data.(neuralBand).(behavior).Comb.gammaTimeVec);
        % average med/mean R2 values
        data.(neuralBand).(behavior).Comb.medR2 = cat(1,data.(neuralBand).(behavior).adjLH.medR2,data.(neuralBand).(behavior).adjRH.medR2);
        data.(neuralBand).(behavior).Comb.meanR2 = cat(1,data.(neuralBand).(behavior).adjLH.meanR2,data.(neuralBand).(behavior).adjRH.meanR2);
        data.(neuralBand).(behavior).avgMedR2 = mean(data.(neuralBand).(behavior).Comb.medR2);
        data.(neuralBand).(behavior).stdMedR2 = std(data.(neuralBand).(behavior).Comb.medR2,0,1);
        data.(neuralBand).(behavior).avgMeanR2 = mean(data.(neuralBand).(behavior).Comb.meanR2);
        data.(neuralBand).(behavior).stdMeanR2 = std(data.(neuralBand).(behavior).Comb.meanR2,0,1);
    end
end
%% Supplemental Figure Panel 12
summaryFigure = figure;
sgtitle('Supplemental Figure Panel 12 - Turner Manuscript 2020')
xInds = ones(1,length(animalIDs)*2);
%% [A] Gamma-band-derived IR analytic
ax1 = subplot(4,2,1);
p1 = plot(data.gammaBandPower.Contra.meanIRtimeVec,data.gammaBandPower.Contra.meanIR,'color',colorE,'LineWidth',2);
hold on
p2 = plot(data.gammaBandPower.Whisk.meanIRtimeVec,data.gammaBandPower.Whisk.meanIR,'color',colorD,'LineWidth',2);
p3 = plot(data.gammaBandPower.Rest.meanIRtimeVec,data.gammaBandPower.Rest.meanIR,'color',colorA,'LineWidth',2);
p4 = plot(data.gammaBandPower.NREM.meanIRtimeVec,data.gammaBandPower.NREM.meanIR,'color',colorB,'LineWidth',2);
p5 = plot(data.gammaBandPower.REM.meanIRtimeVec,data.gammaBandPower.REM.meanIR,'color',colorC,'LineWidth',2);
title({'[A] Gamma-band [30-100 Hz]','Impulse response function'})
xlabel('Time (s)')
ylabel({'Gamma-band [30-100 Hz]';'HRF amplitude (A.U.)'})
legend([p1,p2,p3,p4,p5],'Sensory-evoked','Whisk-evoked','Rest','NREM','REM','Location','NorthEast')
set(gca,'box','off')
xlim([0,5])
ax1.TickLength = [0.03,0.03];
%% [B] Gamma-band-derived gamma function
ax2 = subplot(4,2,2);
plot(data.gammaBandPower.Contra.meanGammaTimeVec,data.gammaBandPower.Contra.meanGamma,'color',colorE,'LineWidth',2)
hold on
plot(data.gammaBandPower.Whisk.meanGammaTimeVec,data.gammaBandPower.Whisk.meanGamma,'color',colorD,'LineWidth',2)
plot(data.gammaBandPower.Rest.meanGammaTimeVec,data.gammaBandPower.Rest.meanGamma,'color',colorA,'LineWidth',2)
plot(data.gammaBandPower.NREM.meanGammaTimeVec,data.gammaBandPower.NREM.meanGamma,'color',colorB,'LineWidth',2)
plot(data.gammaBandPower.REM.meanGammaTimeVec,data.gammaBandPower.REM.meanGamma,'color',colorC,'LineWidth',2)
title({'[B] Gamma-band [30-100 Hz]','Gamma distribution function'})
xlabel('Time (s)')
ylabel({'Gamma-band [30-100 Hz]';'HRF amplitude (A.U.)'})
set(gca,'box','off')
xlim([0,5])
ax2.TickLength = [0.03,0.03];
%% [C] Gamma-band-derived median R2 values
ax3 = subplot(4,2,3);
scatter(xInds*1,data.gammaBandPower.Contra.Comb.medR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorE,'jitter','on','jitterAmount',0.25)
hold on
e1 = errorbar(1,data.gammaBandPower.Contra.avgMedR2,data.gammaBandPower.Contra.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.gammaBandPower.Whisk.Comb.medR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25)
e2 = errorbar(2,data.gammaBandPower.Whisk.avgMedR2,data.gammaBandPower.Whisk.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.gammaBandPower.Rest.Comb.medR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25)
e3 = errorbar(3,data.gammaBandPower.Rest.avgMedR2,data.gammaBandPower.Rest.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.gammaBandPower.NREM.Comb.medR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25)
e4 = errorbar(4,data.gammaBandPower.NREM.avgMedR2,data.gammaBandPower.NREM.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.gammaBandPower.REM.Comb.medR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25)
e5 = errorbar(5,data.gammaBandPower.REM.avgMedR2,data.gammaBandPower.REM.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title({'[C] Gamma-band [30-100 Hz]';'Median R^2 kernel predictions'})
ylabel('R^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [D] Gamma-band-derived mean R2 values
ax4 = subplot(4,2,4);
scatter(xInds*1,data.gammaBandPower.Contra.Comb.meanR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorE,'jitter','on','jitterAmount',0.25)
hold on
e1 = errorbar(1,data.gammaBandPower.Contra.avgMeanR2,data.gammaBandPower.Contra.stdMeanR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.gammaBandPower.Whisk.Comb.meanR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25)
e2 = errorbar(2,data.gammaBandPower.Whisk.avgMeanR2,data.gammaBandPower.Whisk.stdMeanR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.gammaBandPower.Rest.Comb.meanR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25)
e3 = errorbar(3,data.gammaBandPower.Rest.avgMeanR2,data.gammaBandPower.Rest.stdMeanR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.gammaBandPower.NREM.Comb.meanR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25)
e4 = errorbar(4,data.gammaBandPower.NREM.avgMeanR2,data.gammaBandPower.NREM.stdMeanR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.gammaBandPower.REM.Comb.meanR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25)
e5 = errorbar(5,data.gammaBandPower.REM.avgMeanR2,data.gammaBandPower.REM.stdMeanR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title({'[D] Gamma-band [30-100 Hz]';'Mean R^2 kernel predictions'})
ylabel('R^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [E] MUA-derived IR analytic
ax5 = subplot(4,2,5);
plot(data.muaPower.Contra.meanIRtimeVec,data.muaPower.Contra.meanIR,'color',colorE,'LineWidth',2);
hold on
plot(data.muaPower.Whisk.meanIRtimeVec,data.muaPower.Whisk.meanIR,'color',colorD,'LineWidth',2);
plot(data.muaPower.Rest.meanIRtimeVec,data.muaPower.Rest.meanIR,'color',colorA,'LineWidth',2);
plot(data.muaPower.NREM.meanIRtimeVec,data.muaPower.NREM.meanIR,'color',colorB,'LineWidth',2);
plot(data.muaPower.REM.meanIRtimeVec,data.muaPower.REM.meanIR,'color',colorC,'LineWidth',2);
title({'[E] MUA [300-3000 Hz]','Impulse response function'})
xlabel('Time (s)')
ylabel({'MUA [300-3000 Hz]';'HRF amplitude (A.U.)'})
set(gca,'box','off')
xlim([0,5])
ax5.TickLength = [0.03,0.03];
%% [F] MUA-derived gamma function
ax6 = subplot(4,2,6);
plot(data.muaPower.Contra.meanGammaTimeVec,data.muaPower.Contra.meanGamma,'color',colorE,'LineWidth',2)
hold on
plot(data.muaPower.Whisk.meanGammaTimeVec,data.muaPower.Whisk.meanGamma,'color',colorD,'LineWidth',2)
plot(data.muaPower.Rest.meanGammaTimeVec,data.muaPower.Rest.meanGamma,'color',colorA,'LineWidth',2)
plot(data.muaPower.NREM.meanGammaTimeVec,data.muaPower.NREM.meanGamma,'color',colorB,'LineWidth',2)
plot(data.muaPower.REM.meanGammaTimeVec,data.muaPower.REM.meanGamma,'color',colorC,'LineWidth',2)
title({'[F] MUA[300-3000 Hz]','Gamma distribution function'})
xlabel('Time (s)')
ylabel({'MUA [300-3000 Hz]';'HRF amplitude (A.U.)'})
set(gca,'box','off')
xlim([0,5])
ax6.TickLength = [0.03,0.03];
%% [G] Gamma-band-derived median R2 values
ax7 = subplot(4,2,7);
scatter(xInds*1,data.muaPower.Contra.Comb.medR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorE,'jitter','on','jitterAmount',0.25)
hold on
e1 = errorbar(1,data.muaPower.Contra.avgMedR2,data.muaPower.Contra.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.muaPower.Whisk.Comb.medR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25)
e2 = errorbar(2,data.muaPower.Whisk.avgMedR2,data.muaPower.Whisk.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.muaPower.Rest.Comb.medR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25)
e3 = errorbar(3,data.muaPower.Rest.avgMedR2,data.muaPower.Rest.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.muaPower.NREM.Comb.medR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25)
e4 = errorbar(4,data.muaPower.NREM.avgMedR2,data.muaPower.NREM.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.muaPower.REM.Comb.medR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25)
e5 = errorbar(5,data.muaPower.REM.avgMedR2,data.muaPower.REM.stdMedR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title({'[E] MUA [300-3000 Hz]';'Median R^2 kernel predictions'})
ylabel('R^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [H] Gamma-band-derived mean R2 values
ax8 = subplot(4,2,8);
scatter(xInds*1,data.muaPower.Contra.Comb.meanR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorE,'jitter','on','jitterAmount',0.25)
hold on
e1 = errorbar(1,data.muaPower.Contra.avgMeanR2,data.muaPower.Contra.stdMeanR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(xInds*2,data.muaPower.Whisk.Comb.meanR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25)
e2 = errorbar(2,data.muaPower.Whisk.avgMeanR2,data.muaPower.Whisk.stdMeanR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(xInds*3,data.muaPower.Rest.Comb.meanR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25)
e3 = errorbar(3,data.muaPower.Rest.avgMeanR2,data.muaPower.Rest.stdMeanR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(xInds*4,data.muaPower.NREM.Comb.meanR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25)
e4 = errorbar(4,data.muaPower.NREM.avgMeanR2,data.muaPower.NREM.stdMeanR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(xInds*5,data.muaPower.REM.Comb.meanR2,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25)
e5 = errorbar(5,data.muaPower.REM.avgMeanR2,data.muaPower.REM.stdMeanR2,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title({'[H] MUA [300-3000 Hz]';'Mean R^2 kernel predictions'})
ylabel('R^2')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Figure Panel 12']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Supplemental Figure Panel 12'])

end
