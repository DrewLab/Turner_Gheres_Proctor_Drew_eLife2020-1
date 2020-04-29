function [] = FigurePanelFourB_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T115','T116','T117','T118','T125','T126'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
samplingRate = 5;   % Hz
%% Mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
data.AWAKEtoNREM.vesselDiameter = []; data.NREMtoAWAKE.vesselDiameter = [];
data.NREMtoREM.vesselDiameter = []; data.REMtoAWAKE.vesselDiameter = [];
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for b = 1:length(transitions)
        transition = transitions{1,b};
        if isfield(AnalysisResults.(animalID).Transitions,transition)
            data.(transition).vesselDiameter = vertcat(data.(transition).vesselDiameter,AnalysisResults.(animalID).Transitions.(transition).vesselDiameter);
        end
    end
end
% take average for each behavioral transition
for c = 1:length(transitions)
    transition = transitions{1,c};
    data.(transition).meanVD = mean(data.(transition).vesselDiameter,1);
    data.(transition).stdVD = std(data.(transition).vesselDiameter,0,1);
end
T1 = -30 + (1/5):(1/5):30;
T2 = -30 + (1/5):(1/5):70;
%% Figure panel 4
summaryFigure = figure;
sgtitle('Figure panel 4B - Turner Manuscript 2020')
%% [A] Awake to NREM
ax1 = subplot(6,2,1);
% HbT and EMG
p1 = plot(T1,data.AWAKEtoNREM.meanHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanHbT + data.AWAKEtoNREM.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.meanHbT - data.AWAKEtoNREM.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\DeltaHbT (\muM)')
xlim([-30,30])
yyaxis right
p2 = plot(T1,data.AWAKEtoNREM.meanEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
plot(T1,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
plot(T1,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
title('[A] Awake to NREM behavioral transition')
xlabel('Time (s)')
ylabel('EMG log10(pwr)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'HbT','EMG')
ax1.YAxis(1).Color = colors_Manuscript2020('dark candy apple red');
ax1.YAxis(2).Color = colors_Manuscript2020('rich black');
axis tight
ax1.TickLength = [0.03,0.03];
% cort neural
ax2 = subplot(6,2,3);
semilog_imagesc_Manuscript2020(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
axis xy
% c1 = colorbar;
% ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-300,300])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
% hippocampal neural
ax3 = subplot(6,2,5);
semilog_imagesc_Manuscript2020(T2,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
% c2 = colorbar;
% ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-300,300])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [B] NREM to Awake
ax4 = subplot(6,2,2);
% HbT and EMG
plot(T1,data.NREMtoAWAKE.meanHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanHbT + data.NREMtoAWAKE.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.meanHbT - data.NREMtoAWAKE.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\DeltaHbT (\muM)')
xlim([-30,30])
yyaxis right
plot(T1,data.NREMtoAWAKE.meanEMG ,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
plot(T1,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
plot(T1,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
title('[B] NREM to Awake behavioral transition')
xlabel('Time (s)')
ylabel('EMG log10(pwr)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors_Manuscript2020('dark candy apple red');
ax4.YAxis(2).Color = colors_Manuscript2020('rich black');
axis tight
ax4.TickLength = [0.03,0.03];
% cort neural
ax5 = subplot(6,2,4);
semilog_imagesc_Manuscript2020(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-300,300])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
% hippocampal neural
ax6 = subplot(6,2,6);
semilog_imagesc_Manuscript2020(T2,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-300,300])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [C] NREM to REM
ax7 = subplot(6,2,7);
% HbT and EMG
plot(T1,data.NREMtoREM.meanHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.meanHbT + data.NREMtoREM.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.NREMtoREM.meanHbT - data.NREMtoREM.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\DeltaHbT (\muM)')
xlim([-30,30])
yyaxis right
plot(T1,data.NREMtoREM.meanEMG ,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
plot(T1,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
plot(T1,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
title('[C] NREM to REM behavioral transition')
xlabel('Time (s)')
ylabel('EMG log10(pwr)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax7.YAxis(1).Color = colors_Manuscript2020('dark candy apple red');
ax7.YAxis(2).Color = colors_Manuscript2020('rich black');
axis tight
ax7.TickLength = [0.03,0.03];
% cort neural
ax8 = subplot(6,2,9);
semilog_imagesc_Manuscript2020(T2,data.NREMtoREM.F,data.NREMtoREM.meanCort,'y')
axis xy
% c5 = colorbar;
% ylabel(c5,'\DeltaP/P (%)')
caxis([-300,300])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
% hippocampal neural
ax9 = subplot(6,2,11);
semilog_imagesc_Manuscript2020(T2,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
% c6 = colorbar;
% ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-300,300])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [D] REM to Awake
ax10 = subplot(6,2,8);
plot(T1,data.REMtoAWAKE.meanHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.meanHbT + data.REMtoAWAKE.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
plot(T1,data.REMtoAWAKE.meanHbT - data.REMtoAWAKE.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\DeltaHbT (\muM)')
xlim([-30,30])
yyaxis right
plot(T1,data.REMtoAWAKE.meanEMG ,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
plot(T1,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
plot(T1,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
title('[D] REM to Awake behavioral transition')
xlabel('Time (s)')
ylabel('EMG log10(pwr)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax10.YAxis(1).Color = colors_Manuscript2020('dark candy apple red');
ax10.YAxis(2).Color = colors_Manuscript2020('rich black');
axis tight
ax10.TickLength = [0.03,0.03];
% cort neural
ax11 = subplot(6,2,10);
semilog_imagesc_Manuscript2020(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanCort,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-300,300])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
% hippocampal neural
ax12 = subplot(6,2,12);
semilog_imagesc_Manuscript2020(T2,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
c8 = colorbar;
ylabel(c8,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-300,300])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([-30,30])
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% axes positionns
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax8Pos = get(ax8,'position');
ax9Pos = get(ax9,'position');
ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
ax12Pos = get(ax12,'position');
ax2Pos(3:4) = ax1Pos(3:4);
ax3Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax4Pos(3:4);
ax6Pos(3:4) = ax4Pos(3:4);
ax8Pos(3:4) = ax7Pos(3:4);
ax9Pos(3:4) = ax7Pos(3:4);
ax11Pos(3:4) = ax10Pos(3:4);
ax12Pos(3:4) = ax10Pos(3:4);
set(ax2,'position',ax2Pos);
set(ax3,'position',ax3Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax8,'position',ax8Pos);
set(ax9,'position',ax9Pos);
set(ax11,'position',ax11Pos);
set(ax12,'position',ax12Pos);
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Figure Panel 4']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Figure Panel 4'])

end
