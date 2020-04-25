function [] = FigurePanelThree_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

IOSanimalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
samplingRate = 30;   % Hz
%% Mean transitions between each arousal-state
% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,a};
    for b = 1:length(transitions)
        transition = transitions{1,b};
        data.(transition).EMG(a,:) = AnalysisResults.(animalID).Transitions.(transition).EMG;
        data.(transition).Hip(:,:,a) = AnalysisResults.(animalID).Transitions.(transition).Hip;
        data.(transition).T = AnalysisResults.(animalID).Transitions.(transition).T;
        data.(transition).F = AnalysisResults.(animalID).Transitions.(transition).F;
        data.(transition).Cort(:,:,a) = AnalysisResults.(animalID).Transitions.(transition).Cort;
        data.(transition).HbT(a,:) = AnalysisResults.(animalID).Transitions.(transition).HbT;
    end
end
% take average for each behavioral transition
for c = 1:length(transitions)
    transition = transitions{1,c};
    data.(transition).meanEMG = mean(data.(transition).EMG,1)*100;
    data.(transition).stdEMG = std(data.(transition).EMG,0,1)*100;
    data.(transition).meanHbT = mean(data.(transition).HbT,1);
    data.(transition).stdHbT = std(data.(transition).HbT,0,1);
    data.(transition).meanHip = mean(data.(transition).Hip,3)*100;
    data.(transition).meanCort = mean(data.(transition).Cort,3)*100;
end
T = (1/10):(1/10):60;
%% Figure panel three
summaryFigure = figure;
sgtitle('Turner Manuscript 2020 - Figure panel 3')
%% [A] Awake to NREM
ax1 = subplot(6,2,1);
% HbT and EMG
p1 = plot((1:length(data.AWAKEtoNREM.meanHbT))/samplingRate,data.AWAKEtoNREM.meanHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',2);
hold on
plot((1:length(data.AWAKEtoNREM.meanHbT))/samplingRate,data.AWAKEtoNREM.meanHbT + data.AWAKEtoNREM.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
plot((1:length(data.AWAKEtoNREM.meanHbT))/samplingRate,data.AWAKEtoNREM.meanHbT - data.AWAKEtoNREM.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\DeltaHbT (\muM)')
xlim([0,60])
yyaxis right
p2 = plot((1:length(data.AWAKEtoNREM.meanEMG ))/samplingRate,data.AWAKEtoNREM.meanEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
plot((1:length(data.AWAKEtoNREM.meanEMG))/samplingRate,data.AWAKEtoNREM.meanEMG + data.AWAKEtoNREM.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
plot((1:length(data.AWAKEtoNREM.meanEMG))/samplingRate,data.AWAKEtoNREM.meanEMG - data.AWAKEtoNREM.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
title('[A] Awake to NREM behavioral transition')
xlabel('Time (s)')
ylabel('\DeltaEMG power (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
legend([p1,p2],'HbT','EMG')
ax1.YAxis(1).Color = colors_Manuscript2020('dark candy apple red');
ax1.YAxis(2).Color = colors_Manuscript2020('rich black');
axis tight
ax1.TickLength = [0.03,0.03];
% cort neural
ax2 = subplot(6,2,3);
semilog_imagesc_Manuscript2020(T,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
axis xy
% c1 = colorbar;
% ylabel(c1,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-150,150])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
% hippocampal neural
ax3 = subplot(6,2,5);
semilog_imagesc_Manuscript2020(T,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
% c2 = colorbar;
% ylabel(c2,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-150,150])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [B] NREM to Awake
ax4 = subplot(6,2,2);
% HbT and EMG
plot((1:length(data.NREMtoAWAKE.meanHbT))/samplingRate,data.NREMtoAWAKE.meanHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',2);
hold on
plot((1:length(data.NREMtoAWAKE.meanHbT))/samplingRate,data.NREMtoAWAKE.meanHbT + data.NREMtoAWAKE.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
plot((1:length(data.NREMtoAWAKE.meanHbT))/samplingRate,data.NREMtoAWAKE.meanHbT - data.NREMtoAWAKE.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\DeltaHbT (\muM)')
xlim([0,60])
yyaxis right
plot((1:length(data.NREMtoAWAKE.meanEMG ))/samplingRate,data.NREMtoAWAKE.meanEMG ,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
plot((1:length(data.NREMtoAWAKE.meanEMG))/samplingRate,data.NREMtoAWAKE.meanEMG + data.NREMtoAWAKE.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
plot((1:length(data.NREMtoAWAKE.meanEMG))/samplingRate,data.NREMtoAWAKE.meanEMG - data.NREMtoAWAKE.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
title('[B] NREM to Awake behavioral transition')
xlabel('Time (s)')
ylabel('\DeltaEMG power (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax4.YAxis(1).Color = colors_Manuscript2020('dark candy apple red');
ax4.YAxis(2).Color = colors_Manuscript2020('rich black');
axis tight
ax4.TickLength = [0.03,0.03];
% cort neural
ax5 = subplot(6,2,4);
semilog_imagesc_Manuscript2020(T,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-150,150])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
% hippocampal neural
ax6 = subplot(6,2,6);
semilog_imagesc_Manuscript2020(T,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-150,150])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [C] NREM to REM
ax7 = subplot(6,2,7);
% HbT and EMG
plot((1:length(data.NREMtoREM.meanHbT))/samplingRate,data.NREMtoREM.meanHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',2);
hold on
plot((1:length(data.NREMtoREM.meanHbT))/samplingRate,data.NREMtoREM.meanHbT + data.NREMtoREM.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
plot((1:length(data.NREMtoREM.meanHbT))/samplingRate,data.NREMtoREM.meanHbT - data.NREMtoREM.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\DeltaHbT (\muM)')
xlim([0,60])
yyaxis right
plot((1:length(data.NREMtoREM.meanEMG ))/samplingRate,data.NREMtoREM.meanEMG ,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
plot((1:length(data.NREMtoREM.meanEMG))/samplingRate,data.NREMtoREM.meanEMG + data.NREMtoREM.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
plot((1:length(data.NREMtoREM.meanEMG))/samplingRate,data.NREMtoREM.meanEMG - data.NREMtoREM.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
title('[C] NREM to REM behavioral transition')
xlabel('Time (s)')
ylabel('\DeltaEMG power (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax7.YAxis(1).Color = colors_Manuscript2020('dark candy apple red');
ax7.YAxis(2).Color = colors_Manuscript2020('rich black');
axis tight
ax7.TickLength = [0.03,0.03];
% cort neural
ax8 = subplot(6,2,9);
semilog_imagesc_Manuscript2020(T,data.NREMtoREM.F,data.NREMtoREM.meanCort,'y')
axis xy
% c5 = colorbar;
% ylabel(c5,'\DeltaP/P (%)')
caxis([-150,150])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
% hippocampal neural
ax9 = subplot(6,2,11);
semilog_imagesc_Manuscript2020(T,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
% c6 = colorbar;
% ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-150,150])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [D] REM to Awake
ax10 = subplot(6,2,8);
plot((1:length(data.REMtoAWAKE.meanHbT))/samplingRate,data.REMtoAWAKE.meanHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',2);
hold on
plot((1:length(data.REMtoAWAKE.meanHbT))/samplingRate,data.REMtoAWAKE.meanHbT + data.REMtoAWAKE.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
plot((1:length(data.REMtoAWAKE.meanHbT))/samplingRate,data.REMtoAWAKE.meanHbT - data.REMtoAWAKE.stdHbT,'-','color',colors_Manuscript2020('dark candy apple red'),'LineWidth',0.5);
ylabel('\DeltaHbT (\muM)')
xlim([0,60])
yyaxis right
plot((1:length(data.REMtoAWAKE.meanEMG ))/samplingRate,data.REMtoAWAKE.meanEMG ,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
plot((1:length(data.REMtoAWAKE.meanEMG))/samplingRate,data.REMtoAWAKE.meanEMG + data.REMtoAWAKE.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
plot((1:length(data.REMtoAWAKE.meanEMG))/samplingRate,data.REMtoAWAKE.meanEMG - data.REMtoAWAKE.stdEMG,'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
title('[D] REM to Awake behavioral transition')
xlabel('Time (s)')
ylabel('\DeltaEMG power (%)','rotation',-90,'VerticalAlignment','bottom')
set(gca,'box','off')
ax10.YAxis(1).Color = colors_Manuscript2020('dark candy apple red');
ax10.YAxis(2).Color = colors_Manuscript2020('rich black');
axis tight
ax10.TickLength = [0.03,0.03];
% cort neural
ax11 = subplot(6,2,10);
semilog_imagesc_Manuscript2020(T,data.REMtoAWAKE.F,data.REMtoAWAKE.meanCort,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-150,150])
xlabel('Time (s)')
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
% hippocampal neural
ax12 = subplot(6,2,12);
semilog_imagesc_Manuscript2020(T,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
c8 = colorbar;
ylabel(c8,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-150,150])
xlabel('Time (s)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
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
savefig(summaryFigure,[dirpath 'Figure Panel Three']);
print('-painters','-dpdf','-fillpage',[dirpath 'Figure Panel Three'])

end
