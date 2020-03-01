function [] = AvgBehaviorTransitions_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Calculate the average physiological transitions between different behavioral states
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120'};
transitions = {'AWAKEtoNREM','NREMtoAWAKE','NREMtoREM','REMtoAWAKE'};
samplingRate = 30;   % Hz

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
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
    data.(transition).meanEMG = mean(data.(transition).EMG,1);
    data.(transition).meanHbT = mean(data.(transition).HbT,1);
    data.(transition).meanHip = mean(data.(transition).Hip,3)*100;
    data.(transition).meanCort = mean(data.(transition).Cort,3)*100;
end

%% summary figure(s)
spec10 = find(round(data.AWAKEtoNREM.F) == 10);
spec10 = spec10(end);
spec100 = find(round(data.AWAKEtoNREM.F) == 99);
spec100 = spec100(end);
summaryFigure = figure;
sgtitle('Average behavioral transitions')
%% Awake to NREM
% behavior
ax1 = subplot(6,2,1);
p1 = plot((1:length(data.AWAKEtoNREM.meanHbT))/samplingRate,data.AWAKEtoNREM.meanHbT,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
for d = 1:size(data.AWAKEtoNREM.HbT,1)
    plot((1:length(data.AWAKEtoNREM.HbT(d,:)))/samplingRate,data.AWAKEtoNREM.HbT(d,:),'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
end
ylabel('\DeltaHbT (\muM)')
xlim([0,60])
yyaxis right
p2 = plot((1:length(data.AWAKEtoNREM.meanEMG ))/samplingRate,data.AWAKEtoNREM.meanEMG,'-','color',colors_Manuscript2020('deep carrot orange'),'LineWidth',2);
hold on
for e = 1:size(data.AWAKEtoNREM.EMG,1)
    plot((1:length(data.AWAKEtoNREM.EMG(e,:)))/samplingRate,data.AWAKEtoNREM.EMG(e,:),'-','color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
end
ylabel('EMG (Volts^2)')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
legend([p1,p2],'HbT','EMG')
title('Awake to NREM')
% LH cort neural
ax2 = subplot(6,2,3);
semilog_imagesc_Manuscript2020(data.AWAKEtoNREM.T,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanCort,'y')
axis xy
c1 = colorbar;
ylabel(c1,'\DeltaP/P (%)')
caxis([-100,200])
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
% hip
ax3 = subplot(6,2,5);
semilog_imagesc_Manuscript2020(data.AWAKEtoNREM.T,data.AWAKEtoNREM.F,data.AWAKEtoNREM.meanHip,'y')
c2 = colorbar;
ylabel(c2,'\DeltaP/P (%)')
caxis([-100,200])
xlabel('Time (sec)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'TickLength',[0,0])
set(gca,'box','off')

%% NREM to Awake
% behavior
ax4 = subplot(6,2,2);
plot((1:length(data.NREMtoAWAKE.meanHbT))/samplingRate,data.NREMtoAWAKE.meanHbT,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
hold on
for d = 1:size(data.NREMtoAWAKE.HbT,1)
    plot((1:length(data.NREMtoAWAKE.HbT(d,:)))/samplingRate,data.NREMtoAWAKE.HbT(d,:),'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
end
ylabel('\DeltaHbT (\muM)')
xlim([0,60])
yyaxis right
plot((1:length(data.NREMtoAWAKE.meanEMG ))/samplingRate,data.NREMtoAWAKE.meanEMG ,'-','color',colors_Manuscript2020('deep carrot orange'),'LineWidth',2);
hold on
for e = 1:size(data.NREMtoAWAKE.EMG,1)
    plot((1:length(data.NREMtoAWAKE.EMG(e,:)))/samplingRate,data.NREMtoAWAKE.EMG(e,:),'-','color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
end
ylabel('EMG (Volts^2)')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
title('NREM to Awake')
% LH cort neural
ax5 = subplot(6,2,4);
semilog_imagesc_Manuscript2020(data.NREMtoAWAKE.T,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanCort,'y')
axis xy
c3 = colorbar;
ylabel(c3,'\DeltaP/P (%)')
caxis([-100,200])
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
% hip
ax6 = subplot(6,2,6);
semilog_imagesc_Manuscript2020(data.NREMtoAWAKE.T,data.NREMtoAWAKE.F,data.NREMtoAWAKE.meanHip,'y')
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)')
caxis([-100,200])
xlabel('Time (sec)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'TickLength',[0,0])
set(gca,'box','off')

%% NREM to REM
% behavior
ax7 = subplot(6,2,7);
plot((1:length(data.NREMtoREM.meanHbT))/samplingRate,data.NREMtoREM.meanHbT,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
for d = 1:size(data.NREMtoREM.HbT,1)
    plot((1:length(data.NREMtoREM.HbT(d,:)))/samplingRate,data.NREMtoREM.HbT(d,:),'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
end
ylabel('\DeltaHbT (\muM)')
xlim([0,60])
yyaxis right
plot((1:length(data.NREMtoREM.meanEMG ))/samplingRate,data.NREMtoREM.meanEMG ,'-','color',colors_Manuscript2020('deep carrot orange'),'LineWidth',2);
hold on
for e = 1:size(data.NREMtoREM.EMG,1)
    plot((1:length(data.NREMtoREM.EMG(e,:)))/samplingRate,data.NREMtoREM.EMG(e,:),'-','color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
end
ylabel('EMG (Volts^2)')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
title('NREM to REM')
% LH cort neural
ax8 = subplot(6,2,9);
semilog_imagesc_Manuscript2020(data.NREMtoREM.T,data.NREMtoREM.F,data.NREMtoREM.meanCort,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)')
caxis([-100,200])
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
% hip
ax9 = subplot(6,2,11);
semilog_imagesc_Manuscript2020(data.NREMtoREM.T,data.NREMtoREM.F,data.NREMtoREM.meanHip,'y')
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)')
caxis([-100,200])
xlabel('Time (sec)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'TickLength',[0,0])
set(gca,'box','off')

%% REM to Awake
% behavior
ax10 = subplot(6,2,8);
plot((1:length(data.REMtoAWAKE.meanHbT))/samplingRate,data.REMtoAWAKE.meanHbT,'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
for d = 1:size(data.REMtoAWAKE.HbT,1)
    plot((1:length(data.REMtoAWAKE.HbT(d,:)))/samplingRate,data.REMtoAWAKE.HbT(d,:),'-','color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
end
ylabel('\DeltaHbT (\muM)')
xlim([0,60])
yyaxis right
plot((1:length(data.REMtoAWAKE.meanEMG ))/samplingRate,data.REMtoAWAKE.meanEMG ,'-','color',colors_Manuscript2020('deep carrot orange'),'LineWidth',2);
hold on
for e = 1:size(data.REMtoAWAKE.EMG,1)
    plot((1:length(data.REMtoAWAKE.EMG(e,:)))/samplingRate,data.REMtoAWAKE.EMG(e,:),'-','color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
end
ylabel('EMG (Volts^2)')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
title('REM to Awake')
% LH cort neural
ax11 = subplot(6,2,10);
semilog_imagesc_Manuscript2020(data.REMtoAWAKE.T,data.REMtoAWAKE.F,data.REMtoAWAKE.meanCort,'y')
axis xy
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)')
caxis([-100,200])
ylabel({'Cortical LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
% hip
ax12 = subplot(6,2,12);
semilog_imagesc_Manuscript2020(data.REMtoAWAKE.T,data.REMtoAWAKE.F,data.REMtoAWAKE.meanHip,'y')
c8 = colorbar;
ylabel(c8,'\DeltaP/P (%)')
caxis([-100,200])
xlabel('Time (sec)')
ylabel({'Hippocampal LFP';'Frequency (Hz)'})
set(gca,'Yticklabel','10^1')
xlim([0,60])
set(gca,'TickLength',[0,0])
set(gca,'box','off')

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

% save figure(s)
dirpath = [rootFolder '\Analysis Figures\'];
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Behavior Transitions']);

end
