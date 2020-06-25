function [AnalysisResults] = FigS8_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel S8 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

colorRest = [(51/256),(160/256),(44/256)];
colorNREM = [(192/256),(0/256),(256/256)];
colorREM = [(255/256),(140/256),(0/256)];
% colorAwake = [(256/256),(192/256),(0/256)];
% colorSleep = [(0/256),(128/256),(256/256)];
% colorAll = [(184/256),(115/256),(51/256)];
% colorWhisk = [(31/256),(120/256),(180/256)];
% colorStim = [(256/256),(28/256),(207/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
% information and data for first example
animalID = 'T108';
dataLocation = [rootFolder '\' animalID '\Bilateral Imaging\'];
cd(dataLocation)
exampleProcDataFileID = 'T108_190822_13_23_18_ProcData.mat';
load(exampleProcDataFileID,'-mat')
exampleSpecDataFileID = 'T108_190822_13_23_18_SpecDataA.mat';
load(exampleSpecDataFileID,'-mat')
exampleBaselineFileID = 'T108_RestingBaselines.mat';
load(exampleBaselineFileID,'-mat')
[~,fileDate,~] = GetFileInfo_IOS_Manuscript2020(exampleProcDataFileID);
strDay = ConvertDate_IOS_Manuscript2020(fileDate);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1_A,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1_A,k1);
[z2,p2_A,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2_A,k2);
% whisker angle
filtWhiskerAngle = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% force sensor
filtForceSensor = filtfilt(sos1,g1,abs(ProcData.data.forceSensor));
% EMG
EMG = ProcData.data.EMG.emg;
normEMG = EMG - RestingBaselines.manualSelection.EMG.emg.(strDay);
filtEMG = filtfilt(sos1,g1,normEMG);
% heart rate
heartRate = ProcData.data.heartRate;
% CBV data
LH_HbT = ProcData.data.CBV_HbT.adjLH;
filtLH_HbT = filtfilt(sos2,g2,LH_HbT);
RH_HbT = ProcData.data.CBV_HbT.adjRH;
filtRH_HbT = filtfilt(sos2,g2,RH_HbT);
% cortical and hippocampal spectrograms
cortical_LHnormS = SpecData.cortical_LH.normS.*100;
cortical_RHnormS = SpecData.cortical_RH.normS.*100;
hippocampusNormS = SpecData.hippocampus.normS.*100;
T = SpecData.cortical_LH.T;
F = SpecData.cortical_LH.F;
cd(rootFolder)
%% Fig. S8
summaryFigure = figure('Name','FigS8 (a-f)'); %#ok<*NASGU>
sgtitle('Figure Panel S8 (a-f) Turner Manuscript 2020')
%% EMG and force sensor
ax1 = subplot(7,1,1);
p1 = plot((1:length(filtEMG))/ProcData.notes.dsFs,filtEMG,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3])
yyaxis right
p2 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
xlim([300,900]) 
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
%% Whisker angle and heart rate
ax2 = subplot(7,1,2);
p3 = plot((1:length(filtWhiskerAngle))/ProcData.notes.dsFs,-filtWhiskerAngle,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
ylim([-20,60])
yyaxis right
p4 = plot((1:length(heartRate)),heartRate,'color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
legend([p3,p4],'whisker angle','heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
xlim([300,900]) 
ylim([5,15])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors_Manuscript2020('rich black');
ax2.YAxis(2).Color = colors_Manuscript2020('deep carrot orange');
%% CBV and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6 = plot((1:length(filtRH_HbT))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT,'color',colors_Manuscript2020('sapphire'),'LineWidth',1);
hold on
p5 = plot((1:length(filtLH_HbT))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
x1 = xline(300,'color',colorNREM,'LineWidth',2);
x2 = xline(600,'color',colorREM,'LineWidth',2);
x3 = xline(707,'color',colorRest,'LineWidth',2);
ylabel('\DeltaHbT')
legend([p5,p6,x3,x1,x2],'Left hem','Right hem','Awake','NREM','REM')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
axis tight
xlim([300,900]) 
ax34.TickLength = [0.01,0.01];
%% Left cortical electrode spectrogram
ax5 = subplot(7,1,5);
semilog_imagesc_Manuscript2020(T,F,cortical_LHnormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'LH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
xlim([300,900]) 
ax5.TickLength = [0.01,0.01];
%% Right cortical electrode spectrogram
ax6 = subplot(7,1,6);
semilog_imagesc_Manuscript2020(T,F,cortical_RHnormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'RH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
xlim([300,900]) 
ax6.TickLength = [0.01,0.01];
%% Hippocampal electrode spectrogram
ax7 = subplot(7,1,7);
semilog_imagesc_Manuscript2020(T,F,hippocampusNormS,'y')
c7 = colorbar;
ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (min)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([300,900]) 
ax7.TickLength = [0.01,0.01];
%% Axes properties
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax7Pos = get(ax7,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
ax7Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax7,'position',ax7Pos);
%% save figure(s)
% dirpath = [rootFolder '\Summary Figures and Structures\'];
% if ~exist(dirpath,'dir')
%     mkdir(dirpath);
% end
% savefig(summaryFigure,[dirpath 'FigS8']);
% % remove surface subplots because they take forever to render
% cla(ax5);
% set(ax5,'YLim',[1,99]);
% cla(ax6);
% set(ax6,'YLim',[1,99]);
% cla(ax7);
% set(ax7,'YLim',[1,99]);
% set(summaryFigure,'PaperPositionMode','auto');
% print('-painters','-dpdf','-bestfit',[dirpath 'FigS8'])
% close(summaryFigure)
% %% subplot figures
% summaryFigure_imgs = figure;
% % example 5 LH cortical LFP
% subplot(3,1,1);
% semilog_imagesc_Manuscript2020(T,F,cortical_LHnormS,'y')
% caxis([-100,200])
% set(gca,'box','off')
% axis xy
% axis tight
% axis off
% xlim([300,900]) 
% % example 5 RH cortical LFP
% subplot(3,1,2);
% semilog_imagesc_Manuscript2020(T,F,cortical_RHnormS,'y')
% caxis([-100,200])
% set(gca,'box','off')
% axis xy
% axis tight
% axis off
% xlim([300,900]) 
% % example 5 hippocampal LFP
% subplot(3,1,3);
% semilog_imagesc_Manuscript2020(T,F,hippocampusNormS,'y')
% caxis([-100,200])
% set(gca,'box','off')
% axis xy
% axis tight
% axis off
% xlim([300,900]) 
% print('-painters','-dtiffn',[dirpath 'FigS8 subplot images'])
% close(summaryFigure_imgs)
% %% Figure panel S8
% figure('Name','FigS8 (a-f)');
% sgtitle('Figure Panel S8 (a-f) Turner Manuscript 2020')
% %% EMG and force sensor
% ax1 = subplot(7,1,1);
% p1 = plot((1:length(filtEMG))/ProcData.notes.dsFs,filtEMG,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
% ylabel({'EMG','log10(pwr)'})
% ylim([-2.5,3])
% yyaxis right
% p2 = plot((1:length(filtForceSensor))/ProcData.notes.dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
% ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
% legend([p1,p2],'EMG','pressure')
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% xticks([300,360,420,480,540,600,660,720,780,840,900])
% xlim([300,900]) 
% ylim([-0.1,2.5])
% ax1.TickLength = [0.01,0.01];
% ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
% ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% %% Whisker angle and heart rate
% ax2 = subplot(7,1,2);
% p3 = plot((1:length(filtWhiskerAngle))/ProcData.notes.dsFs,-filtWhiskerAngle,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
% ylabel({'Whisker','angle (deg)'})
% ylim([-20,60])
% yyaxis right
% p4 = plot((1:length(heartRate)),heartRate,'color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
% ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
% legend([p3,p4],'whisker angle','heart rate')
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% xticks([300,360,420,480,540,600,660,720,780,840,900])
% xlim([300,900]) 
% ylim([5,15])
% ax2.TickLength = [0.01,0.01];
% ax2.YAxis(1).Color = colors_Manuscript2020('rich black');
% ax2.YAxis(2).Color = colors_Manuscript2020('deep carrot orange');
% %% CBV and behavioral indeces
% ax34 =subplot(7,1,[3,4]);
% p6 = plot((1:length(filtRH_HbT))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT,'color',colors_Manuscript2020('sapphire'),'LineWidth',1);
% hold on
% p5 = plot((1:length(filtLH_HbT))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
% x1 = xline(300,'color',colorNREM,'LineWidth',2);
% x2 = xline(600,'color',colorREM,'LineWidth',2);
% x3 = xline(707,'color',colorRest,'LineWidth',2);
% ylabel('\DeltaHbT')
% legend([p5,p6,x3,x1,x2],'Left hem','Right hem','Awake','NREM','REM')
% set(gca,'TickLength',[0,0])
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% xticks([300,360,420,480,540,600,660,720,780,840,900])
% axis tight
% xlim([300,900]) 
% ax34.TickLength = [0.01,0.01];
% %% Left cortical electrode spectrogram
% ax5 = subplot(7,1,5);
% semilog_imagesc_Manuscript2020(T,F,cortical_LHnormS,'y')
% axis xy
% c5 = colorbar;
% ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,200])
% ylabel({'LH cortical LFP','Freq (Hz)'})
% set(gca,'Yticklabel','10^1')
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% xticks([300,360,420,480,540,600,660,720,780,840,900])
% xlim([300,900]) 
% ax5.TickLength = [0.01,0.01];
% %% Right cortical electrode spectrogram
% ax6 = subplot(7,1,6);
% semilog_imagesc_Manuscript2020(T,F,cortical_RHnormS,'y')
% axis xy
% c6 = colorbar;
% ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,200])
% ylabel({'RH cortical LFP','Freq (Hz)'})
% set(gca,'Yticklabel','10^1')
% set(gca,'Xticklabel',[])
% set(gca,'box','off')
% xticks([300,360,420,480,540,600,660,720,780,840,900])
% xlim([300,900]) 
% ax6.TickLength = [0.01,0.01];
% %% Hippocampal electrode spectrogram
% ax7 = subplot(7,1,7);
% semilog_imagesc_Manuscript2020(T,F,hippocampusNormS,'y')
% c7 = colorbar;
% ylabel(c7,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
% caxis([-100,200])
% xlabel('Time (min)')
% ylabel({'Hippocampal LFP','Freq (Hz)'})
% set(gca,'box','off')
% xticks([300,360,420,480,540,600,660,720,780,840,900])
% xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
% xlim([300,900]) 
% ax7.TickLength = [0.01,0.01];
% %% Axes properties
% ax1Pos = get(ax1,'position');
% ax5Pos = get(ax5,'position');
% ax6Pos = get(ax6,'position');
% ax7Pos = get(ax7,'position');
% ax5Pos(3:4) = ax1Pos(3:4);
% ax6Pos(3:4) = ax1Pos(3:4);
% ax7Pos(3:4) = ax1Pos(3:4);
% set(ax5,'position',ax5Pos);
% set(ax6,'position',ax6Pos);
% set(ax7,'position',ax7Pos);

end
