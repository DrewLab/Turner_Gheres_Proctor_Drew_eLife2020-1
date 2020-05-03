function [] = FigurePanelThree_Manuscript2020(rootFolder)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose:
%________________________________________________________________________________________________________________________

colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
%% information and data for first example
animalID = 'T118';
dataLocation = [rootFolder '\' animalID '\2P Data\'];
cd(dataLocation)
exampleMergedFileID = 'T118_RH_191211_14_44_20_015_A1_MergedData.mat';
load(exampleMergedFileID,'-mat')
exampleSpecFileID = 'T118_RH_191211_14_44_20_015_A1_SpecData.mat';
load(exampleSpecFileID,'-mat')
exampleBaselinesFileID = 'T118_RestingBaselines.mat';
load(exampleBaselinesFileID,'-mat')
[~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P_Manuscript2020(exampleMergedFileID);
strDay = ConvertDate_2P_Manuscript2020(fileDate);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(MergedData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(MergedData.notes.p2Fs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filtWhiskerAngle = filtfilt(sos1,g1,MergedData.data.whiskerAngle);
% pressure sensor
filtForceSensor = filtfilt(sos1,g1,abs(MergedData.data.forceSensorL));
% EMG
EMG = MergedData.data.EMG.data;
normEMG = EMG - RestingBaselines.manualSelection.EMG.data.(strDay);
filtEMG = filtfilt(sos1,g1,normEMG);
% vessel diameter
vesselDiameter = MergedData.data.vesselDiameter.data;
normVesselDiameter = (vesselDiameter - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay))./(RestingBaselines.manualSelection.vesselDiameter.data.(vesselID).(strDay));
filtVesselDiameter = filtfilt(sos2,g2,normVesselDiameter)*100;
% cortical and hippocampal spectrograms
cortNormS = SpecData.corticalNeural.fiveSec.normS.*100;
hipNormS = SpecData.hippocampalNeural.fiveSec.normS.*100;
T = SpecData.corticalNeural.fiveSec.T;
F = SpecData.corticalNeural.fiveSec.F;
cd(rootFolder)
%% Figure Panel Two
summaryFigure = figure;
sgtitle('Figure Panel 3 - Turner Manuscript 2020')
%% [B] second single trial 2P sleep example
% EMG and force sensor
ax1 = subplot(6,1,1);
p1 = plot((1:length(filtEMG))/MergedData.notes.dsFs,filtEMG,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2,2.5])
yyaxis right
p2 = plot((1:length(filtForceSensor))/MergedData.notes.dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
title('[B] 2PLSM sample data')
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
xlim([15,615])
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% whisker angle
ax2 = subplot(6,1,2);
plot((1:length(filtWhiskerAngle))/MergedData.notes.dsFs,-filtWhiskerAngle,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax2.TickLength = [0.01,0.01];
xlim([15,615])
ylim([-10,40])
% vessel diameter
ax34 = subplot(6,1,[3,4]);
p3 = plot((1:length(filtVesselDiameter))/MergedData.notes.p2Fs,filtVesselDiameter,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
hold on
xline(15,'color',colorA,'LineWidth',2)
x2 = xline(55,'color',colorB,'LineWidth',2);
x1 = xline(97,'color',colorA,'LineWidth',2);
xline(105,'color',colorB,'LineWidth',2)
x3 = xline(156,'color',colorC,'LineWidth',2);
xline(224,'color',colorA,'LineWidth',2)
xline(248,'color',colorB,'LineWidth',2)
xline(342,'color',colorA,'LineWidth',2)
xline(360,'color',colorB,'LineWidth',2)
xline(450,'color',colorC,'LineWidth',2)
xline(537,'color',colorA,'LineWidth',2)
ylabel('\DeltaD/D (%)')
legend([p3,x1,x2,x3],'Arteriole diameter','Awake','NREM','REM')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax34.TickLength = [0.01,0.01];
xlim([15,615])
ylim([-30,60])
% cortical LFP
ax5 = subplot(6,1,5);
semilog_imagesc_Manuscript2020(T,F,cortNormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'Cort LFP','Freq (Hz)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax5.TickLength = [0.01,0.01];
xlim([15,615])
% hippocampal LFP
ax6 = subplot(6,1,6);
semilog_imagesc_Manuscript2020(T,F,hipNormS,'y')
axis xy
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax6.TickLength = [0.01,0.01];
xlim([15,615])
% Axes properties
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Figure Panel 3']);
% remove surface subplots because they take forever to render
cla(ax5);
set(ax5,'YLim',[1,99]);
cla(ax6);
set(ax6,'YLim',[1,99]);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Figure Panel 3'])
%% subplot figures
figure;
% example 1 cortical LFP
subplot(2,1,1);
semilog_imagesc_Manuscript2020(T,F,cortNormS,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([15,615])
% example 1 hippocampal LFP
subplot(2,1,2);
semilog_imagesc_Manuscript2020(T,F,hipNormS,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([15,615])
print('-painters','-dtiffn',[dirpath 'Figure Panel 3 subplot images'])

end
