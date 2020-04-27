function [] = FigurePanelTwo_Manuscript2020(rootFolder)
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
animalID_A = 'T116';
dataLocation = [rootFolder '\' animalID_A '\2P Data\'];
cd(dataLocation)
exampleMergedFileID_A = 'T116_RH_191120_11_18_02_009_A2_MergedData.mat';
load(exampleMergedFileID_A,'-mat')
exampleSpecFileID_A = 'T116_RH_191120_11_18_02_009_A2_SpecData.mat';
load(exampleSpecFileID_A,'-mat')
exampleBaselinesFileID_A = 'T116_RestingBaselines.mat';
load(exampleBaselinesFileID_A,'-mat')
[~,~,fileDate_A,~,~,vesselID_A] = GetFileInfo2_2P_Manuscript2020(exampleMergedFileID_A);
strDay_A = ConvertDate_2P_Manuscript2020(fileDate_A);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(MergedData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(MergedData.notes.p2Fs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filtWhiskerAngle_A = filtfilt(sos1,g1,MergedData.data.whiskerAngle);
% vessel diameter
vesselDiameter_A = MergedData.data.vesselDiameter.data;
normVesselDiameter_A = (vesselDiameter_A - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_A).(strDay_A))./(RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_A).(strDay_A));
filtVesselDiameter_A = filtfilt(sos2,g2,normVesselDiameter_A)*100;
% cortical and hippocampal spectrograms
cortNormS_A = SpecData.corticalNeural.fiveSec.normS.*100;
hipNormS_A = SpecData.hippocampalNeural.fiveSec.normS.*100;
cd(rootFolder)
%% information and data for second example
animalID_B = 'T118';
dataLocation = [rootFolder '\' animalID_B '\2P Data\'];
cd(dataLocation)
exampleMergedFileID_B = 'T118_RH_191211_14_44_20_015_A1_MergedData.mat';
load(exampleMergedFileID_B,'-mat')
exampleSpecFileID_B = 'T118_RH_191211_14_44_20_015_A1_SpecData.mat';
load(exampleSpecFileID_B,'-mat')
exampleBaselinesFileID_B = 'T118_RestingBaselines.mat';
load(exampleBaselinesFileID_B,'-mat')
[~,~,fileDate_B,~,~,vesselID_B] = GetFileInfo2_2P_Manuscript2020(exampleMergedFileID_B);
strDay_B = ConvertDate_2P_Manuscript2020(fileDate_B);
% whisker angle
filtWhiskerAngle_B = filtfilt(sos1,g1,MergedData.data.whiskerAngle);
% vessel diameter
vesselDiameter_B = MergedData.data.vesselDiameter.data;
normVesselDiameter_B = (vesselDiameter_B - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_B).(strDay_B))./(RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_B).(strDay_B));
filtVesselDiameter_B = filtfilt(sos2,g2,normVesselDiameter_B)*100;
% cortical and hippocampal spectrograms
cortNormS_B = SpecData.corticalNeural.fiveSec.normS.*100;
hipNormS_B = SpecData.hippocampalNeural.fiveSec.normS.*100;
T = SpecData.corticalNeural.fiveSec.T;
F = SpecData.corticalNeural.fiveSec.F;
cd(rootFolder)
%% Figure Panel Two
summaryFigure = figure;
sgtitle('Turner Manuscript 2020 - Figure Panel Two')
%% [A] first single trial 2P sleep example
% whisker angle
ax1 = subplot(10,1,1);
plot((1:length(filtWhiskerAngle_A))/MergedData.notes.dsFs,-filtWhiskerAngle_A,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([270,330,390,450,510,570,630,690,750,810,870])
ax1.TickLength = [0.03,0.03];
xlim([270,870])
ylim([-10,40])
% vessel diameter
ax23 = subplot(10,1,[2,3]);
p1 = plot((1:length(filtVesselDiameter_A))/MergedData.notes.p2Fs,filtVesselDiameter_A,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
hold on
x1 = xline(416,'color',colorA,'LineWidth',2);
xline(295,'color',colorA,'LineWidth',2)
xline(342,'color',colorA,'LineWidth',2)
xline(511,'color',colorA,'LineWidth',2)
xline(648,'color',colorA,'LineWidth',2)
xline(697,'color',colorA,'LineWidth',2)
xline(771,'color',colorA,'LineWidth',2)
ylabel('\DeltaD/D (%)')
legend([p1,x1],'Arteriole','Brief waking')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([270,330,390,450,510,570,630,690,750,810,870])
ax23.TickLength = [0.03,0.03];
xlim([270,870]) 
ylim([-30,60])
% cortical LFP
ax4 = subplot(10,1,4);
semilog_imagesc_Manuscript2020(T,F,cortNormS_A,'y')
axis xy
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
ylabel({'Cortical LFP','Freq (Hz)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([270,330,390,450,510,570,630,690,750,810,870])
ax4.TickLength = [0.03,0.03];
xlim([270,870]) 
% hippocampal LFP
ax5 = subplot(10,1,5);
semilog_imagesc_Manuscript2020(T,F,hipNormS_A,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
set(gca,'box','off')
axis tight
xticks([270,330,390,450,510,570,630,690,750,810,870])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax5.TickLength = [0.03,0.03];
xlim([270,870]) 
% Axes properties
ax1Pos = get(ax1,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax1Pos(3:4);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
%% [B] second single trial 2P sleep example
% whisker angle
ax6 = subplot(10,1,6);
plot((1:length(filtWhiskerAngle_B))/MergedData.notes.dsFs,-filtWhiskerAngle_B,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax6.TickLength = [0.03,0.03];
xlim([15,615])
ylim([-10,40])
% vessel diameter
ax78 = subplot(10,1,[7,8]);
p2 = plot((1:length(filtVesselDiameter_B))/MergedData.notes.p2Fs,filtVesselDiameter_B,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
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
legend([p2,x1,x2,x3],'Arteriole','Awake onset','NREM onset','REM onset')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax78.TickLength = [0.03,0.03];
xlim([15,615])
ylim([-30,60])
% cortical LFP
ax9 = subplot(10,1,9);
semilog_imagesc_Manuscript2020(T,F,cortNormS_B,'y')
axis xy
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
ylabel({'Cortical LFP','Freq (Hz)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
ax9.TickLength = [0.03,0.03];
xlim([15,615])
% hippocampal LFP
ax10 = subplot(10,1,10);
semilog_imagesc_Manuscript2020(T,F,hipNormS_B,'y')
axis xy
c10 = colorbar;
ylabel(c10,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
set(gca,'box','off')
axis tight
xticks([15,75,135,195,255,315,375,435,495,555,615])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax10.TickLength = [0.03,0.03];
xlim([15,615])
% Axes properties
ax6Pos = get(ax6,'position');
ax9Pos = get(ax9,'position');
ax10Pos = get(ax10,'position');
ax9Pos(3:4) = ax6Pos(3:4);
ax10Pos(3:4) = ax6Pos(3:4);
set(ax9,'position',ax9Pos);
set(ax10,'position',ax10Pos);
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Figure Panel 2']);
% remove surface subplots because they take forever to render
subplot(10,1,4); cla; subplot(10,1,5); cla; subplot(10,1,9); cla; subplot(10,1,10); cla;
print('-painters','-dpdf','-fillpage',[dirpath 'Figure Panel 2'])
%% subplot figures
figure;
% example 1 cortical LFP
subplot(4,1,1);
semilog_imagesc_Manuscript2020(T,F,cortNormS_A,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([270,870]) 
% example 1 hippocampal LFP
subplot(4,1,2);
semilog_imagesc_Manuscript2020(T,F,hipNormS_A,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([270,870]) 
% example 2 cortical LFP
subplot(4,1,3);
semilog_imagesc_Manuscript2020(T,F,cortNormS_B,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([15,615])
% example 2 hippocampal LFP
subplot(4,1,4);
semilog_imagesc_Manuscript2020(T,F,hipNormS_B,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([15,615])
print('-painters','-dtiffn',[dirpath 'Figure Panel 2 subplot images'])

end
