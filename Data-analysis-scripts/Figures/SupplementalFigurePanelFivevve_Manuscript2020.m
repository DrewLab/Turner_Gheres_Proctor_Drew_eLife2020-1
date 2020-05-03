function [] = SupplementalFigurePanelFivevve_Manuscript2020(rootFolder)
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
%% information and data for first example
animalID_A = 'T126';
dataLocation = [rootFolder '\' animalID_A '\2P Data\'];
cd(dataLocation)
exampleMergedFileID_A = 'T126_RH_200310_12_46_04_011_A3_MergedData.mat';
load(exampleMergedFileID_A,'-mat')
exampleSpecFileID_A = 'T126_RH_200310_12_46_04_011_A3_SpecData.mat';
load(exampleSpecFileID_A,'-mat')
exampleBaselinesFileID_A = 'T126_RestingBaselines.mat';
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
% pressure sensor
filtForceSensor_A = filtfilt(sos1,g1,abs(MergedData.data.forceSensorL));
% EMG
EMG_A = MergedData.data.EMG.data;
normEMG_A = EMG_A - RestingBaselines.manualSelection.EMG.data.(strDay_A);
filtEMG_A = filtfilt(sos1,g1,normEMG_A);
% vessel diameter
vesselDiameter_A = MergedData.data.vesselDiameter.data;
normVesselDiameter_A = (vesselDiameter_A - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_A).(strDay_A))./(RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_A).(strDay_A));
filtVesselDiameter_A = filtfilt(sos2,g2,normVesselDiameter_A)*100;
% cortical and hippocampal spectrograms
cortNormS_A = SpecData.corticalNeural.fiveSec.normS.*100;
hipNormS_A = SpecData.hippocampalNeural.fiveSec.normS.*100;
T_A = SpecData.corticalNeural.fiveSec.T;
F_A = SpecData.corticalNeural.fiveSec.F;
cd(rootFolder)
%% information and data for second example
animalID_B = 'T115';
dataLocation = [rootFolder '\' animalID_B '\2P Data\'];
cd(dataLocation)
exampleMergedFileID_B = 'T115_RH_191125_09_22_43_004_P1_MergedData.mat';
load(exampleMergedFileID_B,'-mat')
exampleSpecFileID_B = 'T115_RH_191125_09_22_43_004_P1_SpecData.mat';
load(exampleSpecFileID_B,'-mat')
exampleBaselinesFileID_B = 'T115_RestingBaselines.mat';
load(exampleBaselinesFileID_B,'-mat')
[~,~,fileDate_B,~,~,vesselID_B] = GetFileInfo2_2P_Manuscript2020(exampleMergedFileID_B);
strDay_B = ConvertDate_2P_Manuscript2020(fileDate_B);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(MergedData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(MergedData.notes.p2Fs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filtWhiskerAngle_B = filtfilt(sos1,g1,MergedData.data.whiskerAngle);
% pressure sensor
filtForceSensor_B = filtfilt(sos1,g1,abs(MergedData.data.forceSensorL));
% EMG
EMG_B = MergedData.data.EMG.data;
normEMG_B = EMG_B - RestingBaselines.manualSelection.EMG.data.(strDay_B);
filtEMG_B = filtfilt(sos1,g1,normEMG_B);
% vessel diameter
vesselDiameter_B = MergedData.data.vesselDiameter.data;
normVesselDiameter_B = (vesselDiameter_B - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_B).(strDay_B))./(RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_B).(strDay_B));
filtVesselDiameter_B = filtfilt(sos2,g2,normVesselDiameter_B)*100;
% cortical and hippocampal spectrograms
cortNormS_B = SpecData.corticalNeural.fiveSec.normS.*100;
hipNormS_B = SpecData.hippocampalNeural.fiveSec.normS.*100;
T_B = SpecData.corticalNeural.fiveSec.T;
F_B = SpecData.corticalNeural.fiveSec.F;
cd(rootFolder)
%% information and data for third example
animalID_C = 'T115';
dataLocation = [rootFolder '\' animalID_C '\2P Data\'];
cd(dataLocation)
exampleMergedFileID_C = 'T115_RH_191119_15_47_21_024_A2_MergedData.mat';
load(exampleMergedFileID_C,'-mat')
exampleSpecFileID_C = 'T115_RH_191119_15_47_21_024_A2_SpecData.mat';
load(exampleSpecFileID_C,'-mat')
exampleBaselinesFileID_C = 'T115_RestingBaselines.mat';
load(exampleBaselinesFileID_C,'-mat')
[~,~,fileDate_C,~,~,vesselID_C] = GetFileInfo2_2P_Manuscript2020(exampleMergedFileID_C);
strDay_C = ConvertDate_2P_Manuscript2020(fileDate_C);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(MergedData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(MergedData.notes.p2Fs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filtWhiskerAngle_C = filtfilt(sos1,g1,MergedData.data.whiskerAngle);
% pressure sensor
filtForceSensor_C = filtfilt(sos1,g1,abs(MergedData.data.forceSensorL));
% EMG
EMG_C = MergedData.data.EMG.data;
normEMG_C = EMG_C - RestingBaselines.manualSelection.EMG.data.(strDay_C);
filtEMG_C = filtfilt(sos1,g1,normEMG_C);
% vessel diameter
vesselDiameter_C = MergedData.data.vesselDiameter.data;
normVesselDiameter_C = (vesselDiameter_C - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_C).(strDay_C))./(RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_C).(strDay_C));
filtVesselDiameter_C = filtfilt(sos2,g2,normVesselDiameter_C)*100;
% cortical and hippocampal spectrograms
cortNormS_C = SpecData.corticalNeural.fiveSec.normS.*100;
hipNormS_C = SpecData.hippocampalNeural.fiveSec.normS.*100;
T_C = SpecData.corticalNeural.fiveSec.T;
F_C = SpecData.corticalNeural.fiveSec.F;
cd(rootFolder)
%% information and data for fourth example
animalID_D = 'T116';
dataLocation = [rootFolder '\' animalID_D '\2P Data\'];
cd(dataLocation)
exampleMergedFileID_D = 'T116_RH_191120_11_18_02_009_A2_MergedData.mat';
load(exampleMergedFileID_D,'-mat')
exampleSpecFileID_D = 'T116_RH_191120_11_18_02_009_A2_SpecData.mat';
load(exampleSpecFileID_D,'-mat')
exampleBaselinesFileID_D = 'T116_RestingBaselines.mat';
load(exampleBaselinesFileID_D,'-mat')
[~,~,fileDate_D,~,~,vesselID_D] = GetFileInfo2_2P_Manuscript2020(exampleMergedFileID_D);
strDay_D = ConvertDate_2P_Manuscript2020(fileDate_D);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1,k1] = butter(4,10/(MergedData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1,k1);
[z2,p2,k2] = butter(4,0.5/(MergedData.notes.p2Fs/2),'low');
[sos2,g2] = zp2sos(z2,p2,k2);
% whisker angle
filtWhiskerAngle_D = filtfilt(sos1,g1,MergedData.data.whiskerAngle);
% pressure sensor
filtForceSensor_D = filtfilt(sos1,g1,abs(MergedData.data.forceSensorL));
% EMG
EMG_D = MergedData.data.EMG.data;
normEMG_D = EMG_D - RestingBaselines.manualSelection.EMG.data.(strDay_D);
filtEMG_D = filtfilt(sos1,g1,normEMG_D);
% vessel diameter
vesselDiameter_D = MergedData.data.vesselDiameter.data;
normVesselDiameter_D = (vesselDiameter_D - RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_D).(strDay_D))./(RestingBaselines.manualSelection.vesselDiameter.data.(vesselID_D).(strDay_D));
filtVesselDiameter_D = filtfilt(sos2,g2,normVesselDiameter_D)*100;
% cortical and hippocampal spectrograms
cortNormS_D = SpecData.corticalNeural.fiveSec.normS.*100;
hipNormS_D = SpecData.hippocampalNeural.fiveSec.normS.*100;
T_D = SpecData.corticalNeural.fiveSec.T;
F_D = SpecData.corticalNeural.fiveSec.F;
cd(rootFolder)
%% Supplemental Figure Panel 3 - Example A
summaryFigure_A = figure;
sgtitle({'Supplemental Figure Panel 3 - Turner Manuscript 2020','Example A'})
% EMG and force sensor
ax1 = subplot(6,1,1);
p1_A = plot((1:length(filtEMG_A))/MergedData.notes.dsFs,filtEMG_A,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3]) 
yyaxis right
p2_A = plot((1:length(filtForceSensor_A))/MergedData.notes.dsFs,filtForceSensor_A,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1_A,p2_A],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600])
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% whisker angle
ax2 = subplot(6,1,2);
plot((1:length(filtWhiskerAngle_A))/MergedData.notes.dsFs,-filtWhiskerAngle_A,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
ax2.TickLength = [0.01,0.01];
xlim([0,600])
ylim([-20,60])
% vessel diameter
ax34 = subplot(6,1,[3,4]);
plot((1:length(filtVesselDiameter_A))/MergedData.notes.p2Fs,filtVesselDiameter_A,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
ylabel('\DeltaD/D (%)')
legend('Arteriole diameter')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
ax34.TickLength = [0.01,0.01];
axis tight
xlim([0,600])
% cortical LFP
ax5 = subplot(6,1,5);
semilog_imagesc_Manuscript2020(T_A,F_A,cortNormS_A,'y')
axis xy
c5_A = colorbar;
ylabel(c5_A,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
ylabel({'Cort LFP','Freq (Hz)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
ax5.TickLength = [0.01,0.01];
xlim([0,600])
% hippocampal LFP
ax6 = subplot(6,1,6);
semilog_imagesc_Manuscript2020(T_A,F_A,hipNormS_A,'y')
axis xy
c6_A = colorbar;
ylabel(c6_A,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax6.TickLength = [0.01,0.01];
xlim([0,600])
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
savefig(summaryFigure_A,[dirpath 'Supplemental Figure Panel 3 Example A']);
% remove surface subplots because they take forever to render
cla(ax5);
set(ax5,'YLim',[1,99]);
cla(ax6);
set(ax6,'YLim',[1,99]);
set(summaryFigure_A,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Supplemental Figure Panel 3 Example A'])
close(summaryFigure_A)
%% subplot figures
summaryFigure_Aimgs = figure;
% example 1 cortical LFP
subplot(2,1,1);
semilog_imagesc_Manuscript2020(T_A,F_A,cortNormS_A,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([0,600])
% example 1 hippocampal LFP
subplot(2,1,2);
semilog_imagesc_Manuscript2020(T_A,F_A,hipNormS_A,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([0,600])
print('-painters','-dtiffn',[dirpath 'Supplemental Figure Panel 3 Example A subplot images'])
close(summaryFigure_Aimgs)
%% Supplemental Figure Panel 3 - Example B
summaryFigure_B = figure;
sgtitle({'Supplemental Figure Panel 3 - Turner Manuscript 2020','Example B'})
% EMG and force sensor
ax1 = subplot(6,1,1);
p1_B = plot((1:length(filtEMG_B))/MergedData.notes.dsFs,filtEMG_B,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3]) 
yyaxis right
p2_B = plot((1:length(filtForceSensor_B))/MergedData.notes.dsFs,filtForceSensor_B,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1_B,p2_B],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600])
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% whisker angle
ax2 = subplot(6,1,2);
plot((1:length(filtWhiskerAngle_B))/MergedData.notes.dsFs,-filtWhiskerAngle_B,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
ax2.TickLength = [0.01,0.01];
xlim([0,600])
ylim([-20,60])
% vessel diameter
ax34 = subplot(6,1,[3,4]);
plot((1:length(filtVesselDiameter_B))/MergedData.notes.p2Fs,filtVesselDiameter_B,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
ylabel('\DeltaD/D (%)')
legend('Pen. Arteriole diameter')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
ax34.TickLength = [0.01,0.01];
axis tight
xlim([0,600])
% cortical LFP
ax5 = subplot(6,1,5);
semilog_imagesc_Manuscript2020(T_B,F_B,cortNormS_B,'y')
axis xy
c5_B = colorbar;
ylabel(c5_B,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
ylabel({'Cort LFP','Freq (Hz)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
ax5.TickLength = [0.01,0.01];
xlim([0,600])
% hippocampal LFP
ax6 = subplot(6,1,6);
semilog_imagesc_Manuscript2020(T_B,F_B,hipNormS_B,'y')
axis xy
c6_B = colorbar;
ylabel(c6_B,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax6.TickLength = [0.01,0.01];
xlim([0,600])
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
savefig(summaryFigure_B,[dirpath 'Supplemental Figure Panel 3 Example B']);
% remove surface subplots because they take forever to render
cla(ax5);
set(ax5,'YLim',[1,99]);
cla(ax6);
set(ax6,'YLim',[1,99]);
set(summaryFigure_B,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Supplemental Figure Panel 3 Example B'])
close(summaryFigure_B)
%% subplot figures
summaryFigure_Bimgs = figure;
% example 1 cortical LFP
subplot(2,1,1);
semilog_imagesc_Manuscript2020(T_B,F_B,cortNormS_B,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([0,600])
% example 1 hippocampal LFP
subplot(2,1,2);
semilog_imagesc_Manuscript2020(T_B,F_B,hipNormS_B,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([0,600])
print('-painters','-dtiffn',[dirpath 'Supplemental Figure Panel 3 Example B subplot images'])
close(summaryFigure_Bimgs)
%% Supplemental Figure Panel 3 - Example C
summaryFigure_C = figure;
sgtitle({'Supplemental Figure Panel 3 - Turner Manuscript 2020','Example C'})
% EMG and force sensor
ax1 = subplot(6,1,1);
p1_C = plot((1:length(filtEMG_C))/MergedData.notes.dsFs,filtEMG_C,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3]) 
yyaxis right
p2_C = plot((1:length(filtForceSensor_C))/MergedData.notes.dsFs,filtForceSensor_C,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1_C,p2_C],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([300,360,420,480,540,600,660,720,780,840,900])
xlim([300,900])
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% whisker angle
ax2 = subplot(6,1,2);
plot((1:length(filtWhiskerAngle_C))/MergedData.notes.dsFs,-filtWhiskerAngle_C,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
ax2.TickLength = [0.01,0.01];
xlim([300,900])
ylim([-20,60])
% vessel diameter
ax34 = subplot(6,1,[3,4]);
p3_C = plot((1:length(filtVesselDiameter_C))/MergedData.notes.p2Fs,filtVesselDiameter_C,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
hold on
x1_C = xline(300,'color',colorA,'LineWidth',2);
x2_C = xline(550,'color',colorB,'LineWidth',2);
xline(600,'color',colorA,'LineWidth',2);
xline(645,'color',colorB,'LineWidth',2);
xline(865,'color',colorA,'LineWidth',2);
ylabel('\DeltaD/D (%)')
legend([p3_C,x1_C,x2_C],'Arteriole diameter','Awake','NREM')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([300,360,420,480,540,600,660,720,780,840,900])
ax34.TickLength = [0.01,0.01];
axis tight
xlim([300,900])
% cortical LFP
ax5 = subplot(6,1,5);
semilog_imagesc_Manuscript2020(T_C,F_C,cortNormS_C,'y')
axis xy
c5_C = colorbar;
ylabel(c5_C,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'Cort LFP','Freq (Hz)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([300,360,420,480,540,600,660,720,780,840,900])
ax5.TickLength = [0.01,0.01];
xlim([300,900])
% hippocampal LFP
ax6 = subplot(6,1,6);
semilog_imagesc_Manuscript2020(T_C,F_C,hipNormS_C,'y')
axis xy
c6_C = colorbar;
ylabel(c6_C,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
axis tight
xticks([300,360,420,480,540,600,660,720,780,840,900])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax6.TickLength = [0.01,0.01];
xlim([300,900])
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
savefig(summaryFigure_C,[dirpath 'Supplemental Figure Panel 3 Example C']);
% remove surface subplots because they take forever to render
cla(ax5);
set(ax5,'YLim',[1,99]);
cla(ax6);
set(ax6,'YLim',[1,99]);
set(summaryFigure_C,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Supplemental Figure Panel 3 Example C'])
close(summaryFigure_C)
%% subplot figures
summaryFigure_Cimgs = figure;
% example 1 cortical LFP
subplot(2,1,1);
semilog_imagesc_Manuscript2020(T_C,F_C,cortNormS_C,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([300,900])
% example 1 hippocampal LFP
subplot(2,1,2);
semilog_imagesc_Manuscript2020(T_C,F_C,hipNormS_C,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([300,900])
print('-painters','-dtiffn',[dirpath 'Supplemental Figure Panel 3 Example C subplot images'])
close(summaryFigure_Cimgs)
%% Supplemental Figure Panel 3 - Example D
summaryFigure_D = figure;
sgtitle({'Supplemental Figure Panel 3 - Turner Manuscript 2020','Example D'})
% EMG and force sensor
ax1 = subplot(6,1,1);
p1_D = plot((1:length(filtEMG_D))/MergedData.notes.dsFs,filtEMG_D,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3]) 
yyaxis right
p2_D = plot((1:length(filtForceSensor_D))/MergedData.notes.dsFs,filtForceSensor_D,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1_D,p2_D],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([200,260,320,380,440,500,560,620,680,740,800])
xlim([200,800])
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% whisker angle
ax2 = subplot(6,1,2);
plot((1:length(filtWhiskerAngle_D))/MergedData.notes.dsFs,-filtWhiskerAngle_D,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([200,260,320,380,440,500,560,620,680,740,800])
ax2.TickLength = [0.01,0.01];
xlim([200,800])
ylim([-20,60])
% vessel diameter
ax34 = subplot(6,1,[3,4]);
p3_D = plot((1:length(filtVesselDiameter_D))/MergedData.notes.p2Fs,filtVesselDiameter_D,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
hold on
x1_D = xline(200,'color',colorA,'LineWidth',2);
xline(265,'color',colorA,'LineWidth',2);
xline(300,'color',colorA,'LineWidth',2);
xline(343,'color',colorA,'LineWidth',2);
xline(420,'color',colorA,'LineWidth',2);
xline(515,'color',colorA,'LineWidth',2);
xline(650,'color',colorA,'LineWidth',2);
xline(670,'color',colorA,'LineWidth',2);
xline(706,'color',colorA,'LineWidth',2);
xline(776,'color',colorA,'LineWidth',2);
ylabel('\DeltaD/D (%)')
legend([p3_D,x1_D],'Arteriole diameter','brief Awake')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([200,260,320,380,440,500,560,620,680,740,800])
ax34.TickLength = [0.01,0.01];
axis tight
xlim([200,800])
% cortical LFP
ax5 = subplot(6,1,5);
semilog_imagesc_Manuscript2020(T_D,F_D,cortNormS_D,'y')
axis xy
c5_D = colorbar;
ylabel(c5_D,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'Cort LFP','Freq (Hz)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([200,260,320,380,440,500,560,620,680,740,800])
ax5.TickLength = [0.01,0.01];
xlim([200,800])
% hippocampal LFP
ax6 = subplot(6,1,6);
semilog_imagesc_Manuscript2020(T_D,F_D,hipNormS_D,'y')
axis xy
c6_D = colorbar;
ylabel(c6_D,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hipp LFP','Freq (Hz)'})
set(gca,'box','off')
axis tight
xticks([200,260,320,380,440,500,560,620,680,740,800])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax6.TickLength = [0.01,0.01];
xlim([200,800])
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
savefig(summaryFigure_D,[dirpath 'Supplemental Figure Panel 3 Example D']);
% remove surface subplots because they take forever to render
cla(ax5);
set(ax5,'YLim',[1,99]);
cla(ax6);
set(ax6,'YLim',[1,99]);
set(summaryFigure_D,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Supplemental Figure Panel 3 Example D'])
close(summaryFigure_D)
%% subplot figures
summaryFigure_Dimgs = figure;
% example 1 cortical LFP
subplot(2,1,1);
semilog_imagesc_Manuscript2020(T_D,F_D,cortNormS_D,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([200,800])
% example 1 hippocampal LFP
subplot(2,1,2);
semilog_imagesc_Manuscript2020(T_D,F_D,hipNormS_D,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([200,800])
print('-painters','-dtiffn',[dirpath 'Supplemental Figure Panel 3 Example D subplot images'])
close(summaryFigure_Dimgs)
end
