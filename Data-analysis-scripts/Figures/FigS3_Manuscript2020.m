function [] = FigS3_Manuscript2020(rootFolder,AnalysisResults)
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
colorD = [(31/256),(120/256),(180/256)];  % whisk color
colorE = [(0/256),(256/256),(256/256)];   % Isoflurane color
%% information and data for first example
animalID_A = 'T108';
dataLocation = [rootFolder '\' animalID_A '\Bilateral Imaging\'];
cd(dataLocation)
exampleProcDataFileID_A = 'T108_190822_11_52_51_ProcData.mat';
load(exampleProcDataFileID_A,'-mat')
exampleSpecDataFileID_A = 'T108_190822_11_52_51_SpecDataA.mat';
load(exampleSpecDataFileID_A,'-mat')
exampleBaselineFileID_A = 'T108_RestingBaselines.mat';
load(exampleBaselineFileID_A,'-mat')
[~,fileDate_A,~] = GetFileInfo_IOS_Manuscript2020(exampleProcDataFileID_A);
strDay_A = ConvertDate_IOS_Manuscript2020(fileDate_A);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1_A,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1_A,k1);
[z2,p2_A,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2_A,k2);
% whisker angle
filtWhiskerAngle_A = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% force sensor
filtForceSensor_A = filtfilt(sos1,g1,abs(ProcData.data.forceSensor));
% emg
EMG_A = ProcData.data.EMG.emg;
normEMG_A = EMG_A - RestingBaselines.manualSelection.EMG.emg.(strDay_A);
filtEMG_A = filtfilt(sos1,g1,normEMG_A);
% heart rate
heartRate_A = ProcData.data.heartRate;
% CBV data
LH_HbT_A = ProcData.data.CBV_HbT.adjLH;
filtLH_HbT_A = filtfilt(sos2,g2,LH_HbT_A);
RH_HbT_A = ProcData.data.CBV_HbT.adjRH;
filtRH_HbT_A = filtfilt(sos2,g2,RH_HbT_A);
% cortical and hippocampal spectrograms
cortical_LHnormS_A = SpecData.cortical_LH.normS.*100;
cortical_RHnormS_A = SpecData.cortical_RH.normS.*100;
hippocampusNormS_A = SpecData.hippocampus.normS.*100;
T_A = SpecData.cortical_LH.T;
F_A = SpecData.cortical_LH.F;
% solenoids
LPadSol_A = ProcData.data.solenoids.LPadSol;
RPadSol_A = ProcData.data.solenoids.RPadSol;
AudSol_A = ProcData.data.solenoids.AudSol;
% indeces for solenoids
if max(filtLH_HbT_A) >= max(filtRH_HbT_A)
    LPad_Yvals_A = 1.30*max(filtLH_HbT_A)*ones(size(LPadSol_A));
    RPad_Yvals_A = 1.30*max(filtLH_HbT_A)*ones(size(RPadSol_A));
    Aud_Yvals_A = 1.30*max(filtLH_HbT_A)*ones(size(AudSol_A));
else
    LPad_Yvals_A = 1.30*max(filtRH_HbT_A)*ones(size(LPadSol_A));
    RPad_Yvals_A = 1.30*max(filtRH_HbT_A)*ones(size(RPadSol_A));
    Aud_Yvals_A = 1.30*max(filtRH_HbT_A)*ones(size(AudSol_A));
end
cd(rootFolder)
%% information and data for second example
animalID_B = 'T108';
dataLocation = [rootFolder '\' animalID_B '\Bilateral Imaging\'];
cd(dataLocation)
exampleProcDataFileID_B = 'T122_200218_10_40_55_ProcData.mat';
load(exampleProcDataFileID_B,'-mat')
exampleSpecDataFileID_B = 'T122_200218_10_40_55_SpecDataA.mat';
load(exampleSpecDataFileID_B,'-mat')
exampleBaselineFileID_B = 'T122_RestingBaselines.mat';
load(exampleBaselineFileID_B,'-mat')
[~,fileDate_B,~] = GetFileInfo_IOS_Manuscript2020(exampleProcDataFileID_B);
strDay_B = ConvertDate_IOS_Manuscript2020(fileDate_B);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1_A,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1_A,k1);
[z2,p2_A,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2_A,k2);
% whisker angle
filtWhiskerAngle_B = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% force sensor
filtForceSensor_B = filtfilt(sos1,g1,abs(ProcData.data.forceSensor));
% emg
EMG_B = ProcData.data.EMG.emg;
normEMG_B = EMG_B - RestingBaselines.manualSelection.EMG.emg.(strDay_B);
filtEMG_B = filtfilt(sos1,g1,normEMG_B);
% heart rate
heartRate_B = ProcData.data.heartRate;
% CBV data
LH_HbT_B = ProcData.data.CBV_HbT.adjLH;
filtLH_HbT_B = filtfilt(sos2,g2,LH_HbT_B);
RH_HbT_B = ProcData.data.CBV_HbT.adjRH;
filtRH_HbT_B = filtfilt(sos2,g2,RH_HbT_B);
% cortical and hippocampal spectrograms
cortical_LHnormS_B = SpecData.cortical_LH.normS.*100;
cortical_RHnormS_B = SpecData.cortical_RH.normS.*100;
hippocampusNormS_B = SpecData.hippocampus.normS.*100;
T_B = SpecData.cortical_LH.T;
F_B = SpecData.cortical_LH.F;
% solenoids
LPadSol_B = ProcData.data.solenoids.LPadSol;
RPadSol_B = ProcData.data.solenoids.RPadSol;
AudSol_B = ProcData.data.solenoids.AudSol;
% indeces for solenoids
if max(filtLH_HbT_B) >= max(filtRH_HbT_B)
    LPad_Yvals_B = 1.30*max(filtLH_HbT_B)*ones(size(LPadSol_B));
    RPad_Yvals_B = 1.30*max(filtLH_HbT_B)*ones(size(RPadSol_B));
    Aud_Yvals_B = 1.30*max(filtLH_HbT_B)*ones(size(AudSol_B));
else
    LPad_Yvals_B = 1.30*max(filtRH_HbT_B)*ones(size(LPadSol_B));
    RPad_Yvals_B = 1.30*max(filtRH_HbT_B)*ones(size(RPadSol_B));
    Aud_Yvals_B = 1.30*max(filtRH_HbT_B)*ones(size(AudSol_B));
end
cd(rootFolder)
%% information and data for third example
animalID_C = 'T123';
dataLocation = [rootFolder '\' animalID_C '\Bilateral Imaging\'];
cd(dataLocation)
exampleProcDataFileID_C = 'T123_200224_16_27_59_ProcData.mat';
load(exampleProcDataFileID_C,'-mat')
exampleSpecDataFileID_C = 'T123_200224_16_27_59_SpecDataA.mat';
load(exampleSpecDataFileID_C,'-mat')
exampleBaselineFileID_C = 'T123_RestingBaselines.mat';
load(exampleBaselineFileID_C,'-mat')
[~,fileDate_C,~] = GetFileInfo_IOS_Manuscript2020(exampleProcDataFileID_C);
strDay_C = ConvertDate_IOS_Manuscript2020(fileDate_C);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1_A,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1_A,k1);
[z2,p2_A,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2_A,k2);
% whisker angle
filtWhiskerAngle_C = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% force sensor
filtForceSensor_C = filtfilt(sos1,g1,abs(ProcData.data.forceSensor));
% emg
EMG_C = ProcData.data.EMG.emg;
normEMG_C = EMG_C - RestingBaselines.manualSelection.EMG.emg.(strDay_C);
filtEMG_C = filtfilt(sos1,g1,normEMG_C);
% heart rate
heartRate_C = ProcData.data.heartRate;
% CBV data
LH_HbT_C = ProcData.data.CBV_HbT.adjLH;
filtLH_HbT_C = filtfilt(sos2,g2,LH_HbT_C);
RH_HbT_C = ProcData.data.CBV_HbT.adjRH;
filtRH_HbT_C = filtfilt(sos2,g2,RH_HbT_C);
% cortical and hippocampal spectrograms
cortical_LHnormS_C = SpecData.cortical_LH.normS.*100;
cortical_RHnormS_C = SpecData.cortical_RH.normS.*100;
hippocampusNormS_C = SpecData.hippocampus.normS.*100;
T_C = SpecData.cortical_LH.T;
F_C = SpecData.cortical_LH.F;
cd(rootFolder)
%% information and data for fourth example
animalID_D = 'T108';
dataLocation = [rootFolder '\' animalID_D '\Bilateral Imaging\'];
cd(dataLocation)
exampleProcDataFileID_D = 'T108_190822_13_59_30_ProcData.mat';
load(exampleProcDataFileID_D,'-mat')
exampleSpecDataFileID_D = 'T108_190822_13_59_30_SpecDataA.mat';
load(exampleSpecDataFileID_D,'-mat')
exampleBaselineFileID_D = 'T108_RestingBaselines.mat';
load(exampleBaselineFileID_D,'-mat')
[~,fileDate_D,~] = GetFileInfo_IOS_Manuscript2020(exampleProcDataFileID_D);
strDay_D = ConvertDate_IOS_Manuscript2020(fileDate_D);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1_A,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1_A,k1);
[z2,p2_A,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2_A,k2);
% whisker angle
filtWhiskerAngle_D = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% force sensor
filtForceSensor_D = filtfilt(sos1,g1,abs(ProcData.data.forceSensor));
% emg
EMG_D = ProcData.data.EMG.emg;
normEMG_D = EMG_D - RestingBaselines.manualSelection.EMG.emg.(strDay_D);
filtEMG_D = filtfilt(sos1,g1,normEMG_D);
% heart rate
heartRate_D = ProcData.data.heartRate;
% CBV data
LH_HbT_D = ProcData.data.CBV_HbT.adjLH;
filtLH_HbT_D = filtfilt(sos2,g2,LH_HbT_D);
RH_HbT_D = ProcData.data.CBV_HbT.adjRH;
filtRH_HbT_D = filtfilt(sos2,g2,RH_HbT_D);
% cortical and hippocampal spectrograms
cortical_LHnormS_D = SpecData.cortical_LH.normS.*100;
cortical_RHnormS_D = SpecData.cortical_RH.normS.*100;
hippocampusNormS_D = SpecData.hippocampus.normS.*100;
T_D = SpecData.cortical_LH.T;
F_D = SpecData.cortical_LH.F;
cd(rootFolder)
%% information and data for fifth example
animalID_E = 'T108';
dataLocation = [rootFolder '\' animalID_E '\Bilateral Imaging\'];
cd(dataLocation)
exampleProcDataFileID_E = 'T108_190822_13_23_18_ProcData.mat';
load(exampleProcDataFileID_E,'-mat')
exampleSpecDataFileID_E = 'T108_190822_13_23_18_SpecDataA.mat';
load(exampleSpecDataFileID_E,'-mat')
exampleBaselineFileID_E = 'T108_RestingBaselines.mat';
load(exampleBaselineFileID_E,'-mat')
[~,fileDate_E,~] = GetFileInfo_IOS_Manuscript2020(exampleProcDataFileID_E);
strDay_E = ConvertDate_IOS_Manuscript2020(fileDate_E);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1_A,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1_A,k1);
[z2,p2_A,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2_A,k2);
% whisker angle
filtWhiskerAngle_E = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% force sensor
filtForceSensor_E = filtfilt(sos1,g1,abs(ProcData.data.forceSensor));
% emg
EMG_E = ProcData.data.EMG.emg;
normEMG_E = EMG_E - RestingBaselines.manualSelection.EMG.emg.(strDay_E);
filtEMG_E = filtfilt(sos1,g1,normEMG_E);
% heart rate
heartRate_E = ProcData.data.heartRate;
% CBV data
LH_HbT_E = ProcData.data.CBV_HbT.adjLH;
filtLH_HbT_E = filtfilt(sos2,g2,LH_HbT_E);
RH_HbT_E = ProcData.data.CBV_HbT.adjRH;
filtRH_HbT_E = filtfilt(sos2,g2,RH_HbT_E);
% cortical and hippocampal spectrograms
cortical_LHnormS_E = SpecData.cortical_LH.normS.*100;
cortical_RHnormS_E = SpecData.cortical_RH.normS.*100;
hippocampusNormS_E = SpecData.hippocampus.normS.*100;
T_E = SpecData.cortical_LH.T;
F_E = SpecData.cortical_LH.F;
cd(rootFolder)
%% information and data for sixth example
animalID_F = 'T123';
baselineLocation = [rootFolder '\' animalID_F '\Bilateral Imaging\'];
cd(baselineLocation)
exampleBaselineFileID_F = 'T123_RestingBaselines.mat';
load(exampleBaselineFileID_F,'-mat')
cd(rootFolder)
dataLocation = [rootFolder '\' animalID_F '\Isoflurane Trials\'];
cd(dataLocation)
exampleProcDataFileID_F = 'T123_200304_15_48_29_ProcData.mat';
load(exampleProcDataFileID_F,'-mat')
exampleSpecDataFileID_F = 'T123_200304_15_48_29_SpecDataA.mat';
load(exampleSpecDataFileID_F,'-mat')
[~,fileDate_F,~] = GetFileInfo_IOS_Manuscript2020(exampleProcDataFileID_F);
strDay_F = ConvertDate_IOS_Manuscript2020(fileDate_F);
% setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
[z1,p1_A,k1] = butter(4,10/(ProcData.notes.dsFs/2),'low');
[sos1,g1] = zp2sos(z1,p1_A,k1);
[z2,p2_A,k2] = butter(4,0.5/(ProcData.notes.dsFs/2),'low');
[sos2,g2] = zp2sos(z2,p2_A,k2);
% whisker angle
filtWhiskerAngle_F = filtfilt(sos1,g1,ProcData.data.whiskerAngle);
% force sensor
filtForceSensor_F = filtfilt(sos1,g1,abs(ProcData.data.forceSensor));
% emg
EMG_F = ProcData.data.EMG.emg;
normEMG_F = EMG_F - RestingBaselines.manualSelection.EMG.emg.(strDay_F);
filtEMG_F = filtfilt(sos1,g1,normEMG_F);
% heart rate
heartRate_F = ProcData.data.heartRate;
% CBV data
LH_HbT_F = ProcData.data.CBV_HbT.adjLH;
filtLH_HbT_F = filtfilt(sos2,g2,LH_HbT_F);
RH_HbT_F = ProcData.data.CBV_HbT.adjRH;
filtRH_HbT_F = filtfilt(sos2,g2,RH_HbT_F);
% cortical and hippocampal spectrograms
cortical_LHnormS_F = SpecData.cortical_LH.normS.*100;
cortical_RHnormS_F = SpecData.cortical_RH.normS.*100;
hippocampusNormS_F = SpecData.hippocampus.normS.*100;
T_F = SpecData.cortical_LH.T;
F_F = SpecData.cortical_LH.F;
cd(rootFolder)
%% Supplemental Figure Panel 3 - Example A
summaryFigure_A = figure;
sgtitle({'Supplemental Figure Panel 3 - Turner Manuscript 2020','Example A'})
% EMG and force sensor
ax1 = subplot(7,1,1);
p1_A = plot((1:length(filtEMG_A))/ProcData.notes.dsFs,filtEMG_A,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3])
yyaxis right
p2_A = plot((1:length(filtForceSensor_A))/ProcData.notes.dsFs,filtForceSensor_A,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1_A,p2_A],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600]) 
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% Whisker angle and heart rate
ax2 = subplot(7,1,2);
p3_A = plot((1:length(filtWhiskerAngle_A))/ProcData.notes.dsFs,-filtWhiskerAngle_A,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
xlim([0,600]) 
ylim([-20,60])
yyaxis right
p4_A = plot((1:length(heartRate_A)),heartRate_A,'color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
legend([p3_A,p4_A],'whisker angle','heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600]) 
ylim([5,15])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors_Manuscript2020('rich black');
ax2.YAxis(2).Color = colors_Manuscript2020('deep carrot orange');
% CBV and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6_A = plot((1:length(filtRH_HbT_A))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT_A,'color',colors_Manuscript2020('sapphire'),'LineWidth',1);
hold on
p5_A = plot((1:length(filtLH_HbT_A))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT_A,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
s1_A = scatter(LPadSol_A,LPad_Yvals_A,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
s2_A = scatter(RPadSol_A,RPad_Yvals_A,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
s3_A = scatter(AudSol_A,Aud_Yvals_A,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
ylabel('\DeltaHbT')
legend([p5_A,p6_A,s1_A,s2_A,s3_A],'Left hem','Right hem','LPad sol','RPad sol','Aud sol')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
axis tight
xlim([0,600]) 
ax34.TickLength = [0.01,0.01];
% Left cortical electrode spectrogram
ax5 = subplot(7,1,5);
semilog_imagesc_Manuscript2020(T_A,F_A,cortical_LHnormS_A,'y')
axis xy
c5_A = colorbar;
ylabel(c5_A,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
ylabel({'LH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600]) 
ax5.TickLength = [0.01,0.01];
% Right cortical electrode spectrogram
ax6 = subplot(7,1,6);
semilog_imagesc_Manuscript2020(T_A,F_A,cortical_RHnormS_A,'y')
axis xy
c6_A = colorbar;
ylabel(c6_A,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
ylabel({'RH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600]) 
ax6.TickLength = [0.01,0.01];
% Hippocampal electrode spectrogram
ax7 = subplot(7,1,7);
semilog_imagesc_Manuscript2020(T_A,F_A,hippocampusNormS_A,'y')
c7_A = colorbar;
ylabel(c7_A,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
xlabel('Time (min)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([0,600]) 
ax7.TickLength = [0.01,0.01];
% Axes properties
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
cla(ax7);
set(ax7,'YLim',[1,99]);
set(summaryFigure_A,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Supplemental Figure Panel 3 Example A'])
close(summaryFigure_A)
%% subplot figures
summaryFigure_Aimgs = figure;
% example 1 LH cortical LFP
subplot(3,1,1);
semilog_imagesc_Manuscript2020(T_A,F_A,cortical_LHnormS_A,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([0,600]) 
% example 1 RH cortical LFP
subplot(3,1,2);
semilog_imagesc_Manuscript2020(T_A,F_A,cortical_RHnormS_A,'y')
caxis([-100,100])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([0,600]) 
% example 1 hippocampal LFP
subplot(3,1,3);
semilog_imagesc_Manuscript2020(T_A,F_A,hippocampusNormS_A,'y')
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
ax1 = subplot(7,1,1);
p1_B = plot((1:length(filtEMG_B))/ProcData.notes.dsFs,filtEMG_B,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3])
yyaxis right
p2_B = plot((1:length(filtForceSensor_B))/ProcData.notes.dsFs,filtForceSensor_B,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1_B,p2_B],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([260,320,380,440,500,560,620,680,740,800,860])
xlim([260,860]) 
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% Whisker angle and heart rate
ax2 = subplot(7,1,2);
p3_B = plot((1:length(filtWhiskerAngle_B))/ProcData.notes.dsFs,-filtWhiskerAngle_B,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
ylim([-20,60])
yyaxis right
p4_B = plot((1:length(heartRate_B)),heartRate_B,'color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
legend([p3_B,p4_B],'whisker angle','heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([260,320,380,440,500,560,620,680,740,800,860])
xlim([260,860]) 
ylim([5,15])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors_Manuscript2020('rich black');
ax2.YAxis(2).Color = colors_Manuscript2020('deep carrot orange');
% CBV and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6_B = plot((1:length(filtRH_HbT_B))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT_B,'color',colors_Manuscript2020('sapphire'),'LineWidth',1);
hold on
p5_B = plot((1:length(filtLH_HbT_B))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT_B,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
scatter(LPadSol_B,LPad_Yvals_B,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
scatter(RPadSol_B,RPad_Yvals_B,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
scatter(AudSol_B,Aud_Yvals_B,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
x1_B = xline(260,'color',colorB,'LineWidth',2);
x2_B = xline(650,'color',colorC,'LineWidth',2);
x3_B = xline(808,'color',colorA,'LineWidth',2);
ylabel('\DeltaHbT')
legend([p5_B,p6_B,x3_B,x1_B,x2_B],'Left hem','Right hem','Awake','NREM','REM')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([260,320,380,440,500,560,620,680,740,800,860])
axis tight
xlim([260,860]) 
ax34.TickLength = [0.01,0.01];
% Left cortical electrode spectrogram
ax5 = subplot(7,1,5);
semilog_imagesc_Manuscript2020(T_B,F_B,cortical_LHnormS_B,'y')
axis xy
c5_B = colorbar;
ylabel(c5_B,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'LH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([260,320,380,440,500,560,620,680,740,800,860])
xlim([260,860]) 
ax5.TickLength = [0.01,0.01];
% Right cortical electrode spectrogram
ax6 = subplot(7,1,6);
semilog_imagesc_Manuscript2020(T_B,F_B,cortical_RHnormS_B,'y')
axis xy
c6_B = colorbar;
ylabel(c6_B,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'RH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([260,320,380,440,500,560,620,680,740,800,860])
xlim([260,860]) 
ax6.TickLength = [0.01,0.01];
% Hippocampal electrode spectrogram
ax7 = subplot(7,1,7);
semilog_imagesc_Manuscript2020(T_B,F_B,hippocampusNormS_B,'y')
c7_B = colorbar;
ylabel(c7_B,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (min)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([260,320,380,440,500,560,620,680,740,800,860])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([260,860])
ax7.TickLength = [0.01,0.01];
% Axes properties
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
cla(ax7);
set(ax7,'YLim',[1,99]);
set(summaryFigure_B,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Supplemental Figure Panel 3 Example B'])
close(summaryFigure_B)
%% subplot figures
summaryFigure_Bimgs = figure;
% example 2 LH cortical LFP
subplot(3,1,1);
semilog_imagesc_Manuscript2020(T_B,F_B,cortical_LHnormS_B,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([260,860]) 
% example 2 RH cortical LFP
subplot(3,1,2);
semilog_imagesc_Manuscript2020(T_B,F_B,cortical_RHnormS_B,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([260,860]) 
% example 2 hippocampal LFP
subplot(3,1,3);
semilog_imagesc_Manuscript2020(T_B,F_B,hippocampusNormS_B,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([260,860]) 
print('-painters','-dtiffn',[dirpath 'Supplemental Figure Panel 3 Example B subplot images'])
close(summaryFigure_Bimgs)
%% Supplemental Figure Panel 3 - Example C
summaryFigure_C = figure;
sgtitle({'Supplemental Figure Panel 3 - Turner Manuscript 2020','Example C'})
% EMG and force sensor
ax1 = subplot(7,1,1);
p1_C = plot((1:length(filtEMG_C))/ProcData.notes.dsFs,filtEMG_C,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3])
yyaxis right
p2_C = plot((1:length(filtForceSensor_C))/ProcData.notes.dsFs,filtForceSensor_C,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1_C,p2_C],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600]) 
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% Whisker angle and heart rate
ax2 = subplot(7,1,2);
p3_C = plot((1:length(filtWhiskerAngle_C))/ProcData.notes.dsFs,-filtWhiskerAngle_C,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
ylim([-20,60])
yyaxis right
p4_C = plot((1:length(heartRate_C)),heartRate_C,'color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
legend([p3_C,p4_C],'whisker angle','heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600]) 
ylim([5,15])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors_Manuscript2020('rich black');
ax2.YAxis(2).Color = colors_Manuscript2020('deep carrot orange');
% CBV and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6_C = plot((1:length(filtRH_HbT_C))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT_C,'color',colors_Manuscript2020('sapphire'),'LineWidth',1);
hold on
p5_C = plot((1:length(filtLH_HbT_C))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT_C,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
x1_C = xline(0,'color',colorB,'LineWidth',2);
ylabel('\DeltaHbT')
legend([p5_C,p6_C,x1_C],'Left hem','Right hem','NREM')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
axis tight
xlim([0,600]) 
ax34.TickLength = [0.01,0.01];
% Left cortical electrode spectrogram
ax5 = subplot(7,1,5);
semilog_imagesc_Manuscript2020(T_C,F_C,cortical_LHnormS_C,'y')
axis xy
c5_C = colorbar;
ylabel(c5_C,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'LH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600]) 
ax5.TickLength = [0.01,0.01];
% Right cortical electrode spectrogram
ax6 = subplot(7,1,6);
semilog_imagesc_Manuscript2020(T_C,F_C,cortical_RHnormS_C,'y')
axis xy
c6_C = colorbar;
ylabel(c6_C,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'RH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600]) 
ax6.TickLength = [0.01,0.01];
% Hippocampal electrode spectrogram
ax7 = subplot(7,1,7);
semilog_imagesc_Manuscript2020(T_C,F_C,hippocampusNormS_C,'y')
c7_C = colorbar;
ylabel(c7_C,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (min)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([0,600]) 
ax7.TickLength = [0.01,0.01];
% Axes properties
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
cla(ax7);
set(ax7,'YLim',[1,99]);
set(summaryFigure_C,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Supplemental Figure Panel 3 Example C'])
close(summaryFigure_C)
%% subplot figures
summaryFigure_Cimgs = figure;
% example 3 LH cortical LFP
subplot(3,1,1);
semilog_imagesc_Manuscript2020(T_C,F_C,cortical_LHnormS_C,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([0,600]) 
% example 3 RH cortical LFP
subplot(3,1,2);
semilog_imagesc_Manuscript2020(T_C,F_C,cortical_RHnormS_C,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([0,600]) 
% example 3 hippocampal LFP
subplot(3,1,3);
semilog_imagesc_Manuscript2020(T_C,F_C,hippocampusNormS_C,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([0,600]) 
print('-painters','-dtiffn',[dirpath 'Supplemental Figure Panel 3 Example C subplot images'])
close(summaryFigure_Cimgs)
%% Supplemental Figure Panel 3 - Example D
summaryFigure_D = figure;
sgtitle({'Supplemental Figure Panel 3 - Turner Manuscript 2020','Example D'})
% EMG and force sensor
ax1 = subplot(7,1,1);
p1_D = plot((1:length(filtEMG_D))/ProcData.notes.dsFs,filtEMG_D,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3])
yyaxis right
p2_D = plot((1:length(filtForceSensor_D))/ProcData.notes.dsFs,filtForceSensor_D,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1_D,p2_D],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([35,95,155,215,275,335,395,455,515,575,635])
xlim([35,635]) 
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% Whisker angle and heart rate
ax2 = subplot(7,1,2);
p3_D = plot((1:length(filtWhiskerAngle_D))/ProcData.notes.dsFs,-filtWhiskerAngle_D,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
ylim([-20,60])
yyaxis right
p4_D = plot((1:length(heartRate_D)),heartRate_D,'color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
legend([p3_D,p4_D],'whisker angle','heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([35,95,155,215,275,335,395,455,515,575,635])
xlim([35,635]) 
ylim([5,15])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors_Manuscript2020('rich black');
ax2.YAxis(2).Color = colors_Manuscript2020('deep carrot orange');
% CBV and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6_D = plot((1:length(filtRH_HbT_D))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT_D,'color',colors_Manuscript2020('sapphire'),'LineWidth',1);
hold on
p5_D = plot((1:length(filtLH_HbT_D))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT_D,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
x1_D = xline(35,'color',colorA,'LineWidth',2);
x2_D = xline(175,'color',colorB,'LineWidth',2);
x3_D = xline(445,'color',colorC,'LineWidth',2);
xline(610,'color',colorA,'LineWidth',2);
ylabel('\DeltaHbT')
legend([p5_D,p6_D,x3_D,x1_D,x2_D],'Left hem','Right hem','Awake','NREM','REM')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([35,95,155,215,275,335,395,455,515,575,635])
axis tight
xlim([35,635]) 
ax34.TickLength = [0.01,0.01];
% Left cortical electrode spectrogram
ax5 = subplot(7,1,5);
semilog_imagesc_Manuscript2020(T_D,F_D,cortical_LHnormS_D,'y')
axis xy
c5_D = colorbar;
ylabel(c5_D,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'LH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([35,95,155,215,275,335,395,455,515,575,635])
xlim([35,635]) 
ax5.TickLength = [0.01,0.01];
% Right cortical electrode spectrogram
ax6 = subplot(7,1,6);
semilog_imagesc_Manuscript2020(T_D,F_D,cortical_RHnormS_D,'y')
axis xy
c6_D = colorbar;
ylabel(c6_D,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'RH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([35,95,155,215,275,335,395,455,515,575,635])
xlim([35,635]) 
ax6.TickLength = [0.01,0.01];
% Hippocampal electrode spectrogram
ax7 = subplot(7,1,7);
semilog_imagesc_Manuscript2020(T_D,F_D,hippocampusNormS_D,'y')
c7_D = colorbar;
ylabel(c7_D,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (min)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([35,95,155,215,275,335,395,455,515,575,635])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([35,635]) 
ax7.TickLength = [0.01,0.01];
% Axes properties
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
cla(ax7);
set(ax7,'YLim',[1,99]);
set(summaryFigure_D,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Supplemental Figure Panel 3 Example D'])
close(summaryFigure_D)
%% subplot figures
summaryFigure_Dimgs = figure;
% example 4 LH cortical LFP
subplot(3,1,1);
semilog_imagesc_Manuscript2020(T_D,F_D,cortical_LHnormS_D,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([35,635]) 
% example 4 RH cortical LFP
subplot(3,1,2);
semilog_imagesc_Manuscript2020(T_D,F_D,cortical_RHnormS_D,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([35,635]) 
% example 4 hippocampal LFP
subplot(3,1,3);
semilog_imagesc_Manuscript2020(T_D,F_D,hippocampusNormS_D,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([35,635]) 
print('-painters','-dtiffn',[dirpath 'Supplemental Figure Panel 3 Example D subplot images'])
close(summaryFigure_Dimgs)
%% Supplemental Figure Panel 3 - Example E
summaryFigure_E = figure;
sgtitle({'Supplemental Figure Panel 3 - Turner Manuscript 2020','Example E'})
% EMG and force sensor
ax1 = subplot(7,1,1);
p1_E = plot((1:length(filtEMG_E))/ProcData.notes.dsFs,filtEMG_E,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3])
yyaxis right
p2_E = plot((1:length(filtForceSensor_E))/ProcData.notes.dsFs,filtForceSensor_E,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1_E,p2_E],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
xlim([300,900]) 
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% Whisker angle and heart rate
ax2 = subplot(7,1,2);
p3_E = plot((1:length(filtWhiskerAngle_E))/ProcData.notes.dsFs,-filtWhiskerAngle_E,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
ylim([-20,60])
yyaxis right
p4_E = plot((1:length(heartRate_E)),heartRate_E,'color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
legend([p3_E,p4_E],'whisker angle','heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
xlim([300,900]) 
ylim([5,15])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors_Manuscript2020('rich black');
ax2.YAxis(2).Color = colors_Manuscript2020('deep carrot orange');
% CBV and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6_E = plot((1:length(filtRH_HbT_E))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT_E,'color',colors_Manuscript2020('sapphire'),'LineWidth',1);
hold on
p5_E = plot((1:length(filtLH_HbT_E))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT_E,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
x1_E = xline(300,'color',colorB,'LineWidth',2);
x2_E = xline(600,'color',colorC,'LineWidth',2);
x3_E = xline(707,'color',colorA,'LineWidth',2);
ylabel('\DeltaHbT')
legend([p5_E,p6_E,x3_E,x1_E,x2_E],'Left hem','Right hem','Awake','NREM','REM')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
axis tight
xlim([300,900]) 
ax34.TickLength = [0.01,0.01];
% Left cortical electrode spectrogram
ax5 = subplot(7,1,5);
semilog_imagesc_Manuscript2020(T_E,F_E,cortical_LHnormS_E,'y')
axis xy
c5_E = colorbar;
ylabel(c5_E,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'LH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
xlim([300,900]) 
ax5.TickLength = [0.01,0.01];
% Right cortical electrode spectrogram
ax6 = subplot(7,1,6);
semilog_imagesc_Manuscript2020(T_E,F_E,cortical_RHnormS_E,'y')
axis xy
c6_E = colorbar;
ylabel(c6_E,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'RH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
xlim([300,900]) 
ax6.TickLength = [0.01,0.01];
% Hippocampal electrode spectrogram
ax7 = subplot(7,1,7);
semilog_imagesc_Manuscript2020(T_E,F_E,hippocampusNormS_E,'y')
c7_E = colorbar;
ylabel(c7_E,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (min)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([300,900]) 
ax7.TickLength = [0.01,0.01];
% Axes properties
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
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure_E,[dirpath 'Supplemental Figure Panel 3 Example E']);
% remove surface subplots because they take forever to render
cla(ax5);
set(ax5,'YLim',[1,99]);
cla(ax6);
set(ax6,'YLim',[1,99]);
cla(ax7);
set(ax7,'YLim',[1,99]);
set(summaryFigure_E,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Supplemental Figure Panel 3 Example E'])
close(summaryFigure_E)
%% subplot figures
summaryFigure_Eimgs = figure;
% example 5 LH cortical LFP
subplot(3,1,1);
semilog_imagesc_Manuscript2020(T_E,F_E,cortical_LHnormS_E,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([300,900]) 
% example 5 RH cortical LFP
subplot(3,1,2);
semilog_imagesc_Manuscript2020(T_E,F_E,cortical_RHnormS_E,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([300,900]) 
% example 5 hippocampal LFP
subplot(3,1,3);
semilog_imagesc_Manuscript2020(T_E,F_E,hippocampusNormS_E,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([300,900]) 
print('-painters','-dtiffn',[dirpath 'Supplemental Figure Panel 3 Example E subplot images'])
close(summaryFigure_Eimgs)
%% Supplemental Figure Panel 3 - Example F
summaryFigure_F = figure;
sgtitle({'Supplemental Figure Panel 3 - Turner Manuscript 2020','Example F'})
% EMG and force sensor
ax1 = subplot(7,1,1);
p1_F = plot((1:length(filtEMG_F))/ProcData.notes.dsFs,filtEMG_F,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3])
yyaxis right
p2_F = plot((1:length(filtForceSensor_F))/ProcData.notes.dsFs,filtForceSensor_F,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1_F,p2_F],'EMG','pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
% Whisker angle and heart rate
ax2 = subplot(7,1,2);
p3_F = plot((1:length(filtWhiskerAngle_F))/ProcData.notes.dsFs,-filtWhiskerAngle_F,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'Whisker','angle (deg)'})
ylim([-20,60])
yyaxis right
p4_F = plot((1:length(heartRate_F)),heartRate_F,'color',colors_Manuscript2020('deep carrot orange'),'LineWidth',0.5);
ylabel('Heart rate (Hz)','rotation',-90,'VerticalAlignment','bottom')
legend([p3_F,p4_F],'whisker angle','heart rate')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
ylim([5,15])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors_Manuscript2020('rich black');
ax2.YAxis(2).Color = colors_Manuscript2020('deep carrot orange');
% CBV and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6_F = plot((1:length(filtRH_HbT_F))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT_F,'color',colors_Manuscript2020('sapphire'),'LineWidth',1);
hold on
p5_F = plot((1:length(filtLH_HbT_F))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT_F,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
x1_F = xline(130,'color',colorB,'LineWidth',2);
x2_F = xline(360,'color',colorA,'LineWidth',2);
x3_F = xline(465,'color','k','LineWidth',2);
ylabel('\DeltaHbT')
legend([p5_F,p6_F,x2_F,x1_F,x3_F],'Left hem','Right hem','Awake','NREM','Isoflurane')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
axis tight
xlim([130,730]) 
ax34.TickLength = [0.01,0.01];
% Left cortical electrode spectrogram
ax5 = subplot(7,1,5);
semilog_imagesc_Manuscript2020(T_F,F_F,cortical_LHnormS_F,'y')
axis xy
c5_F = colorbar;
ylabel(c5_F,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'LH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
ax5.TickLength = [0.01,0.01];
% Right cortical electrode spectrogram
ax6 = subplot(7,1,6);
semilog_imagesc_Manuscript2020(T_F,F_F,cortical_RHnormS_F,'y')
axis xy
c6_F = colorbar;
ylabel(c6_F,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
ylabel({'RH cortical LFP','Freq (Hz)'})
set(gca,'Yticklabel','10^1')
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
ax6.TickLength = [0.01,0.01];
% Hippocampal electrode spectrogram
ax7 = subplot(7,1,7);
semilog_imagesc_Manuscript2020(T_F,F_F,hippocampusNormS_F,'y')
c7_F = colorbar;
ylabel(c7_F,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,200])
xlabel('Time (min)')
ylabel({'Hippocampal LFP','Freq (Hz)'})
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([130,730]) 
ax7.TickLength = [0.01,0.01];
% Axes properties
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
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure_F,[dirpath 'Supplemental Figure Panel 3 Example F']);
% remove surface subplots because they take forever to render
cla(ax5);
set(ax5,'YLim',[1,99]);
cla(ax6);
set(ax6,'YLim',[1,99]);
cla(ax7);
set(ax7,'YLim',[1,99]);
set(summaryFigure_F,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Supplemental Figure Panel 3 Example F'])
close(summaryFigure_F)
%% subplot figures
summaryFigure_Fimgs = figure;
% example 6 LH cortical LFP
subplot(3,1,1);
semilog_imagesc_Manuscript2020(T_F,F_F,cortical_LHnormS_F,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([130,730]) 
% example 6 RH cortical LFP
subplot(3,1,2);
semilog_imagesc_Manuscript2020(T_F,F_F,cortical_RHnormS_F,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([130,730]) 
% example 6 hippocampal LFP
subplot(3,1,3);
semilog_imagesc_Manuscript2020(T_F,F_F,hippocampusNormS_F,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([130,730]) 
print('-painters','-dtiffn',[dirpath 'Supplemental Figure Panel 3 Example F subplot images'])
close(summaryFigure_Fimgs)
%% Isoflurane figure
IOSanimalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
isoAnimalIDs = {'T108','T109','T110','T111','T119','T120','T121','T122','T123'};
modelType = 'Forest';
numComparisons = 3;
%% Mean HbT and heart rate comparison between behaviors
% cd through each animal's directory and extract the appropriate analysis results
behavFields = {'Rest','Whisk','NREM','REM'};
for a = 1:length(IOSanimalIDs)
    animalID = IOSanimalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
            data.(behavField).CBV_HbT.meanLH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjLH);
            data.(behavField).CBV_HbT.meanRH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjRH);
            data.(behavField).CBV_HbT.allLH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjLH;
            data.(behavField).CBV_HbT.allRH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.adjRH;
        elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            data.(behavField).CBV_HbT.meanLH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjLH);
            data.(behavField).CBV_HbT.meanRH(a,1) = mean(AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjRH);
            data.(behavField).CBV_HbT.allLH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjLH;
            data.(behavField).CBV_HbT.allRH{a,1} = AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.adjRH;
        end
        data.(behavField).CBV_HbT.animalID{a,1} = animalID;
        data.(behavField).CBV_HbT.behavior{a,1} = behavField;
        data.(behavField).CBV_HbT.LH{a,1} = 'LH';
        data.(behavField).CBV_HbT.RH{a,1} = 'RH';
    end
end
% pull out isoflurane data
for a = 1:length(isoAnimalIDs)
    isoAnimalID = isoAnimalIDs{1,a};
    data.Iso.CBV_HbT.meanLH(a,1) = AnalysisResults.(isoAnimalID).MeanCBV.Iso.CBV_HbT.adjLH;
    data.Iso.CBV_HbT.meanRH(a,1) = AnalysisResults.(isoAnimalID).MeanCBV.Iso.CBV_HbT.adjRH;
    data.Iso.CBV_HbT.animalID{a,1} = isoAnimalID;
    data.Iso.CBV_HbT.behavior{a,1} = 'Iso';
    data.Iso.CBV_HbT.LH{a,1} = 'LH';
    data.Iso.CBV_HbT.RH{a,1} = 'RH';
end
% take the mean and standard deviation of each set of signals
behavFields = {'Rest','Whisk','NREM','REM','Iso'};
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    data.(behavField).CBV_HbT.Comb = cat(1,data.(behavField).CBV_HbT.meanLH,data.(behavField).CBV_HbT.meanRH);
    data.(behavField).CBV_HbT.catAllLH = [];
    data.(behavField).CBV_HbT.catAllRH = [];
    if strcmp(behavField,'Iso') == false
        for h = 1:length(data.(behavField).CBV_HbT.allLH)
            data.(behavField).CBV_HbT.catAllLH = cat(1,data.(behavField).CBV_HbT.catAllLH,data.(behavField).CBV_HbT.allLH{h,1});
            data.(behavField).CBV_HbT.catAllRH = cat(1,data.(behavField).CBV_HbT.catAllRH,data.(behavField).CBV_HbT.allRH{h,1});
        end
        data.(behavField).CBV_HbT.allComb = cat(1,data.(behavField).CBV_HbT.catAllLH,data.(behavField).CBV_HbT.catAllRH);
    end
    data.(behavField).CBV_HbT.meanCBV = mean(data.(behavField).CBV_HbT.Comb);
    data.(behavField).CBV_HbT.stdCBV = std(data.(behavField).CBV_HbT.Comb,0,1);
end
% statistics - linear mixed effects model
% HbT
HbT_alphaConf = 0.005;
HbTtableSize = cat(1,data.Rest.CBV_HbT.meanLH,data.Rest.CBV_HbT.meanRH,data.Whisk.CBV_HbT.meanLH,data.Whisk.CBV_HbT.meanRH,...
    data.NREM.CBV_HbT.meanLH,data.NREM.CBV_HbT.meanRH,data.REM.CBV_HbT.meanLH,data.REM.CBV_HbT.meanRH,data.Iso.CBV_HbT.meanLH,data.Iso.CBV_HbT.meanRH);
HbTTable = table('Size',[size(HbTtableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','HbT','Behavior','Hemisphere'});
HbTTable.Mouse = cat(1,data.Rest.CBV_HbT.animalID,data.Rest.CBV_HbT.animalID,data.Whisk.CBV_HbT.animalID,data.Whisk.CBV_HbT.animalID,...
    data.NREM.CBV_HbT.animalID,data.NREM.CBV_HbT.animalID,data.REM.CBV_HbT.animalID,data.REM.CBV_HbT.animalID,data.Iso.CBV_HbT.animalID,data.Iso.CBV_HbT.animalID);
HbTTable.HbT = cat(1,data.Rest.CBV_HbT.meanLH,data.Rest.CBV_HbT.meanRH,data.Whisk.CBV_HbT.meanLH,data.Whisk.CBV_HbT.meanRH,...
    data.NREM.CBV_HbT.meanLH,data.NREM.CBV_HbT.meanRH,data.REM.CBV_HbT.meanLH,data.REM.CBV_HbT.meanRH,data.Iso.CBV_HbT.meanLH,data.Iso.CBV_HbT.meanRH);
HbTTable.Behavior = cat(1,data.Rest.CBV_HbT.behavior,data.Rest.CBV_HbT.behavior,data.Whisk.CBV_HbT.behavior,data.Whisk.CBV_HbT.behavior,...
    data.NREM.CBV_HbT.behavior,data.NREM.CBV_HbT.behavior,data.REM.CBV_HbT.behavior,data.REM.CBV_HbT.behavior,data.Iso.CBV_HbT.behavior,data.Iso.CBV_HbT.behavior);
HbTTable.Hemisphere = cat(1,data.Rest.CBV_HbT.LH,data.Rest.CBV_HbT.RH,data.Whisk.CBV_HbT.LH,data.Whisk.CBV_HbT.RH,...
    data.NREM.CBV_HbT.LH,data.NREM.CBV_HbT.RH,data.REM.CBV_HbT.LH,data.REM.CBV_HbT.RH,data.Iso.CBV_HbT.LH,data.Iso.CBV_HbT.RH);
HbTFitFormula = 'HbT ~ 1 + Behavior + (1|Mouse) + (1|Hemisphere)';
HbTStats = fitglme(HbTTable,HbTFitFormula);
HbTCI = coefCI(HbTStats,'Alpha',(HbT_alphaConf/numComparisons));
%% [A] Mean HbT during different behaviors
summaryFigure = figure;
HbT_xInds = ones(1,length(IOSanimalIDs)*2);
Iso_xInds = ones(1,length(isoAnimalIDs)*2);
s1 = scatter(HbT_xInds*1,data.Rest.CBV_HbT.Comb,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Rest.CBV_HbT.meanCBV,data.Rest.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(HbT_xInds*2,data.Whisk.CBV_HbT.Comb,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Whisk.CBV_HbT.meanCBV,data.Whisk.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(HbT_xInds*3,data.NREM.CBV_HbT.Comb,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.NREM.CBV_HbT.meanCBV,data.NREM.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(HbT_xInds*4,data.REM.CBV_HbT.Comb,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.REM.CBV_HbT.meanCBV,data.REM.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(Iso_xInds*5,data.Iso.CBV_HbT.Comb,75,'MarkerEdgeColor','k','MarkerFaceColor',colorE,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,data.Iso.CBV_HbT.meanCBV,data.Iso.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
title({'[A] Mean \DeltaHbT (\muM)','during arousal-states',''})
ylabel('\DeltaHbT (\muM)')
legend([s1,s2,s3,s4,s5],'Awake Rest','Whisking','NREM','REM','Isoflurane','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
ylim([-10,250])
set(gca,'box','off')
set(gca,'TickLength',[0.03,0.03]);
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
set(summaryFigure,'PaperPositionMode','auto');
savefig(summaryFigure,[dirpath 'Supplemental Figure Panel 3 - isoflurane']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Supplemental Figure Panel 3 - isoflurane'])
%% statistical diary
diaryFile = [dirpath 'SupplementalFigurePanel3_Statistics.txt'];
if exist(diaryFile,'file') == true
    delete(diaryFile)
end
diary(diaryFile)
diary on
% HbT statistical diary
disp('======================================================================================================================')
disp('[A] Generalized linear mixed-effects model statistics for mean HbT during Rest, Whisking, NREM, and REM, Isoflurane')
disp('======================================================================================================================')
disp(HbTStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.005 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HbTCI(1,:))])
disp(['Whisk: ' num2str(HbTCI(2,:))])
disp(['NREM: ' num2str(HbTCI(3,:))])
disp(['REM: ' num2str(HbTCI(4,:))])
disp(['Isoflurane: ' num2str(HbTCI(5,:))])
diary off
end
