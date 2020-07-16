function [AnalysisResults] = FigS13_Manuscript2020(rootFolder,saveFigs,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel S13 for Turner_Kederasetti_Gheres_Proctorostanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

% colorBlack = [(0/256),(0/256),(0/256)];
% colorGrey = [(209/256),(211/256),(212/256)];
colorRfcAwake = [(0/256),(64/256),(64/256)];
colorRfcNREM = [(0/256),(174/256),(239/256)];
% colorRfcREM = [(190/256),(30/256),(45/256)];
% colorRest = [(0/256),(166/256),(81/256)];
% colorWhisk = [(31/256),(120/256),(179/256)];
% colorStim = [(255/256),(28/256),(206/256)];
% colorNREM = [(191/256),(0/256),(255/256)];
% colorREM = [(254/256),(139/256),(0/256)];
% colorAlert = [(255/256),(191/256),(0/256)];
% colorAsleep = [(0/256),(128/256),(255/256)];
% colorAll = [(183/256),(115/256),(51/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% set-up and process data
% information and data for first example
animalID = 'T115';
dataLocation = [rootFolder '\' animalID '\2P Data\'];
cd(dataLocation)
exampleMergedFileID = 'T115_RH_191119_15_47_21_024_A2_MergedData.mat';
load(exampleMergedFileID,'-mat')
exampleSpecFileID = 'T115_RH_191119_15_47_21_024_A2_SpecData.mat';
load(exampleSpecFileID,'-mat')
exampleBaselinesFileID = 'T115_RestingBaselines.mat';
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
%% Fig. S13
summaryFigure = figure('Name','FigS13 (a-e)'); %#ok<*NASGU>
sgtitle('Figure Panel S13 (a-e) Turner Manuscript 2020')
%% EMG and force sensor
ax1 = subplot(6,1,1);
p1 = plot((1:length(filtEMG))/MergedData.notes.dsFs,filtEMG,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3])
yyaxis right
p2 = plot((1:length(filtForceSensor))/MergedData.notes.dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'EMG','Pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([300,360,420,480,540,600,660,720,780,840,900])
xlim([300,900])
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
%% whisker angle
ax2 = subplot(6,1,2);
plot((1:length(filtWhiskerAngle))/MergedData.notes.dsFs,-filtWhiskerAngle,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([300,360,420,480,540,600,660,720,780,840,900])
ax2.TickLength = [0.01,0.01];
xlim([300,900])
ylim([-20,60])
%% vessel diameter
ax34 = subplot(6,1,[3,4]);
p3 = plot((1:length(filtVesselDiameter))/MergedData.notes.p2Fs,filtVesselDiameter,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
hold on
x1 = xline(300,'color',colorRfcAwake,'LineWidth',2);
x2 = xline(550,'color',colorRfcNREM,'LineWidth',2);
xline(600,'color',colorRfcAwake,'LineWidth',2);
xline(645,'color',colorRfcNREM,'LineWidth',2);
xline(865,'color',colorRfcAwake,'LineWidth',2);
ylabel('\DeltaD/D (%)')
legend([p3,x1,x2],'Arteriole diameter','Awake','NREM')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([300,360,420,480,540,600,660,720,780,840,900])
ax34.TickLength = [0.01,0.01];
axis tight
xlim([300,900])
%% cortical LFP
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
xticks([300,360,420,480,540,600,660,720,780,840,900])
ax5.TickLength = [0.01,0.01];
xlim([300,900])
%% hippocampal LFP
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
xticks([300,360,420,480,540,600,660,720,780,840,900])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax6.TickLength = [0.01,0.01];
xlim([300,900])
%% Axes properties
ax1Pos = get(ax1,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax5Pos(3:4) = ax1Pos(3:4);
ax6Pos(3:4) = ax1Pos(3:4);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder '\Summary Figures and Structures\MATLAB Analysis Figures\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'FigS13']);
    % remove surface subplots because they take forever to render
    cla(ax5);
    set(ax5,'YLim',[1,99]);
    cla(ax6);
    set(ax6,'YLim',[1,99]);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'FigS13'])
    close(summaryFigure)
    %% subplot figures
    summaryFigure_imgs = figure;
    % example 1 cortical LFP
    subplot(2,1,1);
    semilog_imagesc_Manuscript2020(T,F,cortNormS,'y')
    caxis([-100,200])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([300,900])
    % example 1 hippocampal LFP
    subplot(2,1,2);
    semilog_imagesc_Manuscript2020(T,F,hipNormS,'y')
    caxis([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([300,900])
    print('-painters','-dtiffn',[dirpath 'FigS13 subplot images'])
    close(summaryFigure_imgs)
    %% Fig. S13
    figure('Name','FigS13 (a-e)');
    sgtitle('Figure Panel S13 (a-e) Turner Manuscript 2020')
    %% EMG and force sensor
    ax1 = subplot(6,1,1);
    p1 = plot((1:length(filtEMG))/MergedData.notes.dsFs,filtEMG,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
    ylabel({'EMG','log10(pwr)'})
    ylim([-2.5,3])
    yyaxis right
    p2 = plot((1:length(filtForceSensor))/MergedData.notes.dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
    ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
    legend([p1,p2],'EMG','Pressure')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    xticks([300,360,420,480,540,600,660,720,780,840,900])
    xlim([300,900])
    ylim([-0.1,2.5])
    ax1.TickLength = [0.01,0.01];
    ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
    ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
    %% whisker angle
    ax2 = subplot(6,1,2);
    plot((1:length(filtWhiskerAngle))/MergedData.notes.dsFs,-filtWhiskerAngle,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
    ylabel({'Whisker','angle (deg)'})
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([300,360,420,480,540,600,660,720,780,840,900])
    ax2.TickLength = [0.01,0.01];
    xlim([300,900])
    ylim([-20,60])
    %% vessel diameter
    ax34 = subplot(6,1,[3,4]);
    p3 = plot((1:length(filtVesselDiameter))/MergedData.notes.p2Fs,filtVesselDiameter,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
    hold on
    x1 = xline(300,'color',colorRfcAwake,'LineWidth',2);
    x2 = xline(550,'color',colorRfcNREM,'LineWidth',2);
    xline(600,'color',colorRfcAwake,'LineWidth',2);
    xline(645,'color',colorRfcNREM,'LineWidth',2);
    xline(865,'color',colorRfcAwake,'LineWidth',2);
    ylabel('\DeltaD/D (%)')
    legend([p3,x1,x2],'Arteriole diameter','Awake','NREM')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    xticks([300,360,420,480,540,600,660,720,780,840,900])
    ax34.TickLength = [0.01,0.01];
    axis tight
    xlim([300,900])
    %% cortical LFP
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
    xticks([300,360,420,480,540,600,660,720,780,840,900])
    ax5.TickLength = [0.01,0.01];
    xlim([300,900])
    %% hippocampal LFP
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
    xticks([300,360,420,480,540,600,660,720,780,840,900])
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    ax6.TickLength = [0.01,0.01];
    xlim([300,900])
    %% Axes properties
    ax1Pos = get(ax1,'position');
    ax5Pos = get(ax5,'position');
    ax6Pos = get(ax6,'position');
    ax5Pos(3:4) = ax1Pos(3:4);
    ax6Pos(3:4) = ax1Pos(3:4);
    set(ax5,'position',ax5Pos);
    set(ax6,'position',ax6Pos);
end

end