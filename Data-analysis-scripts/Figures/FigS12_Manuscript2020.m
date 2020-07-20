function [AnalysisResults] = FigS12_Manuscript2020(rootFolder,saveFigs,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel S12 for Turner_Kederasetti_Gheres_Proctorostanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% set-up and process data
% information and data for first example
if isfield(AnalysisResults,'ExampleTrials') == false
    AnalysisResults.ExampleTrials = [];
end
if isfield(AnalysisResults.ExampleTrials,'T115B') == true
    dsFs = AnalysisResults.ExampleTrials.T115B.dsFs;
    p2Fs = AnalysisResults.ExampleTrials.T115B.p2Fs;
    filtEMG = AnalysisResults.ExampleTrials.T115B.filtEMG;
    filtForceSensor = AnalysisResults.ExampleTrials.T115B.filtForceSensor;
    filtWhiskerAngle = AnalysisResults.ExampleTrials.T115B.filtWhiskerAngle;
    filtVesselDiameter = AnalysisResults.ExampleTrials.T115B.filtVesselDiameter;
    T = AnalysisResults.ExampleTrials.T115B.T;
    F = AnalysisResults.ExampleTrials.T115B.F;
    cortNormS = AnalysisResults.ExampleTrials.T115B.cortNormS;
    hipNormS = AnalysisResults.ExampleTrials.T115B.hipNormS;
else
    animalID = 'T115';
    dataLocation = [rootFolder '\' animalID '\2P Data\'];
    cd(dataLocation)
    exampleMergedFileID = 'T115_RH_191125_09_22_43_004_P1_MergedData.mat';
    load(exampleMergedFileID,'-mat')
    exampleSpecFileID = 'T115_RH_191125_09_22_43_004_P1_SpecData.mat';
    load(exampleSpecFileID,'-mat')
    exampleBaselinesFileID = 'T115_RestingBaselines.mat';
    load(exampleBaselinesFileID,'-mat')
    [~,~,fileDate,~,~,vesselID] = GetFileInfo2_2P_Manuscript2020(exampleMergedFileID);
    strDay = ConvertDate_2P_Manuscript2020(fileDate);
    dsFs = MergedData.notes.dsFs;
    p2Fs = MergedData.notes.p2Fs;
    % setup butterworth filter coefficients for a 1 Hz and 10 Hz lowpass based on the sampling rate
    [z1,p1,k1] = butter(4,10/(dsFs/2),'low');
    [sos1,g1] = zp2sos(z1,p1,k1);
    [z2,p2,k2] = butter(4,0.5/(p2Fs/2),'low');
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
    % update analysis structure
    AnalysisResults.ExampleTrials.T115B.dsFs = dsFs;
    AnalysisResults.ExampleTrials.T115B.p2Fs = p2Fs;
    AnalysisResults.ExampleTrials.T115B.filtEMG = filtEMG;
    AnalysisResults.ExampleTrials.T115B.filtForceSensor = filtForceSensor;
    AnalysisResults.ExampleTrials.T115B.filtWhiskerAngle = filtWhiskerAngle;
    AnalysisResults.ExampleTrials.T115B.filtVesselDiameter = filtVesselDiameter;
    AnalysisResults.ExampleTrials.T115B.T = T;
    AnalysisResults.ExampleTrials.T115B.F = F;
    AnalysisResults.ExampleTrials.T115B.cortNormS = cortNormS;
    AnalysisResults.ExampleTrials.T115B.hipNormS = hipNormS;
    % save results
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end
%% Fig. S12
summaryFigure = figure('Name','FigS12 (a-e)');
sgtitle('Figure Panel S12 (a-e) Turner Manuscript 2020')
%% EMG and force sensor
ax1 = subplot(6,1,1);
p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
ylabel({'EMG','log10(pwr)'})
ylim([-2.5,3])
yyaxis right
p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
legend([p1,p2],'EMG','Pressure')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
xlim([0,600])
ylim([-0.1,2.5])
ax1.TickLength = [0.01,0.01];
ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
%% whisker angle
ax2 = subplot(6,1,2);
plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
ylabel({'Whisker','angle (deg)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([0,60,120,180,240,300,360,420,480,540,600])
ax2.TickLength = [0.01,0.01];
xlim([0,600])
ylim([-20,60])
%% vessel diameter
ax34 = subplot(6,1,[3,4]);
plot((1:length(filtVesselDiameter))/p2Fs,filtVesselDiameter,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
ylabel('\DeltaD/D (%)')
legend('Pen. arteriole diameter')
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
ax34.TickLength = [0.01,0.01];
axis tight
xlim([0,600])
%% cortical LFP
ax5 = subplot(6,1,5);
semilog_imagesc_Manuscript2020(T,F,cortNormS,'y')
axis xy
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-100,100])
ylabel({'Cort LFP','Freq (Hz)'})
set(gca,'Xticklabel',[])
set(gca,'box','off')
axis tight
xticks([0,60,120,180,240,300,360,420,480,540,600])
ax5.TickLength = [0.01,0.01];
xlim([0,600])
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
xticks([0,60,120,180,240,300,360,420,480,540,600])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
ax6.TickLength = [0.01,0.01];
xlim([0,600])
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
    savefig(summaryFigure,[dirpath 'FigS12']);
    % remove surface subplots because they take forever to render
    cla(ax5);
    set(ax5,'YLim',[1,99]);
    cla(ax6);
    set(ax6,'YLim',[1,99]);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'FigS12'])
    close(summaryFigure)
    %% subplot figures
    summaryFigure_imgs = figure;
    % example 1 cortical LFP
    subplot(2,1,1);
    semilog_imagesc_Manuscript2020(T,F,cortNormS,'y')
    caxis([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([0,600])
    % example 1 hippocampal LFP
    subplot(2,1,2);
    semilog_imagesc_Manuscript2020(T,F,hipNormS,'y')
    caxis([-100,100])
    set(gca,'box','off')
    axis xy
    axis tight
    axis off
    xlim([0,600])
    print('-painters','-dtiffn',[dirpath 'FigS12 subplot images'])
    close(summaryFigure_imgs)
    %% Fig. S12
    figure('Name','FigS12 (a-e)');
    sgtitle('Figure Panel S12 (a-e) Turner Manuscript 2020')
    %% EMG and force sensor
    ax1 = subplot(6,1,1);
    p1 = plot((1:length(filtEMG))/dsFs,filtEMG,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5);
    ylabel({'EMG','log10(pwr)'})
    ylim([-2.5,3])
    yyaxis right
    p2 = plot((1:length(filtForceSensor))/dsFs,filtForceSensor,'color',[(256/256),(28/256),(207/256)],'LineWidth',0.5);
    ylabel({'Pressure','(a.u.)'},'rotation',-90,'VerticalAlignment','bottom')
    legend([p1,p2],'EMG','Pressure')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    xlim([0,600])
    ylim([-0.1,2.5])
    ax1.TickLength = [0.01,0.01];
    ax1.YAxis(1).Color = colors_Manuscript2020('rich black');
    ax1.YAxis(2).Color = [(256/256),(28/256),(207/256)];
    %% whisker angle
    ax2 = subplot(6,1,2);
    plot((1:length(filtWhiskerAngle))/dsFs,-filtWhiskerAngle,'color',colors_Manuscript2020('rich black'),'LineWidth',0.5)
    ylabel({'Whisker','angle (deg)'})
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    ax2.TickLength = [0.01,0.01];
    xlim([0,600])
    ylim([-20,60])
    %% vessel diameter
    ax34 = subplot(6,1,[3,4]);
    plot((1:length(filtVesselDiameter))/p2Fs,filtVesselDiameter,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
    ylabel('\DeltaD/D (%)')
    legend('Pen. arteriole diameter')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    ax34.TickLength = [0.01,0.01];
    axis tight
    xlim([0,600])
    %% cortical LFP
    ax5 = subplot(6,1,5);
    semilog_imagesc_Manuscript2020(T,F,cortNormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
    caxis([-100,100])
    ylabel({'Cort LFP','Freq (Hz)'})
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    ax5.TickLength = [0.01,0.01];
    xlim([0,600])
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
    xticks([0,60,120,180,240,300,360,420,480,540,600])
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    ax6.TickLength = [0.01,0.01];
    xlim([0,600])
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
