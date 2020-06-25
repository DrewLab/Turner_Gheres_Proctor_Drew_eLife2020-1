function [] = FigS16_Manuscript2020_fin(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel S16 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% Set-up and process data for Fig S8 (a-f)
IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
Iso_AnimalIDs = {'T108','T109','T110','T111','T119','T120','T121','T122','T123'};
modelType = 'Forest';
colorA = [(0/256),(166/256),(81/256)];    % rest color
colorB = [(191/256),(0/256),(255/256)];   % NREM color
colorC = [(254/256),(139/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(179/256)];  % whisk color
colorE = [(256/256),(28/256),(207/256)];  % stim color
colorF = [(0/256),(256/256),(256/256)];   % isoflurane color
% information and data for first example
animalID = 'T123';
baselineLocation = [rootFolder '\' animalID '\Bilateral Imaging\'];
cd(baselineLocation)
exampleBaselineFileID = 'T123_RestingBaselines.mat';
load(exampleBaselineFileID,'-mat')
cd(rootFolder)
dataLocation = [rootFolder '\' animalID '\Isoflurane Trials\'];
cd(dataLocation)
exampleProcDataFileID = 'T123_200304_15_48_29_ProcData.mat';
load(exampleProcDataFileID,'-mat')
exampleSpecDataFileID = 'T123_200304_15_48_29_SpecDataA.mat';
load(exampleSpecDataFileID,'-mat')
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
% emg
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
%% Figure panel S16
summaryFigure = figure('Name','FigS16 (a-f)');
sgtitle('Figure Panel S16 (a-f) Turner Manuscript 2020')
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
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
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
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
ylim([5,15])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors_Manuscript2020('rich black');
ax2.YAxis(2).Color = colors_Manuscript2020('deep carrot orange');
%% CBV and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6 = plot((1:length(filtRH_HbT))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT,'color',colors_Manuscript2020('sapphire'),'LineWidth',1);
hold on
p5 = plot((1:length(filtLH_HbT))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
x1 = xline(130,'color',colorB,'LineWidth',2);
x2 = xline(360,'color',colorA,'LineWidth',2);
x3 = xline(465,'color','k','LineWidth',2);
ylabel('\DeltaHbT')
legend([p5,p6,x2,x1,x3],'Left hem','Right hem','Awake','NREM','Isoflurane')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
axis tight
xlim([130,730]) 
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
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
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
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
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
xticks([130,190,250,310,370,430,490,550,610,670,730])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([130,730]) 
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
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'FigS16_A']);
% remove surface subplots because they take forever to render
cla(ax5);
set(ax5,'YLim',[1,99]);
cla(ax6);
set(ax6,'YLim',[1,99]);
cla(ax7);
set(ax7,'YLim',[1,99]);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'FigS16_A'])
close(summaryFigure)
%% subplot figures
summaryFigure_imgs = figure;
% example 6 LH cortical LFP
subplot(3,1,1);
semilog_imagesc_Manuscript2020(T,F,cortical_LHnormS,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([130,730]) 
% example 6 RH cortical LFP
subplot(3,1,2);
semilog_imagesc_Manuscript2020(T,F,cortical_RHnormS,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([130,730]) 
% example 6 hippocampal LFP
subplot(3,1,3);
semilog_imagesc_Manuscript2020(T,F,hippocampusNormS,'y')
caxis([-100,200])
set(gca,'box','off')
axis xy
axis tight
axis off
xlim([130,730]) 
print('-painters','-dtiffn',[dirpath 'FigS16 subplot images'])
close(summaryFigure_imgs)
%% Figure panel S16
figure('Name','FigS16 (a-f)');
sgtitle('Figure Panel S16 (a-f) Turner Manuscript 2020')
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
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
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
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
ylim([5,15])
ax2.TickLength = [0.01,0.01];
ax2.YAxis(1).Color = colors_Manuscript2020('rich black');
ax2.YAxis(2).Color = colors_Manuscript2020('deep carrot orange');
%% CBV and behavioral indeces
ax34 =subplot(7,1,[3,4]);
p6 = plot((1:length(filtRH_HbT))/ProcData.notes.CBVCamSamplingRate,filtRH_HbT,'color',colors_Manuscript2020('sapphire'),'LineWidth',1);
hold on
p5 = plot((1:length(filtLH_HbT))/ProcData.notes.CBVCamSamplingRate,filtLH_HbT,'color',colors_Manuscript2020('dark candy apple red'),'LineWidth',1);
x1 = xline(130,'color',colorB,'LineWidth',2);
x2 = xline(360,'color',colorA,'LineWidth',2);
x3 = xline(465,'color','k','LineWidth',2);
ylabel('\DeltaHbT')
legend([p5,p6,x2,x1,x3],'Left hem','Right hem','Awake','NREM','Isoflurane')
set(gca,'TickLength',[0,0])
set(gca,'Xticklabel',[])
set(gca,'box','off')
xticks([130,190,250,310,370,430,490,550,610,670,730])
axis tight
xlim([130,730]) 
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
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
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
xticks([130,190,250,310,370,430,490,550,610,670,730])
xlim([130,730]) 
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
xticks([130,190,250,310,370,430,490,550,610,670,730])
xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
xlim([130,730]) 
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
%% Mean HbT comparison between behaviors
% pre-allocate the date for each day
IOS_behavFields = {'Rest','Whisk','Stim','NREM','REM'};
for aa = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,aa};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.FileIDs);
    whiskFileDates = [];
    % identify the unique days present for each animal using the whisking field.
    for bb = 1:length(whiskFileIDs)
        whiskFileDates{bb,1} = ConvertDate_IOS_Manuscript2020(whiskFileIDs{bb,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % put pre-allocate each date
    for dd = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,dd};
        for ee = 1:length(uniqueWhiskFileDates)
            fileDate = uniqueWhiskFileDates{ee,1};
            data.HbT.(animalID).(behavField).(fileDate).MeanLH = [];
            data.HbT.(animalID).(behavField).(fileDate).MeanRH = [];
            data.HbT.(animalID).(behavField).(fileDate).IndLH = {};
            data.HbT.(animalID).(behavField).(fileDate).IndRH = {};
        end
        procData.HbT.(behavField).animalID{aa,1} = animalID;
        procData.HbT.(behavField).behavior{aa,1} = behavField;
        procData.HbT.(behavField).LH{aa,1} = 'LH';
        procData.HbT.(behavField).RH{aa,1} = 'RH';
    end
end
% put data into cell for each unique date
for ff = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,ff};
    for gg = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,gg};
        % data is structured slightly differently depending on class
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'Whisk') == true
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS_Manuscript2020(fileIDs{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjLH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjRH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).IndLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).IndRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{hh,1});
            end
        elseif strcmp(behavField,'Stim') == true
            % Left hem stims
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.LH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS_Manuscript2020(fileIDs{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjLH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).IndLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjLH{hh,1});
            end
            % Right hem stims
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.RH_FileIDs;
            for hh = 1:length(fileIDs)
                fileDate = ConvertDate_IOS_Manuscript2020(fileIDs{hh,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.MeanAdjRH(hh,1));
                data.HbT.(animalID).(behavField).(fileDate).IndRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanCBV.(behavField).CBV_HbT.IndAdjRH{hh,1});
            end
        else
            fileIDs = AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.FileIDs;
            for ii = 1:length(fileIDs)
                fileDate = ConvertDate_IOS_Manuscript2020(fileIDs{ii,1});
                data.HbT.(animalID).(behavField).(fileDate).MeanLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanLH,AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.MeanAdjLH(ii,1));
                data.HbT.(animalID).(behavField).(fileDate).MeanRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).MeanRH,AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.MeanAdjRH(ii,1));
                data.HbT.(animalID).(behavField).(fileDate).IndLH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndLH,AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.IndAdjLH{ii,1});
                data.HbT.(animalID).(behavField).(fileDate).IndRH = cat(1,data.HbT.(animalID).(behavField).(fileDate).IndRH,AnalysisResults.(animalID).MeanCBV.(behavField).(modelType).CBV_HbT.IndAdjRH{ii,1});
            end
        end
    end
end
% find the mean of the 10-second resting periods from each day to determine a baseline
for jj = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,jj};
    whiskFileIDs = unique(AnalysisResults.(animalID).MeanCBV.Whisk.CBV_HbT.FileIDs);
    whiskFileDates = [];
    for kk = 1:length(whiskFileIDs)
        whiskFileDates{kk,1} = ConvertDate_IOS_Manuscript2020(whiskFileIDs{kk,1}); %#ok<*AGROW>
    end
    uniqueWhiskFileDates = unique(whiskFileDates);
    % take mean from each day. Days with no data will show up as NaN and be excluded
    for ll = 1:length(uniqueWhiskFileDates)
        fileDate = uniqueWhiskFileDates{ll,1};
        data.HbT.(animalID).Rest.(fileDate).baselineLH = mean(data.HbT.(animalID).Rest.(fileDate).MeanLH);
        data.HbT.(animalID).Rest.(fileDate).baselineRH = mean(data.HbT.(animalID).Rest.(fileDate).MeanRH);
    end
end
% Subtract the 10-second resting baseline for each day from the other data types. If the day doesn't have resting data,
% exclude it from analysis
for mm = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,mm};
    for nn = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,nn};
        % Subtract each day's 10-second baseline from each behavior field
        fileDates = fieldnames(data.HbT.(animalID).(behavField));
        for oo = 1:length(fileDates)
            fileDate = fileDates{oo,1};
            if strcmp(behavField,'Stim') == true || strcmp(behavField,'Whisk') == true
                data.HbT.(animalID).(behavField).(fileDate).CorrMeanLH = data.HbT.(animalID).(behavField).(fileDate).MeanLH;% - data.HbT.(animalID).Rest.(fileDate).baselineLH;
                data.HbT.(animalID).(behavField).(fileDate).CorrMeanRH = data.HbT.(animalID).(behavField).(fileDate).MeanRH;%; - data.HbT.(animalID).Rest.(fileDate).baselineRH;
                for pp = 1:length(data.HbT.(animalID).(behavField).(fileDate).IndLH)
                    data.HbT.(animalID).(behavField).(fileDate).CorrIndLH{pp,1} = data.HbT.(animalID).(behavField).(fileDate).IndLH{pp,1};% - data.HbT.(animalID).Rest.(fileDate).baselineLH;
                end
                for pp = 1:length(data.HbT.(animalID).(behavField).(fileDate).IndRH)
                    data.HbT.(animalID).(behavField).(fileDate).CorrIndRH{pp,1} = data.HbT.(animalID).(behavField).(fileDate).IndRH{pp,1};% - data.HbT.(animalID).Rest.(fileDate).baselineRH;
                end
            else
                data.HbT.(animalID).(behavField).(fileDate).CorrMeanLH = data.HbT.(animalID).(behavField).(fileDate).MeanLH - data.HbT.(animalID).Rest.(fileDate).baselineLH;
                data.HbT.(animalID).(behavField).(fileDate).CorrMeanRH = data.HbT.(animalID).(behavField).(fileDate).MeanRH - data.HbT.(animalID).Rest.(fileDate).baselineRH;
                for pp = 1:length(data.HbT.(animalID).(behavField).(fileDate).IndLH)
                    data.HbT.(animalID).(behavField).(fileDate).CorrIndLH{pp,1} = data.HbT.(animalID).(behavField).(fileDate).IndLH{pp,1} - data.HbT.(animalID).Rest.(fileDate).baselineLH;
                end
                for pp = 1:length(data.HbT.(animalID).(behavField).(fileDate).IndRH)
                    data.HbT.(animalID).(behavField).(fileDate).CorrIndRH{pp,1} = data.HbT.(animalID).(behavField).(fileDate).IndRH{pp,1} - data.HbT.(animalID).Rest.(fileDate).baselineRH;
                end
            end
        end
    end
end
% Take the mean of the corrected data from each unique day
for qq = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,qq};
    for rr = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,rr};
        fileDates = fieldnames(data.HbT.(animalID).(behavField));
        for ss = 1:length(fileDates)
            fileDate = fileDates{ss,1};
            data.HbT.(animalID).(behavField).(fileDate).DayMeanLH = mean(data.HbT.(animalID).(behavField).(fileDate).CorrMeanLH);
            data.HbT.(animalID).(behavField).(fileDate).DayMeanRH = mean(data.HbT.(animalID).(behavField).(fileDate).CorrMeanRH);
            data.HbT.(animalID).(behavField).(fileDate).DayAllMeanLH = [];
            data.HbT.(animalID).(behavField).(fileDate).DayAllMeanRH = [];
            data.HbT.(animalID).(behavField).(fileDate).DayIndLH = [];
            data.HbT.(animalID).(behavField).(fileDate).DayIndRH = [];
            % concatenate individual trials into a single array for each unique day
            if isfield(data.HbT.(animalID).(behavField).(fileDate),'CorrIndLH') == true
                % LH means - diff loop is necessary as STIM field has diff number of events
                for tt = 1:length(data.HbT.(animalID).(behavField).(fileDate).CorrMeanLH)
                    data.HbT.(animalID).(behavField).(fileDate).DayAllMeanLH = cat(2,data.HbT.(animalID).(behavField).(fileDate).DayIndLH,data.HbT.(animalID).(behavField).(fileDate).CorrIndLH{tt,1});
                end
                % RH means
                for tt = 1:length(data.HbT.(animalID).(behavField).(fileDate).CorrMeanRH)
                    data.HbT.(animalID).(behavField).(fileDate).DayAllMeanRH = cat(2,data.HbT.(animalID).(behavField).(fileDate).DayIndRH,data.HbT.(animalID).(behavField).(fileDate).CorrIndRH{tt,1});
                end
                % LH individual data pts
                for tt = 1:length(data.HbT.(animalID).(behavField).(fileDate).CorrIndLH)
                    data.HbT.(animalID).(behavField).(fileDate).DayIndLH = cat(2,data.HbT.(animalID).(behavField).(fileDate).DayIndLH,data.HbT.(animalID).(behavField).(fileDate).CorrIndLH{tt,1});
                end
                % RH individual data pts
                for tt = 1:length(data.HbT.(animalID).(behavField).(fileDate).CorrIndRH)
                    data.HbT.(animalID).(behavField).(fileDate).DayIndRH = cat(2,data.HbT.(animalID).(behavField).(fileDate).DayIndRH,data.HbT.(animalID).(behavField).(fileDate).CorrIndRH{tt,1});
                end
            end
        end
    end
end
% Put all the corrected means from each unique day into a single vector
nans = 1;
for uu = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,uu};
    for vv = 1:length(IOS_behavFields)
        behavField = IOS_behavFields{1,vv};
        fileDates = fieldnames(data.HbT.(animalID).(behavField));
        procData.HbT.(animalID).(behavField).DayMeansLH = [];
        procData.HbT.(animalID).(behavField).DayMeansRH = [];
        procData.HbT.(animalID).(behavField).CatIndLH = [];
        procData.HbT.(animalID).(behavField).CatIndRH = [];
        for ww = 1:length(fileDates)
            fileDate = fileDates{ww,1};
            if isnan(data.HbT.(animalID).(behavField).(fileDate).DayMeanLH) == false
                procData.HbT.(animalID).(behavField).DayMeansLH = cat(1,procData.HbT.(animalID).(behavField).DayMeansLH,data.HbT.(animalID).(behavField).(fileDate).DayMeanLH);
                procData.HbT.(animalID).(behavField).DayMeansRH = cat(1,procData.HbT.(animalID).(behavField).DayMeansRH,data.HbT.(animalID).(behavField).(fileDate).DayMeanRH);
                procData.HbT.(animalID).(behavField).CatIndLH = cat(2,procData.HbT.(animalID).(behavField).CatIndLH,data.HbT.(animalID).(behavField).(fileDate).DayIndLH);
                procData.HbT.(animalID).(behavField).CatIndRH = cat(2,procData.HbT.(animalID).(behavField).CatIndRH,data.HbT.(animalID).(behavField).(fileDate).DayIndRH);
            else
                nans = nans + 1;
            end
        end
    end
end
% Put all the means (of the corrected means) from each unique day into a single vector
for yy = 1:length(IOS_behavFields)
    behavField = IOS_behavFields{1,yy};
    procData.HbT.(behavField).IndMeanCBV = [];
    procData.HbT.(behavField).CatCBV = [];
    procData.HbT.(behavField).meanLH = [];
    procData.HbT.(behavField).meanRH = [];
    for zz = 1:length(IOS_animalIDs)
        animalID = IOS_animalIDs{1,zz};
        procData.HbT.(behavField).IndMeanCBV = cat(1,procData.HbT.(behavField).IndMeanCBV,mean(procData.HbT.(animalID).(behavField).DayMeansLH),mean(procData.HbT.(animalID).(behavField).DayMeansRH));
        procData.HbT.(behavField).meanLH = cat(1,procData.HbT.(behavField).meanLH,mean(procData.HbT.(animalID).(behavField).DayMeansLH));
        procData.HbT.(behavField).meanRH = cat(1,procData.HbT.(behavField).meanRH,mean(procData.HbT.(animalID).(behavField).DayMeansRH));
        procData.HbT.(behavField).CatCBV = cat(2,procData.HbT.(behavField).CatCBV,procData.HbT.(animalID).(behavField).CatIndLH,procData.HbT.(animalID).(behavField).CatIndRH);
    end
end
% Take the mean and stdev across animals
for aaa = 1:length(IOS_behavFields)
    behavField = IOS_behavFields{1,aaa};
    procData.HbT.(behavField).MeanCBV = mean(procData.HbT.(behavField).IndMeanCBV,1);
    procData.HbT.(behavField).StdMeanCBV = std(procData.HbT.(behavField).IndMeanCBV,0,1);
end
% % pull out isoflurane data
for a = 1:length(Iso_AnimalIDs)
    isoAnimalID = Iso_AnimalIDs{1,a};
    data.Iso.CBV_HbT.meanLH(a,1) = AnalysisResults.(isoAnimalID).MeanCBV.Iso.CBV_HbT.adjLH;
    data.Iso.CBV_HbT.meanRH(a,1) = AnalysisResults.(isoAnimalID).MeanCBV.Iso.CBV_HbT.adjRH;
    data.Iso.CBV_HbT.animalID{a,1} = isoAnimalID;
    data.Iso.CBV_HbT.behavior{a,1} = 'Iso';
    data.Iso.CBV_HbT.LH{a,1} = 'LH';
    data.Iso.CBV_HbT.RH{a,1} = 'RH';
end
% take the mean and standard deviation of each set of signals
behavFields = {'Iso'};
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
%% statistics - linear mixed effects model
HbTtableSize = cat(1,procData.HbT.Rest.meanLH,procData.HbT.Rest.meanRH,procData.HbT.Whisk.meanLH,procData.HbT.Whisk.meanRH,...
    procData.HbT.Stim.meanLH,procData.HbT.Stim.meanRH,procData.HbT.NREM.meanLH,procData.HbT.NREM.meanRH,procData.HbT.REM.meanLH,procData.HbT.REM.meanRH,data.Iso.CBV_HbT.meanLH,data.Iso.CBV_HbT.meanRH);
HbTTable = table('Size',[size(HbTtableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','HbT','Behavior','Hemisphere'});
HbTTable.Mouse = cat(1,procData.HbT.Rest.animalID,procData.HbT.Rest.animalID,procData.HbT.Whisk.animalID,procData.HbT.Whisk.animalID,...
    procData.HbT.Stim.animalID,procData.HbT.Stim.animalID,procData.HbT.NREM.animalID,procData.HbT.NREM.animalID,procData.HbT.REM.animalID,procData.HbT.REM.animalID,data.Iso.CBV_HbT.animalID,data.Iso.CBV_HbT.animalID);
HbTTable.HbT = cat(1,procData.HbT.Rest.meanLH,procData.HbT.Rest.meanRH,procData.HbT.Whisk.meanLH,procData.HbT.Whisk.meanRH,...
    procData.HbT.Stim.meanLH,procData.HbT.Stim.meanRH,procData.HbT.NREM.meanLH,procData.HbT.NREM.meanRH,procData.HbT.REM.meanLH,procData.HbT.REM.meanRH,data.Iso.CBV_HbT.meanLH,data.Iso.CBV_HbT.meanRH);
HbTTable.Behavior = cat(1,procData.HbT.Rest.behavior,procData.HbT.Rest.behavior,procData.HbT.Whisk.behavior,procData.HbT.Whisk.behavior,...
    procData.HbT.Stim.behavior,procData.HbT.Stim.behavior,procData.HbT.NREM.behavior,procData.HbT.NREM.behavior,procData.HbT.REM.behavior,procData.HbT.REM.behavior,data.Iso.CBV_HbT.behavior,data.Iso.CBV_HbT.behavior);
HbTTable.Hemisphere = cat(1,procData.HbT.Rest.LH,procData.HbT.Rest.RH,procData.HbT.Whisk.LH,procData.HbT.Whisk.RH,...
    procData.HbT.Stim.LH,procData.HbT.Stim.RH,procData.HbT.NREM.LH,procData.HbT.NREM.RH,procData.HbT.REM.LH,procData.HbT.REM.RH,data.Iso.CBV_HbT.LH,data.Iso.CBV_HbT.RH);
HbTFitFormula = 'HbT ~ 1 + Behavior + (1|Mouse) + (1|Mouse:Hemisphere)';
HbTStats = fitglme(HbTTable,HbTFitFormula);
%% Pixel panel S16
summaryFigure = figure('Name','FigS16 (g)');
sgtitle('Figure panel S16 (g) Turner Manuscript 2020')
%% [S16g] Mean HbT during different behaviors
HbT_xInds = ones(1,length(IOS_animalIDs)*2);
Iso_xInds = ones(1,length(Iso_AnimalIDs)*2);
s1 = scatter(HbT_xInds*1,procData.HbT.Rest.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,procData.HbT.Rest.MeanCBV,procData.HbT.Rest.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(HbT_xInds*2,procData.HbT.Whisk.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,procData.HbT.Whisk.MeanCBV,procData.HbT.Whisk.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(HbT_xInds*3,procData.HbT.Stim.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorE,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,procData.HbT.Stim.MeanCBV,procData.HbT.Stim.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(HbT_xInds*4,procData.HbT.NREM.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,procData.HbT.NREM.MeanCBV,procData.HbT.NREM.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(HbT_xInds*5,procData.HbT.REM.IndMeanCBV,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e5 = errorbar(5,procData.HbT.REM.MeanCBV,procData.HbT.REM.StdMeanCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
s6 = scatter(Iso_xInds*6,data.Iso.CBV_HbT.Comb,75,'MarkerEdgeColor','k','MarkerFaceColor',colorF,'jitter','on','jitterAmount',0.25);
e6 = errorbar(6,data.Iso.CBV_HbT.meanCBV,data.Iso.CBV_HbT.stdCBV,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
ylabel('\DeltaHbT (\muM)')
legend([s1,s2,s3,s4,s5,s6],'Awake Rest','Whisk','Stim','NREM','REM','Iso','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(IOS_behavFields) + 2])
ylim([-10,250])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
set(summaryFigure,'PaperPositionMode','auto');
savefig(summaryFigure,[dirpath 'FigS16_B']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'FigS16_B'])
%% statistical diary
diaryFile = [dirpath 'FigS16_Statistics.txt'];
if exist(diaryFile,'file') == true
    delete(diaryFile)
end
diary(diaryFile)
diary on
% HbT statistical diary
disp('======================================================================================================================')
disp('[S16a] Generalized linear mixed-effects model statistics for mean HbT during Rest, Whisk, Stim, NREM, REM, Isoflurane')
disp('======================================================================================================================')
disp(HbTStats)
disp('----------------------------------------------------------------------------------------------------------------------')
diary off

end
