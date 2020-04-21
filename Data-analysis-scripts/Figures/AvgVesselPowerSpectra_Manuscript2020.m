function [SimulationData] = AvgVesselPowerSpectra_Manuscript2020(rootFolder,AnalysisResults,SimulationData)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T115','T116','T117','T118','T125','T126'};
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(180/256)];  % whisk color

%% cd through each animal's directory and extract the appropriate analysis results
data.Whisk.S = [];
data.Rest.S = [];
data.NREM.S = [];
data.REM.S = [];
data.CombRestWhisk.S = [];
data.AllData.S = [];
data.Whisk.f = [];
data.Rest.f = [];
data.NREM.f = [];
data.REM.f = [];
data.CombRestWhisk.f = [];
data.AllData.f = [];
data.whiskingPerc = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    behavFields = fieldnames(AnalysisResults.(animalID).PowerSpectra);
    for bb = 1:length(behavFields)
        behavField = behavFields{bb,1};
        vesselIDs = fieldnames(AnalysisResults.(animalID).PowerSpectra.(behavField));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            data.(behavField).S = horzcat(data.(behavField).S,AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).S);
            data.(behavField).f = vertcat(data.(behavField).f,AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).f);
            if strcmp(behavField,'AllData') == true
                data.whiskingPerc = horzcat(data.whiskingPerc,AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).whiskingPerc);
            end
        end
    end
end
% take the average of the vessels for each behavior
behavFields = {'Whisk','Rest','CombRestWhisk','AllData','NREM','REM'};
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    data.(behavField).meanS = mean(data.(behavField).S,2);
    data.(behavField).StDS = std(data.(behavField).S,0,2);
    data.(behavField).meanf = mean(data.(behavField).f,1);
end
meanWhiskingPerc = round(mean(data.whiskingPerc),1);
disp(['Mean whisking percent for AllData: ' num2str(meanWhiskingPerc) '%']);

%% save data for simulations
behavFields = {'AllData','NREM','REM'};
for ee = 1:length(behavFields)
    behavField = behavFields{1,ee};
    SimulationData.VesselPowerSpectra.(behavField).indVesselS = data.(behavField).S;
    SimulationData.VesselPowerSpectra.(behavField).frequencyVector = data.(behavField).meanf;
    SimulationData.VesselPowerSpectra.(behavField).meanS = data.(behavField).meanS;
    SimulationData.VesselPowerSpectra.(behavField).StD = data.(behavField).StDS;
    if strcmp(behavField,'AllData') == true
        SimulationData.VesselPowerSpectra.(behavField).percentTimeWhisking = meanWhiskingPerc;
    end
end

%% summary figure(s)
summaryFigure = figure;
ax1 = subplot(1,2,1);
L1 = loglog(data.Whisk.meanf,data.Whisk.meanS,'color',colorD,'LineWidth',3);
hold on
L2 = loglog(data.Rest.meanf,data.Rest.meanS,'color',colorA,'LineWidth',3);
L3 = loglog(data.CombRestWhisk.meanf,data.CombRestWhisk.meanS,'color','k','LineWidth',3);
L4 = loglog(data.AllData.meanf,data.AllData.meanS,'color','r','LineWidth',3);
L5 = loglog(data.NREM.meanf,data.NREM.meanS,'color',colorB,'LineWidth',3);
L6 = loglog(data.REM.meanf,data.REM.meanS,'color',colorC,'LineWidth',3);
title('\DeltaD/D power spectrum')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([L1,L2,L3,L4,L5,L6],'Whisk','Rest','CombRestWhisk','All Awake','NREM','REM','Location','NorthEast')
axis square
xlim([0.1,0.5])
set(gca,'box','off')

ax2 = subplot(1,2,2);
L1 = loglog(data.AllData.meanf,data.AllData.meanS,'color','r','LineWidth',3);
hold on
L2 = loglog(data.NREM.meanf,data.NREM.meanS,'color',colorB,'LineWidth',3);
L3 = loglog(data.REM.meanf,data.REM.meanS,'color',colorC,'LineWidth',3);
title('\DeltaD/D power spectrum')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([L1,L2,L3],'All Awake','NREM','REM','Location','SouthWest')
axis square
xlim([0.033,0.5])
set(gca,'box','off')

linkaxes([ax1,ax2],'y')

% individual figures
indSummaryFigure = figure;
sgtitle('\DeltaD/D behavior-dependent power spectrum')
ax1 = subplot(2,3,1);
for aa = 1:size(data.Whisk.S,2)
    loglog(data.Whisk.f(aa,:),data.Whisk.S(:,aa),'color',colors_Manuscript2020('battleship grey'),'LineWidth',1);
    hold on
end
loglog(data.Whisk.meanf,data.Whisk.meanS,'color',colorD,'LineWidth',3);
title('Whisking')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
set(gca,'box','off')

ax2 = subplot(2,3,2);
for aa = 1:size(data.Rest.S,2)
    loglog(data.Rest.f(aa,:),data.Rest.S(:,aa),'color',colors_Manuscript2020('battleship grey'),'LineWidth',1);
    hold on
end
loglog(data.Rest.meanf,data.Rest.meanS,'color',colorA,'LineWidth',3);
title('Rest')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
set(gca,'box','off')

ax3 = subplot(2,3,3);
for aa = 1:size(data.CombRestWhisk.S,2)
    loglog(data.CombRestWhisk.f(aa,:),data.CombRestWhisk.S(:,aa),'color',colors_Manuscript2020('battleship grey'),'LineWidth',1);
    hold on
end
loglog(data.CombRestWhisk.meanf,data.CombRestWhisk.meanS,'color','k','LineWidth',3);
title('Combined rest and whisking')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
set(gca,'box','off')

ax4 = subplot(2,3,4);
for aa = 1:size(data.AllData.S,2)
    loglog(data.AllData.f(aa,:),data.AllData.S(:,aa),'color',colors_Manuscript2020('battleship grey'),'LineWidth',1);
    hold on
end
loglog(data.AllData.meanf,data.AllData.meanS,'color','r','LineWidth',3);
title('All awake data')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
set(gca,'box','off')

ax5 = subplot(2,3,5);
for aa = 1:size(data.NREM.S,2)
    loglog(data.NREM.f(aa,:),data.NREM.S(:,aa),'color',colors_Manuscript2020('battleship grey'),'LineWidth',1);
    hold on
end
loglog(data.NREM.meanf,data.NREM.meanS,'color',colorB,'LineWidth',3);
title('NREM')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
set(gca,'box','off')

ax6 = subplot(2,3,6);
for aa = 1:size(data.REM.S,2)
    loglog(data.REM.f(aa,:),data.REM.S(:,aa),'color',colors_Manuscript2020('battleship grey'),'LineWidth',1);
    hold on
end
loglog(data.REM.meanf,data.REM.meanS,'color',colorC,'LineWidth',3);
title('REM')
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
set(gca,'box','off')

linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'y')

% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Vessel Power']);
savefig(indSummaryFigure,[dirpath 'Summary Figure - Ind Vessel Power']);

end
