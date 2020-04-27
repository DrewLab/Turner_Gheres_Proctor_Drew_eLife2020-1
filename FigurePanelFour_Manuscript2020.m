function [] = FigurePanelFour_Manuscript2020(rootFolder,AnalysisResults)
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
%% REM dilation and trnasition REM to Awake
% cd through each animal's directory and extract the appropriate analysis results
data.Whisk.Evoked.data = []; data.REM.Evoked.data = []; data.REMtoAwake.Evoked.data = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    evokedBehavFields = fieldnames(AnalysisResults.(animalID).EvokedAvgs);
    for bb = 1:length(evokedBehavFields)
        behavField = evokedBehavFields{bb,1};
        vesselIDs = fieldnames(AnalysisResults.(animalID).EvokedAvgs.(behavField));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            data.(behavField).Evoked.data = vertcat(data.(behavField).Evoked.data,AnalysisResults.(animalID).EvokedAvgs.(behavField).(vesselID).mean);
            data.(behavField).Evoked.timeVector = AnalysisResults.(animalID).EvokedAvgs.(behavField).(vesselID).timeVector;
        end
    end
end
% take the average of the vessels for each behavior
evokedBehavFields = {'Whisk','REM','REMtoAwake'};
for dd = 1:length(evokedBehavFields)
    behavField = evokedBehavFields{1,dd};
    data.(behavField).Evoked.mean = mean(data.(behavField).Evoked.data,1);
    data.(behavField).Evoked.StD = std(data.(behavField).Evoked.data,0,1);
end
%% Power spectra of different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Rest.PowerSpec.S = []; data.Whisk.PowerSpec.S = []; data.NREM.PowerSpec.S = []; data.REM.PowerSpec.S = []; data.CombRestWhisk.PowerSpec.S = []; data.AllData.PowerSpec.S = [];
data.Rest.PowerSpec.f = []; data.Whisk.PowerSpec.f = []; data.NREM.PowerSpec.f = []; data.REM.PowerSpec.f = []; data.CombRestWhisk.PowerSpec.f = []; data.AllData.PowerSpec.f = [];
data.PowerSpec.whiskingPerc = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    powerSpecBehavFields = fieldnames(AnalysisResults.(animalID).PowerSpectra);
    for bb = 1:length(powerSpecBehavFields)
        behavField = powerSpecBehavFields{bb,1};
        vesselIDs = fieldnames(AnalysisResults.(animalID).PowerSpectra.(behavField));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            data.(behavField).PowerSpec.S = horzcat(data.(behavField).PowerSpec.S,AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).S);
            data.(behavField).PowerSpec.f = vertcat(data.(behavField).PowerSpec.f,AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).f);
            if strcmp(behavField,'AllData') == true
                data.PowerSpec.whiskingPerc = horzcat(data.PowerSpec.whiskingPerc,AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).whiskingPerc);
            end
        end
    end
end
% take the average power of the vessels for each behavior before normalizing
powerSpecBehavFields = {'Whisk','Rest','CombRestWhisk','AllData','NREM','REM'};
for dd = 1:length(powerSpecBehavFields)
    behavField = powerSpecBehavFields{1,dd};
    data.(behavField).PowerSpec.preMeanS = mean(data.(behavField).PowerSpec.S,2);
end
% normalize the vessel power by the peak average power during the restng trace
restScaleFactor = 1/max(data.Rest.PowerSpec.preMeanS);
for dd = 1:length(powerSpecBehavFields)
    behavField = powerSpecBehavFields{1,dd};
    for ee = 1:size(data.(behavField).PowerSpec.S,2)
        data.(behavField).PowerSpec.normS(:,ee) = (data.(behavField).PowerSpec.S(:,ee))*restScaleFactor;% - restBaseline)./restBaseline;
    end
end
% take the average power of the vessels for each behavior
for dd = 1:length(powerSpecBehavFields)
    behavField = powerSpecBehavFields{1,dd};
    data.(behavField).PowerSpec.meanS = mean(data.(behavField).PowerSpec.normS,2);
    data.(behavField).PowerSpec.StDS = std(data.(behavField).PowerSpec.normS,0,2);
    data.(behavField).PowerSpec.meanf = mean(data.(behavField).PowerSpec.f,1);
end
%% Figure panel four
summaryFigure = figure;
sgtitle('Turner Manuscript 2020 - Figure panel 4')
%% [A] Slow arteriole dilation during REM
ax1 = subplot(2,4,[1,2]);
plot(data.REM.Evoked.timeVector,data.REM.Evoked.mean,'color',colors_Manuscript2020('rich black'),'LineWidth',2)
hold on;
plot(data.REM.Evoked.timeVector,data.REM.Evoked.mean + data.REM.Evoked.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.REM.Evoked.timeVector,data.REM.Evoked.mean - data.REM.Evoked.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[A] Slow arteriole dilation during REM')
xlabel('Time after REM onset (s)')
ylabel('\DeltaD/D (%)')
axis tight
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [B] REM to Awake transition
ax2 = subplot(2,4,[3,4]);
plot(data.REMtoAwake.Evoked.timeVector,data.REMtoAwake.Evoked.mean,'color',colors_Manuscript2020('rich black'),'LineWidth',2)
hold on;
plot(data.REMtoAwake.Evoked.timeVector,data.REMtoAwake.Evoked.mean + data.REMtoAwake.Evoked.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.REMtoAwake.Evoked.timeVector,data.REMtoAwake.Evoked.mean - data.REMtoAwake.Evoked.StD,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[B] REM to Awake behavioral transition')
xlabel('Periwaking time (s)')
ylabel('\DeltaD/D (%)')
axis tight
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [C] Arteriole power spectra during awake rest
ax3 = subplot(2,4,5);
for aa = 1:size(data.Rest.PowerSpec.S,2)
    loglog(data.Rest.PowerSpec.f(aa,:),data.Rest.PowerSpec.normS(:,aa),'color',colors_Manuscript2020('battleship grey'),'LineWidth',1);
    hold on
end
loglog(data.Rest.PowerSpec.meanf,data.Rest.PowerSpec.meanS,'color',colorA,'LineWidth',3);
title({'[C] Rest arteriole','power spectra',''})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [D] Arteriole power spectra during whisking
ax4 = subplot(2,4,6);
for aa = 1:size(data.Whisk.PowerSpec.S,2)
    loglog(data.Whisk.PowerSpec.f(aa,:),data.Whisk.PowerSpec.normS(:,aa),'color',colors_Manuscript2020('battleship grey'),'LineWidth',1);
    hold on
end
loglog(data.Whisk.PowerSpec.meanf,data.Whisk.PowerSpec.meanS,'color',colorD,'LineWidth',3);
title({'[D] Whisking arteriole','power spectra',''})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.1,0.5])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [E] Arteriole power spectra during NREM
ax5 = subplot(2,4,7);
for aa = 1:size(data.NREM.PowerSpec.S,2)
    loglog(data.NREM.PowerSpec.f(aa,:),data.NREM.PowerSpec.normS(:,aa),'color',colors_Manuscript2020('battleship grey'),'LineWidth',1);
    hold on
end
loglog(data.NREM.PowerSpec.meanf,data.NREM.PowerSpec.meanS,'color',colorB,'LineWidth',3);
title({'[E] NREM arteriole','power spectra',''})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([1/30,0.5])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [F] Arterole power spectra during REM
ax6 = subplot(2,4,8);
for aa = 1:size(data.REM.PowerSpec.S,2)
    loglog(data.REM.PowerSpec.f(aa,:),data.REM.PowerSpec.normS(:,aa),'color',colors_Manuscript2020('battleship grey'),'LineWidth',1);
    hold on
end
loglog(data.REM.PowerSpec.meanf,data.REM.PowerSpec.meanS,'color',colorC,'LineWidth',3);
title({'[F] REM arteriole','power spectra',''})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([1/60,0.5])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
linkaxes([ax3,ax4,ax5,ax6],'y')
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Figure Panel 4']);
print('-painters','-dpdf','-bestfit',[dirpath 'Figure Panel 4'])

end
