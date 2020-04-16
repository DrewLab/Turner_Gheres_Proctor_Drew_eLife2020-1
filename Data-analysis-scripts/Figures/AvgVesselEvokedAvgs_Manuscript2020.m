function [] = AvgVesselEvokedAvgs_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T115','T116','T117','T118','T125','T126'};

%% cd through each animal's directory and extract the appropriate analysis results
data.Whisk.data = [];
data.REM.data = [];
data.REMtoAwake.data = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    behavFields = fieldnames(AnalysisResults.(animalID).EvokedAvgs);
    for bb = 1:length(behavFields)
        behavField = behavFields{bb,1};
        vesselIDs = fieldnames(AnalysisResults.(animalID).EvokedAvgs.(behavField));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            data.(behavField).data = vertcat(data.(behavField).data,AnalysisResults.(animalID).EvokedAvgs.(behavField).(vesselID).mean);
            data.(behavField).timeVector = AnalysisResults.(animalID).EvokedAvgs.(behavField).(vesselID).timeVector;
        end
    end
end
% take the average of the vessels for each behavior
behavFields = {'Whisk','REM','REMtoAwake'};
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    data.(behavField).mean = mean(data.(behavField).data,1);
    data.(behavField).StD = std(data.(behavField).data,0,1);
end

%% summary figure(s)
summaryFigure = figure;
sgtitle('Event-triggered vessel diameter')
subplot(1,3,1)
plot(data.Whisk.timeVector,data.Whisk.mean,'k','LineWidth',2)
hold on;
plot(data.Whisk.timeVector,data.Whisk.mean + data.Whisk.StD,'color',colors_Manuscript2020('battleship grey'))
plot(data.Whisk.timeVector,data.Whisk.mean - data.Whisk.StD,'color',colors_Manuscript2020('battleship grey'))
title('Whisking (1-5 sec)')
xlabel('Periwhisking time (s)')
ylabel('\DeltaD/D (%)')
axis square
axis tight
set(gca,'box','off')
subplot(1,3,2)
plot(data.REM.timeVector,data.REM.mean,'k','LineWidth',2)
hold on;
plot(data.REM.timeVector,data.REM.mean + data.REM.StD,'color',colors_Manuscript2020('battleship grey'))
plot(data.REM.timeVector,data.REM.mean - data.REM.StD,'color',colors_Manuscript2020('battleship grey'))
title('REM')
xlabel('Time from REM onset (s)')
ylabel('\DeltaD/D (%)')
axis square
axis tight
set(gca,'box','off')
subplot(1,3,3)
plot(data.REMtoAwake.timeVector,data.REMtoAwake.mean,'k','LineWidth',2)
hold on;
plot(data.REMtoAwake.timeVector,data.REMtoAwake.mean + data.REMtoAwake.StD,'color',colors_Manuscript2020('battleship grey'))
plot(data.REMtoAwake.timeVector,data.REMtoAwake.mean - data.REMtoAwake.StD,'color',colors_Manuscript2020('battleship grey'))
title('REM to Awake transition')
xlabel('Periwaking time (s)')
ylabel('\DeltaD/D (%)')
axis square
axis tight
set(gca,'box','off')
% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Vessel Evoked Responses']);

end

