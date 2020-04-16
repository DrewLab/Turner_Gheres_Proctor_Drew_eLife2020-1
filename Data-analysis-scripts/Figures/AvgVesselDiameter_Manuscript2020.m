function [] = AvgVesselDiameter_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: Calculate the average vessel diameter during different behavioral states
%________________________________________________________________________________________________________________________

animalIDs = {'T115','T116','T117','T118','T125','T126'};
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(180/256)];  % whisk color

%% cd through each animal's directory and extract the appropriate analysis results
data.Whisk.data = []; data.Whisk.animalID = {}; data.Whisk.behavior = {}; data.Whisk.vID = {};
data.Rest.data = []; data.Rest.animalID = {}; data.Rest.behavior = {}; data.Rest.vID = {};
data.NREM.data = []; data.NREM.animalID = {}; data.NREM.behavior = {}; data.NREM.vID = {};
data.REM.data = []; data.REM.animalID = {}; data.REM.behavior = {}; data.REM.vID = {};
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    behavFields = fieldnames(AnalysisResults.(animalID).MeanVesselDiameter);
    for bb = 1:length(behavFields)
        behavField = behavFields{bb,1};
        vesselIDs = fieldnames(AnalysisResults.(animalID).MeanVesselDiameter.(behavField));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            if strcmp(vesselID(1),'V') == false
                data.(behavField).data = vertcat(data.(behavField).data,AnalysisResults.(animalID).MeanVesselDiameter.(behavField).(vesselID));
                data.(behavField).animalID = vertcat(data.(behavField).animalID,animalID);
                data.(behavField).behavior = vertcat(data.(behavField).behavior,behavField);
                data.(behavField).vID = vertcat(data.(behavField).vID,vesselID);
            end
        end
    end
end
% take the average of the vessels for each behavior
behavFields = {'Whisk','Rest','NREM','REM'};
for dd = 1:length(behavFields)
    behavField = behavFields{1,dd};
    data.(behavField).mean = mean(data.(behavField).data);
    data.(behavField).StD = std(data.(behavField).data,0,1);
end

%% statistics - linear mixed effects model
alphaConf = 0.001;
numComparisons = 3;
tableSize = cat(1,data.Rest.animalID,data.Whisk.animalID,data.NREM.animalID,data.REM.animalID);
vesselDiameterTable = table('Size',[size(tableSize,1),4],'VariableTypes',{'string','double','string','string'},'VariableNames',{'Mouse','Diameter','Behavior','Vessel'});
vesselDiameterTable.Mouse = cat(1,data.Rest.animalID,data.Whisk.animalID,data.NREM.animalID,data.REM.animalID);
vesselDiameterTable.Diameter = cat(1,data.Rest.data,data.Whisk.data,data.NREM.data,data.REM.data);
vesselDiameterTable.Behavior = cat(1,data.Rest.behavior,data.Whisk.behavior,data.NREM.behavior,data.REM.behavior);
vesselDiameterTable.Vessel = cat(1,data.Rest.vID,data.Whisk.vID,data.NREM.vID,data.REM.vID);
vesselFitFormula = 'Diameter ~ 1 + Behavior + (1|Mouse) + (1|Vessel)';
vesselStats = fitglme(vesselDiameterTable,vesselFitFormula);
vesselCI = coefCI(vesselStats,'Alpha',(alphaConf/numComparisons));

%% summary figure(s)
summaryFigure = figure;
xIndsWhisk = ones(1,length(data.Whisk.data));
xIndsRest = ones(1,length(data.Rest.data));
xIndsNREM = ones(1,length(data.NREM.data));
xIndsREM = ones(1,length(data.REM.data));
% scatter plot of mean vessel diameter per behavior
s1 = scatter(xIndsWhisk*1,data.Whisk.data,100,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.mean,data.Whisk.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 15;
e1.CapSize = 15;
s2 = scatter(xIndsRest*2,data.Rest.data,100,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Rest.mean,data.Rest.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 15;
e2.CapSize = 15;
s3 = scatter(xIndsNREM*3,data.NREM.data,100,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.NREM.mean,data.NREM.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 15;
e3.CapSize = 15;
s4 = scatter(xIndsREM*4,data.REM.data,100,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.REM.mean,data.REM.StD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 15;
e4.CapSize = 15;
title('Maximum vessel diameter (%)')
ylabel('\DeltaD/D (%)')
legend([s1,s2,s3,s4],'Whisking','Awake Rest','NREM','REM','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')
% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Vessel Diameter']);
% statistical diary
diary([dirpath 'Behavior_MaxVesselDiameter_Stats.txt'])
diary on
disp('Generalized linear mixed-effects model statistics for mean (max) vessel diameter during Rest, Whisking, NREM, and REM')
disp('======================================================================================================================')
disp(vesselStats)
disp('======================================================================================================================')
disp('Alpha = 0.001 confidence intervals with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(vesselCI(1,:))])
disp(['Whisk: ' num2str(vesselCI(2,:))])
disp(['NREM: ' num2str(vesselCI(3,:))])
disp(['REM: ' num2str(vesselCI(4,:))])
diary off

end

