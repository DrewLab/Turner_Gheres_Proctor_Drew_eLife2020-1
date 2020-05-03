%% Power spectra of different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Rest.VesselPowerSpec.S = []; data.Whisk.VesselPowerSpec.S = []; data.NREM.VesselPowerSpec.S = []; data.REM.VesselPowerSpec.S = []; data.CombRestWhisk.VesselPowerSpec.S = []; data.Awake.VesselPowerSpec.S = [];
data.Rest.VesselPowerSpec.f = []; data.Whisk.VesselPowerSpec.f = []; data.NREM.VesselPowerSpec.f = []; data.REM.VesselPowerSpec.f = []; data.CombRestWhisk.VesselPowerSpec.f = []; data.Awake.VesselPowerSpec.f = [];
data.VesselPowerSpec.whiskingPerc = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    powerSpecBehavFields = fieldnames(AnalysisResults.(animalID).VesselPowerSpectra);
    for bb = 1:length(powerSpecBehavFields)
        behavField = powerSpecBehavFields{bb,1};
        vesselIDs = fieldnames(AnalysisResults.(animalID).VesselPowerSpectra.(behavField));
        for cc = 1:length(vesselIDs)
            vesselID = vesselIDs{cc,1};
            data.(behavField).VesselPowerSpec.S = horzcat(data.(behavField).VesselPowerSpec.S,AnalysisResults.(animalID).VesselPowerSpectra.(behavField).(vesselID).S);
            data.(behavField).VesselPowerSpec.f = vertcat(data.(behavField).VesselPowerSpec.f,AnalysisResults.(animalID).VesselPowerSpectra.(behavField).(vesselID).f);
            if strcmp(behavField,'Awake') == true
                data.VesselPowerSpec.whiskingPerc = horzcat(data.VesselPowerSpec.whiskingPerc,AnalysisResults.(animalID).VesselPowerSpectra.(behavField).(vesselID).whiskingPerc);
            end
        end
    end
end
% take the average power of the vessels for each behavior before normalizing
powerSpecBehavFields = {'Whisk','Rest','CombRestWhisk','Awake','NREM','REM'};
for dd = 1:length(powerSpecBehavFields)
    behavField = powerSpecBehavFields{1,dd};
    data.(behavField).VesselPowerSpec.preMeanS = mean(data.(behavField).VesselPowerSpec.S,2);
end
% normalize the vessel power by the peak average power during the restng trace
restScaleFactor = 1/max(data.Rest.VesselPowerSpec.preMeanS);
for dd = 1:length(powerSpecBehavFields)
    behavField = powerSpecBehavFields{1,dd};
    for ee = 1:size(data.(behavField).VesselPowerSpec.S,2)
        data.(behavField).VesselPowerSpec.normS(:,ee) = (data.(behavField).VesselPowerSpec.S(:,ee))*restScaleFactor;% - restBaseline)./restBaseline;
    end
end
% take the average power of the vessels for each behavior
for dd = 1:length(powerSpecBehavFields)
    behavField = powerSpecBehavFields{1,dd};
    data.(behavField).VesselPowerSpec.meanS = mean(data.(behavField).VesselPowerSpec.normS,2);
    data.(behavField).VesselPowerSpec.StDS = std(data.(behavField).VesselPowerSpec.normS,0,2);
    data.(behavField).VesselPowerSpec.meanf = mean(data.(behavField).VesselPowerSpec.f,1);
end