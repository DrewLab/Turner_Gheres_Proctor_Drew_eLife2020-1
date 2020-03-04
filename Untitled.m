
mergedDirectory = dir('*_MergedData.mat');
mergedDataFiles = {mergedDirectory.name}';
mergedDataFileIDs = char(mergedDataFiles);
for a = 1:17
    disp(['Combining the data from LabVIEW and MScan file(s) number ' num2str(a) ' of ' num2str(size(mergedDataFileIDs, 1)) '...']); disp(' ');
    mergedDataFileID = mergedDataFileIDs(a,:);
    load(mergedDataFileID);

    M2.notes = MergedData.notes;
    
    % Pull the notes and data from LabVIEW
    M2.data.whiskerAngle = MergedData.data.whiskerAngle;
    M2.data.rawWhiskerAngle = MergedData.data.rawWhiskerAngle;
    M2.data.binWhiskerAngle = MergedData.data.binWhiskerAngle;
    M2.data.forceSensorL = MergedData.data.forceSensorL;
    M2.data.binForceSensorL = MergedData.data.binForceSensorL;
    
    % Save solenoid times (in seconds). Identify the solenoids by amplitude.
    M2.data.solenoids = MergedData.data.solenoids.LPadSol;
    M2.data.solenoids.LPadSol = MergedData.data.solenoids.RPadSol;
    M2.data.solenoids.LPadSol = MergedData.data.solenoids.AudSol;
    
    % Pull the notes and data from MScan
    M2.data.rawCorticalNeural = MergedData.data.corticalNeural;
    M2.data.rawHippocampalNeural = MergedData.data.hippocampalNeural;
    M2.data.corticalNeural.muaPower = MergedData.data.corticalMUAPower;
    M2.data.corticalNeural.gammaBandPower = MergedData.data.corticalGammaBandPower;
    M2.data.corticalNeural.betaBandPower = MergedData.data.corticalBetaBandPower;
    M2.data.corticalNeural.alphaBandPower  = MergedData.data.corticalAlphaBandPower;
    M2.data.corticalNeural.thetaBandPower = MergedData.data.corticalThetaBandPower;
    M2.data.corticalNeural.deltaBandPower = MergedData.data.corticalDeltaBandPower;
    M2.data.hippocampalNeural.muaPower = MergedData.data.hippocampalMUAPower;
    M2.data.hippocampalNeural.gammaBandPower = MergedData.data.hippocampalGammaBandPower;
    M2.data.hippocampalNeural.betaBandPower = MergedData.data.hippocampalBetaBandPower;
    M2.data.hippocampalNeural.alphaBandPower = MergedData.data.hippocampalAlphaBandPower;
    M2.data.hippocampalNeural.thetaBandPower = MergedData.data.hippocampalThetaBandPower;
    M2.data.hippocampalNeural.deltaBandPower  = MergedData.data.hippocampalDeltaBandPower;
    M2.data.forceSensorM = MergedData.data.forceSensorM;
    M2.data.EMG.data = MergedData.data.EMG;
    M2.data.binForceSensorM = MergedData.data.binForceSensorM;
    M2.data.vesselDiameter.data = MergedData.data.vesselDiameter;
%     M2.flags = MergedData.flags;
    MergedData = [];
    MergedData = M2;
    save(mergedDataFileID,'MergedData')
end
