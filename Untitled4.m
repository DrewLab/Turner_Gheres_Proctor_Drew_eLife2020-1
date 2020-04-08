function [] = Untitled4()
specDirectory = dir('*_SpecData.mat');
specDataFiles = {specDirectory.name}';
specDataFileIDs = char(specDataFiles);

for aa = 1:size(specDataFileIDs)
    disp([num2str(aa) '/' num2str(size(specDataFileIDs,1))]);
    specDataFileID = specDataFileIDs(aa,:);
    load(specDataFileID)
    [animalID,~,fileID] = GetFileInfo_IOS_Manuscript2020(specDataFileID);
    specDataFileIDA = [animalID '_' fileID '_SpecDataA.mat'];
%     specDataFileIDB = [animalID '_' fileID '_SpecDataB.mat'];
    specDataFileIDC = [animalID '_' fileID '_SpecDataC.mat'];
    TempSpec = SpecData;
    SpecData = [];
    SpecData.cortical_LH = TempSpec.cortical_LH.fiveSec;
    SpecData.cortical_RH = TempSpec.cortical_RH.fiveSec;
    SpecData.hippocampus = TempSpec.hippocampus.fiveSec;
    save(specDataFileIDA,'SpecData');
%     SpecData = [];
%     SpecData.cortical_LH = TempSpec.cortical_LH.oneSecB;
%     SpecData.cortical_RH = TempSpec.cortical_RH.oneSecB;
%     SpecData.hippocampus = TempSpec.hippocampus.oneSecB;
%     save(specDataFileIDB,'SpecData');
    SpecData = [];
    SpecData.cortical_LH = TempSpec.cortical_LH.oneSec;
    SpecData.cortical_RH = TempSpec.cortical_RH.oneSec;
    SpecData.hippocampus = TempSpec.hippocampus.oneSec;
    save(specDataFileIDC,'SpecData');
end