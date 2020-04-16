animalIDs = {'T110','T111'};
rootFolder = cd;
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    animalDir = [rootFolder '\' animalID '\Single Hemisphere\'];
    cd(animalDir)
    StageTwoProcessing2_IOS_Manuscript2020(rootFolder)
end