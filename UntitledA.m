animalIDs = {'T99','T101','T102'};
rootFolder = cd;
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    animalDir = [rootFolder '\' animalID '\Bilateral Imaging\'];
    cd(animalDir)
    StageTwoProcessing2_IOS_Manuscript2020()
    cd(rootFolder)
end