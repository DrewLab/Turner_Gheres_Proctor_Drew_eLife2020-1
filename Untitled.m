animalIDs = {'T108','T122','T123'};
rootDir = cd;
for aa = 1:length(animalIDs)
    clearvars -except animalIDs rootDir aa
    animalID = animalIDs{1,aa};
    baseDir = [rootDir '\' animalID '\Bilateral Imaging\'];
    cd(baseDir)
    StageThreeProcessing_IOS_Manuscript2020
    cd(rootDir)
end

