animalIDs = {'T108','T109','T110','T111','T119','T120','T121','T122','T123'};
startingDir = cd;
for aa = 1:size(animalIDs,2)
    clearvars -except animalIDs startingDir aa
    animalID = animalIDs{1,aa};
    animalDir2 = [startingDir '\' animalID '\Single Hemisphere\'];
    cd(animalDir2)
    StageThreeProcessing_IOS_Manuscript2020
    cd(startingDir)
end
