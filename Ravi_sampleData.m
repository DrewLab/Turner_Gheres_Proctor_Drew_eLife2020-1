% T115
596:815   % Whisk events
455:485   % rest event one
% T125
420:520  % NREM
520:680 % REM
770:796 % rest event two
SampleData.animalIDs = {'T115','T115','T125','T125','T125'}';
SampleData.fileIDs = {'T115_RH_191122_16_39_22_011_P1_MergedData.mat','T115_RH_191122_16_39_22_011_P1_MergedData.mat','T125_RH_200308_11_39_57_006_P1_MergedData.mat','T125_RH_200308_11_39_57_006_P1_MergedData.mat','T125_RH_200308_11_39_57_006_P1_MergedData.mat'}';
SampleData.startTimes = [596,455,420,520,770]';
SampleData.endTimes = [815,485,520,680,796]';
SampleData.vFs = 5;
SampleData.wFs = 30;
SampleData.behavior = {'Whisk','Rest','NREM','REM','Rest'}';
% T115
vesselA = (MergedData.data.vesselDiameter.data(596*SampleData.vFs:815*SampleData.vFs) - RestingBaselines.manualSelection.vesselDiameter.data.P1.Nov22)./RestingBaselines.manualSelection.vesselDiameter.data.P1.Nov22;
whiskA = MergedData.data.whiskerAngle(596*SampleData.wFs:815*SampleData.wFs);
vesselB = (MergedData.data.vesselDiameter.data(455*SampleData.vFs:485*SampleData.vFs) - RestingBaselines.manualSelection.vesselDiameter.data.P1.Nov22)./RestingBaselines.manualSelection.vesselDiameter.data.P1.Nov22;
whiskB = MergedData.data.whiskerAngle(455*SampleData.wFs:485*SampleData.wFs);
SampleData.T115_BaselineDiameter = RestingBaselines.manualSelection.vesselDiameter.data.P1.Nov22;
% T125
vesselC = (MergedData.data.vesselDiameter.data(420*SampleData.vFs:520*SampleData.vFs) - RestingBaselines.manualSelection.vesselDiameter.data.P1.Mar08)./RestingBaselines.manualSelection.vesselDiameter.data.P1.Mar08;
whiskC = MergedData.data.whiskerAngle(420*SampleData.wFs:520*SampleData.wFs);
vesselD = (MergedData.data.vesselDiameter.data(520*SampleData.vFs:680*SampleData.vFs) - RestingBaselines.manualSelection.vesselDiameter.data.P1.Mar08)./RestingBaselines.manualSelection.vesselDiameter.data.P1.Mar08;
whiskD = MergedData.data.whiskerAngle(520*SampleData.wFs:680*SampleData.wFs);
vesselE = (MergedData.data.vesselDiameter.data(770*SampleData.vFs:796*SampleData.vFs) - RestingBaselines.manualSelection.vesselDiameter.data.P1.Mar08)./RestingBaselines.manualSelection.vesselDiameter.data.P1.Mar08;
whiskE = MergedData.data.whiskerAngle(770*SampleData.wFs:796*SampleData.wFs);
SampleData.vesselDiameters = {vesselA,vesselB,vesselC,vesselD,vesselE}';
SampleData.whiskerAngles = {whiskA,whiskB,whiskC,whiskD,whiskE}';
SampleData.T125_BaselineDiameter = RestingBaselines.manualSelection.vesselDiameter.data.P1.Mar08;
% figures
figure;
subplot(2,1,1)
plot(whiskA)
subplot(2,1,2)
plot(vesselA)

fileID = fopen('T125_RestEvent.txt','w');
fprintf(fileID,'%f\n',vesselE);
fclose(fileID);
