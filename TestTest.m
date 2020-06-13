function [] = TestTest(rootFolder)

colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
% extract data and mean of each event
nremCatLH_HbT = []; nremMeanLH_HbT = []; remCatLH_HbT = []; remMeanLH_HbT = [];
nremCatRH_HbT = []; nremMeanRH_HbT = []; remCatRH_HbT = []; remMeanRH_HbT = [];
nremCatLH_Gam = []; nremMeanLH_Gam = []; remCatLH_Gam = []; remMeanLH_Gam = [];
nremCatRH_Gam = []; nremMeanRH_Gam = []; remCatRH_Gam = []; remMeanRH_Gam = [];
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    dataLoc = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLoc)
    % add this animal's scoring labels with the other animal'ss
    sleepStruct = [animalID '_SleepData.mat'];
    load(sleepStruct,'-mat')
    whiskStruct = [animalID '_EventData.mat'];
    load(whiskStruct,'-mat')
    restStruct = [animalID '_RestData.mat'];
    load(restStruct,'-mat')
    % nrem
    for bb = 1:length(SleepData.Forest.NREM.data.CBV_HbT.LH)
        % HbT
        nremCatLH_HbT = cat(2,nremCatLH_HbT,SleepData.Forest.NREM.data.CBV_HbT.LH{bb,1});
        nremMeanLH_HbT = cat(2,nremMeanLH_HbT,mean(SleepData.Forest.NREM.data.CBV_HbT.RH{bb,1}));
        nremCatRH_HbT = cat(2,nremCatRH_HbT,SleepData.Forest.NREM.data.CBV_HbT.LH{bb,1});
        nremMeanRH_HbT = cat(2,nremMeanRH_HbT,mean(SleepData.Forest.NREM.data.CBV_HbT.LH{bb,1}));
        % Gamma
        nremCatLH_Gam = cat(2,nremCatLH_Gam,SleepData.Forest.NREM.data.cortical_LH.gammaBandPower{bb,1});
        nremMeanLH_Gam = cat(2,nremMeanLH_Gam,mean(SleepData.Forest.NREM.data.cortical_LH.gammaBandPower{bb,1}));
        nremCatRH_Gam = cat(2,nremCatRH_Gam,SleepData.Forest.NREM.data.cortical_RH.gammaBandPower{bb,1});
        nremMeanRH_Gam = cat(2,nremMeanRH_Gam,mean(SleepData.Forest.NREM.data.cortical_RH.gammaBandPower{bb,1}));
    end
    % rem
    for cc = 1:length(SleepData.Forest.REM.data.CBV_HbT.LH)
        % HbT
        remCatLH_HbT = cat(2,remCatLH_HbT,SleepData.Forest.REM.data.CBV_HbT.LH{cc,1});
        remMeanLH_HbT = cat(2,remMeanLH_HbT,mean(SleepData.Forest.REM.data.CBV_HbT.RH{cc,1}));
        remCatRH_HbT = cat(2,remCatRH_HbT,SleepData.Forest.REM.data.CBV_HbT.LH{cc,1});
        remMeanRH_HbT = cat(2,remMeanRH_HbT,mean(SleepData.Forest.REM.data.CBV_HbT.LH{cc,1}));
        % Gamma
        remCatLH_Gam = cat(2,remCatLH_Gam,SleepData.Forest.REM.data.cortical_LH.gammaBandPower{cc,1});
        remMeanLH_Gam = cat(2,remMeanLH_Gam,mean(SleepData.Forest.REM.data.cortical_LH.gammaBandPower{cc,1}));
        remCatRH_Gam = cat(2,remCatRH_Gam,SleepData.Forest.REM.data.cortical_RH.gammaBandPower{cc,1});
        remMeanRH_Gam = cat(2,remMeanRH_Gam,mean(SleepData.Forest.REM.data.cortical_RH.gammaBandPower{cc,1}));
    end
end
% awakeGamma = cat(
% awakeHbT = cat(
nremGamma = cat(2,nremMeanLH_Gam,nremMeanRH_Gam);
nremHbT = cat(2,nremMeanLH_HbT,nremMeanRH_HbT);
remGamma = cat(2,remMeanLH_Gam,remMeanRH_Gam);
remHbT = cat(2,remMeanLH_HbT,remMeanRH_HbT);


nremGamma = cat(2,nremCatLH_Gam,nremCatRH_Gam);
nremHbT = cat(2,nremCatLH_HbT,nremCatRH_HbT);
remGamma = cat(2,remCatLH_Gam,remCatRH_Gam);
remHbT = cat(2,remCatLH_HbT,remCatRH_HbT);

%
figure;
% scatter(awakeGamma,awakeHbT,'.','MarkerFaceColor',colorA)
hold on
scatter(nremGamma,nremHbT,'MarkerFaceColor',colorB)
hold on
scatter(remGamma,remHbT,'MarkerFaceColor',colorC)

figure; 
% h1 = histogram2(awakeGamma,awakeHbT,'Normalization','probability','BinWidth',[0.25,10],'XBinLimits',[-.5,2],'YBinLimits',[-50,150]);
hold on
h2 = histogram2(nremGamma,nremHbT,'Normalization','probability','BinWidth',[0.25,10],'XBinLimits',[-.5,2],'YBinLimits',[-50,150]);
hold on
h3 = histogram2(remGamma,remHbT,'Normalization','probability','BinWidth',[0.25,10],'XBinLimits',[-.5,2],'YBinLimits',[-50,150]);

end