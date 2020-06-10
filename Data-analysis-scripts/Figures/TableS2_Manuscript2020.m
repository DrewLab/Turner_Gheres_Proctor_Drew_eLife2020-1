function [] = TableS2_Manuscript2020(rootFolder)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate Table S2 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% Set-up and process data for Table S2
animalIDs = {'T115','T116','T117','T118','T125','T126'};
% extract data from each animal's sleep scoring results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    dataLoc = [rootFolder '/' animalID '/2P Data/'];
    cd(dataLoc)
    % find and load RestingBaselines.mat strut
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % character list of all merged data file IDs
    mergedDirectory = dir('*_MergedData.mat');
    mergedDataFiles = {mergedDirectory.name}';
    mergedDataFileIDs = char(mergedDataFiles);
    % take out logicals from each vessel to determine how much data it has
    data.(animalID) = [];
    for bb = 1:size(mergedDataFileIDs,1)
        [~,~,~,~,~,vID] = GetFileInfo2_2P_Manuscript2020(mergedDataFileIDs(bb,:));
        if strcmp(vID(1),'V') == false
            load(mergedDataFileIDs(bb,:),'-mat')
            if isfield(data.(animalID),vID) == false
                data.(animalID).(vID).awakeData = sum(MergedData.sleep.logicals.Manual.awakeLogical);
                data.(animalID).(vID).nremData = sum(MergedData.sleep.logicals.Manual.nremLogical);
                data.(animalID).(vID).remData = sum(MergedData.sleep.logicals.Manual.remLogical);
            else
                data.(animalID).(vID).awakeData = data.(animalID).(vID).awakeData + sum(MergedData.sleep.logicals.Manual.awakeLogical);
                data.(animalID).(vID).nremData = data.(animalID).(vID).nremData + sum(MergedData.sleep.logicals.Manual.nremLogical);
                data.(animalID).(vID).remData =  data.(animalID).(vID).remData + sum(MergedData.sleep.logicals.Manual.remLogical);
            end
        end
    end
    % find baseline of each vessel ID
    uniqueVIDs = fieldnames(data.(animalID));
    for cc = 1:length(uniqueVIDs)
        diamBaseline = [];
        strDays = fieldnames(RestingBaselines.manualSelection.vesselDiameter.data.(uniqueVIDs{cc,1}));
        for dd = 1:length(strDays)
            diamBaseline(dd,1) = RestingBaselines.manualSelection.vesselDiameter.data.(uniqueVIDs{cc,1}).(strDays{dd,1}); %#ok<*AGROW>
        end
        data.(animalID).(uniqueVIDs{cc,1}).baseline = mean(diamBaseline);
    end
    cd(rootFolder)
end
% go through and put each animal/vessel ID into an table
gg = 1;
labelTime = 5;   % seconds
for ee = 1:length(animalIDs)
    animalID = animalIDs{1,ee};
    uniqueVIDs = fieldnames(data.(animalID));
    TwoP_indTotalTimeAwake = [];
    TwoP_indTotalTimeNREM = [];
    TwoP_indTotalTimeREM = [];
    for ff = 1:length(uniqueVIDs)
        TwoP_animalIDs{gg,1} = [animalID '_' uniqueVIDs{ff,1}];
        TwoP_baselineDiams(gg,1) = round(data.(animalID).(uniqueVIDs{ff,1}).baseline,1);
        TwoP_totalTimeAwake(gg,1) = round(((data.(animalID).(uniqueVIDs{ff,1}).awakeData*labelTime)/60),1);   % sec -> min
        TwoP_totalTimeNREM(gg,1) = round(((data.(animalID).(uniqueVIDs{ff,1}).nremData*labelTime)/60),1);   % sec -> min -> hrs
        TwoP_totalTimeREM(gg,1) = round(((data.(animalID).(uniqueVIDs{ff,1}).remData*labelTime)/60),1);   % sec -> min -> hrs
        TwoP_totalTimeMins(gg,1) = round(TwoP_totalTimeAwake(gg,1) + TwoP_totalTimeNREM(gg,1) + TwoP_totalTimeREM(gg,1),1); 
        gg = gg + 1;
        % per animal
        TwoP_indTotalTimeAwake(ff,1) = ((data.(animalID).(uniqueVIDs{ff,1}).awakeData*labelTime)/60)/60;   % sec -> min -> hrs
        TwoP_indTotalTimeNREM(ff,1) = ((data.(animalID).(uniqueVIDs{ff,1}).nremData*labelTime)/60)/60;   % sec -> min -> hrs
        TwoP_indTotalTimeREM(ff,1) = ((data.(animalID).(uniqueVIDs{ff,1}).remData*labelTime)/60)/60;   % sec -> min -> hrs    
    end
    TwoP_totalTimePerAnimal(ee,1) = sum(TwoP_indTotalTimeAwake) + sum(TwoP_indTotalTimeNREM) + sum(TwoP_indTotalTimeREM);
end
TwoP_allTimeHours = round(sum(TwoP_totalTimeMins)/60);
TwoP_meanTimeHours = round(mean(TwoP_totalTimePerAnimal,1),1);
TwoP_stdTimeHours = round(std(TwoP_totalTimePerAnimal,0,1),1);

%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
%% Table S2 
summaryTable = figure('Name','TableS2');
sgtitle('Table S2 Turner Manuscript 2020')
variableNames = {'TotalTimeMins','BaseDiamUm','AwakeTimeMins','NREMTimeMins','REMTimeMins'};
T = table(TwoP_totalTimeMins,TwoP_baselineDiams,TwoP_totalTimeAwake,TwoP_totalTimeNREM,TwoP_totalTimeREM,'RowNames',TwoP_animalIDs,'VariableNames',variableNames);
T2 = sortrows(T,'RowNames');
uitable('Data',T2{:,:},'ColumnName',T2.Properties.VariableNames,'RowName',T2.Properties.RowNames,'Units','Normalized','Position',[0,0,1,1]);
uicontrol('Style','text','Position',[700,600,100,150],'String',{'Total Time (Hrs): ' num2str(TwoP_allTimeHours),'Mean time per animal (Hrs): ' num2str(TwoP_meanTimeHours) ' +/- ' num2str(TwoP_stdTimeHours)});
savefig(summaryTable,[dirpath 'TableS2']);

end
