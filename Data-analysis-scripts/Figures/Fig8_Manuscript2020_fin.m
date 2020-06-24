function [] = Fig8_Manuscript2020_fin(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
modelType = 'Forest';
dataTypes = {'gammaBandPower'};
colorA = [(51/256),(160/256),(44/256)];   % Rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(180/256)];  % Whisk color
% colorE = [(0/256),(256/256),(256/256)];  % Isoflurane color
colorF = [(256/256),(192/256),(0/256)];   % Awake color
colorG = [(0/256),(128/256),(256/256)];   % Sleep color
colorH = [(184/256),(115/256),(51/256)];  % All color
LH_allCatLabels = [];
RH_allCatLabels = [];
LH_HbTallCatMeans = [];
RH_HbTallCatMeans = [];
LH_gamAllCatMeans = [];
RH_gamAllCatMeans = [];
%% extract data from each animal's sleep scoring results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    dataLoc = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLoc)
    % add this animal's scoring labels with the other animal'ss
    scoringResults = 'Forest_ScoringResults.mat';
    load(scoringResults,'-mat')
    baselinesFile = [animalID '_RestingBaselines.mat'];
    load(baselinesFile)
    LH_allCatLabels = cat(1,LH_allCatLabels,ScoringResults.alllabels);
    RH_allCatLabels = cat(1,RH_allCatLabels,ScoringResults.alllabels);
    % take the mean of each 5 second bin
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    binSize = 5;   % seconds
    numBins = 180;   % 15 minutes with 5 sec bins
    samplingRate = 30;   % Hz
    samplesPerBin = binSize*samplingRate;
    for bb = 1:size(procDataFileIDs,1)
        procDataFileID = procDataFileIDs(bb,:);
        load(procDataFileID,'-mat')
        [~,fileDate,~] = GetFileInfo_IOS_Manuscript2020(procDataFileID);
        strDay = ConvertDate_IOS_Manuscript2020(fileDate);
        for cc = 1:numBins
            if cc == 1
                LH_HbTbinSamples = ProcData.data.CBV_HbT.adjLH(1:samplesPerBin);
                RH_HbTbinSamples = ProcData.data.CBV_HbT.adjRH(1:samplesPerBin);               
                LH_gamBinSamples = (ProcData.data.cortical_LH.gammaBandPower(1:samplesPerBin) - RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay))./RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay);
                RH_gamBinSamples = (ProcData.data.cortical_RH.gammaBandPower(1:samplesPerBin) - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay))./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay);
            else
                LH_HbTbinSamples = ProcData.data.CBV_HbT.adjLH((cc - 1)*samplesPerBin + 1:cc*samplesPerBin);
                RH_HbTbinSamples = ProcData.data.CBV_HbT.adjRH((cc - 1)*samplesPerBin + 1:cc*samplesPerBin);
                LH_gamBinSamples = (ProcData.data.cortical_LH.gammaBandPower((cc - 1)*samplesPerBin + 1:cc*samplesPerBin) - RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay))./RestingBaselines.manualSelection.cortical_LH.gammaBandPower.(strDay);
                RH_gamBinSamples = (ProcData.data.cortical_RH.gammaBandPower((cc - 1)*samplesPerBin + 1:cc*samplesPerBin) - RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay))./RestingBaselines.manualSelection.cortical_RH.gammaBandPower.(strDay);
            end
            LH_HbTallCatMeans = cat(1,LH_HbTallCatMeans,mean(LH_HbTbinSamples));
            RH_HbTallCatMeans = cat(1,RH_HbTallCatMeans,mean(RH_HbTbinSamples));
            LH_gamAllCatMeans = cat(1,LH_gamAllCatMeans,mean(LH_gamBinSamples));
            RH_gamAllCatMeans = cat(1,RH_gamAllCatMeans,mean(RH_gamBinSamples));
        end
    end
end
% concatenage LH/RH labels and mean values.
allCatLabels = cat(1,LH_allCatLabels,RH_allCatLabels);
HbTallCatMeans = cat(1,LH_HbTallCatMeans,RH_HbTallCatMeans);
gamAllCatMeans = cat(1,LH_gamAllCatMeans,RH_gamAllCatMeans);
%% 
zAwake = 1; zNrem = 1; zRem = 1;
for zz = 1:length(allCatLabels)
    if strcmp(allCatLabels{zz,1},'Not Sleep') == true
        awakeHbT(zAwake,1) = HbTallCatMeans(zz,1); %#ok<*NASGU>
        awakeGamma(zAwake,1) = gamAllCatMeans(zz,1);
        zAwake = zAwake + 1;
    elseif strcmp(allCatLabels{zz,1},'NREM Sleep') == true
        nremHbT(zNrem,1) = HbTallCatMeans(zz,1);
        nremGamma(zNrem,1) = gamAllCatMeans(zz,1);
        zNrem = zNrem + 1;
    elseif strcmp(allCatLabels{zz,1},'REM Sleep') == true
        remHbT(zRem,1) = HbTallCatMeans(zz,1);
        remGamma(zRem,1) = gamAllCatMeans(zz,1);
        zRem = zRem + 1;
    end
end
% 
% figure;
% scatter(awakeGamma,awakeHbT,'.','MarkerFaceColor',colorA)
% hold on
% scatter(nremGamma,nremHbT,'.','MarkerFaceColor',colorB)
% scatter(remGamma,remHbT,'.','MarkerFaceColor',colorC)
% xlim([0,1])
% 
% figure; 
% h1 = histogram2(awakeGamma,awakeHbT,'Normalization','probability','BinWidth',[0.25,10],'XBinLimits',[-.5,2],'YBinLimits',[-50,150]);
% hold on
% h2 = histogram2(nremGamma,nremHbT,'Normalization','probability','BinWidth',[0.25,10],'XBinLimits',[-.5,2],'YBinLimits',[-50,150]);
% h3 = histogram2(remGamma,remHbT,'Normalization','probability','BinWidth',[0.25,10],'XBinLimits',[-.5,2],'YBinLimits',[-50,150]);

%% put each mean and scoring label into a cell
minHbT = floor(min(HbTallCatMeans));
maxHbT = ceil(max(HbTallCatMeans));
awakeBins = minHbT:1:maxHbT;
cutDown = abs(minHbT - (-35));
cutUp = maxHbT - 115;
probBinLabels = cell(length(minHbT:1:maxHbT),1);
probBinMeans = cell(length(minHbT:1:maxHbT),1);
discBins = discretize(HbTallCatMeans,awakeBins);
for dd = 1:length(discBins)
    probBinLabels{discBins(dd),1} = cat(1,probBinLabels{discBins(dd),1},{allCatLabels(dd,1)});
    probBinMeans{discBins(dd),1} = cat(1,probBinMeans{discBins(dd),1},{HbTallCatMeans(dd,1)});
end
% condense the left edges of the histogram bins to -35:1:120
cutDownLabels = [];
cutDownMeans = [];
for ee = 1:cutDown
    cutDownLabels = cat(1,cutDownLabels,probBinLabels{ee,1});
    cutDownMeans = cat(1,cutDownMeans,probBinMeans{ee,1});
end
% condense the right edges of the histogram to -35:1:120
cutUpLabels = [];
cutUpMeans = [];
for ff = 1:cutUp + 1
    cutUpLabels = cat(1,cutUpLabels,probBinLabels{end - ff,1});
    cutUpMeans = cat(1,cutUpMeans,probBinMeans{end - ff,1});
end
% reconstruct array of labels based on new edges
finCatLabels = cat(1,{cutDownLabels},probBinLabels(cutDown + 1:end - (cutUp + 2)),{cutUpLabels});
% strcmp the bins and if the bin is asleep (NREM/REM) set to 0, else set 1
for gg = 1:length(finCatLabels)
    for hh = 1:length(finCatLabels{gg,1})
        if strcmp(finCatLabels{gg,1}{hh,1}{1,1},'Not Sleep') == true
            awakeProbEvents{gg,1}(hh,1) = 1; %#ok<*AGROW>
        else
            awakeProbEvents{gg,1}(hh,1) = 0;
        end
    end
end
% strcmp the bins and if the bin is not in NREM (Awake/REM) set to 0, else set 1
for gg = 1:length(finCatLabels)
    for hh = 1:length(finCatLabels{gg,1})
        if strcmp(finCatLabels{gg,1}{hh,1}{1,1},'NREM Sleep') == true
            nremProbEvents{gg,1}(hh,1) = 1;
        else
            nremProbEvents{gg,1}(hh,1) = 0;
        end
    end
end
% strcmp the bins and if the bin is not in REM (Awake/NREM) set to 0, else set 1
for gg = 1:length(finCatLabels)
    for hh = 1:length(finCatLabels{gg,1})
        if strcmp(finCatLabels{gg,1}{hh,1}{1,1},'REM Sleep') == true
            remProbEvents{gg,1}(hh,1) = 1;
        else
            remProbEvents{gg,1}(hh,1) = 0;
        end
    end
end
% take probability of each bin
for ii = 1:length(awakeProbEvents)
    awakeProbPerc(ii,1) = sum(awakeProbEvents{ii,1})/length(awakeProbEvents{ii,1})*100;
    nremProbPerc(ii,1) = sum(nremProbEvents{ii,1})/length(nremProbEvents{ii,1})*100;
    remProbPerc(ii,1) = sum(remProbEvents{ii,1})/length(remProbEvents{ii,1})*100;
end
%% Average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Coherr = [];
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.Coherr,behavField) == false
            data.Coherr.(behavField) = [];
        end
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).Coherence.(behavField).(dataType).C) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.Coherr.(behavField),dataType) == false
                    data.Coherr.(behavField).(dataType).LH_C = [];
                    data.Coherr.(behavField).(dataType).RH_C = [];
                    data.Coherr.(behavField).(dataType).f = [];
                    data.Coherr.(behavField).(dataType).confC = [];
                end
                % concatenate C/f for existing data - exclude any empty sets
                data.Coherr.(behavField).(dataType).LH_C = cat(2,data.Coherr.(behavField).(dataType).C,AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                data.Coherr.(behavField).(dataType).RH_C = cat(2,data.Coherr.(behavField).(dataType).C,AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                data.Coherr.(behavField).(dataType).f = cat(1,data.Coherr.(behavField).(dataType).f,AnalysisResults.(animalID).Coherence.(behavField).(dataType).f);
                data.Coherr.(behavField).(dataType).confC = cat(1,data.Coherr.(behavField).(dataType).confC,AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC);
            end
        end
    end
end
% take mean/StD of C/f and determine confC line
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        data.Coherr.(behavField).(dataType).meanC = mean(data.Coherr.(behavField).(dataType).C,2);
        data.Coherr.(behavField).(dataType).stdC = std(data.Coherr.(behavField).(dataType).C,0,2);
        data.Coherr.(behavField).(dataType).meanf = mean(data.Coherr.(behavField).(dataType).f,1);
        data.Coherr.(behavField).(dataType).maxConfC = geomean(data.Coherr.(behavField).(dataType).confC);
        data.Coherr.(behavField).(dataType).maxConfC_Y = ones(length(data.Coherr.(behavField).(dataType).meanf),1)*data.Coherr.(behavField).(dataType).maxConfC;
    end
end
%% Figure Panel 8
summaryFigure = figure('Name','Fig8 (a,c)');
sgtitle('Figure Panel 8 (a,c) Turner Manuscript 2020')
%% [8a] HbT vs. arousal state probability
ax1 = subplot(1,2,1);
edges = -35:1:115;
yyaxis right
h1 = histogram(HbTallCatMeans,edges,'Normalization','probability','EdgeColor','k','FaceColor',colors_Manuscript2020('dark candy apple red'));
ylabel({'5-sec Mean \DeltaHbT','Probability distribution'},'rotation',-90,'VerticalAlignment','bottom')
yyaxis left
p1 = plot(edges,sgolayfilt(medfilt1(awakeProbPerc,10,'truncate'),3,17),'-','color',colors_Manuscript2020('rich black'),'LineWidth',2);
hold on
p2 = plot(edges,sgolayfilt(medfilt1(nremProbPerc,10,'truncate'),3,17),'-','color',colorB,'LineWidth',2);
p3 = plot(edges,sgolayfilt(medfilt1(remProbPerc,10,'truncate'),3,17),'-','color',colorC,'LineWidth',2);
ylabel({'Arousal-state probability (%)'})
xlim([-35,115])
ylim([0,85])
legend([p1,p2,p3,h1],'Awake','NREM','REM','\DeltaHbT','Location','NorthEast')
title({'[8a] 5-sec mean \DeltaHbT (\muM)','vs. arousal state probability',''})
xlabel({'\DeltaHbT (\muM)','1 \muM bins'})
axis square
set(gca,'box','off')
set(gca,'TickLength',[0.03,0.03]);
ylim([0,90])
xlim([-35,115])
set(h1,'facealpha',0.2);
ax1.TickLength = [0.03,0.03];
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = colors_Manuscript2020('dark candy apple red');
%% [8c] Coherence between HbT and gamma-band power during different arousal-states
ax2 = subplot(1,2,2);
s1 = semilogx(data.Coherr.Awake.gammaBandPower.meanf,data.Coherr.Awake.gammaBandPower.meanC,'color',colorF,'LineWidth',2);
hold on
s2 = semilogx(data.Coherr.Rest.gammaBandPower.meanf,data.Coherr.Rest.gammaBandPower.meanC,'color',colorA,'LineWidth',2);
s3 = semilogx(data.Coherr.NREM.(modelType).gammaBandPower.meanf,data.Coherr.NREM.(modelType).gammaBandPower.meanC,'color',colorB,'LineWidth',2);
s4 = semilogx(data.Coherr.REM.(modelType).gammaBandPower.meanf,data.Coherr.REM.(modelType).gammaBandPower.meanC,'color',colorC,'LineWidth',2);
s5 = semilogx(data.Coherr.All.gammaBandPower.meanf,data.Coherr.All.gammaBandPower.meanC,'color','r','LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[8c] Neural-hemo coherence','Gamma-band power and \DeltaHbT \muM (%)',''})
legend([s1,s2,s3,s4,s5],'Awake','Rest','NREM','REM','All','Location','SouthEast')
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Fig8']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Fig8'])

end
