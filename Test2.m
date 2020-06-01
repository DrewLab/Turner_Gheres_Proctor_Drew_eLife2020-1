function [] = Test2(rootFolder)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
LH_allCatLabels = [];
RH_allCatLabels = [];
LH_allCatMeans = [];
RH_allCatMeans = [];
% extract data from each animal's sleep scoring results
for aa = 1:length(animalIDs)
    animalID = animalIDs{1,aa};
    dataLoc = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLoc)
    % add this animal's scoring labels with the other animal'ss
    scoringResults = 'Forest_ScoringResults.mat';
    load(scoringResults,'-mat')
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
        for cc = 1:numBins
            if cc == 1
                LH_binSamples = ProcData.data.CBV_HbT.adjLH(1:samplesPerBin);
                RH_binSamples = ProcData.data.CBV_HbT.adjRH(1:samplesPerBin);
            else
                LH_binSamples = ProcData.data.CBV_HbT.adjLH((cc - 1)*samplesPerBin + 1:cc*samplesPerBin);
                RH_binSamples = ProcData.data.CBV_HbT.adjRH((cc - 1)*samplesPerBin + 1:cc*samplesPerBin);
            end
            LH_allCatMeans = cat(1,LH_allCatMeans,mean(LH_binSamples));
            RH_allCatMeans = cat(1,RH_allCatMeans,mean(RH_binSamples));
        end
    end
end
% concatenage LH/RH labels and mean values.
allCatLabels = cat(1,LH_allCatLabels,RH_allCatLabels);
allCatMeans = cat(1,LH_allCatMeans,RH_allCatMeans);
% put each mean and scoring label into a cell
minHbT = floor(min(allCatMeans));
maxHbT = ceil(max(allCatMeans));
awakeBins = minHbT:1:maxHbT;
cutDown = abs(minHbT - (-35));
cutUp = maxHbT - 115;
probBinLabels = cell(length(minHbT:1:maxHbT),1);
probBinMeans = cell(length(minHbT:1:maxHbT),1);
discBins = discretize(allCatMeans,awakeBins);
for dd = 1:length(discBins)
    probBinLabels{discBins(dd),1} = cat(1,probBinLabels{discBins(dd),1},{allCatLabels(dd,1)});
    probBinMeans{discBins(dd),1} = cat(1,probBinMeans{discBins(dd),1},{allCatMeans(dd,1)});
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
finCatMeans = cat(1,{cutDownMeans},probBinMeans(cutDown + 1:end - (cutUp + 2)),{cutUpMeans});
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
%% summary figure
figure;
ax1 = subplot(1,1,1);
edges = -35:1:115;
yyaxis right
h1 = histogram(allCatMeans,edges,'Normalization','probability','EdgeColor','k','FaceColor',colors_Manuscript2020('dark candy apple red'));
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
title({'5-sec mean \DeltaHbT (\muM)','vs. arousal state probability',''})
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

end
