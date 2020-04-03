function [] = AvgSleepProbability_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Calculate the average coherence of different behavioral states
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
bins = {'five','ten','fifteen','twenty','twentyfive','thirty','thirtyfive','forty','fortyfive','fifty','fiftyfive','sixty','sixtyplus'};
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color

%% cd through each animal's directory and extract the appropriate analysis results
c = 1;
e = 1;
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    % rest event probability
    for b = 1:length(bins)
        bin = bins{1,b};
        data.(bin).ind{c,1} = AnalysisResults.(animalID).SleepProbability.(bin).awakeLogical;
    end
    % awake hypnogram probability
    c = c + 1;
    strDays = fields(AnalysisResults.(animalID).SleepProbability.Hypnogram);
    q = 1;
    for d = 1:length(strDays)
        strDay = strDays{d,1};
        data.hypAwakeProb.all{e,1} = AnalysisResults.(animalID).SleepProbability.Hypnogram.(strDay).AwakeProb_inds;
        data.hypNREMProb.all{e,1} = AnalysisResults.(animalID).SleepProbability.Hypnogram.(strDay).NREMProb_inds;
        data.hypREMProb.all{e,1} = AnalysisResults.(animalID).SleepProbability.Hypnogram.(strDay).REMProb_inds;
        data.hypAwakeProb.(animalID).ind{q,1} = AnalysisResults.(animalID).SleepProbability.Hypnogram.(strDay).AwakeProb_inds;
        data.hypNREMProb.((animalID)).ind{q,1} = AnalysisResults.(animalID).SleepProbability.Hypnogram.(strDay).NREMProb_inds;
        data.hypREMProb.((animalID)).ind{q,1} = AnalysisResults.(animalID).SleepProbability.Hypnogram.(strDay).REMProb_inds;
        e = e + 1;
        q = q + 1;
    end
end
% concatenate and organize the data
for f = 1:length(bins)
    bin = bins{1,f};
    % rest event probability
    data.(bin).all = [];
    for g = 1:length(data.(bin).ind)
        data.(bin).all = cat(1,data.(bin).all,data.(bin).ind{g,1});
    end
    for q = 1:length(data.(bin).ind)
        data.(bin).indProb{q,1} = sum(data.(bin).ind{q,1})/length(data.(bin).ind{q,1});
    end
    % awake hypnogram probability
    for h = 1:length(data.hypAwakeProb.all)
        awakeHypLength(h,1) = length(data.hypAwakeProb.all{h,1}); %#ok<*AGROW>
    end
    maxHypLength = max(awakeHypLength);
    for i = 1:length(data.hypAwakeProb.all)
        indHypLength = length(data.hypAwakeProb.all{i,1});
        lenDiff = maxHypLength - indHypLength;
        nanPad = NaN(1,lenDiff);
        awakePadHypData = cat(2,data.hypAwakeProb.all{i,1},nanPad);
        awakeAllHypData(i,:) = awakePadHypData;
        nremPadHypData = cat(2,data.hypNREMProb.all{i,1},nanPad);
        nremAllHypData(i,:) = nremPadHypData;
        remPadHypData = cat(2,data.hypREMProb.all{i,1},nanPad);
        remAllHypData(i,:) = remPadHypData;
    end
end
%
for x = 1:length(animalIDs)
    animalID = animalIDs{1,x};
    indRestEventProbability = [];
    for f = 1:length(bins)
        bin = bins{1,f};
        indRestEventProbability = cat(1,indRestEventProbability,data.(bin).indProb{x,1});
    end
    finalRestEventProb{x,1} = indRestEventProbability;
    %
    for z = 1:length(data.hypAwakeProb.(animalID).ind)
        indHypLength = length(data.hypAwakeProb.(animalID).ind{z,1});
        lenDiff = maxHypLength - indHypLength;
        nanPad = NaN(1,lenDiff);
        indAwakePadHypData = cat(2,data.hypAwakeProb.(animalID).ind{z,1},nanPad);
        allIndAwakeAllHypData(z,:) = indAwakePadHypData;
        indNremPadHypData = cat(2,data.hypNREMProb.(animalID).ind{z,1},nanPad);
        allIndNremAllHypData(z,:) = indNremPadHypData;
        indRemPadHypData = cat(2,data.hypREMProb.(animalID).ind{z,1},nanPad);
        allIndRemAllHypData(z,:) = indRemPadHypData;
    end
    data.hypAwakeProb.(animalID).awakeProb = allIndAwakeAllHypData;
    data.hypNREMProb.(animalID).nremProb = allIndNremAllHypData;
    data.hypREMProb.(animalID).remProb = allIndRemAllHypData;
end
% calculate rest event awakeness probability over time
for j = 1:length(bins)
    bin = bins{1,j};
    restEventProbability(j,1) = sum(data.(bin).all)/length(data.(bin).all);
end
% calculate awake probabilty over time
awakeProbability = nansum(awakeAllHypData)./(length(data.hypAwakeProb.all) - sum(isnan(awakeAllHypData)));
awakeNaNInds = ~isnan(awakeProbability);
awakeDiffs = cumsum(awakeNaNInds - diff([1,awakeNaNInds])/2);
patchedAwakeProbability = interp1(1:nnz(awakeNaNInds),awakeProbability(awakeNaNInds),awakeDiffs);
% NREM
nremProbability = nansum(nremAllHypData)./(length(data.hypNREMProb.all) - sum(isnan(nremAllHypData)));
nremNaNInds = ~isnan(nremProbability);
nremDiffs = cumsum(nremNaNInds - diff([1,nremNaNInds])/2);
patchedNREMProbability = interp1(1:nnz(nremNaNInds),nremProbability(nremNaNInds),nremDiffs);
% REM
remProbability = nansum(remAllHypData)./(length(data.hypREMProb.all) - sum(isnan(remAllHypData)));
remNaNInds = ~isnan(remProbability);
remDiffs = cumsum(remNaNInds - diff([1,remNaNInds])/2);
patchedREMProbability = interp1(1:nnz(remNaNInds),remProbability(remNaNInds),remDiffs);
%
dataLength = (3*60*60)/5;   % 3 hrs - 60 minutes - 60 seconds - 5 second bins
patchedAwakeProbability = patchedAwakeProbability(1:dataLength);
patchedNREMProbability = patchedNREMProbability(1:dataLength);
patchedREMProbability = patchedREMProbability(1:dataLength);
%
for qx = 1:length(animalIDs)
    animalID = animalIDs{1,qx};
    % calculate awake probabilty over time
    indAwakeProbability{qx,1} = nansum(data.hypAwakeProb.(animalID).awakeProb)./(size(data.hypAwakeProb.(animalID).awakeProb,1) - sum(isnan(data.hypAwakeProb.(animalID).awakeProb)));
    indAwakeNaNInds{qx,1} = ~isnan(indAwakeProbability{qx,1});
    indAwakeDiffs{qx,1} = cumsum(indAwakeNaNInds{qx,1} - diff([1,indAwakeNaNInds{qx,1}])/2);
    indPatchedAwakeProbability{qx,1} = interp1(1:nnz(indAwakeNaNInds{qx,1}),indAwakeProbability{qx,1}(indAwakeNaNInds{qx,1}),indAwakeDiffs{qx,1});
    finalPatchedAwakeProbability{qx,1} = indPatchedAwakeProbability{qx,1}(1:dataLength);
    % NREM
    indNREMProbability{qx,1} = nansum(data.hypNREMProb.(animalID).nremProb)./(size(data.hypNREMProb.(animalID).nremProb,1) - sum(isnan(data.hypNREMProb.(animalID).nremProb)));
    indNREMNaNInds{qx,1} = ~isnan(indNREMProbability{qx,1});
    indNREMDiffs{qx,1} = cumsum(indNREMNaNInds{qx,1} - diff([1,indNREMNaNInds{qx,1}])/2);
    indPatchedNREMProbability{qx,1} = interp1(1:nnz(indNREMNaNInds{qx,1}),indNREMProbability{qx,1}(indNREMNaNInds{qx,1}),indNREMDiffs{qx,1});
    finalPatchedNREMProbability{qx,1} = indPatchedNREMProbability{qx,1}(1:dataLength);
    % REM
    indREMProbability{qx,1} = nansum(data.hypREMProb.(animalID).remProb)./(size(data.hypREMProb.(animalID).remProb,1) - sum(isnan(data.hypREMProb.(animalID).remProb)));
    indREMNaNInds{qx,1} = ~isnan(indREMProbability{qx,1});
    indREMDiffs{qx,1} = cumsum(indREMNaNInds{qx,1} - diff([1,indREMNaNInds{qx,1}])/2);
    indPatchedREMProbability{qx,1} = interp1(1:nnz(indREMNaNInds{qx,1}),indREMProbability{qx,1}(indREMNaNInds{qx,1}),indREMDiffs{qx,1});
    finalPatchedREMProbability{qx,1} = indPatchedREMProbability{qx,1}(1:dataLength);
end
%
binSize = 60/5;   % 60 sec divided by 5 sec bins
numBins = length(patchedAwakeProbability)/binSize;
for k = 1:numBins
    if k == 1
        binnedAwakeProbability(1,k) = mean(patchedAwakeProbability(1:binSize));
        binnedNREMProbability(1,k) = mean(patchedNREMProbability(1:binSize));
        binnedREMProbability(1,k) = mean(patchedNREMProbability(1:binSize));
    else
        binnedAwakeProbability(1,k) = mean(patchedAwakeProbability((k - 1)*binSize + 1:k*binSize));
        binnedNREMProbability(1,k) = mean(patchedNREMProbability((k - 1)*binSize + 1:k*binSize));
        binnedREMProbability(1,k) = mean(patchedREMProbability((k - 1)*binSize + 1:k*binSize));
    end
end
%
for qx = 1:length(animalIDs)
    for k = 1:numBins
        if k == 1
            binnedIndFinalPatchedAwakeProbability{qx,1}(1,k) = mean(finalPatchedAwakeProbability{qx,1}(1:binSize));
            binnedIndFinalPatchedNREMProbability{qx,1}(1,k) = mean(finalPatchedNREMProbability{qx,1}(1:binSize));
            binnedIndFinalPatchedREMProbability{qx,1}(1,k) = mean(finalPatchedREMProbability{qx,1}(1:binSize));
        else
            binnedIndFinalPatchedAwakeProbability{qx,1}(1,k) =  mean(finalPatchedAwakeProbability{qx,1}((k - 1)*binSize + 1:k*binSize));
            binnedIndFinalPatchedNREMProbability{qx,1}(1,k) =  mean(finalPatchedNREMProbability{qx,1}((k - 1)*binSize + 1:k*binSize));
            binnedIndFinalPatchedREMProbability{qx,1}(1,k) =  mean(finalPatchedREMProbability{qx,1}((k - 1)*binSize + 1:k*binSize));
        end
    end
end

%% summary figure(s)
summaryFigure = figure;
%% Hypnogram probability
subplot(2,2,1);
xinds1 = (1:numBins)/(numBins/3);
% awake
s1 = scatter(xinds1,binnedAwakeProbability,'MarkerEdgeColor','k','MarkerFaceColor',colorA);
[awakeHypExpCurve,awakeHypGOF] = fit(xinds1',binnedAwakeProbability','exp1');
awakeHypExpFit = awakeHypExpCurve(xinds1);
hold on
p1 = plot(xinds1,awakeHypExpFit,'color',colorA,'LineWidth',2);
% nrem
s2 = scatter(xinds1,binnedNREMProbability,'MarkerEdgeColor','k','MarkerFaceColor',colorB);
[nremHypExpCurve,nremHypGOF] = fit(xinds1',binnedNREMProbability','exp1');
nremHypExpFit = nremHypExpCurve(xinds1);
p2 = plot(xinds1,nremHypExpFit,'color',colorB,'LineWidth',2);
% rem
s3 = scatter(xinds1,binnedREMProbability,'MarkerEdgeColor','k','MarkerFaceColor',colorC);
[remHypExpCurve,remHypGOF] = fit(xinds1',binnedREMProbability','exp1');
remHypExpFit = remHypExpCurve(xinds1);
p3 = plot(xinds1,remHypExpFit,'color',colorC,'LineWidth',2);
% legend and figure labels
legend([s1,p1,s2,p2,s3,p3],'Awake bin prob',['(Exp2) adjR^2 = ' num2str(round(awakeHypGOF.adjrsquare,2))],'NREM bin prob',['(Exp2) adjR^2 = ' num2str(round(nremHypGOF.adjrsquare,2))],'REM bin prob',['(Exp2) adjR^2 = ' num2str(round(remHypGOF.adjrsquare,2))])
title('Trial duration awake probability')
xlabel('Trial duration (Hr)')
ylabel('Probability of wakefulness')
ylim([0,1])
set(gca,'box','off')
axis square

%% Rest event probability
subplot(2,2,2);
xinds2 = 0:length(bins) - 1;
xinds3 = 0:0.01:length(bins) - 1;
s2 = scatter(xinds2,restEventProbability,'MarkerEdgeColor','k','MarkerFaceColor',colorA);
[restExpCurve,restGOF] = fit(xinds2',restEventProbability,'exp1');
restExpFit = restExpCurve(xinds3);
hold on
p2 = plot(xinds3,restExpFit,'k','LineWidth',2);
legend([s2,p2],'Time point probability',['(Exp2) adjR^2 = ' num2str(round(restGOF.adjrsquare,3))])
title('Rest event awake probability')
xticks([1,3,5,7,9,11,12])
xticklabels({'10','20','30','40','50','60','60+'})
xlabel('Rest duration (s)')
ylabel('Probability of wakefulness')
ylim([0,1])
set(gca,'box','off')
axis square

subplot(2,2,3)
plot(xinds1,awakeHypExpFit,'-','color','k','LineWidth',4);
hold on
for a = 1:length(animalIDs)
    s = scatter(xinds1,binnedIndFinalPatchedAwakeProbability{a,1},'MarkerEdgeColor','k','MarkerFaceColor',colorA);
    [indAwakeHypExpCurve{a,1},indAwakeHypGOF{a,1}] = fit(xinds1',binnedIndFinalPatchedAwakeProbability{a,1}','exp1');
    indAwakeHypExpFit = indAwakeHypExpCurve{a,1}(xinds1);
    plot(xinds1,indAwakeHypExpFit,'-','color',colorA,'LineWidth',0.5);
    delete(s)
end
% nrem
plot(xinds1,nremHypExpFit,'-','color','k','LineWidth',4);
hold on
for a = 1:length(animalIDs)
    s = scatter(xinds1,binnedIndFinalPatchedNREMProbability{a,1},'MarkerEdgeColor','k','MarkerFaceColor',colorB);
    [indNREMHypExpCurve{a,1},indNREMHypGOF{a,1}] = fit(xinds1',binnedIndFinalPatchedNREMProbability{a,1}','exp1');
    indNREMHypExpFit{a,1} = indNREMHypExpCurve{a,1}(xinds1);
    plot(xinds1,indNREMHypExpFit{a,1},'-','color',colorB,'LineWidth',0.5);
    delete(s)
end
% rem
plot(xinds1,remHypExpFit,'-','color','k','LineWidth',4);
hold on
for a = 1:length(animalIDs)
    s = scatter(xinds1,binnedIndFinalPatchedREMProbability{a,1},'MarkerEdgeColor','k','MarkerFaceColor',colorC);
    [indREMHypExpCurve{a,1},indREMHypGOF{a,1}] = fit(xinds1',binnedIndFinalPatchedREMProbability{a,1}','exp1');
    indREMHypExpFit{a,1} = indREMHypExpCurve{a,1}(xinds1);
    plot(xinds1,indREMHypExpFit{a,1},'-','color',colorC,'LineWidth',0.5);
    delete(s)
end
title('Trial duration awake probability')
xlabel('Trial duration (Hr)')
ylabel('Probability of wakefulness')
ylim([0,1])
set(gca,'box','off')
axis square

subplot(2,2,4)
plot(xinds3,restExpFit,'k','LineWidth',4);
hold on
for a = 1:length(animalIDs)
    s = scatter(xinds2,finalRestEventProb{a,1},'MarkerEdgeColor','k','MarkerFaceColor',colorA);
    [indRestExpCurve{a,1},indRestGOF{a,1}] = fit(xinds2',finalRestEventProb{a,1},'exp1');
    indRestExpFit{a,1} = indRestExpCurve{a,1}(xinds3);
    plot(xinds3,indRestExpFit{a,1},'-','color',colorA,'LineWidth',0.5);
    delete(s)
end
title('Rest event awake probability')
xticks([1,3,5,7,9,11,12])
xticklabels({'10','20','30','40','50','60','60+'})
xlabel('Rest duration (s)')
ylabel('Probability of wakefulness')
ylim([0,1])
set(gca,'box','off')
axis square

%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Awake Probability']);

end

