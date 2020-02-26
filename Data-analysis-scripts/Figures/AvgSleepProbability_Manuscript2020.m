function [] = AvgSleepProbability_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Calculate the average coherence of different behavioral states
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120'};
bins = {'five','ten','fifteen','twenty','twentyfive','thirty','thirtyfive','forty','fortyfive','fifty','fiftyfive','sixty','sixtyplus'};

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
    for d = 1:length(strDays)
        strDay = strDays{d,1};
        data.hypAwakeProb.ind{e,1} = AnalysisResults.(animalID).SleepProbability.Hypnogram.(strDay).AwakeProb_inds;
        e = e + 1;
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
    % awake hypnogram probability
    for h = 1:length(data.hypAwakeProb.ind)
        hypLength(h,1) = length(data.hypAwakeProb.ind{h,1}); %#ok<*AGROW>
    end
    maxHypLength = max(hypLength);
    for i = 1:length(data.hypAwakeProb.ind)
        indHypLength = length(data.hypAwakeProb.ind{i,1});
        lenDiff = maxHypLength - indHypLength;
        nanPad = NaN(1,lenDiff);
        padHypData = cat(2,data.hypAwakeProb.ind{i,1},nanPad);
        allHypData(i,:) = padHypData;
    end
end
% calculate rest event awakeness probability over time
for j = 1:length(bins)
    bin = bins{1,j};
    restEventProbability(j,1) = sum(data.(bin).all)/length(data.(bin).all);
end
% calculate awake probabilty over time
awakeProbability = nansum(allHypData)./(length(data.hypAwakeProb.ind) - sum(isnan(allHypData)));
nanInds = ~isnan(awakeProbability);
diffs = cumsum(nanInds - diff([1,nanInds])/2);
patchedAwakeProbability = interp1(1:nnz(nanInds),awakeProbability(nanInds),diffs);
dataLength = (3*60*60)/5;   % 4 hrs - 60 minutes - 60 seconds - 5 second bins
patchedAwakeProbability = patchedAwakeProbability(1:dataLength);

%% summary figure(s)
summaryFigure = figure;
%% Hypnogram probability
subplot(1,2,1);
xinds1 = (1:dataLength)/(dataLength/3);
p1 = plot(xinds1,patchedAwakeProbability,'color',colors_Manuscript2020('ash grey'));
[hypExpCurve,hypGOF] = fit(xinds1',patchedAwakeProbability','exp2');
hypExpFit = hypExpCurve(xinds1);
hold on
p2 = plot(xinds1,hypExpFit,'k','LineWidth',2);
legend([p1,p2],'Individial bin probability',['(Exp2) adjR^2 = ' num2str(round(hypGOF.adjrsquare,2))])
title('Trial duration awake probability')
xlabel('Trial duration (Hr)')
ylabel('Probability of wakefulness')
ylim([0,1])
set(gca,'box','off')
axis square

%% Rest event probability
subplot(1,2,2);
xinds2 = 0:length(bins) - 1;
xinds3 = 0:0.01:length(bins) - 1; 
s1 = scatter(xinds2,restEventProbability,'MarkerEdgeColor','r');
[restExpCurve,restGOF] = fit(xinds2',restEventProbability,'exp2');
restExpFit = restExpCurve(xinds3);
hold on
p3 = plot(xinds3,restExpFit,'k','LineWidth',2);
legend([s1,p3],'Time point probability',['(Exp2) adjR^2 = ' num2str(round(restGOF.adjrsquare,3))])
title('Rest event awake probability')
xticks([1,3,5,7,9,11,12])
xticklabels({'10','20','30','40','50','60','60+'})
xlabel('Rest duration (s)')
ylabel('Probability of wakefulness')
ylim([0,1])
set(gca,'box','off')
axis square

%% save figure(s)
dirpath = [rootFolder '\Analysis Figures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Awake Probability']);

end

