function [AnalysisResults] = Fig8_Manuscript2020(rootFolder,saveFigs,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Generate figure panel 8 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
dataTypes = {'gammaBandPower'};
colorRest = [(51/256),(160/256),(44/256)];
colorNREM = [(192/256),(0/256),(256/256)];
colorREM = [(255/256),(140/256),(0/256)];
colorAwake = [(256/256),(192/256),(0/256)];
colorSleep = [(0/256),(128/256),(256/256)];
colorAll = [(184/256),(115/256),(51/256)];
% colorWhisk = [(31/256),(120/256),(180/256)];
% colorStim = [(256/256),(28/256),(207/256)];
% colorIso = [(0/256),(256/256),(256/256)];
%% extract data from each animal's sleep scoring results
HbTallCatMeans = AnalysisResults.HbTSleepProbability.HbTCatMeans;
awakeProbPerc = AnalysisResults.HbTSleepProbability.awakeProbPerc;
nremProbPerc = AnalysisResults.HbTSleepProbability.nremProbPerc;
remProbPerc = AnalysisResults.HbTSleepProbability.remProbPerc;
%% average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.NeuralHemoCoherence = [];
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.NeuralHemoCoherence,behavField) == false
            data.NeuralHemoCoherence.(behavField) = [];
        end
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.NeuralHemoCoherence.(behavField),dataType) == false
                    data.NeuralHemoCoherence.(behavField).(dataType).C = [];
                    data.NeuralHemoCoherence.(behavField).(dataType).f = [];
                    data.NeuralHemoCoherence.(behavField).(dataType).confC = [];
                end
                % concatenate C/f for existing data - exclude any empty sets
                data.NeuralHemoCoherence.(behavField).(dataType).C = cat(2,data.NeuralHemoCoherence.(behavField).(dataType).C,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.C);
                data.NeuralHemoCoherence.(behavField).(dataType).f = cat(1,data.NeuralHemoCoherence.(behavField).(dataType).f,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.f,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.f);
                data.NeuralHemoCoherence.(behavField).(dataType).confC = cat(1,data.NeuralHemoCoherence.(behavField).(dataType).confC,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.confC,AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.confC);
            end
        end
    end
end
% take mean/StD of C/f and determine confC line
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        data.NeuralHemoCoherence.(behavField).(dataType).meanC = mean(data.NeuralHemoCoherence.(behavField).(dataType).C,2);
        data.NeuralHemoCoherence.(behavField).(dataType).stdC = std(data.NeuralHemoCoherence.(behavField).(dataType).C,0,2);
        data.NeuralHemoCoherence.(behavField).(dataType).meanf = mean(data.NeuralHemoCoherence.(behavField).(dataType).f,1);
        data.NeuralHemoCoherence.(behavField).(dataType).maxConfC = geomean(data.NeuralHemoCoherence.(behavField).(dataType).confC);
        data.NeuralHemoCoherence.(behavField).(dataType).maxConfC_Y = ones(length(data.NeuralHemoCoherence.(behavField).(dataType).meanf),1)*data.NeuralHemoCoherence.(behavField).(dataType).maxConfC;
    end
end
%% Fig. 8
summaryFigure = figure('Name','Fig8 (a,c)'); %#ok<*NASGU>
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
p2 = plot(edges,sgolayfilt(medfilt1(nremProbPerc,10,'truncate'),3,17),'-','color',colorNREM,'LineWidth',2);
p3 = plot(edges,sgolayfilt(medfilt1(remProbPerc,10,'truncate'),3,17),'-','color',colorREM,'LineWidth',2);
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
s1 = semilogx(data.NeuralHemoCoherence.Rest.gammaBandPower.meanf,data.NeuralHemoCoherence.Rest.gammaBandPower.meanC,'color',colorRest,'LineWidth',2);
hold on
s2 = semilogx(data.NeuralHemoCoherence.NREM.gammaBandPower.meanf,data.NeuralHemoCoherence.NREM.gammaBandPower.meanC,'color',colorNREM,'LineWidth',2);
s3 = semilogx(data.NeuralHemoCoherence.REM.gammaBandPower.meanf,data.NeuralHemoCoherence.REM.gammaBandPower.meanC,'color',colorREM,'LineWidth',2);
s4 = semilogx(data.NeuralHemoCoherence.Awake.gammaBandPower.meanf,data.NeuralHemoCoherence.Awake.gammaBandPower.meanC,'color',colorAwake,'LineWidth',2);
s5 = semilogx(data.NeuralHemoCoherence.Sleep.gammaBandPower.meanf,data.NeuralHemoCoherence.Sleep.gammaBandPower.meanC,'color',colorSleep,'LineWidth',2);
s6 = semilogx(data.NeuralHemoCoherence.All.gammaBandPower.meanf,data.NeuralHemoCoherence.All.gammaBandPower.meanC,'color',colorAll,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[8c] Neural-hemo coherence','Gamma-band power and \DeltaHbT \muM (%)',''})
legend([s1,s2,s3,s4,s5,s6],'Rest','NREM','REM','Awake','Sleep','All','Location','SouthEast')
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder '\Summary Figures and Structures\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'Fig8']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-fillpage',[dirpath 'Fig8'])
end

end
