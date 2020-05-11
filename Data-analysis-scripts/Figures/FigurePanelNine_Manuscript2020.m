function [] = FigurePanelNine_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose:
%________________________________________________________________________________________________________________________


animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','Awake','NREM','REM','All'};
modelType = 'Forest';
dataTypes = {'gammaBandPower'};
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
colorF = [(197/256),(179/256),(90/256)];  % Awake color
%% Average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
z = 1;
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'All') == true
            for c = 1:length(dataTypes)
                dataType = dataTypes{1,c};
                data.Coherr.(behavField).(dataType).LH_C(:,a) = (AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C);
                data.Coherr.(behavField).(dataType).RH_C(:,a) = (AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.C);
                data.Coherr.(behavField).(dataType).f(:,a) = AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.f;
                data.Coherr.(behavField).(dataType).confC(:,a) = AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.confC;
            end
        elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            for d = 1:length(dataTypes)
                dataType = dataTypes{1,d};
                data.Coherr.(behavField).(modelType).(dataType).LH_C(:,a) = (AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(modelType).(dataType).adjLH.C);
                data.Coherr.(behavField).(modelType).(dataType).RH_C(:,a) = (AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(modelType).(dataType).adjRH.C);
                data.Coherr.(behavField).(modelType).(dataType).f(:,a) = AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(modelType).(dataType).adjLH.f;
                data.Coherr.(behavField).(modelType).(dataType).confC(:,a) = AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(modelType).(dataType).adjLH.confC;
            end
        elseif strcmp(behavField,'Awake') == true
            if isfield(AnalysisResults.(animalID).Coherence,'Awake') == true
                for e = 1:length(dataTypes)
                    dataType = dataTypes{1,e};
                    data.Coherr.(behavField).(dataType).LH_C(:,z) = (AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.C);
                    data.Coherr.(behavField).(dataType).RH_C(:,z) = (AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjRH.C);
                    data.Coherr.(behavField).(dataType).f(:,z) = AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.f;
                    data.Coherr.(behavField).(dataType).confC(:,z) = AnalysisResults.(animalID).NeuralHemoCoherence.(behavField).(dataType).adjLH.confC;
                end
                z = z + 1;
            end
        end
    end
end
% concatenate left and right signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    if strcmp(behavField,'Rest') == true || strcmp(behavField,'Awake') == true || strcmp(behavField,'All') == true
        for f = 1:length(dataTypes)
            dataType = dataTypes{1,f};
            data.Coherr.(behavField).(dataType).catC = cat(2,data.Coherr.(behavField).(dataType).LH_C,data.Coherr.(behavField).(dataType).RH_C);
        end
    elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
        for h = 1:length(dataTypes)
            dataType = dataTypes{1,h};
            data.Coherr.(behavField).(modelType).(dataType).catC = cat(2,data.Coherr.(behavField).(modelType).(dataType).LH_C,data.Coherr.(behavField).(modelType).(dataType).RH_C);
        end
    end
end
% take the mean and standard deviation of each set of signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    if strcmp(behavField,'Rest') == true || strcmp(behavField,'Awake') == true || strcmp(behavField,'All') == true
        for f = 1:length(dataTypes)
            dataType = dataTypes{1,f};
            data.Coherr.(behavField).(dataType).meanC = mean(data.Coherr.(behavField).(dataType).catC,2);
            data.Coherr.(behavField).(dataType).stdC = std(data.Coherr.(behavField).(dataType).catC,0,2);
            data.Coherr.(behavField).(dataType).meanf = mean(data.Coherr.(behavField).(dataType).f,2);
            data.Coherr.(behavField).(dataType).maxConfC = max(data.Coherr.(behavField).(dataType).confC);
            data.Coherr.(behavField).(dataType).maxConfC_Y = ones(length(data.Coherr.(behavField).(dataType).meanf),1)*data.Coherr.(behavField).(dataType).maxConfC;
        end
    elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
        for h = 1:length(dataTypes)
            dataType = dataTypes{1,h};
            data.Coherr.(behavField).(modelType).(dataType).meanC = mean(data.Coherr.(behavField).(modelType).(dataType).catC,2);
            data.Coherr.(behavField).(modelType).(dataType).stdC = std(data.Coherr.(behavField).(modelType).(dataType).catC,0,2);
            data.Coherr.(behavField).(modelType).(dataType).meanf = mean(data.Coherr.(behavField).(modelType).(dataType).f,2);
            data.Coherr.(behavField).(modelType).(dataType).maxConfC = max(data.Coherr.(behavField).(modelType).(dataType).confC);
            data.Coherr.(behavField).(modelType).(dataType).maxConfC_Y = ones(length(data.Coherr.(behavField).(modelType).(dataType).meanf),1)*data.Coherr.(behavField).(modelType).(dataType).maxConfC;
        end
    end
end
%% Figure Panel 9
summaryFigure = figure;
sgtitle('Figure Panel 9 - Turner Manuscript 2020')
%% [A] Coherence between HbT and gamma-band power during different arousal-states
ax1 = subplot(1,1,1);
s1 = semilogx(data.Coherr.Awake.gammaBandPower.meanf,data.Coherr.Awake.gammaBandPower.meanC,'color',colorF,'LineWidth',2);
hold on
s2 = semilogx(data.Coherr.Rest.gammaBandPower.meanf,data.Coherr.Rest.gammaBandPower.meanC,'color',colorA,'LineWidth',2);
s3 = semilogx(data.Coherr.NREM.(modelType).gammaBandPower.meanf,data.Coherr.NREM.(modelType).gammaBandPower.meanC,'color',colorB,'LineWidth',2);
s4 = semilogx(data.Coherr.REM.(modelType).gammaBandPower.meanf,data.Coherr.REM.(modelType).gammaBandPower.meanC,'color',colorC,'LineWidth',2);
s5 = semilogx(data.Coherr.All.gammaBandPower.meanf,data.Coherr.All.gammaBandPower.meanC,'color','r','LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('\DeltaHbT (\muM)')
ylabel('Coherence')
xlabel('Freq (Hz)')
title({'[A] Neural-hemo coherence','Gamma-band power and \DeltaHbT \muM (%)',''})
legend([s1,s2,s3,s4,s5],'Awake','Rest','NREM','REM','All','Location','SouthEast')
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Figure Panel 9']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Figure Panel 9'])

end
