function [] = FigurePanelSeven_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

coherrAnimalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Rest','NREM','REM'};
modelType = 'Forest';
coherr_dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
%% Average coherence during different behaviors 
% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(coherrAnimalIDs)
    animalID = coherrAnimalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        if strcmp(behavField,'Rest') == true
            for c = 1:length(coherr_dataTypes)
                coherr_dataType = coherr_dataTypes{1,c};
                data.Coherr.(behavField).(coherr_dataType).C(:,a) = (AnalysisResults.(animalID).Coherence.(behavField).(coherr_dataType).C);
                data.Coherr.(behavField).(coherr_dataType).f(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(coherr_dataType).f;
                data.Coherr.(behavField).(coherr_dataType).confC(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(coherr_dataType).confC;
            end
        elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            for d = 1:length(coherr_dataTypes)
                coherr_dataType = coherr_dataTypes{1,d};
                data.Coherr.(behavField).(modelType).(coherr_dataType).C(:,a) = (AnalysisResults.(animalID).Coherence.(behavField).(modelType).(coherr_dataType).C);
                data.Coherr.(behavField).(modelType).(coherr_dataType).f(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(modelType).(coherr_dataType).f;
                data.Coherr.(behavField).(modelType).(coherr_dataType).confC(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(modelType).(coherr_dataType).confC;
            end
        end
    end
end
% take the mean and standard deviation of each set of signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    if strcmp(behavField,'Rest') == true
        for f = 1:length(coherr_dataTypes)
            coherr_dataType = coherr_dataTypes{1,f};
            data.Coherr.(behavField).(coherr_dataType).meanC = mean(data.Coherr.(behavField).(coherr_dataType).C,2);
            data.Coherr.(behavField).(coherr_dataType).stdC = std(data.Coherr.(behavField).(coherr_dataType).C,0,2);
            data.Coherr.(behavField).(coherr_dataType).meanf = mean(data.Coherr.(behavField).(coherr_dataType).f,2);
            data.Coherr.(behavField).(coherr_dataType).maxConfC = max(data.Coherr.(behavField).(coherr_dataType).confC);
            data.Coherr.(behavField).(coherr_dataType).maxConfC_Y = ones(length(data.Coherr.(behavField).(coherr_dataType).meanf),1)*data.Coherr.(behavField).(coherr_dataType).maxConfC;
        end
    elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
        for h = 1:length(coherr_dataTypes)
            coherr_dataType = coherr_dataTypes{1,h};
            data.Coherr.(behavField).(modelType).(coherr_dataType).meanC = mean(data.Coherr.(behavField).(modelType).(coherr_dataType).C,2);
            data.Coherr.(behavField).(modelType).(coherr_dataType).stdC = std(data.Coherr.(behavField).(modelType).(coherr_dataType).C,0,2);
            data.Coherr.(behavField).(modelType).(coherr_dataType).meanf = mean(data.Coherr.(behavField).(modelType).(coherr_dataType).f,2);
            data.Coherr.(behavField).(modelType).(coherr_dataType).maxConfC = max(data.Coherr.(behavField).(modelType).(coherr_dataType).confC);
            data.Coherr.(behavField).(modelType).(coherr_dataType).maxConfC_Y = ones(length(data.Coherr.(behavField).(modelType).(coherr_dataType).meanf),1)*data.Coherr.(behavField).(modelType).(coherr_dataType).maxConfC;
        end
    end
end
%%

%% summary figure(s)
summaryFigure = figure;
sgtitle('Spectral Coherence Between Bilateral Signals')
%% CBV_HbT
subplot(2,3,1);
semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.meanC + data.Coherr.Rest.CBV_HbT.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.meanC - data.Coherr.Rest.CBV_HbT.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).CBV_HbT.meanf,data.Coherr.NREM.(modelType).CBV_HbT.meanC,'color',colorB,'LineWidth',3);
semilogx(data.Coherr.NREM.(modelType).CBV_HbT.meanf,data.Coherr.NREM.(modelType).CBV_HbT.meanC + data.Coherr.NREM.(modelType).CBV_HbT.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).CBV_HbT.meanf,data.Coherr.NREM.(modelType).CBV_HbT.meanC - data.Coherr.NREM.(modelType).CBV_HbT.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).CBV_HbT.meanf,data.Coherr.REM.(modelType).CBV_HbT.meanC,'color',colorC,'LineWidth',3);
semilogx(data.Coherr.REM.(modelType).CBV_HbT.meanf,data.Coherr.REM.(modelType).CBV_HbT.meanC + data.Coherr.REM.(modelType).CBV_HbT.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).CBV_HbT.meanf,data.Coherr.REM.(modelType).CBV_HbT.meanC - data.Coherr.REM.(modelType).CBV_HbT.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.(modelType).CBV_HbT.meanf,data.Coherr.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.(modelType).CBV_HbT.meanf,data.Coherr.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('\DeltaHbT (\muM)')
ylabel('Coherence')
xlabel('Freq (Hz)')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

%% Delta-band power
subplot(2,3,2);
semilogx(data.Coherr.Rest.deltaBandPower.meanf,data.Coherr.Rest.deltaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Coherr.Rest.deltaBandPower.meanf,data.Coherr.Rest.deltaBandPower.meanC + data.Coherr.Rest.deltaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.Rest.deltaBandPower.meanf,data.Coherr.Rest.deltaBandPower.meanC - data.Coherr.Rest.deltaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).deltaBandPower.meanf,data.Coherr.NREM.(modelType).deltaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.Coherr.NREM.(modelType).deltaBandPower.meanf,data.Coherr.NREM.(modelType).deltaBandPower.meanC + data.Coherr.NREM.(modelType).deltaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).deltaBandPower.meanf,data.Coherr.NREM.(modelType).deltaBandPower.meanC - data.Coherr.NREM.(modelType).deltaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).deltaBandPower.meanf,data.Coherr.REM.(modelType).deltaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.Coherr.REM.(modelType).deltaBandPower.meanf,data.Coherr.REM.(modelType).deltaBandPower.meanC + data.Coherr.REM.(modelType).deltaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).deltaBandPower.meanf,data.Coherr.REM.(modelType).deltaBandPower.meanC - data.Coherr.REM.(modelType).deltaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.(modelType).CBV_HbT.meanf,data.Coherr.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.(modelType).CBV_HbT.meanf,data.Coherr.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('Delta-band [1-4 Hz]')
ylabel('Coherence')
xlabel('Freq (Hz)')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

%% Theta-band power
subplot(2,3,3)
L1 = semilogx(data.Coherr.Rest.thetaBandPower.meanf,data.Coherr.Rest.thetaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Coherr.Rest.thetaBandPower.meanf,data.Coherr.Rest.thetaBandPower.meanC + data.Coherr.Rest.thetaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.Rest.thetaBandPower.meanf,data.Coherr.Rest.thetaBandPower.meanC - data.Coherr.Rest.thetaBandPower.stdC,'color',colorA,'LineWidth',1)
L2 = semilogx(data.Coherr.NREM.(modelType).thetaBandPower.meanf,data.Coherr.NREM.(modelType).thetaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.Coherr.NREM.(modelType).thetaBandPower.meanf,data.Coherr.NREM.(modelType).thetaBandPower.meanC + data.Coherr.NREM.(modelType).thetaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).thetaBandPower.meanf,data.Coherr.NREM.(modelType).thetaBandPower.meanC - data.Coherr.NREM.(modelType).thetaBandPower.stdC,'color',colorB,'LineWidth',1)
L3 = semilogx(data.Coherr.REM.(modelType).thetaBandPower.meanf,data.Coherr.REM.(modelType).thetaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.Coherr.REM.(modelType).thetaBandPower.meanf,data.Coherr.REM.(modelType).thetaBandPower.meanC + data.Coherr.REM.(modelType).thetaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).thetaBandPower.meanf,data.Coherr.REM.(modelType).thetaBandPower.meanC - data.Coherr.REM.(modelType).thetaBandPower.stdC,'color',colorC,'LineWidth',1)
L4 = semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
L5 = semilogx(data.Coherr.NREM.(modelType).CBV_HbT.meanf,data.Coherr.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
L6 = semilogx(data.Coherr.REM.(modelType).CBV_HbT.meanf,data.Coherr.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('Theta-band [4-10 Hz]')
ylabel('Coherence')
xlabel('Freq (Hz)')
legend([L1,L2,L3,L4,L5,L6],'Rest','NREM','REM','Rest 95% conf','NREM 95% conf','REM 95% conf')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

%% Alpha-band power
subplot(2,3,4);
semilogx(data.Coherr.Rest.alphaBandPower.meanf,data.Coherr.Rest.alphaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Coherr.Rest.alphaBandPower.meanf,data.Coherr.Rest.alphaBandPower.meanC + data.Coherr.Rest.alphaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.Rest.alphaBandPower.meanf,data.Coherr.Rest.alphaBandPower.meanC - data.Coherr.Rest.alphaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).alphaBandPower.meanf,data.Coherr.NREM.(modelType).alphaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.Coherr.NREM.(modelType).alphaBandPower.meanf,data.Coherr.NREM.(modelType).alphaBandPower.meanC + data.Coherr.NREM.(modelType).alphaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).alphaBandPower.meanf,data.Coherr.NREM.(modelType).alphaBandPower.meanC - data.Coherr.NREM.(modelType).alphaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).alphaBandPower.meanf,data.Coherr.REM.(modelType).alphaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.Coherr.REM.(modelType).alphaBandPower.meanf,data.Coherr.REM.(modelType).alphaBandPower.meanC + data.Coherr.REM.(modelType).alphaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).alphaBandPower.meanf,data.Coherr.REM.(modelType).alphaBandPower.meanC - data.Coherr.REM.(modelType).alphaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.(modelType).CBV_HbT.meanf,data.Coherr.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.(modelType).CBV_HbT.meanf,data.Coherr.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('Alpha-band [10-13 Hz]')
ylabel('Coherence')
xlabel('Freq (Hz)')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

%% Beta-band power
subplot(2,3,5);
semilogx(data.Coherr.Rest.betaBandPower.meanf,data.Coherr.Rest.betaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Coherr.Rest.betaBandPower.meanf,data.Coherr.Rest.betaBandPower.meanC + data.Coherr.Rest.betaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.Rest.betaBandPower.meanf,data.Coherr.Rest.betaBandPower.meanC - data.Coherr.Rest.betaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).betaBandPower.meanf,data.Coherr.NREM.(modelType).betaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.Coherr.NREM.(modelType).betaBandPower.meanf,data.Coherr.NREM.(modelType).betaBandPower.meanC + data.Coherr.NREM.(modelType).betaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).betaBandPower.meanf,data.Coherr.NREM.(modelType).betaBandPower.meanC - data.Coherr.NREM.(modelType).betaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).betaBandPower.meanf,data.Coherr.REM.(modelType).betaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.Coherr.REM.(modelType).betaBandPower.meanf,data.Coherr.REM.(modelType).betaBandPower.meanC + data.Coherr.REM.(modelType).betaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).betaBandPower.meanf,data.Coherr.REM.(modelType).betaBandPower.meanC - data.Coherr.REM.(modelType).betaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.(modelType).CBV_HbT.meanf,data.Coherr.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.(modelType).CBV_HbT.meanf,data.Coherr.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('Beta-band [13-30 Hz]')
ylabel('Coherence')
xlabel('Freq (Hz)')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

%% Gamma-band power
subplot(2,3,6);
semilogx(data.Coherr.Rest.gammaBandPower.meanf,data.Coherr.Rest.gammaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Coherr.Rest.gammaBandPower.meanf,data.Coherr.Rest.gammaBandPower.meanC + data.Coherr.Rest.gammaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.Rest.gammaBandPower.meanf,data.Coherr.Rest.gammaBandPower.meanC - data.Coherr.Rest.gammaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).gammaBandPower.meanf,data.Coherr.NREM.(modelType).gammaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.Coherr.NREM.(modelType).gammaBandPower.meanf,data.Coherr.NREM.(modelType).gammaBandPower.meanC + data.Coherr.NREM.(modelType).gammaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.NREM.(modelType).gammaBandPower.meanf,data.Coherr.NREM.(modelType).gammaBandPower.meanC - data.Coherr.NREM.(modelType).gammaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).gammaBandPower.meanf,data.Coherr.REM.(modelType).gammaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.Coherr.REM.(modelType).gammaBandPower.meanf,data.Coherr.REM.(modelType).gammaBandPower.meanC + data.Coherr.REM.(modelType).gammaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.REM.(modelType).gammaBandPower.meanf,data.Coherr.REM.(modelType).gammaBandPower.meanC - data.Coherr.REM.(modelType).gammaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.(modelType).CBV_HbT.meanf,data.Coherr.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.(modelType).CBV_HbT.meanf,data.Coherr.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('Gamma-band [30-100 Hz]')
ylabel('Coherence')
xlabel('Freq (Hz)')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Coherence']);

end

