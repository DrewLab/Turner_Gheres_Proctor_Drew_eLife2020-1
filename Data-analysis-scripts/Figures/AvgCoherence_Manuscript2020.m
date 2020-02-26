function [] = AvgAwakeProbability_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Calculate the average coherence of different behavioral states
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120'};
bins = {'Rest','NREM','REM'};
modelType = 'SVM';
coherr_dataTypes = {'CBV_HbT','deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        if strcmp(behavField,'Rest') == true
            for c = 1:length(coherr_dataTypes)
                coherr_dataType = coherr_dataTypes{1,c};
                data.(behavField).(coherr_dataType).C(:,a) = (AnalysisResults.(animalID).Coherence.(behavField).(coherr_dataType).C).^2;
                data.(behavField).(coherr_dataType).f(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(coherr_dataType).f;
                data.(behavField).(coherr_dataType).confC(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(coherr_dataType).confC;
            end
        elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            for d = 1:length(coherr_dataTypes)
                coherr_dataType = coherr_dataTypes{1,d};
                data.(behavField).(modelType).(coherr_dataType).C(:,a) = (AnalysisResults.(animalID).Coherence.(behavField).(modelType).(coherr_dataType).C).^2;
                data.(behavField).(modelType).(coherr_dataType).f(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(modelType).(coherr_dataType).f;
                data.(behavField).(modelType).(coherr_dataType).confC(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(modelType).(coherr_dataType).confC;
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
            data.(behavField).(coherr_dataType).meanC = mean(data.(behavField).(coherr_dataType).C,2);
            data.(behavField).(coherr_dataType).stdC = std(data.(behavField).(coherr_dataType).C,0,2);
            data.(behavField).(coherr_dataType).meanf = mean(data.(behavField).(coherr_dataType).f,2);
            data.(behavField).(coherr_dataType).maxConfC = max(data.(behavField).(coherr_dataType).confC);
            data.(behavField).(coherr_dataType).maxConfC_Y = ones(length(data.(behavField).(coherr_dataType).meanf),1)*data.(behavField).(coherr_dataType).maxConfC;
        end
    elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
        for h = 1:length(coherr_dataTypes)
            coherr_dataType = coherr_dataTypes{1,h};
            data.(behavField).(modelType).(coherr_dataType).meanC = mean(data.(behavField).(modelType).(coherr_dataType).C,2);
            data.(behavField).(modelType).(coherr_dataType).stdC = std(data.(behavField).(modelType).(coherr_dataType).C,0,2);
            data.(behavField).(modelType).(coherr_dataType).meanf = mean(data.(behavField).(modelType).(coherr_dataType).f,2);
            data.(behavField).(modelType).(coherr_dataType).maxConfC = max(data.(behavField).(modelType).(coherr_dataType).confC);
            data.(behavField).(modelType).(coherr_dataType).maxConfC_Y = ones(length(data.(behavField).(modelType).(coherr_dataType).meanf),1)*data.(behavField).(modelType).(coherr_dataType).maxConfC;
        end
    end
end

%% summary figure(s)
summaryFigure = figure;
sgtitle('Spectral Coherence Between Bilateral Signals')
%% CBV_HbT
subplot(2,3,1);
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.meanC + data.Rest.CBV_HbT.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.meanC - data.Rest.CBV_HbT.stdC,'color',colorA,'LineWidth',1)
semilogx(data.NREM.(modelType).CBV_HbT.meanf,data.NREM.(modelType).CBV_HbT.meanC,'color',colorB,'LineWidth',3);
semilogx(data.NREM.(modelType).CBV_HbT.meanf,data.NREM.(modelType).CBV_HbT.meanC + data.NREM.(modelType).CBV_HbT.stdC,'color',colorB,'LineWidth',1)
semilogx(data.NREM.(modelType).CBV_HbT.meanf,data.NREM.(modelType).CBV_HbT.meanC - data.NREM.(modelType).CBV_HbT.stdC,'color',colorB,'LineWidth',1)
semilogx(data.REM.(modelType).CBV_HbT.meanf,data.REM.(modelType).CBV_HbT.meanC,'color',colorC,'LineWidth',3);
semilogx(data.REM.(modelType).CBV_HbT.meanf,data.REM.(modelType).CBV_HbT.meanC + data.REM.(modelType).CBV_HbT.stdC,'color',colorC,'LineWidth',1)
semilogx(data.REM.(modelType).CBV_HbT.meanf,data.REM.(modelType).CBV_HbT.meanC - data.REM.(modelType).CBV_HbT.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
semilogx(data.NREM.(modelType).CBV_HbT.meanf,data.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
semilogx(data.REM.(modelType).CBV_HbT.meanf,data.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('\DeltaHbT (\muM)')
ylabel('Coherence^2')
xlabel('Frequency (Hz)')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

%% Delta-band power
subplot(2,3,2);
semilogx(data.Rest.deltaBandPower.meanf,data.Rest.deltaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Rest.deltaBandPower.meanf,data.Rest.deltaBandPower.meanC + data.Rest.deltaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Rest.deltaBandPower.meanf,data.Rest.deltaBandPower.meanC - data.Rest.deltaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.NREM.(modelType).deltaBandPower.meanf,data.NREM.(modelType).deltaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.NREM.(modelType).deltaBandPower.meanf,data.NREM.(modelType).deltaBandPower.meanC + data.NREM.(modelType).deltaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.NREM.(modelType).deltaBandPower.meanf,data.NREM.(modelType).deltaBandPower.meanC - data.NREM.(modelType).deltaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.REM.(modelType).deltaBandPower.meanf,data.REM.(modelType).deltaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.REM.(modelType).deltaBandPower.meanf,data.REM.(modelType).deltaBandPower.meanC + data.REM.(modelType).deltaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.REM.(modelType).deltaBandPower.meanf,data.REM.(modelType).deltaBandPower.meanC - data.REM.(modelType).deltaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
semilogx(data.NREM.(modelType).CBV_HbT.meanf,data.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
semilogx(data.REM.(modelType).CBV_HbT.meanf,data.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('Delta-band [1-4 Hz]')
ylabel('Coherence^2')
xlabel('Frequency (Hz)')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

%% Theta-band power
subplot(2,3,3)
L1 = semilogx(data.Rest.thetaBandPower.meanf,data.Rest.thetaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Rest.thetaBandPower.meanf,data.Rest.thetaBandPower.meanC + data.Rest.thetaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Rest.thetaBandPower.meanf,data.Rest.thetaBandPower.meanC - data.Rest.thetaBandPower.stdC,'color',colorA,'LineWidth',1)
L2 = semilogx(data.NREM.(modelType).thetaBandPower.meanf,data.NREM.(modelType).thetaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.NREM.(modelType).thetaBandPower.meanf,data.NREM.(modelType).thetaBandPower.meanC + data.NREM.(modelType).thetaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.NREM.(modelType).thetaBandPower.meanf,data.NREM.(modelType).thetaBandPower.meanC - data.NREM.(modelType).thetaBandPower.stdC,'color',colorB,'LineWidth',1)
L3 = semilogx(data.REM.(modelType).thetaBandPower.meanf,data.REM.(modelType).thetaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.REM.(modelType).thetaBandPower.meanf,data.REM.(modelType).thetaBandPower.meanC + data.REM.(modelType).thetaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.REM.(modelType).thetaBandPower.meanf,data.REM.(modelType).thetaBandPower.meanC - data.REM.(modelType).thetaBandPower.stdC,'color',colorC,'LineWidth',1)
L4 = semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
L5 = semilogx(data.NREM.(modelType).CBV_HbT.meanf,data.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
L6 = semilogx(data.REM.(modelType).CBV_HbT.meanf,data.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('Theta-band [4-10 Hz]')
ylabel('Coherence^2')
xlabel('Frequency (Hz)')
legend([L1,L2,L3,L4,L5,L6],'Rest','NREM','REM','Rest 95% conf','NREM 95% conf','REM 95% conf')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

%% Alpha-band power
subplot(2,3,4);
semilogx(data.Rest.alphaBandPower.meanf,data.Rest.alphaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Rest.alphaBandPower.meanf,data.Rest.alphaBandPower.meanC + data.Rest.alphaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Rest.alphaBandPower.meanf,data.Rest.alphaBandPower.meanC - data.Rest.alphaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.NREM.(modelType).alphaBandPower.meanf,data.NREM.(modelType).alphaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.NREM.(modelType).alphaBandPower.meanf,data.NREM.(modelType).alphaBandPower.meanC + data.NREM.(modelType).alphaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.NREM.(modelType).alphaBandPower.meanf,data.NREM.(modelType).alphaBandPower.meanC - data.NREM.(modelType).alphaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.REM.(modelType).alphaBandPower.meanf,data.REM.(modelType).alphaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.REM.(modelType).alphaBandPower.meanf,data.REM.(modelType).alphaBandPower.meanC + data.REM.(modelType).alphaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.REM.(modelType).alphaBandPower.meanf,data.REM.(modelType).alphaBandPower.meanC - data.REM.(modelType).alphaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
semilogx(data.NREM.(modelType).CBV_HbT.meanf,data.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
semilogx(data.REM.(modelType).CBV_HbT.meanf,data.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('Alpha-band [10-13 Hz]')
ylabel('Coherence^2')
xlabel('Frequency (Hz)')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

%% Beta-band power
subplot(2,3,5);
semilogx(data.Rest.betaBandPower.meanf,data.Rest.betaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Rest.betaBandPower.meanf,data.Rest.betaBandPower.meanC + data.Rest.betaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Rest.betaBandPower.meanf,data.Rest.betaBandPower.meanC - data.Rest.betaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.NREM.(modelType).betaBandPower.meanf,data.NREM.(modelType).betaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.NREM.(modelType).betaBandPower.meanf,data.NREM.(modelType).betaBandPower.meanC + data.NREM.(modelType).betaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.NREM.(modelType).betaBandPower.meanf,data.NREM.(modelType).betaBandPower.meanC - data.NREM.(modelType).betaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.REM.(modelType).betaBandPower.meanf,data.REM.(modelType).betaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.REM.(modelType).betaBandPower.meanf,data.REM.(modelType).betaBandPower.meanC + data.REM.(modelType).betaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.REM.(modelType).betaBandPower.meanf,data.REM.(modelType).betaBandPower.meanC - data.REM.(modelType).betaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
semilogx(data.NREM.(modelType).CBV_HbT.meanf,data.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
semilogx(data.REM.(modelType).CBV_HbT.meanf,data.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('Beta-band [13-30 Hz]')
ylabel('Coherence^2')
xlabel('Frequency (Hz)')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

%% Gamma-band power
subplot(2,3,6);
semilogx(data.Rest.gammaBandPower.meanf,data.Rest.gammaBandPower.meanC,'color',colorA,'LineWidth',3);
hold on
semilogx(data.Rest.gammaBandPower.meanf,data.Rest.gammaBandPower.meanC + data.Rest.gammaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.Rest.gammaBandPower.meanf,data.Rest.gammaBandPower.meanC - data.Rest.gammaBandPower.stdC,'color',colorA,'LineWidth',1)
semilogx(data.NREM.(modelType).gammaBandPower.meanf,data.NREM.(modelType).gammaBandPower.meanC,'color',colorB,'LineWidth',3);
semilogx(data.NREM.(modelType).gammaBandPower.meanf,data.NREM.(modelType).gammaBandPower.meanC + data.NREM.(modelType).gammaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.NREM.(modelType).gammaBandPower.meanf,data.NREM.(modelType).gammaBandPower.meanC - data.NREM.(modelType).gammaBandPower.stdC,'color',colorB,'LineWidth',1)
semilogx(data.REM.(modelType).gammaBandPower.meanf,data.REM.(modelType).gammaBandPower.meanC,'color',colorC,'LineWidth',3);
semilogx(data.REM.(modelType).gammaBandPower.meanf,data.REM.(modelType).gammaBandPower.meanC + data.REM.(modelType).gammaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.REM.(modelType).gammaBandPower.meanf,data.REM.(modelType).gammaBandPower.meanC - data.REM.(modelType).gammaBandPower.stdC,'color',colorC,'LineWidth',1)
semilogx(data.Rest.CBV_HbT.meanf,data.Rest.CBV_HbT.maxConfC_Y,'--','color',colorA,'LineWidth',1);
semilogx(data.NREM.(modelType).CBV_HbT.meanf,data.NREM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorB,'LineWidth',1);
semilogx(data.REM.(modelType).CBV_HbT.meanf,data.REM.(modelType).CBV_HbT.maxConfC_Y,'--','color',colorC,'LineWidth',1);
title('Gamma-band [30-100 Hz]')
ylabel('Coherence^2')
xlabel('Frequency (Hz)')
axis square
ylim([0,1])
xlim([0.1,0.5])
set(gca,'box','off')

% save figure(s)
dirpath = [rootFolder '\Analysis Figures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Coherence']);

end

