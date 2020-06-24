function [] = Fig7_Manuscript2020_fin(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

IOS_animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
TwoP_animalIDs = {'T115','T116','T117','T118','T125','T126'};
behavFields = {'Rest','NREM','REM','Awake','Sleep','All'};
behavFields2 = {'Rest','NREM','REM','Awake','All'};
behavFields3 = {'Rest','Whisk','NREM','REM','Awake','Sleep','All'};
dataTypes = {'CBV_HbT','gammaBandPower'};
colorA = [(51/256),(160/256),(44/256)];   % Rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(180/256)];  % Whisk color
% colorE = [(0/256),(256/256),(256/256)];  % Isoflurane color
colorF = [(256/256),(192/256),(0/256)];   % Awake color
colorG = [(0/256),(128/256),(256/256)];   % Sleep color
colorH = [(184/256),(115/256),(51/256)];  % All color
%% Average coherence during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.Coherr = [];
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
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
                    data.Coherr.(behavField).(dataType).C = [];
                    data.Coherr.(behavField).(dataType).f = [];
                    data.Coherr.(behavField).(dataType).confC = [];
                end
                % concatenate C/f for existing data - exclude any empty sets
                data.Coherr.(behavField).(dataType).C = cat(2,data.Coherr.(behavField).(dataType).C,AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
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
%% Power spectra during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            data.PowerSpec.(behavField).(dataType).adjLH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.S;
            data.PowerSpec.(behavField).(dataType).adjLH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.f;
            data.PowerSpec.(behavField).(dataType).adjRH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.S;
            data.PowerSpec.(behavField).(dataType).adjRH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.f;
        end
    end
end
% find the peak of the resting PSD for each animal/hemisphere
for a = 1:length(IOS_animalIDs)
    for c = 1:length(dataTypes)
        dataType = dataTypes{1,c};
        data.PowerSpec.baseline.(dataType).LH{a,1} = max(data.PowerSpec.Rest.(dataType).adjLH.S{a,1});
        data.PowerSpec.baseline.(dataType).RH{a,1} = max(data.PowerSpec.Rest.(dataType).adjRH.S{a,1});
    end
end
% DC-shift each animal/hemisphere/behavior PSD with respect to the resting peak
for a = 1:length(IOS_animalIDs)
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        for j = 1:length(dataTypes)
            dataType = dataTypes{1,j};
            for ee = 1:size(data.PowerSpec.(behavField).(dataType).adjLH.S,2)
                data.PowerSpec.(behavField).(dataType).normLH{a,1} = (data.PowerSpec.(behavField).(dataType).adjLH.S{a,1})*(1/(data.PowerSpec.baseline.(dataType).LH{a,1}));
                data.PowerSpec.(behavField).(dataType).normRH{a,1} = (data.PowerSpec.(behavField).(dataType).adjRH.S{a,1})*(1/(data.PowerSpec.baseline.(dataType).RH{a,1}));
            end
        end
    end
end
% concatenate the data from the left and right hemispheres - removes any empty data
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        data.PowerSpec.(behavField).(dataType).cat_S = [];
        data.PowerSpec.(behavField).(dataType).cat_f = [];
        for z = 1:length(data.PowerSpec.(behavField).(dataType).normLH)
            data.PowerSpec.(behavField).(dataType).cat_S = cat(2,data.PowerSpec.(behavField).(dataType).cat_S,data.PowerSpec.(behavField).(dataType).normLH{z,1},data.PowerSpec.(behavField).(dataType).normRH{z,1});
            data.PowerSpec.(behavField).(dataType).cat_f = cat(1,data.PowerSpec.(behavField).(dataType).cat_f,data.PowerSpec.(behavField).(dataType).adjLH.f{z,1},data.PowerSpec.(behavField).(dataType).adjRH.f{z,1});
        end
    end
end
% take mean/StD of S/f
for h = 1:length(behavFields)
    behavField = behavFields{1,h};
    for j = 1:length(dataTypes)
        dataType = dataTypes{1,j};
        data.PowerSpec.(behavField).(dataType).meanCortS = mean(data.PowerSpec.(behavField).(dataType).cat_S,2);
        data.PowerSpec.(behavField).(dataType).stdCortS = std(data.PowerSpec.(behavField).(dataType).cat_S,0,2);
        data.PowerSpec.(behavField).(dataType).meanCortf = mean(data.PowerSpec.(behavField).(dataType).cat_f,1);
    end
end
%% Vessel Power spectra of different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.VesselPowerSpec = [];
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    for bb = 1:length(behavFields2)
        behavField = behavFields2{1,bb};
        if isfield(AnalysisResults.(animalID).PowerSpectra,behavField) == true
            vesselIDs = fieldnames(AnalysisResults.(animalID).PowerSpectra.(behavField));
            for cc = 1:length(vesselIDs)
                vesselID = vesselIDs{cc,1};
                data.VesselPowerSpec.(animalID).(behavField).(vesselID).S = AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).S;
                data.VesselPowerSpec.(animalID).(behavField).(vesselID).f = AnalysisResults.(animalID).PowerSpectra.(behavField).(vesselID).f;
            end
        end
    end
end
% find the peak of the resting PSD for each arteriole
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    vesselIDs = fieldnames(data.VesselPowerSpec.(animalID).Rest);
    for cc = 1:length(vesselIDs)
        vesselID = vesselIDs{cc,1};
        data.VesselPowerSpec.(animalID).baseline.(vesselID) = max(data.VesselPowerSpec.(animalID).Rest.(vesselID).S);
    end
end
% DC-shift each arteriole's PSD with respect to the resting peak
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    for bb = 1:length(behavFields2)
        behavField = behavFields2{1,bb};
        if isfield(data.VesselPowerSpec.(animalID),behavField) == true
            vesselIDs = fieldnames(data.VesselPowerSpec.(animalID).(behavField));
            for cc = 1:length(vesselIDs)
                vesselID = vesselIDs{cc,1};
                data.VesselPowerSpec.(animalID).(behavField).(vesselID).normS = data.VesselPowerSpec.(animalID).(behavField).(vesselID).S*(1/data.VesselPowerSpec.(animalID).baseline.(vesselID));
            end
        end
    end
end
% concatenate the data from all arterioles - removes any empty data
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    for bb = 1:length(behavFields2)
        behavField = behavFields2{1,bb};
        if isfield(data.VesselPowerSpec,behavField) == false
            data.VesselPowerSpec.(behavField).S = [];
            data.VesselPowerSpec.(behavField).f = [];
        end
        if isfield(data.VesselPowerSpec.(animalID),behavField) == true
            vesselIDs = fieldnames(data.VesselPowerSpec.(animalID).(behavField));
            for cc = 1:length(vesselIDs)
                vesselID = vesselIDs{cc,1};
                data.VesselPowerSpec.(behavField).S = cat(2,data.VesselPowerSpec.(behavField).S,data.VesselPowerSpec.(animalID).(behavField).(vesselID).normS);
                data.VesselPowerSpec.(behavField).f = cat(1,data.VesselPowerSpec.(behavField).f,data.VesselPowerSpec.(animalID).(behavField).(vesselID).f);
            end
        end
    end
end
% take mean/StD of S/f
for dd = 1:length(behavFields2)
    behavField = behavFields2{1,dd};
    data.VesselPowerSpec.(behavField).meanS = mean(data.VesselPowerSpec.(behavField).S,2);
    data.VesselPowerSpec.(behavField).StDS = std(data.VesselPowerSpec.(behavField).S,0,2);
    data.VesselPowerSpec.(behavField).meanf = mean(data.VesselPowerSpec.(behavField).f,1);
end
%% Pearson's correlations during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
data.CorrCoef = [];
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    for b = 1:length(behavFields3)
        behavField = behavFields3{1,b};
        % create the behavior folder for the first iteration of the loop
        if isfield(data.CorrCoef,behavField) == false
            data.CorrCoef.(behavField) = [];
        end
        for c = 1:length(dataTypes)
            dataType = dataTypes{1,c};
            % don't concatenate empty arrays where there was no data for this behavior
            if isempty(AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR) == false
                % create the data type folder for the first iteration of the loop
                if isfield(data.CorrCoef.(behavField),dataType) == false
                    data.CorrCoef.(behavField).(dataType).meanRs = [];
                    data.CorrCoef.(behavField).(dataType).animalID = {};
                    data.CorrCoef.(behavField).(dataType).behavior = {};
                end
                % concatenate mean R and the animalID/behavior for statistics table
                data.CorrCoef.(behavField).(dataType).meanRs = cat(1,data.CorrCoef.(behavField).(dataType).meanRs,AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR);
                data.CorrCoef.(behavField).(dataType).animalID = cat(1,data.CorrCoef.(behavField).(dataType).animalID,animalID);
                data.CorrCoef.(behavField).(dataType).behavior = cat(1,data.CorrCoef.(behavField).(dataType).behavior,behavField);
            end
        end
    end
end
% take mean/STD of R
for e = 1:length(behavFields3)
    behavField = behavFields3{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        data.CorrCoef.(behavField).(dataType).meanR = mean(data.CorrCoef.(behavField).(dataType).meanRs,1);
        data.CorrCoef.(behavField).(dataType).stdR = std(data.CorrCoef.(behavField).(dataType).meanRs,0,1);
    end
end
%% statistics - linear mixed effects model
% HbT
HbTtableSize = cat(1,data.CorrCoef.Rest.CBV_HbT.meanRs,data.CorrCoef.Whisk.CBV_HbT.meanRs,data.CorrCoef.NREM.CBV_HbT.meanRs,data.CorrCoef.REM.CBV_HbT.meanRs,...
    data.CorrCoef.Awake.CBV_HbT.meanRs,data.CorrCoef.Sleep.CBV_HbT.meanRs,data.CorrCoef.All.CBV_HbT.meanRs);
HbTTable = table('Size',[size(HbTtableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
HbTTable.Mouse = cat(1,data.CorrCoef.Rest.CBV_HbT.animalID,data.CorrCoef.Whisk.CBV_HbT.animalID,data.CorrCoef.NREM.CBV_HbT.animalID,data.CorrCoef.REM.CBV_HbT.animalID,...
    data.CorrCoef.Awake.CBV_HbT.animalID,data.CorrCoef.Sleep.CBV_HbT.animalID,data.CorrCoef.All.CBV_HbT.animalID);
HbTTable.CorrCoef = cat(1,data.CorrCoef.Rest.CBV_HbT.meanRs,data.CorrCoef.Whisk.CBV_HbT.meanRs,data.CorrCoef.NREM.CBV_HbT.meanRs,data.CorrCoef.REM.CBV_HbT.meanRs,...
    data.CorrCoef.Awake.CBV_HbT.meanRs,data.CorrCoef.Sleep.CBV_HbT.meanRs,data.CorrCoef.All.CBV_HbT.meanRs);
HbTTable.Behavior = cat(1,data.CorrCoef.Rest.CBV_HbT.behavior,data.CorrCoef.Whisk.CBV_HbT.behavior,data.CorrCoef.NREM.CBV_HbT.behavior,data.CorrCoef.REM.CBV_HbT.behavior,...
    data.CorrCoef.Awake.CBV_HbT.behavior,data.CorrCoef.Sleep.CBV_HbT.behavior,data.CorrCoef.All.CBV_HbT.behavior);
HbTFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
HbTStats = fitglme(HbTTable,HbTFitFormula);
% gamma-band power
gammatableSize = cat(1,data.CorrCoef.Rest.gammaBandPower.meanRs,data.CorrCoef.Whisk.gammaBandPower.meanRs,data.CorrCoef.NREM.gammaBandPower.meanRs,data.CorrCoef.REM.gammaBandPower.meanRs,...
    data.CorrCoef.Awake.gammaBandPower.meanRs,data.CorrCoef.Sleep.gammaBandPower.meanRs,data.CorrCoef.All.gammaBandPower.meanRs);
gammaTable = table('Size',[size(gammatableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
gammaTable.Mouse = cat(1,data.CorrCoef.Rest.gammaBandPower.animalID,data.CorrCoef.Whisk.gammaBandPower.animalID,data.CorrCoef.NREM.gammaBandPower.animalID,data.CorrCoef.REM.gammaBandPower.animalID,...
    data.CorrCoef.Awake.gammaBandPower.animalID,data.CorrCoef.Sleep.gammaBandPower.animalID,data.CorrCoef.All.gammaBandPower.animalID);
gammaTable.CorrCoef = cat(1,data.CorrCoef.Rest.gammaBandPower.meanRs,data.CorrCoef.Whisk.gammaBandPower.meanRs,data.CorrCoef.NREM.gammaBandPower.meanRs,data.CorrCoef.REM.gammaBandPower.meanRs,...
    data.CorrCoef.Awake.gammaBandPower.meanRs,data.CorrCoef.Sleep.gammaBandPower.meanRs,data.CorrCoef.All.gammaBandPower.meanRs);
gammaTable.Behavior = cat(1,data.CorrCoef.Rest.gammaBandPower.behavior,data.CorrCoef.Whisk.gammaBandPower.behavior,data.CorrCoef.NREM.gammaBandPower.behavior,data.CorrCoef.REM.gammaBandPower.behavior,...
    data.CorrCoef.Awake.gammaBandPower.behavior,data.CorrCoef.Sleep.gammaBandPower.behavior,data.CorrCoef.All.gammaBandPower.behavior);
gammaFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
gammaStats = fitglme(gammaTable,gammaFitFormula);
%% Figure Panel 7
summaryFigure = figure('Name','Fig7 (a-g)');
sgtitle('Figure Panel 7 (a-g) Turner Manuscript 2020')
CC_xInds = ones(1,length(IOS_animalIDs));
CC_xInds2 = ones(1,length(data.CorrCoef.Awake.CBV_HbT.animalID));
CC_xInds3 = ones(1,length(data.CorrCoef.Sleep.CBV_HbT.animalID));
%% [7a] Power spectra of gamma-band power during different arousal-states
ax1 = subplot(3,3,1);
L1 = loglog(data.PowerSpec.Rest.gammaBandPower.meanCortf,data.PowerSpec.Rest.gammaBandPower.meanCortS,'color',colorA,'LineWidth',2);
hold on
L2 = loglog(data.PowerSpec.NREM.gammaBandPower.meanCortf,data.PowerSpec.NREM.gammaBandPower.meanCortS,'color',colorB,'LineWidth',2);
L3 = loglog(data.PowerSpec.REM.gammaBandPower.meanCortf,data.PowerSpec.REM.gammaBandPower.meanCortS,'color',colorC,'LineWidth',2);
L4 = loglog(data.PowerSpec.Awake.gammaBandPower.meanCortf,data.PowerSpec.Awake.gammaBandPower.meanCortS,'color',colorF,'LineWidth',2);
L5 = loglog(data.PowerSpec.Sleep.gammaBandPower.meanCortf,data.PowerSpec.Sleep.gammaBandPower.meanCortS,'color',colorG,'LineWidth',2);
L6 = loglog(data.PowerSpec.All.gammaBandPower.meanCortf,data.PowerSpec.All.gammaBandPower.meanCortS,'color',colorH,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title({'[7a] Cortical power','Gamma-band [30-100 Hz]',''})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
legend([L1,L2,L3,L4,L5,L6],'Rest','NREM','REM','Awake','Sleep','All','Location','SouthEast')
axis square
xlim([0.003,0.5])
y5 = ylim(ax1);
ylim([y5(1)/2,y5(2)*2])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [7b] Coherence between bilateral gamma-band power during different arousal-states
ax2 = subplot(3,3,2);
semilogx(data.Coherr.Rest.gammaBandPower.meanf,data.Coherr.Rest.gammaBandPower.meanC.^2,'color',colorA,'LineWidth',2);
hold on
semilogx(data.Coherr.NREM.gammaBandPower.meanf,data.Coherr.NREM.gammaBandPower.meanC.^2,'color',colorB,'LineWidth',2);
semilogx(data.Coherr.REM.gammaBandPower.meanf,data.Coherr.REM.gammaBandPower.meanC.^2,'color',colorC,'LineWidth',2);
semilogx(data.Coherr.Awake.gammaBandPower.meanf,data.Coherr.Awake.gammaBandPower.meanC.^2,'color',colorF,'LineWidth',2);
semilogx(data.Coherr.Sleep.gammaBandPower.meanf,data.Coherr.Sleep.gammaBandPower.meanC.^2,'color',colorG,'LineWidth',2);
semilogx(data.Coherr.All.gammaBandPower.meanf,data.Coherr.All.gammaBandPower.meanC.^2,'color',colorH,'LineWidth',2);
% confidence lines
semilogx(data.Coherr.Rest.gammaBandPower.meanf,data.Coherr.Rest.gammaBandPower.maxConfC_Y.^2,'-','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.gammaBandPower.meanf,data.Coherr.NREM.gammaBandPower.maxConfC_Y.^2,'-','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.gammaBandPower.meanf,data.Coherr.REM.gammaBandPower.maxConfC_Y.^2,'-','color',colorC,'LineWidth',1);
semilogx(data.Coherr.Awake.gammaBandPower.meanf,data.Coherr.Awake.gammaBandPower.maxConfC_Y.^2,'-','color',colorF,'LineWidth',1);
semilogx(data.Coherr.Sleep.gammaBandPower.meanf,data.Coherr.Sleep.gammaBandPower.maxConfC_Y.^2,'-','color',colorG,'LineWidth',1);
semilogx(data.Coherr.All.gammaBandPower.meanf,data.Coherr.All.gammaBandPower.maxConfC_Y.^2,'-','color',colorH,'LineWidth',1);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('Gamma-band [30-100 Hz]')
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[7b] Bilateral coherence','Gamma-band [30-100 Hz]',''})
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [7c] Pearson's correlations between bilateral gamma-band power during different arousal-states
ax3 = subplot(3,3,3);
s1 = scatter(CC_xInds*1,data.CorrCoef.Rest.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.CorrCoef.Rest.gammaBandPower.meanR,data.CorrCoef.Rest.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(CC_xInds*2,data.CorrCoef.Whisk.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.CorrCoef.Whisk.gammaBandPower.meanR,data.CorrCoef.Whisk.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(CC_xInds*3,data.CorrCoef.NREM.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.CorrCoef.NREM.gammaBandPower.meanR,data.CorrCoef.NREM.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(CC_xInds*4,data.CorrCoef.REM.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.CorrCoef.REM.gammaBandPower.meanR,data.CorrCoef.REM.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(CC_xInds2*5,data.CorrCoef.Awake.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.CorrCoef.Awake.gammaBandPower.meanR,data.CorrCoef.Awake.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
s6 = scatter(CC_xInds3*6,data.CorrCoef.Sleep.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.CorrCoef.Sleep.gammaBandPower.meanR,data.CorrCoef.Sleep.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
s7 = scatter(CC_xInds*7,data.CorrCoef.All.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e7 = errorbar(7,data.CorrCoef.All.gammaBandPower.meanR,data.CorrCoef.All.gammaBandPower.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
title({'[7c] Cortical Pearson''s corr. coef','Gamma-band [30-100 Hz]',''})
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
legend([s1,s2,s3,s4,s5,s6,s7],'Rest','Whisk','NREM','REM','Awake','Sleep','All')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields3) + 1])
ylim([-0.1,1])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [7d] Power spectra of HbT during different arousal-states
ax4 = subplot(3,3,4);
loglog(data.PowerSpec.Rest.CBV_HbT.meanCortf,data.PowerSpec.Rest.CBV_HbT.meanCortS,'color',colorA,'LineWidth',2);
hold on
loglog(data.PowerSpec.NREM.CBV_HbT.meanCortf,data.PowerSpec.NREM.CBV_HbT.meanCortS,'color',colorB,'LineWidth',2);
loglog(data.PowerSpec.REM.CBV_HbT.meanCortf,data.PowerSpec.REM.CBV_HbT.meanCortS,'color',colorC,'LineWidth',2);
loglog(data.PowerSpec.Awake.CBV_HbT.meanCortf,data.PowerSpec.Awake.CBV_HbT.meanCortS,'color',colorF,'LineWidth',2);
loglog(data.PowerSpec.Sleep.CBV_HbT.meanCortf,data.PowerSpec.Sleep.CBV_HbT.meanCortS,'color',colorG,'LineWidth',2);
loglog(data.PowerSpec.All.CBV_HbT.meanCortf,data.PowerSpec.All.CBV_HbT.meanCortS,'color',colorH,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title({'[7d] Cortical power','\DeltaHbT (\muM)',''})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.003,0.5]);
y2 = ylim(ax4);
ylim([y2(1)/2,y2(2)*2])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [7b] Coherence between bilateral HbT during different arousal-states
ax5 = subplot(3,3,5);
semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.meanC.^2,'color',colorA,'LineWidth',2);
hold on
semilogx(data.Coherr.NREM.CBV_HbT.meanf,data.Coherr.NREM.CBV_HbT.meanC.^2,'color',colorB,'LineWidth',2);
semilogx(data.Coherr.REM.CBV_HbT.meanf,data.Coherr.REM.CBV_HbT.meanC.^2,'color',colorC,'LineWidth',2);
semilogx(data.Coherr.Awake.CBV_HbT.meanf,data.Coherr.Awake.CBV_HbT.meanC.^2,'color',colorF,'LineWidth',2);
semilogx(data.Coherr.Sleep.CBV_HbT.meanf,data.Coherr.Sleep.CBV_HbT.meanC.^2,'color',colorG,'LineWidth',2);
semilogx(data.Coherr.All.CBV_HbT.meanf,data.Coherr.All.CBV_HbT.meanC.^2,'color',colorH,'LineWidth',2);
% confidence lines
semilogx(data.Coherr.Rest.CBV_HbT.meanf,data.Coherr.Rest.CBV_HbT.maxConfC_Y.^2,'-','color',colorA,'LineWidth',1);
semilogx(data.Coherr.NREM.CBV_HbT.meanf,data.Coherr.NREM.CBV_HbT.maxConfC_Y.^2,'-','color',colorB,'LineWidth',1);
semilogx(data.Coherr.REM.CBV_HbT.meanf,data.Coherr.REM.CBV_HbT.maxConfC_Y.^2,'-','color',colorC,'LineWidth',1);
semilogx(data.Coherr.Awake.CBV_HbT.meanf,data.Coherr.Awake.CBV_HbT.maxConfC_Y.^2,'-','color',colorF,'LineWidth',1);
semilogx(data.Coherr.Sleep.CBV_HbT.meanf,data.Coherr.Sleep.CBV_HbT.maxConfC_Y.^2,'-','color',colorG,'LineWidth',1);
semilogx(data.Coherr.All.CBV_HbT.meanf,data.Coherr.All.CBV_HbT.maxConfC_Y.^2,'-','color',colorH,'LineWidth',1);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title('\DeltaHbT (\muM)')
ylabel('Coherence^2')
xlabel('Freq (Hz)')
title({'[7b] Bilateral coherence','\DeltaHbT \muM (%)',''})
axis square
xlim([0.003,0.5])
ylim([0,1])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [7e] Pearson's correlations between bilateral HbT during different arousal-states
ax6 = subplot(3,3,6);
scatter(CC_xInds*1,data.CorrCoef.Rest.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.CorrCoef.Rest.CBV_HbT.meanR,data.CorrCoef.Rest.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(CC_xInds*2,data.CorrCoef.Whisk.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.CorrCoef.Whisk.CBV_HbT.meanR,data.CorrCoef.Whisk.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(CC_xInds*3,data.CorrCoef.NREM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.CorrCoef.NREM.CBV_HbT.meanR,data.CorrCoef.NREM.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(CC_xInds*4,data.CorrCoef.REM.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.CorrCoef.REM.CBV_HbT.meanR,data.CorrCoef.REM.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(CC_xInds2*5,data.CorrCoef.Awake.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.CorrCoef.Awake.CBV_HbT.meanR,data.CorrCoef.Awake.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(CC_xInds3*6,data.CorrCoef.Sleep.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.CorrCoef.Sleep.CBV_HbT.meanR,data.CorrCoef.Sleep.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(CC_xInds*7,data.CorrCoef.All.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e7 = errorbar(7,data.CorrCoef.All.CBV_HbT.meanR,data.CorrCoef.All.CBV_HbT.stdR,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
title({'[7e] Cortical Pearson''s corr. coef','\DeltaHbT \muM (%)',''})
ylabel({'Corr. Coefficient';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields3) + 1])
ylim([0,1])
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [7g] Arteriole power spectra of HbT during different arousal-states
ax7 = subplot(3,3,8);
loglog(data.VesselPowerSpec.Rest.meanf,data.VesselPowerSpec.Rest.meanS,'color',colorA,'LineWidth',2);
hold on
loglog(data.VesselPowerSpec.NREM.meanf,data.VesselPowerSpec.NREM.meanS,'color',colorB,'LineWidth',2);
loglog(data.VesselPowerSpec.REM.meanf,data.VesselPowerSpec.REM.meanS,'color',colorC,'LineWidth',2);
loglog(data.VesselPowerSpec.Awake.meanf,data.VesselPowerSpec.Awake.meanS,'color',colorF,'LineWidth',2);
loglog(data.VesselPowerSpec.All.meanf,data.VesselPowerSpec.All.meanS,'color',colorH,'LineWidth',2);
xline(1/10,'color','k');
xline(1/30,'color','k');
xline(1/60,'color','k');
title({'[7g] Arteriole power','\DeltaD/D (%)',''})
ylabel('Power (a.u.)')
xlabel('Freq (Hz)')
axis square
xlim([0.003,0.5]);
y2 = ylim(ax7);
ylim([y2(1)/2,y2(2)*2])
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Fig7']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-fillpage',[dirpath 'Fig7'])
%% statistical diary
diaryFile = [dirpath 'Fig7_Statistics.txt'];
if exist(diaryFile,'file') == 2
    delete(diaryFile)
end
diary(diaryFile)
diary on
% HbT statistical diary
disp('======================================================================================================================')
disp('[7c] Generalized linear mixed-effects model statistics for mean HbT corr. coef during Rest, Whisk, NREM, and REM')
disp('======================================================================================================================')
disp(gammaStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp(['Rest  gamma P/P R: ' num2str(round(data.CorrCoef.Rest.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Rest.gammaBandPower.stdR,2))]); disp(' ')
disp(['Whisk gamma P/P R: ' num2str(round(data.CorrCoef.Whisk.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Whisk.gammaBandPower.stdR,2))]); disp(' ')
disp(['NREM  gamma P/P R: ' num2str(round(data.CorrCoef.NREM.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.NREM.gammaBandPower.stdR,2))]); disp(' ')
disp(['REM   gamma P/P R: ' num2str(round(data.CorrCoef.REM.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.REM.gammaBandPower.stdR,2))]); disp(' ')
disp(['Awake gamma P/P R: ' num2str(round(data.CorrCoef.Awake.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Awake.gammaBandPower.stdR,2))]); disp(' ')
disp(['Sleep gamma P/P R: ' num2str(round(data.CorrCoef.Sleep.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Sleep.gammaBandPower.stdR,2))]); disp(' ')
disp(['All   gamma P/P R: ' num2str(round(data.CorrCoef.All.gammaBandPower.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.All.gammaBandPower.stdR,2))]); disp(' ')
disp('----------------------------------------------------------------------------------------------------------------------')
% gamma statistical diary
disp('======================================================================================================================')
disp('[7f] Generalized linear mixed-effects model statistics for mean gamma-band corr. coef during Rest, Whisk, NREM, and REM')
disp('======================================================================================================================')
disp(HbTStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp(['Rest  [HbT] (uM) R: ' num2str(round(data.CorrCoef.Rest.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Rest.CBV_HbT.stdR,2))]); disp(' ')
disp(['Whisk [HbT] (uM) R: ' num2str(round(data.CorrCoef.Whisk.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Whisk.CBV_HbT.stdR,2))]); disp(' ')
disp(['NREM  [HbT] (uM) R: ' num2str(round(data.CorrCoef.NREM.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.NREM.CBV_HbT.stdR,2))]); disp(' ')
disp(['REM   [HbT] (uM) R: ' num2str(round(data.CorrCoef.REM.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.REM.CBV_HbT.stdR,2))]); disp(' ')
disp(['Awake [HbT] (uM) R: ' num2str(round(data.CorrCoef.Awake.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Awake.CBV_HbT.stdR,2))]); disp(' ')
disp(['Sleep [HbT] (uM) R: ' num2str(round(data.CorrCoef.Sleep.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.Sleep.CBV_HbT.stdR,2))]); disp(' ')
disp(['All   [HbT] (uM) R: ' num2str(round(data.CorrCoef.All.CBV_HbT.meanR,2)) ' +/- ' num2str(round(data.CorrCoef.All.CBV_HbT.stdR,2))]); disp(' ')
disp('----------------------------------------------------------------------------------------------------------------------')
diary off

end
