function [] = Fig7_Manuscript2020(rootFolder,AnalysisResults)
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
z1 = 1; z2 = 1; z3 = 1;
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            for c = 1:length(dataTypes)
                dataType = dataTypes{1,c};
                data.Coherr.(behavField).(dataType).C(:,a) = (AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                data.Coherr.(behavField).(dataType).f(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(dataType).f;
                data.Coherr.(behavField).(dataType).confC(:,a) = AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC;
            end
        elseif strcmp(behavField,'Awake') == true
            if isfield(AnalysisResults.(animalID).Coherence,behavField) == true
                for e = 1:length(dataTypes)
                    dataType = dataTypes{1,e};
                    data.Coherr.(behavField).(dataType).C(:,z1) = (AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                    data.Coherr.(behavField).(dataType).f(:,z1) = AnalysisResults.(animalID).Coherence.(behavField).(dataType).f;
                    data.Coherr.(behavField).(dataType).confC(:,z1) = AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC;
                end
                z1 = z1 + 1;
            end
        elseif strcmp(behavField,'All') == true
            if isfield(AnalysisResults.(animalID).Coherence,behavField) == true
                for e = 1:length(dataTypes)
                    dataType = dataTypes{1,e};
                    data.Coherr.(behavField).(dataType).C(:,z2) = (AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                    data.Coherr.(behavField).(dataType).f(:,z2) = AnalysisResults.(animalID).Coherence.(behavField).(dataType).f;
                    data.Coherr.(behavField).(dataType).confC(:,z2) = AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC;
                end
                z2 = z2 + 1;
            end
        elseif strcmp(behavField,'Sleep') == true
            if isfield(AnalysisResults.(animalID).Coherence,behavField) == true
                for e = 1:length(dataTypes)
                    dataType = dataTypes{1,e};
                    data.Coherr.(behavField).(dataType).C(:,z3) = (AnalysisResults.(animalID).Coherence.(behavField).(dataType).C);
                    data.Coherr.(behavField).(dataType).f(:,z3) = AnalysisResults.(animalID).Coherence.(behavField).(dataType).f;
                    data.Coherr.(behavField).(dataType).confC(:,z3) = AnalysisResults.(animalID).Coherence.(behavField).(dataType).confC;
                end
                z3 = z3 + 1;
            end
        end
    end
end
% take the mean and standard deviation of each set of signals
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        data.Coherr.(behavField).(dataType).meanC = mean(data.Coherr.(behavField).(dataType).C,2);
        data.Coherr.(behavField).(dataType).stdC = std(data.Coherr.(behavField).(dataType).C,0,2);
        data.Coherr.(behavField).(dataType).meanf = mean(data.Coherr.(behavField).(dataType).f,2);
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
        if strcmp(behavField,'Awake') == true
            if isfield(AnalysisResults.(animalID).PowerSpectra,'Awake') == true
                for c = 1:length(dataTypes)
                    dataType = dataTypes{1,c};
                    data.PowerSpec.(behavField).(dataType).adjLH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.S;
                    data.PowerSpec.(behavField).(dataType).adjLH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.f;
                    data.PowerSpec.(behavField).(dataType).adjRH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.S;
                    data.PowerSpec.(behavField).(dataType).adjRH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.f;
                end
            else
                for c = 1:length(dataTypes)
                    dataType = dataTypes{1,c};
                    data.PowerSpec.(behavField).(dataType).adjLH.S{a,1} = [];
                    data.PowerSpec.(behavField).(dataType).adjLH.f{a,1} = [];
                    data.PowerSpec.(behavField).(dataType).adjRH.S{a,1} = [];
                    data.PowerSpec.(behavField).(dataType).adjRH.f{a,1} = [];
                end
            end
        elseif strcmp(behavField,'All') == true
            if isfield(AnalysisResults.(animalID).PowerSpectra,'All') == true
                for c = 1:length(dataTypes)
                    dataType = dataTypes{1,c};
                    data.PowerSpec.(behavField).(dataType).adjLH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.S;
                    data.PowerSpec.(behavField).(dataType).adjLH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.f;
                    data.PowerSpec.(behavField).(dataType).adjRH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.S;
                    data.PowerSpec.(behavField).(dataType).adjRH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.f;
                end
            end
        elseif strcmp(behavField,'Sleep') == true
            if isfield(AnalysisResults.(animalID).PowerSpectra,'Sleep') == true
                for c = 1:length(dataTypes)
                    dataType = dataTypes{1,c};
                    data.PowerSpec.(behavField).(dataType).adjLH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.S;
                    data.PowerSpec.(behavField).(dataType).adjLH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.f;
                    data.PowerSpec.(behavField).(dataType).adjRH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.S;
                    data.PowerSpec.(behavField).(dataType).adjRH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.f;
                end
            end
        else
            for c = 1:length(dataTypes)
                dataType = dataTypes{1,c};
                data.PowerSpec.(behavField).(dataType).adjLH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.S;
                data.PowerSpec.(behavField).(dataType).adjLH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjLH.f;
                data.PowerSpec.(behavField).(dataType).adjRH.S{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.S;
                data.PowerSpec.(behavField).(dataType).adjRH.f{a,1} = AnalysisResults.(animalID).PowerSpectra.(behavField).(dataType).adjRH.f;
            end
        end
    end
end
%
for a = 1:length(IOS_animalIDs)
    for c = 1:length(dataTypes)
        dataType = dataTypes{1,c};
        data.PowerSpec.baseline.(dataType).LH{a,1} = max(data.PowerSpec.Rest.(dataType).adjLH.S{a,1});
        data.PowerSpec.baseline.(dataType).RH{a,1} = max(data.PowerSpec.Rest.(dataType).adjRH.S{a,1});
    end
end
%%
for a = 1:length(IOS_animalIDs)
    for dd = 1:length(behavFields)
        behavField = behavFields{1,dd};
        for j = 1:length(dataTypes)
            dataType = dataTypes{1,j};
            if isempty(data.PowerSpec.(behavField).(dataType).adjLH.S{a,1}) == false
                for ee = 1:size(data.PowerSpec.(behavField).(dataType).adjLH.S,2)
                    data.PowerSpec.(behavField).(dataType).normLH{a,1} = (data.PowerSpec.(behavField).(dataType).adjLH.S{a,1})*(1/(data.PowerSpec.baseline.(dataType).LH{a,1}));
                    data.PowerSpec.(behavField).(dataType).normRH{a,1} = (data.PowerSpec.(behavField).(dataType).adjRH.S{a,1})*(1/(data.PowerSpec.baseline.(dataType).RH{a,1}));
                end
            end
        end
    end
end
% concatenate the data from the left and right hemispheres
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
% take the mean and standard deviation of each set of signals
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
data.VesselPowerSpec.Rest.S = []; data.VesselPowerSpec.NREM.S = []; data.VesselPowerSpec.REM.S = []; data.VesselPowerSpec.Awake.S = []; data.VesselPowerSpec.All.S = [];
data.VesselPowerSpec.Rest.f = []; data.VesselPowerSpec.NREM.f = []; data.VesselPowerSpec.REM.f = []; data.VesselPowerSpec.Awake.f = []; data.VesselPowerSpec.All.f = [];
data.VesselPowerSpec.whiskingPerc = [];
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
%
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    vesselIDs = fieldnames(data.VesselPowerSpec.(animalID).Rest);
    for cc = 1:length(vesselIDs)
        vesselID = vesselIDs{cc,1};
        data.VesselPowerSpec.(animalID).baseline.(vesselID) = max(data.VesselPowerSpec.(animalID).Rest.(vesselID).S);
    end
end
%
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
%
for aa = 1:length(TwoP_animalIDs)
    animalID = TwoP_animalIDs{1,aa};
    for bb = 1:length(behavFields2)
        behavField = behavFields2{1,bb};
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
% take the average power of the vessels for each behavior
for dd = 1:length(behavFields2)
    behavField = behavFields2{1,dd};
    data.VesselPowerSpec.(behavField).meanS = mean(data.VesselPowerSpec.(behavField).S,2);
    data.VesselPowerSpec.(behavField).StDS = std(data.VesselPowerSpec.(behavField).S,0,2);
    data.VesselPowerSpec.(behavField).meanf = mean(data.VesselPowerSpec.(behavField).f,1);
end
%% Pearson's correlations during different behaviors
% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(IOS_animalIDs)
    animalID = IOS_animalIDs{1,a};
    for b = 1:length(behavFields3)
        behavField = behavFields3{1,b};
        if isfield(AnalysisResults.(animalID).CorrCoeff,behavField) == true
            for c = 1:length(dataTypes)
                dataType = dataTypes{1,c};
                data.CorrCoef.(behavField).(dataType).R{a,1} = AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).R;
                data.CorrCoef.(behavField).(dataType).meanRs(a,1) = AnalysisResults.(animalID).CorrCoeff.(behavField).(dataType).meanR;
                data.CorrCoef.(behavField).animalID{a,1} = animalID;
                data.CorrCoef.(behavField).behavior{a,1} = behavField;
            end
        end
    end
end
% concatenate the data and take mean/STD
for e = 1:length(behavFields3)
    behavField = behavFields3{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        data.CorrCoef.(behavField).(dataType).catR = [];
        for z = 1:length(data.CorrCoef.(behavField).(dataType).R)
            data.CorrCoef.(behavField).(dataType).catR = cat(1,data.CorrCoef.(behavField).(dataType).catR,data.CorrCoef.(behavField).(dataType).R{z,1});
        end
        data.CorrCoef.(behavField).(dataType).meanR = mean(data.CorrCoef.(behavField).(dataType).meanRs,1);
        data.CorrCoef.(behavField).(dataType).stdR = std(data.CorrCoef.(behavField).(dataType).meanRs,0,1);
    end
end
% %% statistics - linear mixed effects model
% numComparisons = 3;
% % HbT
% CCHbT_alphaConf = [0.05,0.01,0.001];
% HbTtableSize = cat(1,data.CorrCoef.Rest.CBV_HbT.meanRs,data.CorrCoef.Whisk.CBV_HbT.meanRs,data.CorrCoef.NREM.CBV_HbT.meanRs,data.CorrCoef.REM.CBV_HbT.meanRs);
% HbTTable = table('Size',[size(HbTtableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
% HbTTable.Mouse = cat(1,data.CorrCoef.Rest.animalID,data.CorrCoef.Whisk.animalID,data.CorrCoef.NREM.animalID,data.CorrCoef.REM.animalID);
% HbTTable.CorrCoef = cat(1,data.CorrCoef.Rest.CBV_HbT.meanRs,data.CorrCoef.Whisk.CBV_HbT.meanRs,data.CorrCoef.NREM.CBV_HbT.meanRs,data.CorrCoef.REM.CBV_HbT.meanRs);
% HbTTable.Behavior = cat(1,data.CorrCoef.Rest.behavior,data.CorrCoef.Whisk.behavior,data.CorrCoef.NREM.behavior,data.CorrCoef.REM.behavior);
% HbTFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
% HbTStats = fitglme(HbTTable,HbTFitFormula);
% for z = 1:length(CCHbT_alphaConf)
%     HbTCI{z,1} = coefCI(HbTStats,'Alpha',(CCHbT_alphaConf(z)/numComparisons)); %#ok<*AGROW>
% end
% % gamma-band power
% CCGamma_alphaConf = [0.05,0.01,0.001];
% gammaTableSize = cat(1,data.CorrCoef.Rest.gammaBandPower.meanRs,data.CorrCoef.Whisk.gammaBandPower.meanRs,data.CorrCoef.NREM.gammaBandPower.meanRs,data.CorrCoef.REM.gammaBandPower.meanRs);
% gammaTable = table('Size',[size(gammaTableSize,1),3],'VariableTypes',{'string','double','string'},'VariableNames',{'Mouse','CorrCoef','Behavior'});
% gammaTable.Mouse = cat(1,data.CorrCoef.Rest.animalID,data.CorrCoef.Whisk.animalID,data.CorrCoef.NREM.animalID,data.CorrCoef.REM.animalID);
% gammaTable.CorrCoef = cat(1,data.CorrCoef.Rest.gammaBandPower.meanRs,data.CorrCoef.Whisk.gammaBandPower.meanRs,data.CorrCoef.NREM.gammaBandPower.meanRs,data.CorrCoef.REM.gammaBandPower.meanRs);
% gammaTable.Behavior = cat(1,data.CorrCoef.Rest.behavior,data.CorrCoef.Whisk.behavior,data.CorrCoef.NREM.behavior,data.CorrCoef.REM.behavior);
% gammaFitFormula = 'CorrCoef ~ 1 + Behavior + (1|Mouse)';
% gammaStats = fitglme(gammaTable,gammaFitFormula);
% for z = 1:length(CCGamma_alphaConf)
%     gammaCI{z,1} = coefCI(gammaStats,'Alpha',(CCGamma_alphaConf(z)/numComparisons));
% end
%% Figure Panel 7
summaryFigure = figure('Name','Fig7 (a-g)');
sgtitle('Figure Panel 7 (a-g) Turner Manuscript 2020')
CC_xInds = ones(1,length(IOS_animalIDs));
CC_xInds2 = ones(1,length(data.CorrCoef.Awake.animalID));
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
s6 = scatter(CC_xInds*6,data.CorrCoef.Sleep.gammaBandPower.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
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
scatter(CC_xInds*6,data.CorrCoef.Sleep.CBV_HbT.meanRs,75,'MarkerEdgeColor','k','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
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
disp(HbTStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.05 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HbTCI{1,1}(1,:))])
disp(['Whisk: ' num2str(HbTCI{1,1}(2,:))])
disp(['NREM: ' num2str(HbTCI{1,1}(3,:))])
disp(['REM: ' num2str(HbTCI{1,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.01 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HbTCI{2,1}(1,:))])
disp(['Whisk: ' num2str(HbTCI{2,1}(2,:))])
disp(['NREM: ' num2str(HbTCI{2,1}(3,:))])
disp(['REM: ' num2str(HbTCI{2,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.001 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(HbTCI{3,1}(1,:))])
disp(['Whisk: ' num2str(HbTCI{3,1}(2,:))])
disp(['NREM: ' num2str(HbTCI{3,1}(3,:))])
disp(['REM: ' num2str(HbTCI{3,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
% gamma statistical diary
disp('======================================================================================================================')
disp('[7f] Generalized linear mixed-effects model statistics for mean gamma-band corr. coef during Rest, Whisk, NREM, and REM')
disp('======================================================================================================================')
disp(gammaStats)
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.05 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(gammaCI{1,1}(1,:))])
disp(['Whisk: ' num2str(gammaCI{1,1}(2,:))])
disp(['NREM: ' num2str(gammaCI{1,1}(3,:))])
disp(['REM: ' num2str(gammaCI{1,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.01 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(gammaCI{2,1}(1,:))])
disp(['Whisk: ' num2str(gammaCI{2,1}(2,:))])
disp(['NREM: ' num2str(gammaCI{2,1}(3,:))])
disp(['REM: ' num2str(gammaCI{2,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
disp('Alpha = 0.001 confidence interval with 3 comparisons to ''Rest'' (Intercept): ')
disp(['Rest: ' num2str(gammaCI{3,1}(1,:))])
disp(['Whisk: ' num2str(gammaCI{3,1}(2,:))])
disp(['NREM: ' num2str(gammaCI{3,1}(3,:))])
disp(['REM: ' num2str(gammaCI{3,1}(4,:))])
disp('----------------------------------------------------------------------------------------------------------------------')
diary off

end
