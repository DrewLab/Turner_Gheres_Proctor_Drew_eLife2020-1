function [] = Fig7b_Manuscript2020(rootFolder,AnalysisResults)
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
% find 0.1/0.01 Hz peaks in coherence
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        for g = 1:size(data.Coherr.(behavField).(dataType).C,2)
            if strcmp(behavField,'Rest') == true
                f_short = data.Coherr.(behavField).(dataType).f(g,:);
                C = data.Coherr.(behavField).(dataType).C(:,g);
                f_long = 0:0.01:0.5;
                C_long = interp1(f_short,C,f_long);
                index01 = find(f_long == 0.1);
                data.Coherr.(behavField).(dataType).C01(g,1) = C_long(index01).^2; %#ok<FNDSB>
            elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
                F = round(data.Coherr.(behavField).(dataType).f(g,:),2);
                C = data.Coherr.(behavField).(dataType).C(:,g);
                index01 = find(F == 0.1);
                data.Coherr.(behavField).(dataType).C01(g,1) = C(index01(1)).^2;
            else
                F = round(data.Coherr.(behavField).(dataType).f(g,:),3);
                C = data.Coherr.(behavField).(dataType).C(:,g);
                index01 = find(F == 0.1);
                index001 = find(F == 0.01);
                data.Coherr.(behavField).(dataType).C01(g,1) = C(index01(1)).^2;
                data.Coherr.(behavField).(dataType).C001(g,1) = C(index001(1)).^2;
            end
        end
    end
end
% take mean/StD of peak C
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            data.Coherr.(behavField).(dataType).meanC01 = mean(data.Coherr.(behavField).(dataType).C01,1);
            data.Coherr.(behavField).(dataType).stdC01 = std(data.Coherr.(behavField).(dataType).C01,0,1);
        else
            data.Coherr.(behavField).(dataType).meanC01 = mean(data.Coherr.(behavField).(dataType).C01,1);
            data.Coherr.(behavField).(dataType).stdC01 = std(data.Coherr.(behavField).(dataType).C01,0,1);
            data.Coherr.(behavField).(dataType).meanC001 = mean(data.Coherr.(behavField).(dataType).C001,1);
            data.Coherr.(behavField).(dataType).stdC001 = std(data.Coherr.(behavField).(dataType).C001,0,1);
        end
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
% find 0.1/0.01 Hz peaks in PSD
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        for g = 1:size(data.PowerSpec.(behavField).(dataType).cat_S,2)
            if strcmp(behavField,'Rest') == true
                f_short = data.PowerSpec.(behavField).(dataType).cat_f(g,:);
                S = data.PowerSpec.(behavField).(dataType).cat_S(:,g);
                f_long = 0:0.01:0.5;
                S_long = interp1(f_short,S,f_long);
                index01 = find(f_long == 0.1);
                data.PowerSpec.(behavField).(dataType).S01(g,1) = S_long(index01); %#ok<FNDSB>
            elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
                F = round(data.PowerSpec.(behavField).(dataType).cat_f(g,:),2);
                S = data.PowerSpec.(behavField).(dataType).cat_S(:,g);
                index01 = find(F == 0.1);
                data.PowerSpec.(behavField).(dataType).S01(g,1) = S(index01(1));
            else
                F = round(data.PowerSpec.(behavField).(dataType).cat_f(g,:),3);
                S = data.PowerSpec.(behavField).(dataType).cat_S(:,g);
                index01 = find(F == 0.1);
                index001 = find(F == 0.01);
                data.PowerSpec.(behavField).(dataType).S01(g,1) = S(index01(1));
                data.PowerSpec.(behavField).(dataType).S001(g,1) = S(index001(1));
            end
        end
    end
end
% take mean/StD of peak S
for e = 1:length(behavFields)
    behavField = behavFields{1,e};
    for f = 1:length(dataTypes)
        dataType = dataTypes{1,f};
        if strcmp(behavField,'Rest') == true || strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            data.PowerSpec.(behavField).(dataType).meanS01 = mean(data.PowerSpec.(behavField).(dataType).S01,1);
            data.PowerSpec.(behavField).(dataType).stdS01 = std(data.PowerSpec.(behavField).(dataType).S01,0,1);
        else
            data.PowerSpec.(behavField).(dataType).meanS01 = mean(data.PowerSpec.(behavField).(dataType).S01,1);
            data.PowerSpec.(behavField).(dataType).stdS01 = std(data.PowerSpec.(behavField).(dataType).S01,0,1);
            data.PowerSpec.(behavField).(dataType).meanS001 = mean(data.PowerSpec.(behavField).(dataType).S001,1);
            data.PowerSpec.(behavField).(dataType).stdS001 = std(data.PowerSpec.(behavField).(dataType).S001,0,1);
        end
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
% find 0.1/0.01 Hz peaks in PSD
for e = 1:length(behavFields2)
    behavField = behavFields2{1,e};
    for g = 1:size(data.VesselPowerSpec.(behavField).S,2)
        if strcmp(behavField,'Rest') == true
            f_short = data.VesselPowerSpec.(behavField).f(g,:);
            S = data.VesselPowerSpec.(behavField).S(:,g);
            f_long = 0:0.01:0.5;
            S_long = interp1(f_short,S,f_long);
            index01 = find(f_long == 0.1);
            data.VesselPowerSpec.(behavField).S01(g,1) = S_long(index01); %#ok<FNDSB>
        elseif strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
            F = round(data.VesselPowerSpec.(behavField).f(g,:),2);
            S = data.VesselPowerSpec.(behavField).S(:,g);
            index01 = find(F == 0.1);
            data.VesselPowerSpec.(behavField).S01(g,1) = S(index01(1));
        else
            F = round(data.VesselPowerSpec.(behavField).f(g,:),3);
            S = data.VesselPowerSpec.(behavField).S(:,g);
            index01 = find(F == 0.1);
            index001 = find(F == 0.01);
            data.VesselPowerSpec.(behavField).S01(g,1) = S(index01(1));
            data.VesselPowerSpec.(behavField).S001(g,1) = S(index001(1));
        end
    end
end
% take mean/StD of peak S
for e = 1:length(behavFields2)
    behavField = behavFields2{1,e};
    if strcmp(behavField,'Rest') == true || strcmp(behavField,'NREM') == true || strcmp(behavField,'REM') == true
        data.VesselPowerSpec.(behavField).meanS01 = mean(data.VesselPowerSpec.(behavField).S01,1);
        data.VesselPowerSpec.(behavField).stdS01 = std(data.VesselPowerSpec.(behavField).S01,0,1);
    else
        data.VesselPowerSpec.(behavField).meanS01 = mean(data.VesselPowerSpec.(behavField).S01,1);
        data.VesselPowerSpec.(behavField).stdS01 = std(data.VesselPowerSpec.(behavField).S01,0,1);
        data.VesselPowerSpec.(behavField).meanS001 = mean(data.VesselPowerSpec.(behavField).S001,1);
        data.VesselPowerSpec.(behavField).stdS001 = std(data.VesselPowerSpec.(behavField).S001,0,1);
    end
end
%% statistics - linear mixed effects model
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
%% [7d] HbT PSD
ax1 = subplot(2,3,1);
scatter(ones(1,length(data.PowerSpec.Rest.gammaBandPower.S01))*1,data.PowerSpec.Rest.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Rest.gammaBandPower.meanS01,data.PowerSpec.Rest.gammaBandPower.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.NREM.gammaBandPower.S01))*2,data.PowerSpec.NREM.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.NREM.gammaBandPower.meanS01,data.PowerSpec.NREM.gammaBandPower.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.REM.gammaBandPower.S01))*3,data.PowerSpec.REM.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.REM.gammaBandPower.meanS01,data.PowerSpec.REM.gammaBandPower.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Awake.gammaBandPower.S01))*4,data.PowerSpec.Awake.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.PowerSpec.Awake.gammaBandPower.meanS01,data.PowerSpec.Awake.gammaBandPower.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Awake.gammaBandPower.S001))*5,data.PowerSpec.Awake.gammaBandPower.S001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.PowerSpec.Awake.gammaBandPower.meanS001,data.PowerSpec.Awake.gammaBandPower.stdS001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.gammaBandPower.S01))*6,data.PowerSpec.Sleep.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.PowerSpec.Sleep.gammaBandPower.meanS01,data.PowerSpec.Sleep.gammaBandPower.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.gammaBandPower.S001))*7,data.PowerSpec.Sleep.gammaBandPower.S001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
e7 = errorbar(7,data.PowerSpec.Sleep.gammaBandPower.meanS001,data.PowerSpec.Sleep.gammaBandPower.stdS001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.gammaBandPower.S01))*8,data.PowerSpec.All.gammaBandPower.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e8 = errorbar(8,data.PowerSpec.All.gammaBandPower.meanS01,data.PowerSpec.All.gammaBandPower.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.gammaBandPower.S001))*9,data.PowerSpec.All.gammaBandPower.S001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e9 = errorbar(9,data.PowerSpec.All.gammaBandPower.meanS001,data.PowerSpec.All.gammaBandPower.stdS001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
title({'[7a] PSD @ 0.1,0.01 Hz','Gamma-band [30-100 Hz]',''})
ylabel('Power (a.u.) @ 1,0.01 Hz')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,10])
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [7b] Gamma-band Coherence
ax2 = subplot(2,3,2);
s1 = scatter(ones(1,length(data.Coherr.Rest.gammaBandPower.C01))*1,data.Coherr.Rest.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Rest.gammaBandPower.meanC01,data.Coherr.Rest.gammaBandPower.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
s2 = scatter(ones(1,length(data.Coherr.NREM.gammaBandPower.C01))*2,data.Coherr.NREM.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.NREM.gammaBandPower.meanC01,data.Coherr.NREM.gammaBandPower.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
s3 = scatter(ones(1,length(data.Coherr.REM.gammaBandPower.C01))*3,data.Coherr.REM.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.REM.gammaBandPower.meanC01,data.Coherr.REM.gammaBandPower.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
s4 = scatter(ones(1,length(data.Coherr.Awake.gammaBandPower.C01))*4,data.Coherr.Awake.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.Coherr.Awake.gammaBandPower.meanC01,data.Coherr.Awake.gammaBandPower.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
s5 = scatter(ones(1,length(data.Coherr.Awake.gammaBandPower.C001))*5,data.Coherr.Awake.gammaBandPower.C001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.Coherr.Awake.gammaBandPower.meanC001,data.Coherr.Awake.gammaBandPower.stdC001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
s6 = scatter(ones(1,length(data.Coherr.Sleep.gammaBandPower.C01))*6,data.Coherr.Sleep.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.Coherr.Sleep.gammaBandPower.meanC01,data.Coherr.Sleep.gammaBandPower.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
s7 = scatter(ones(1,length(data.Coherr.Sleep.gammaBandPower.C001))*7,data.Coherr.Sleep.gammaBandPower.C001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
e7 = errorbar(7,data.Coherr.Sleep.gammaBandPower.meanC001,data.Coherr.Sleep.gammaBandPower.stdC001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
s8 = scatter(ones(1,length(data.Coherr.All.gammaBandPower.C01))*8,data.Coherr.All.gammaBandPower.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e8 = errorbar(8,data.Coherr.All.gammaBandPower.meanC01,data.Coherr.All.gammaBandPower.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
s9 = scatter(ones(1,length(data.Coherr.All.gammaBandPower.C001))*9,data.Coherr.All.gammaBandPower.C001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e9 = errorbar(9,data.Coherr.All.gammaBandPower.meanC001,data.Coherr.All.gammaBandPower.stdC001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
title({'[7b] Coherence^2 @ 0.1,0.01 Hz','Gamma-band [30-100 Hz]',''})
ylabel({'Coherence^2 @ 0.1,0.01 Hz';'Left hem vs. Right hem'})
legend([s1,s2,s3,s4,s5,s6,s7,s8,s9],'Rest','NREM','REM','Awake','Awake 0.01','Sleep','Sleep 0.01','All','All 0.01')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,10])
ylim([0,1])
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [7d] HbT PSD
ax3 = subplot(2,3,4);
scatter(ones(1,length(data.PowerSpec.Rest.CBV_HbT.S01))*1,data.PowerSpec.Rest.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.PowerSpec.Rest.CBV_HbT.meanS01,data.PowerSpec.Rest.CBV_HbT.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.NREM.CBV_HbT.S01))*2,data.PowerSpec.NREM.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.PowerSpec.NREM.CBV_HbT.meanS01,data.PowerSpec.NREM.CBV_HbT.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.REM.CBV_HbT.S01))*3,data.PowerSpec.REM.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.PowerSpec.REM.CBV_HbT.meanS01,data.PowerSpec.REM.CBV_HbT.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Awake.CBV_HbT.S01))*4,data.PowerSpec.Awake.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.PowerSpec.Awake.CBV_HbT.meanS01,data.PowerSpec.Awake.CBV_HbT.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Awake.CBV_HbT.S001))*5,data.PowerSpec.Awake.CBV_HbT.S001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.PowerSpec.Awake.CBV_HbT.meanS001,data.PowerSpec.Awake.CBV_HbT.stdS001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.CBV_HbT.S01))*6,data.PowerSpec.Sleep.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.PowerSpec.Sleep.CBV_HbT.meanS01,data.PowerSpec.Sleep.CBV_HbT.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.Sleep.CBV_HbT.S001))*7,data.PowerSpec.Sleep.CBV_HbT.S001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
e7 = errorbar(7,data.PowerSpec.Sleep.CBV_HbT.meanS001,data.PowerSpec.Sleep.CBV_HbT.stdS001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.CBV_HbT.S01))*8,data.PowerSpec.All.CBV_HbT.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e8 = errorbar(8,data.PowerSpec.All.CBV_HbT.meanS01,data.PowerSpec.All.CBV_HbT.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(ones(1,length(data.PowerSpec.All.CBV_HbT.S001))*9,data.PowerSpec.All.CBV_HbT.S001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e9 = errorbar(9,data.PowerSpec.All.CBV_HbT.meanS001,data.PowerSpec.All.CBV_HbT.stdS001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
title({'[7d] PSD @ 0.1,0.01 Hz','\DeltaHbT',''})
ylabel('Power (a.u.) @ 1,0.01 Hz')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,10])
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [7e] HbT Coherence
ax4 = subplot(2,3,5);
scatter(ones(1,length(data.Coherr.Rest.CBV_HbT.C01))*1,data.Coherr.Rest.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Coherr.Rest.CBV_HbT.meanC01,data.Coherr.Rest.CBV_HbT.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.Coherr.NREM.CBV_HbT.C01))*2,data.Coherr.NREM.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.Coherr.NREM.CBV_HbT.meanC01,data.Coherr.NREM.CBV_HbT.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.Coherr.REM.CBV_HbT.C01))*3,data.Coherr.REM.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.Coherr.REM.CBV_HbT.meanC01,data.Coherr.REM.CBV_HbT.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.Coherr.Awake.CBV_HbT.C01))*4,data.Coherr.Awake.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.Coherr.Awake.CBV_HbT.meanC01,data.Coherr.Awake.CBV_HbT.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.Coherr.Awake.CBV_HbT.C001))*5,data.Coherr.Awake.CBV_HbT.C001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.Coherr.Awake.CBV_HbT.meanC001,data.Coherr.Awake.CBV_HbT.stdC001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.CBV_HbT.C01))*6,data.Coherr.Sleep.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.Coherr.Sleep.CBV_HbT.meanC01,data.Coherr.Sleep.CBV_HbT.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(ones(1,length(data.Coherr.Sleep.CBV_HbT.C001))*7,data.Coherr.Sleep.CBV_HbT.C001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorG,'jitter','on', 'jitterAmount',0.25);
e7 = errorbar(7,data.Coherr.Sleep.CBV_HbT.meanC001,data.Coherr.Sleep.CBV_HbT.stdC001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.CBV_HbT.C01))*8,data.Coherr.All.CBV_HbT.C01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e8 = errorbar(8,data.Coherr.All.CBV_HbT.meanC01,data.Coherr.All.CBV_HbT.stdC01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e8.Color = 'black';
e8.MarkerSize = 10;
e8.CapSize = 10;
scatter(ones(1,length(data.Coherr.All.CBV_HbT.C001))*9,data.Coherr.All.CBV_HbT.C001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e9 = errorbar(9,data.Coherr.All.CBV_HbT.meanC001,data.Coherr.All.CBV_HbT.stdC001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e9.Color = 'black';
e9.MarkerSize = 10;
e9.CapSize = 10;
title({'[7e] Coherence^2 @ 0.1,0.01 Hz','\DeltaHbT',''})
ylabel({'Coherence^2 @ 0.1,0.01 Hz';'Left hem vs. Right hem'})
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,10])
ylim([0,1])
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [7g] Vessel PSD
ax5 = subplot(2,3,6);
scatter(ones(1,length(data.VesselPowerSpec.Rest.S01))*1,data.VesselPowerSpec.Rest.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on', 'jitterAmount',0.25);
hold on
e1 = errorbar(1,data.VesselPowerSpec.Rest.meanS01,data.VesselPowerSpec.Rest.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
e1.MarkerSize = 10;
e1.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.NREM.S01))*2,data.VesselPowerSpec.NREM.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on', 'jitterAmount',0.25);
e2 = errorbar(2,data.VesselPowerSpec.NREM.meanS01,data.VesselPowerSpec.NREM.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
e2.MarkerSize = 10;
e2.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.REM.S01))*3,data.VesselPowerSpec.REM.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on', 'jitterAmount',0.25);
e3 = errorbar(3,data.VesselPowerSpec.REM.meanS01,data.VesselPowerSpec.REM.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
e3.MarkerSize = 10;
e3.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.Awake.S01))*4,data.VesselPowerSpec.Awake.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e4 = errorbar(4,data.VesselPowerSpec.Awake.meanS01,data.VesselPowerSpec.Awake.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
e4.MarkerSize = 10;
e4.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.Awake.S001))*5,data.VesselPowerSpec.Awake.S001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorF,'jitter','on', 'jitterAmount',0.25);
e5 = errorbar(5,data.VesselPowerSpec.Awake.meanS001,data.VesselPowerSpec.Awake.stdS001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e5.Color = 'black';
e5.MarkerSize = 10;
e5.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.All.S01))*6,data.VesselPowerSpec.All.S01,75,'MarkerEdgeColor','k','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e6 = errorbar(6,data.VesselPowerSpec.All.meanS01,data.VesselPowerSpec.All.stdS01,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e6.Color = 'black';
e6.MarkerSize = 10;
e6.CapSize = 10;
scatter(ones(1,length(data.VesselPowerSpec.All.S001))*7,data.VesselPowerSpec.All.S001,75,'MarkerEdgeColor','w','MarkerFaceColor',colorH,'jitter','on', 'jitterAmount',0.25);
e7 = errorbar(7,data.VesselPowerSpec.All.meanS001,data.VesselPowerSpec.All.stdS001,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e7.Color = 'black';
e7.MarkerSize = 10;
e7.CapSize = 10;
title({'[7g] PSD @ 0.1,0.01 Hz','\DeltaD/D',''})
ylabel('Power (a.u.) @ 1,0.01 Hz')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yscale','log')
axis square
xlim([0,10])
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
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
