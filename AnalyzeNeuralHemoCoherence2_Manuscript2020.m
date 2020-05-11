function [AnalysisResults] = AnalyzeNeuralHemoCoherence2_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Analyze the spectral coherence between bilateral hemodynamic and neural signals.
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
dataTypes = {'deltaBandPower','thetaBandPower','alphaBandPower','betaBandPower','gammaBandPower'};
hemDataTypes = {'adjLH','adjRH'};
params.minTime.Rest = 10;   % seconds
params.minTime.NREM = 30;   % seconds
params.minTime.REM = 60;   % seconds

%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    % character list of all ProcData file IDs
    procDataFileStruct = dir('*_ProcData.mat');
    procDataFiles = {procDataFileStruct.name}';
    procDataFileIDs = char(procDataFiles);
    % find and load RestingBaselines.mat struct
    baselineDataFileStruct = dir('*_RestingBaselines.mat');
    baselineDataFile = {baselineDataFileStruct.name}';
    baselineDataFileID = char(baselineDataFile);
    load(baselineDataFileID)
    % identify animal's ID and pull important infortmat
    samplingRate = 30;
    % lowpass filter and detrend each segment
    [z,p,k] = butter(4,1/(samplingRate/2),'low');
    [sos,g] = zp2sos(z,p,k);
    % go through each valid data type for behavior-based coherence analysis
    for zzz = 1:length(hemDataTypes)
        hemDataType = hemDataTypes{1,zzz};
        for aa = 1:length(dataTypes)
            dataType = dataTypes{1,aa};
            %% Analyze coherence during awake periods with no sleep scores
            zz = 1;
            clear HbT_AllUnstimData Gamma_AllUnstimData HbT_ProcAllUnstimData Gamma_ProcAllUnstimData
            HbT_AllUnstimData = [];
            for bb = 1:size(procDataFileIDs,1)
                procDataFileID = procDataFileIDs(bb,:);
                [~,allDataFileDate,~] = GetFileInfo_IOS_Manuscript2020(procDataFileID);
                strDay = ConvertDate_IOS_Manuscript2020(allDataFileDate);
                load(procDataFileID)
                puffs = ProcData.data.solenoids.LPadSol;
                if isempty(puffs) == true
                    HbT_AllUnstimData{zz,1} = ProcData.data.CBV_HbT.(hemDataType); %#ok<*AGROW>
                    Gamma_AllUnstimData{zz,1} = (ProcData.data.(['cortical_' hemDataType(4:5)]).(dataType) - RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay))./RestingBaselines.manualSelection.(['cortical_' hemDataType(4:5)]).(dataType).(strDay);
                    zz = zz + 1;
                end
            end
            % process
            if isempty(HbT_AllUnstimData) == false
                for bb = 1:length(HbT_AllUnstimData)
                    HbT_ProcAllUnstimData{bb,1} = filtfilt(sos,g,detrend(HbT_AllUnstimData{bb,1},'constant'));
                    Gamma_ProcAllUnstimData{bb,1} = filtfilt(sos,g,detrend(Gamma_AllUnstimData{bb,1},'constant'));
                end
                % input data as time(1st dimension, vertical) by trials (2nd dimension, horizontunstimy)
                HbT_awakeData = zeros(length(HbT_ProcAllUnstimData{1,1}),length(HbT_ProcAllUnstimData));
                Gamma_awakeData = zeros(length(Gamma_ProcAllUnstimData{1,1}),length(Gamma_ProcAllUnstimData));
                for cc = 1:length(HbT_ProcAllUnstimData)
                    HbT_awakeData(:,cc) = HbT_ProcAllUnstimData{cc,1};
                    Gamma_awakeData(:,cc) = Gamma_ProcAllUnstimData{cc,1};
                end
                % parameters for coherencyc - information available in function
                params.tapers = [5,9];   % Tapers [n, 2n - 1]
                params.pad = 1;
                params.Fs = samplingRate;   % Sampling Rate
                params.fpass = [0,0.5];   % Pass band [0, nyquist]
                params.trialave = 1;
                params.err = [2,0.05];
                % calculate the coherence between desired signals
                [C_AllUnstimData,~,~,~,~,f_AllUnstimData,confC_AllUnstimData,~,cErr_AllUnstimData] = coherencyc_Manuscript2020(HbT_awakeData,Gamma_awakeData,params);
                % save data and figures
                AnalysisResults.(animalID).NeuralHemoCoherence.All.(dataType).(hemDataType).C = C_AllUnstimData;
                AnalysisResults.(animalID).NeuralHemoCoherence.All.(dataType).(hemDataType).f = f_AllUnstimData;
                AnalysisResults.(animalID).NeuralHemoCoherence.All.(dataType).(hemDataType).confC = confC_AllUnstimData;
                AnalysisResults.(animalID).NeuralHemoCoherence.All.(dataType).(hemDataType).cErr = cErr_AllUnstimData;
                % save figures if desired
                if strcmp(saveFigs,'y') == true
                    awakeCoherence = figure;
                    plot(f_AllUnstimData,C_AllUnstimData,'k')
                    hold on;
                    plot(f_AllUnstimData,cErr_AllUnstimData,'color',colors_Manuscript2020('battleship grey'))
                    xlabel('Freq (Hz)');
                    ylabel('Coherence');
                    title([animalID  ' ' dataType ' coherence for awake data']);
                    set(gca,'Ticklength',[0,0]);
                    legend('Coherence','Jackknife Lower','Jackknife Upper','Location','Southeast');
                    set(legend,'FontSize',6);
                    ylim([0,1])
                    xlim([0.1,0.5])
                    axis square
                    set(gca,'box','off')
                    [pathstr,~,~] = fileparts(cd);
                    dirpath = [pathstr '/Figures/Coherence/'];
                    if ~exist(dirpath,'dir')
                        mkdir(dirpath);
                    end
                    savefig(awakeCoherence,[dirpath animalID '_All_' dataType '_' hemDataType '_Coherence']);
                    close(awakeCoherence)
                end
            end
        end
    end
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
