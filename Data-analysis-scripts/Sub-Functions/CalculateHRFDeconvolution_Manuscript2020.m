function [AnalysisResults] = CalculateHRFDeconvolution_Manuscript2020(animalID,neuralBand,hemisphere,behavior,saveFigs,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner - adapted from code written by Aaron T. Winder
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: calculate the hemodynamic response function from neural data
%________________________________________________________________________________________________________________________

%% setup parameters and load initial data structures
HRFLims = [0,5];
HRFParams.offset = 3;
HRFParams.dur = 10;
EventInds.FitStart = 1;
EventInds.TestStart = 2;
EventInds.Increment = 2;
modelType = 'Forest';
samplingRate = 30;
% load the Resting baselines structure
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% load Manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID)

%% load the data structure relevent to the analysis
if strcmp(behavior,'Rest')
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    BehData = RestData;
elseif strcmp(behavior,'Whisk') || strcmp(behavior,'Contra')
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
    BehData = EventData;
elseif strcmp(behavior,'NREM') == true || strcmp(behavior,'REM')
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
end

%% get the valid neural and hemodynamic data from the data structures
if strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true || strcmp(behavior,'Rest') == true
    % extract neural and hemodynamic data from event structure
    [NeuralDataStruct,NeuralFiltArray] = SelectConvolutionBehavioralEvents_IOS_Manuscript2020(BehData.(['cortical_' hemisphere(4:end)]).(neuralBand),behavior,hemisphere);
    [HemoDataStruct,HemoFiltArray] = SelectConvolutionBehavioralEvents_IOS_Manuscript2020(BehData.CBV_HbT.(hemisphere),behavior,hemisphere);
    % remove events that don't meet criteria
    [NeuralData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(NeuralDataStruct.NormData(NeuralFiltArray,:),NeuralDataStruct.fileIDs(NeuralFiltArray,:),NeuralDataStruct.duration(NeuralFiltArray,:),NeuralDataStruct.eventTime(NeuralFiltArray,:),ManualDecisions);
    [HemoData,~,~,~] = RemoveInvalidData_IOS_Manuscript2020(HemoDataStruct.data(HemoFiltArray,:),HemoDataStruct.fileIDs(HemoFiltArray,:),HemoDataStruct.duration(HemoFiltArray,:),HemoDataStruct.eventTime(HemoFiltArray,:),ManualDecisions);
elseif strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
    % extract neural and hemodynamic data from event structure
    NeuralData = SleepData.(modelType).(behavior).data.(['cortical_' hemisphere(4:end)]).(neuralBand);
    HemoData = SleepData.(modelType).(behavior).data.CBV_HbT.(hemisphere(4:end));
end

%% setup the neural data based on behavior
% Insert padding of zeros with size equal to the HRF between individual
fitIndex = EventInds.FitStart:EventInds.Increment:size(NeuralData,1);
zpadIndex = EventInds.TestStart:EventInds.Increment:size(NeuralData,1);
zpad = zeros(1,HRFParams.dur*samplingRate);
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true 
    % set every other event to zeros for padding
    NeuralData(zpadIndex) = {zpad};
    % process each event individually since they are all different lengths
    processedNeuralFitData = cell(size(NeuralData));
    for aa = 1:length(NeuralData)
        template = zeros(size(NeuralData{aa}));
        strt = 2*samplingRate;
        stp = size(template,2);
        template(:,strt:stp) = detrend(NeuralData{aa}(:,strt:stp) - mean(NeuralData{aa}(:,strt:stp)));
        processedNeuralFitData{aa} = template;
    end
    ProcNeuralFitData = [processedNeuralFitData{:}];
elseif strcmp(behavior,'Whisk') == true
    neuralDataEnd = 5;   % seconds
    strt = round((NeuralDataStruct.epoch.offset - 1)*samplingRate);
    stp = strt + round(neuralDataEnd*samplingRate);
    template = zeros(size(NeuralData(fitIndex,:)));
    neuralOffset = mean(NeuralData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    template(:,strt:stp) = NeuralData(fitIndex,strt:stp) - neuralOffset;
    % padding between events
    neuralDataPad = [template,ones(length(fitIndex),1)*zpad];
    ProcNeuralFitData = reshape(neuralDataPad',1,numel(neuralDataPad));
elseif strcmp(behavior,'Contra') == true
    neuralDataEnd = 1.5;   % seconds
    strt = round((NeuralDataStruct.epoch.offset)*NeuralDataStruct.samplingRate);
    stp = strt + (neuralDataEnd*samplingRate);
    template = zeros(size(NeuralData(fitIndex,:)));
    neuralOffset = mean(NeuralData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    template(:,strt:stp) = NeuralData(fitIndex,strt:stp) - neuralOffset;
    % padding between events
    neuralDataPad = [template,ones(length(fitIndex),1)*zpad];
    ProcNeuralFitData = reshape(neuralDataPad',1,numel(neuralDataPad));
end

%% setup the hemodynamic data based on behavior
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true 
    % set every other event to zeros for padding
    HemoData(zpadIndex) = {zpad};
    % process each event individually since they are all different lengths
    processedNeuralFitData = cell(size(HemoData));
    for bb = 1:length(HemoData)
        template = zeros(size(HemoData{bb}));
        strt = 2*samplingRate;
        stp = size(template,2);
        hemoOffset = mean(HemoData{bb})*ones(1,stp - strt + 1);
        template(:,strt:stp) = detrend(HemoData{bb}(:,strt:stp) - hemoOffset);
        processedNeuralFitData{bb} = template;
    end
    ProcHemoFitData = [processedNeuralFitData{:}];
elseif strcmp(behavior,'Whisk')
    hemoDataEnd = 7;   % seconds
    strt = round((HemoDataStruct.epoch.offset - 1)*samplingRate);
    stp = strt + round(hemoDataEnd*samplingRate);
    template = zeros(size(HemoData(fitIndex,:)));
    hemoOffset = mean(HemoData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    template(:,strt:stp) = HemoData(fitIndex,strt:stp) - hemoOffset;
    % padding between events
    hemoDataPad = [template,ones(length(fitIndex),1)*zpad];
    ProcHemoFitData = reshape(hemoDataPad',1,numel(hemoDataPad));
elseif strcmp(behavior,'Contra')
    hemoDataEnd = 3;   % seconds
    strt = round(HemoDataStruct.epoch.offset*samplingRate);
    stp = strt + (hemoDataEnd*samplingRate);
    template = zeros(size(HemoData(fitIndex,:)));
    hemoOffset = mean(HemoData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    template(:,strt:stp) = HemoData(fitIndex,strt:stp) - hemoOffset;
    % padding between events
    hemoDataPad = [template,ones(length(fitIndex),1)*zpad];
    ProcHemoFitData = reshape(hemoDataPad',1,numel(hemoDataPad));
end

%% calculate HRF based on deconvolution
IR_est = IR_analytic_IOS_Manuscript2020(ProcNeuralFitData',ProcHemoFitData',HRFParams.offset*samplingRate,HRFParams.dur*samplingRate);
HRF = sgolayfilt(IR_est.IR',3,samplingRate + 1);
timevec = (1:length(HRF))/samplingRate - HRFParams.offset;
numEvents = length(fitIndex);
timeLims = timevec >= HRFLims(1) & timevec <= HRFLims(2);
timeLimHRF = HRF(timeLims);
timeLimVec = timevec(timeLims);
% save results
AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).IR = timeLimHRF;
AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).IRtimeVec = timeLimVec;
AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).HRFParams = HRFParams;
AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).num_calc_events = numEvents;
AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).Event_Inds = EventInds;

%% calculate the gamma HRF
options = optimset('MaxFunEvals',2e4,'MaxIter',2e4,'TolFun',1e-7,'TolX',1e-7);
initvals = [1e-1,1,1];
HRFDur = 5; % seconds
[GamParams,~,~] = fminsearch(@(x)gammaconvolve_Manuscript2020(x,ProcNeuralFitData,ProcHemoFitData,samplingRate,HRFDur),initvals,options);



% save results
AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).gammaFunc = gamma;
AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).gammaAlpha = a;
AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).gammaBeta = beta;
AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).gammaTimeVec = t;
AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).GamParams = GamParams;

%% Save figures if desired
if strcmp(saveFigs,'y') == true
    kernelFig = figure;
    sgtitle([animalID ' ' hemisphere ' ' neuralBand ' during ' behavior])
    subplot(1,2,1)
    plot(AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).IRtimeVec,AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).IR,'k')
    title('Impulse response function')
    ylabel('A.U')
    xlabel('Time (s)')
    axis square
    axis tight
    subplot(1,2,2)
    plot(AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).gammaTimeVec,AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).gammaFunc,'k')
    title('Gamma function')
    ylabel('A.U')
    xlabel('Time (s)')
    axis square
    axis tight
    % save figures
    [pathstr,~,~] = fileparts(cd);
    dirpath = [pathstr '/Figures/HRF Kernels/'];
    if ~exist(dirpath, 'dir')
        mkdir(dirpath);
    end
    savefig(kernelFig,[dirpath animalID '_' hemisphere '_' neuralBand '_' behavior '_HRFs']);
    close(kernelFig)
end

end
