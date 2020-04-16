function [AnalysisResults] = EvaluateCBVPredictionAccuracy_Manuscript2020(animalID,neuralBand,hemisphere,behavior,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner - adapted from code written by Aaron T. Winder
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose: convolve HRF with neural activity and estimate accuracy of predictions via R2
%________________________________________________________________________________________________________________________

%% setup parameters and load initial data structures
EventIndex.FitStart = 1;
EventIndex.TestStart = 2;
EventIndex.Increment = 2;
modelType = 'Forest';
samplingRate = 30;
% filter characteristics
[z,p,k] = butter(4,1/(samplingRate/2),'low');
[sos,g] = zp2sos(z,p,k);
% load the Resting baselines structure
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% find and load Manual baseline event information
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
fitIndex = EventIndex.FitStart:EventIndex.Increment:size(NeuralData,1);
testIndex = EventIndex.TestStart:EventIndex.Increment:size(NeuralData,1);
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true 
    % fit data
    neuralFitData = NeuralData(fitIndex);
    processedNeuralFitData = cell(size(neuralFitData));
    for aa = 1:length(neuralFitData)
        template = zeros(size(neuralFitData{aa}));
        strt = 2*samplingRate;
        stp = size(template,2);
        template(:,strt:stp) = neuralFitData{aa}(:,strt:stp) - mean(neuralFitData{aa}(:,strt:stp));
        processedNeuralFitData{aa} = template;
    end
    ProcNeuralFitData = processedNeuralFitData;
    % test data
    neuralTestData = NeuralData(testIndex);
    processedNeuralTestData = cell(size(neuralTestData));
    for bb = 1:length(neuralTestData)
        template = zeros(size(neuralTestData{bb}));
        strt = 2*samplingRate;
        stp = size(template,2);
        template(:,strt:stp) = neuralTestData{bb}(:,strt:stp) - mean(neuralTestData{bb}(:,strt:stp));
        processedNeuralTestData{bb} = template;
    end
    ProcNeuralTestData = processedNeuralTestData;
elseif strcmp(behavior,'Whisk') == true
    neuralDataEnd = 5;   % seconds
    strt = round((NeuralDataStruct.epoch.offset - 1)*samplingRate);
    stp = strt + round(neuralDataEnd*samplingRate);
    % fit data
    ProcNeuralFitData = zeros(size(NeuralData(fitIndex,:)));
    fitNeuralOffset = mean(NeuralData(fitIndex,1:strt),2)*ones(1,stp - strt+1);
    ProcNeuralFitData(:,strt:stp) = NeuralData(fitIndex,strt:stp) - fitNeuralOffset;
    % test data
    ProcNeuralTestData = zeros(size(NeuralData(testIndex,:)));
    testNeuralOffset = mean(NeuralData(testIndex,1:strt),2)*ones(1,stp - strt+1);
    ProcNeuralTestData(:,strt:stp) = NeuralData(testIndex,strt:stp) - testNeuralOffset;
elseif strcmp(behavior,'Contra') == true
    neuralDataEnd = 1.5;   % seconds 
    strt = round((NeuralDataStruct.epoch.offset)*samplingRate); 
    stp = strt + (neuralDataEnd*samplingRate);
    % fit data
    ProcNeuralFitData = zeros(size(NeuralData(fitIndex,:)));
    fitNeuralOffset = mean(NeuralData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    ProcNeuralFitData(:,strt:stp) = NeuralData(fitIndex,strt:stp) - fitNeuralOffset;
    % test data
    ProcNeuralTestData = zeros(size(NeuralData(testIndex,:)));
    testNeuralOffset = mean(NeuralData(testIndex,1:strt),2)*ones(1,stp - strt + 1);
    ProcNeuralTestData(:,strt:stp) = NeuralData(testIndex,strt:stp) - testNeuralOffset;
end

%% setup the hemodynamic data based on behavior
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true 
    % fit data
    fitHemoData = HemoData(fitIndex);
    processedHemoFitData = cell(size(fitHemoData));
    for cc = 1:length(fitHemoData)
        template = zeros(size(fitHemoData{cc}));
        strt = 2*samplingRate;
        stp = size(template,2);
        offset = mean(fitHemoData{cc})*ones(1,stp - strt+1);
        template(:,strt:stp) = detrend(fitHemoData{cc}(:,strt:stp) - offset);
        processedHemoFitData{cc} = template;
    end
    ProcHemoFitData = processedHemoFitData;
    % test data
    testHemoData = HemoData(testIndex);
    processedHemoTestData = cell(size(testHemoData));
    for dd = 1:length(testHemoData)
        template = zeros(size(testHemoData{dd}));
        strt = 2*samplingRate;
        stp = size(template,2);
        offset = mean(testHemoData{dd})*ones(1,stp - strt+1);
        template(:,strt:stp) = detrend(testHemoData{dd}(:,strt:stp) - offset);
        processedHemoTestData{dd} = template;
    end
    ProcHemoTestData = processedHemoTestData;
elseif strcmp(behavior,'Whisk')
    hemoDataEnd = 7;   % seconds
    strt = round((HemoDataStruct.epoch.offset - 1)*samplingRate);
    stp = strt + round(hemoDataEnd*samplingRate);
    % fit data
    ProcHemoFitData = zeros(size(HemoData(fitIndex,:)));
    fitHemoOffset = mean(HemoData(fitIndex,1:strt),2)*ones(1,stp - strt+1);
    ProcHemoFitData(:,strt:stp) = HemoData(fitIndex,strt:stp) - fitHemoOffset;
    % test data
    ProcHemoTestData = zeros(size(HemoData(testIndex,:)));
    testHemoOffset = mean(HemoData(testIndex,1:strt),2)*ones(1,stp - strt+1);
    ProcHemoTestData(:,strt:stp) = HemoData(testIndex,strt:stp) - testHemoOffset;
elseif strcmp(behavior,'Contra')
    hemoDataEnd = 3;   % seconds
    strt = round(HemoDataStruct.epoch.offset*samplingRate);
    stp = strt + (hemoDataEnd*samplingRate);
    % fit data
    ProcHemoFitData = zeros(size(HemoData(fitIndex,:)));
    fitHemoOffset = mean(HemoData(fitIndex,1:strt),2)*ones(1,stp - strt + 1);
    ProcHemoFitData(:,strt:stp) = HemoData(fitIndex,strt:stp) - fitHemoOffset;
    % test data
    ProcHemoTestData = zeros(size(HemoData(testIndex,:)));
    testHemoOffset = mean(HemoData(testIndex,1:strt),2)*ones(1,stp - strt + 1);
    ProcHemoTestData(:,strt:stp) = HemoData(testIndex,strt:stp) - testHemoOffset;
end

%% Calculate R-squared on individual data
IndR2_A = NaN*ones(1,size(ProcHemoTestData,1));
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
    % fit data
    for ee = 1:length(ProcHemoFitData)
        strt = 2*samplingRate;
        stp = length(ProcHemoFitData{ee});
        [Act,Pred] = ConvolveHRF_IOS_Manuscript2020(AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).gammaFunc,detrend(ProcNeuralFitData{ee}),detrend(ProcHemoFitData{ee}),0);
        mPred = Pred(strt:stp) - mean(Pred(strt:stp));
        mAct = Act(strt:stp) - mean(Act(strt:stp));
        FitData.mPred{ee} = mPred;
        FitData.mAct{ee} = mAct;
    end
    % test data
    for ff = 1:length(ProcHemoTestData)
        strt = 2*samplingRate;
        stp = length(ProcHemoTestData{ff});
        [Act,Pred] = ConvolveHRF_IOS_Manuscript2020(AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).gammaFunc,detrend(ProcNeuralTestData{ff}),detrend(ProcHemoTestData{ff}),0);
        mPred = Pred(strt:stp) - mean(Pred(strt:stp));
        mAct = Act(strt:stp) - mean(Act(strt:stp));
        TestData.mPred{ff} = mPred;
        TestData.mAct{ff} = mAct;
    end
    % R2 without correction
    for gg = 1:length(TestData.mPred)
        IndR2_A(gg) = CalculateRsquared_IOS_Manuscript2020(filtfilt(sos,g,TestData.mPred{1,gg}),filtfilt(sos,g,TestData.mAct{1,gg}));
    end
%     uncorrectedMedian = median(IndR2_A);
%     disp(['Uncorrected ' neuralBand ' ' behavior ' Median R2: ' num2str(uncorrectedMedian)]);
%     uncorrectedMean = mean(IndR2_A);
%     disp(['Uncorrected ' neuralBand ' ' behavior ' Mean R2: ' num2str(uncorrectedMean)]);
%     % correct test data using fit data
%     [CorrectedTestData,FitParams] = CorrectKernelError_IOS_Manuscript2020(FitData,TestData,behavior);
%     for gg = 1:length(CorrectedTestData)
%         IndR2_B(gg) = CalculateRsquared_IOS_Manuscript2020(filtfilt(sos,g,CorrectedTestData{gg,1}),filtfilt(sos,g,TestData.mAct{1,gg}));
%     end
%     correctedMedian = median(IndR2_B);
%     disp(['Corrected ' neuralBand ' ' behavior ' Median R2: ' num2str(correctedMedian)]);
%     correctedMean = mean(IndR2_B);
%     disp(['Corrected ' neuralBand ' ' behavior ' Mean R2: ' num2str(correctedMean)]);
    % save results
    AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).meanR2 = mean(IndR2_A);
    AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).medR2 = median(IndR2_A);
%     AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).FitParams = FitParams;
elseif strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true
    % fit data
    for hh = 1:size(ProcHemoFitData,1)
        [Act,Pred] = ConvolveHRF_IOS_Manuscript2020(AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).gammaFunc,ProcNeuralFitData(hh,:),ProcHemoFitData(hh,:),0);
        mPred = Pred(strt:stp) - mean(Pred(strt:stp));
        mAct = Act(strt:stp) - mean(Act(strt:stp));
        FitData.mPred{hh,1} = mPred;
        FitData.mAct{hh,1} = mAct;
    end
    % test data
    for ii = 1:size(ProcHemoTestData,1)
        [Act,Pred] = ConvolveHRF_IOS_Manuscript2020(AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).gammaFunc,ProcNeuralTestData(ii,:),ProcHemoTestData(ii,:),0);
        mPred = Pred(strt:stp) - mean(Pred(strt:stp));
        mAct = Act(strt:stp) - mean(Act(strt:stp));
        TestData.mPred{ii,1} = mPred;
        TestData.mAct{ii,1} = mAct;
    end
    % R2 without correction
    for jj = 1:length(TestData.mPred)
        IndR2_A(jj) = CalculateRsquared_IOS_Manuscript2020(filtfilt(sos,g,TestData.mPred{jj,1}),filtfilt(sos,g,TestData.mAct{jj,1}));
    end
%     uncorrectedMedian = median(IndR2_A);
%     disp(['Uncorrected ' neuralBand ' ' behavior ' Median R2: ' num2str(uncorrectedMedian)]);
%     uncorrectedMean = mean(IndR2_A);
%     disp(['Uncorrected ' neuralBand ' ' behavior ' Mean R2: ' num2str(uncorrectedMean)]);
%     % correct test data using fit data
%     [CorrectedTestData,FitParams] = CorrectKernelError_IOS_Manuscript2020(FitData,TestData,behavior);
%     for jj = 1:size(CorrectedTestData,1)
%         IndR2_B(jj) = CalculateRsquared_IOS_Manuscript2020(filtfilt(sos,g,CorrectedTestData(jj,:)),filtfilt(sos,g,TestData.mAct{jj,1}));
%     end
%     correctedMedian = median(IndR2_B);
%     disp(['Corrected ' neuralBand ' ' behavior ' Median R2: ' num2str(correctedMedian)]);
%     correctedMean = mean(IndR2_B);
%     disp(['Corrected ' neuralBand ' ' behavior ' Mean R2: ' num2str(correctedMean)]);
    % save results
    AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).meanR2 = mean(IndR2_A);
    AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).medR2 = median(IndR2_A);
%     AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).FitParams = FitParams;
end

end
