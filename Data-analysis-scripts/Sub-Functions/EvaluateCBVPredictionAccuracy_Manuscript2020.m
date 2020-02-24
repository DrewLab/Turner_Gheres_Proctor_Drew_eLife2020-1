function [AnalysisResults] = EvaluateCBVPredictionAccuracy_Manuscript2020(animalID,neuralBand,hemisphere,behavior,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: 
%________________________________________________________________________________________________________________________

%% Setup
Event_Inds.CalcStart = 1;
Event_Inds.TestStart = 2;
Event_Inds.Increment = 2;
modelType = 'SVM';
baselineDataFileStruct = dir('*_RestingBaselines.mat');
baselineDataFile = {baselineDataFileStruct.name}';
baselineDataFileID = char(baselineDataFile);
load(baselineDataFileID)
% find and load Manual baseline event information
manualBaselineFileStruct = dir('*_ManualBaselineFileList.mat');
manualBaselineFile = {manualBaselineFileStruct.name}';
manualBaselineFileID = char(manualBaselineFile);
load(manualBaselineFileID)

if strcmp(behavior,'Rest')
    restDataFileStruct = dir('*_RestData.mat');
    restDataFile = {restDataFileStruct.name}';
    restDataFileID = char(restDataFile);
    load(restDataFileID)
    BehData = RestData;
    clear RestData;
elseif strcmp(behavior,'Whisk') || strcmp(behavior,'Contra')
    eventDataFileStruct = dir('*_EventData.mat');
    eventDataFile = {eventDataFileStruct.name}';
    eventDataFileID = char(eventDataFile);
    load(eventDataFileID)
    BehData = EventData;
    clear EventData;
elseif strcmp(behavior,'NREM') == true || strcmp(behavior,'REM')
    sleepDataFileStruct = dir('*_SleepData.mat');
    sleepDataFile = {sleepDataFileStruct.name}';
    sleepDataFileID = char(sleepDataFile);
    load(sleepDataFileID)
end

%% Get the arrays for the calculation
if strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true || strcmp(behavior,'Rest') == true
    [NeuralDataStruct,NeuralFiltArray] = SelectConvolutionBehavioralEvents_IOS_Manuscript2020(BehData.(['cortical_' hemisphere(4:end)]).(neuralBand),behavior,hemisphere);
    [HemoDataStruct,HemoFiltArray] = SelectConvolutionBehavioralEvents_IOS_Manuscript2020(BehData.CBV_HbT.(hemisphere),behavior,hemisphere);
     % remove events that don't meet criteria
    [NormData1,~,~,~] = DecimateRestData_Manuscript2020(NeuralDataStruct.NormData(NeuralFiltArray,:),NeuralDataStruct.fileIDs(NeuralFiltArray,:),NeuralDataStruct.duration(NeuralFiltArray,:),NeuralDataStruct.eventTime(NeuralFiltArray,:),ManualDecisions);
    [NormData2,~,~,~] = DecimateRestData_Manuscript2020(HemoDataStruct.data(HemoFiltArray,:),HemoDataStruct.fileIDs(HemoFiltArray,:),HemoDataStruct.duration(HemoFiltArray,:),HemoDataStruct.eventTime(HemoFiltArray,:),ManualDecisions);
elseif strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true
    NormData1 = SleepData.(modelType).(behavior).data.(['cortical_' hemisphere(4:end)]).(neuralBand);
    NormData2 = SleepData.(modelType).(behavior).data.CBV_HbT.(hemisphere(4:end));
    NeuralDataStruct.samplingRate = 30;
    HemoDataStruct.samplingRate = 30;
end

%% Setup the data
test_inds = Event_Inds.TestStart:Event_Inds.Increment:size(NormData1,1);
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true 
    NormData1 = NormData1(test_inds);
    % Mean subtract the data
    Processed1 = cell(size(NormData1));
    for c = 1:length(NormData1)
        template = zeros(size(NormData1{c}));
        strt = 2*NeuralDataStruct.samplingRate;
        stp = size(template,2);
        template(:,strt:stp) = NormData1{c}(:,strt:stp) - mean(NormData1{c}(:,strt:stp));
        Processed1{c} = template;
    end
    Data1 = Processed1;
    clear Processed1;
elseif strcmp(behavior,'Whisk')
    Data1_end = 5; 
    strt = round((NeuralDataStruct.epoch.offset - 1)*NeuralDataStruct.samplingRate);
    stp = strt + round(Data1_end*NeuralDataStruct.samplingRate);
    Data1 = zeros(size(NormData1(test_inds,:)));
    offset1 = mean(NormData1(test_inds,1:strt),2)*ones(1,stp - strt+1);
    Data1(:,strt:stp) = NormData1(test_inds,strt:stp) - offset1;
elseif strcmp(behavior,'Contra')
    Data1_end = 1.5; 
    strt = round((NeuralDataStruct.epoch.offset)*NeuralDataStruct.samplingRate); 
    stp = strt + (Data1_end*NeuralDataStruct.samplingRate);
    Data1 = zeros(size(NormData1(test_inds,:)));
    offset1 = mean(NormData1(test_inds,1:strt),2)*ones(1,stp-strt+1);
    Data1(:,strt:stp) = NormData1(test_inds,strt:stp) - offset1;
end

%%
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true 
    NormData2 = NormData2(test_inds);
    % Mean subtract the data
    Processed2 = cell(size(NormData2));
    for c = 1:length(NormData2)
        template = zeros(size(NormData2{c}));
        strt = 2*HemoDataStruct.samplingRate;
        stp = size(template,2);
        offset = mean(NormData2{c})*ones(1,stp - strt+1);
        template(:,strt:stp) = detrend(NormData2{c}(:,strt:stp) - offset);
        Processed2{c} = template;
    end
    Data2 = Processed2;
    clear Processed2
elseif strcmp(behavior,'Whisk')
    Data2_end = 7;
    strt = round((HemoDataStruct.epoch.offset - 1)*HemoDataStruct.samplingRate);
    stp = strt + round(Data2_end*HemoDataStruct.samplingRate);
    Data2 = zeros(size(NormData2(test_inds,:)));
    offset2 = mean(NormData2(test_inds,1:strt),2)*ones(1,stp - strt+1);
    Data2(:,strt:stp) = NormData2(test_inds,strt:stp) - offset2;
elseif strcmp(behavior,'Contra')
    Data2_end = 3;
    strt = round(HemoDataStruct.epoch.offset*HemoDataStruct.samplingRate);
    stp = strt + (Data2_end*HemoDataStruct.samplingRate);
    Data2 = zeros(size(NormData2(test_inds,:)));
    offset2 = mean(NormData2(test_inds,1:strt),2)*ones(1,stp - strt+1);
    Data2(:,strt:stp) = NormData2(test_inds,strt:stp) - offset2;
end

%% Calculate R-squared on average data
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true 
    AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).AveR2 = NaN;
else
    [Act,Pred] = ConvolveHRF_IOS_Manuscript2020(AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).Contra.gammaFunc,mean(Data1),mean(Data2),0);
    mPred = Pred(strt:stp) - mean(Pred(strt:stp));
    mAct = Act(strt:stp) - mean(Act(strt:stp));
    AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).AveR2 = CalculateRsquared_IOS_Manuscript2020(mPred,mAct);
end

%% Calculate R-squared on individual data
IndR2 = NaN*ones(1,size(Data2,1));
if strcmp(behavior,'Rest') == true || strcmp(behavior,'NREM') == true || strcmp(behavior,'REM') == true 
    for tc = 1:length(Data2)
        strt = 2*HemoDataStruct.samplingRate;
        stp = length(Data2{tc});
        % convolve the contra HRF with all behaviors
        [Act,Pred] = ConvolveHRF_IOS_Manuscript2020(AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).Contra.gammaFunc,detrend(Data1{tc}),detrend(Data2{tc}),0);
        mPred = Pred(strt:stp) - mean(Pred(strt:stp));
        mAct = Act(strt:stp) - mean(Act(strt:stp));
        IndR2(tc) = CalculateRsquared_IOS_Manuscript2020(mPred,mAct);
    end
    AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).Mean_IndR2 = mean(IndR2);
    AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).Med_IndR2 = median(IndR2);
elseif strcmp(behavior,'Contra') == true || strcmp(behavior,'Whisk') == true
    for tc = 1:size(Data2,1)
        [Act,Pred] = ConvolveHRF_IOS_Manuscript2020(AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).Contra.gammaFunc,Data1(tc,:),Data2(tc,:),0);
        mPred = Pred(strt:stp) - mean(Pred(strt:stp));
        mAct = Act(strt:stp) - mean(Act(strt:stp));
        IndR2(tc) = CalculateRsquared_IOS_Manuscript2020(mPred,mAct);
    end
    AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).Mean_IndR2 = mean(IndR2);
    AnalysisResults.(animalID).HRFs.(neuralBand).(hemisphere).(behavior).Med_IndR2 = median(IndR2);
end

end
