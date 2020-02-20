function [] = TrainSleepModels_IOS_Manuscript2020(animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Train several machine learning techniques on manually scored sleep data, and evaluate each model's accuracy
%________________________________________________________________________________________________________________________

%% temp
for a = 1:length(animalIDs)
    startingDirectory = cd;
    modelDirectory = [animalIDs{1,a} '\Bilateral Imaging\'];
    cd(modelDirectory)
    % character list of all training files
    modelDataFileStruct = dir('*_ModelData.mat');
    modelDataFiles = {modelDataFileStruct.name}';
    modelDataFileIDs = char(modelDataFiles);
    % Load each updated training set and concatenate the data into table
    for b = 1:size(modelDataFileIDs,1)
        modelTableFileID = modelDataFileIDs(b,:);
        if a == 1 && b == 1
            load(modelTableFileID)
            joinedTable = paramsTable;
        else
            load(modelTableFileID)
            joinedTable = vertcat(joinedTable,paramsTable);
        end
    end
    cd(startingDirectory)
end

%% load in all the data to create a table of values
for a = 1:length(animalIDs)
    startingDirectory = cd;
    trainingDirectory = [animalIDs{1,a} '\Training Data\'];
    cd(trainingDirectory)
    % character list of all training files
    trainingDataFileStruct = dir('*_TrainingData.mat');
    trainingDataFiles = {trainingDataFileStruct.name}';
    trainingDataFileIDs = char(trainingDataFiles);
    % Load each updated training set and concatenate the data into table
    for b = 1:size(trainingDataFileIDs,1)
        trainingTableFileID = trainingDataFileIDs(b,:);
        if a == 1 && b == 1
            load(trainingTableFileID)
            dataLength = size(trainingTable,1);
            joinedTableOdd = trainingTable;
        elseif a == 1 && b == 2
            load(trainingTableFileID)
            joinedTableEven = trainingTable;
        elseif rem(b,2) == 1
            load(trainingTableFileID)
            joinedTableOdd = vertcat(joinedTableOdd,trainingTable);
        elseif rem(b,2) == 0
            load(trainingTableFileID)
            joinedTableEven = vertcat(joinedTableEven,trainingTable);
        end
    end
    cd(startingDirectory)
end
% train on odd data
Xodd = joinedTableOdd(:,1:end - 1);
Yodd = joinedTableOdd(:,end);
% test on even data
Xeven = joinedTableEven(:,1:end - 1);
Yeven = joinedTableEven(:,end);

%% Train Support Vector Machine (SVM) classifier
t = templateSVM('Standardize',true,'KernelFunction','gaussian');
disp('Training Support Vector Machine...'); disp(' ')
SVM_MDL = fitcecoc(Xodd,Yodd,'Learners',t,'FitPosterior',true,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'},'Verbose',2);
% determine k-fold loss of the model
disp('Cross-validating (10-fold) the support vector machine classifier...'); disp(' ')
CV_SVM_MDL = crossval(SVM_MDL);
loss = kfoldLoss(CV_SVM_MDL);
disp(['k-fold loss classification error: ' num2str(loss*100) '%']); disp(' ')
% save model and cross validation in desired location
saveLoc = [startingDirectory '\Support Files\'];
save([saveLoc '\IOS_SVM_SleepScoringModel.mat\'],'SVM_MDL')
save([saveLoc '\IOS_SVM_ModelCrossValidation.mat\'],'CV_SVM_MDL')
% use the model to generate a set of scores for the even set of data
[XevenLabels,~] = predict(SVM_MDL,Xeven);
% apply a logical patch on the REM events
REMindex = strcmp(XevenLabels,'REM Sleep');
numFiles = length(XevenLabels)/dataLength;
reshapedREMindex = reshape(REMindex,dataLength,numFiles);
patchedREMindex = [];
% patch missing REM indeces due to theta band falling off
for c = 1:size(reshapedREMindex,2)
    remArray = reshapedREMindex(:,c);
    patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
    patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
end
% change labels for each event
for d = 1:length(XevenLabels)
    if patchedREMindex(d,1) == 1
        XevenLabels{d,1} = 'REM Sleep';
    end
end
% confusion matrix
confMat = figure;
cm = confusionchart(Yeven.behavState,XevenLabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
cm.Title = 'Support Vector Machine Classifier Confusion Matrix';
% Pull data out of confusion matrix
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
SVM_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
disp(['SVM model prediction accuracy: ' num2str(SVM_accuracy) '%']); disp(' ')
% save figure
dirpath = [startingDirectory '\Support Files\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(confMat,[dirpath 'IOS_SVM_ConfusionMatrix']);
close(confMat)

%% Ensemble classification - AdaBoostM2, Subspace, Bag, LPBoost,RUSBoost, TotalBoost
disp('Training Ensemble Classifier...'); disp(' ')
t = templateTree('Reproducible',true);
EC_MDL = fitcensemble(Xodd,Yodd,'OptimizeHyperparameters','auto','Learners',t,'HyperparameterOptimizationOptions',...
    struct('AcquisitionFunctionName','expected-improvement-plus'),'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
% use the model to generate a set of scores for the even set of data
[XevenLabels,~] = predict(EC_MDL,Xeven);
% apply a logical patch on the REM events
REMindex = strcmp(XevenLabels,'REM Sleep');
reshapedREMindex = reshape(REMindex,dataLength,numFiles);
patchedREMindex = [];
% patch missing REM indeces due to theta band falling off
for c = 1:size(reshapedREMindex,2)
    remArray = reshapedREMindex(:,c);
    patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
    patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
end
% change labels for each event
for d = 1:length(XevenLabels)
    if patchedREMindex(d,1) == 1
        XevenLabels{d,1} = 'REM Sleep';
    end
end
% confusion matrix
EC_confMat = figure;
cm = confusionchart(Yeven.behavState,XevenLabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
cm.Title = 'Ensemble Classifier Confusion Matrix';
% Pull data out of confusion matrix
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
EC_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
disp(['Ensemble model prediction accuracy: ' num2str(EC_accuracy) '%']); disp(' ')
% save model and figure
save([saveLoc '\IOS_EC_SleepScoringModel.mat\'],'EC_MDL')
savefig(EC_confMat,[dirpath 'IOS_EC_ConfusionMatrix']);
close(EC_confMat)

%% Decision Tree classification 
disp('Training Decision Tree Classifier...'); disp(' ')
DT_MDL = fitctree(Xodd,Yodd,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
% use the model to generate a set of scores for the even set of data
[XevenLabels,~] = predict(DT_MDL,Xeven);
% apply a logical patch on the REM events
REMindex = strcmp(XevenLabels,'REM Sleep');
reshapedREMindex = reshape(REMindex,dataLength,numFiles);
patchedREMindex = [];
% patch missing REM indeces due to theta band falling off
for c = 1:size(reshapedREMindex,2)
    remArray = reshapedREMindex(:,c);
    patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
    patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
end
% change labels for each event
for d = 1:length(XevenLabels)
    if patchedREMindex(d,1) == 1
        XevenLabels{d,1} = 'REM Sleep';
    end
end
% confusion matrix
DT_confMat = figure;
cm = confusionchart(Yeven.behavState,XevenLabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
cm.Title = 'Decision Tree Classifier Confusion Matrix';
% Pull data out of confusion matrix
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
DT_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
disp(['Decision Tree model prediction accuracy: ' num2str(DT_accuracy) '%']); disp(' ')
% save model and figure
save([saveLoc '\IOS_DT_SleepScoringModel.mat\'],'DT_MDL')
savefig(DT_confMat,[dirpath 'IOS_DT_ConfusionMatrix']);
close(DT_confMat)

%% Random forest
disp('Training Random Forest Classifier...'); disp(' ')
numTrees = 128;
RF_MDL = TreeBagger(numTrees,Xodd,Yodd,'Method','Classification','Surrogate','all','ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
% use the model to generate a set of scores for the even set of data
[XevenLabels,~] = predict(RF_MDL,Xeven);
% apply a logical patch on the REM events
REMindex = strcmp(XevenLabels,'REM Sleep');
reshapedREMindex = reshape(REMindex,dataLength,numFiles);
patchedREMindex = [];
% patch missing REM indeces due to theta band falling off
for c = 1:size(reshapedREMindex,2)
    remArray = reshapedREMindex(:,c);
    patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
    patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
end
% change labels for each event
for d = 1:length(XevenLabels)
    if patchedREMindex(d,1) == 1
        XevenLabels{d,1} = 'REM Sleep';
    end
end
% confusion matrix
RF_confMat = figure;
cm = confusionchart(Yeven.behavState,XevenLabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
cm.Title = 'Random Forest Classifier Confusion Matrix';
% Pull data out of confusion matrix
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
RF_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
disp(['Random Forest model prediction accuracy: ' num2str(RF_accuracy) '%']); disp(' ')
% save model and figure
save([saveLoc '\IOS_RF_SleepScoringModel.mat\'],'RF_MDL')
savefig(RF_confMat,[dirpath 'IOS_RF_ConfusionMatrix']);
close(RF_confMat)

%% k-nearest neighbor classifier
disp('Training k-nearest neighbor Classifier...'); disp(' ')
t = templateKNN('NumNeighbors',5,'Standardize',1);
KNN_MDL = fitcecoc(Xodd,Yodd,'Learners',t);
% use the model to generate a set of scores for the even set of data
[XevenLabels,~] = predict(KNN_MDL,Xeven);
% apply a logical patch on the REM events
REMindex = strcmp(XevenLabels,'REM Sleep');
reshapedREMindex = reshape(REMindex,dataLength,numFiles);
patchedREMindex = [];
% patch missing REM indeces due to theta band falling off
for c = 1:size(reshapedREMindex,2)
    remArray = reshapedREMindex(:,c);
    patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
    patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
end
% change labels for each event
for d = 1:length(XevenLabels)
    if patchedREMindex(d,1) == 1
        XevenLabels{d,1} = 'REM Sleep';
    end
end
% confusion matrix
KNN_confMat = figure;
cm = confusionchart(Yeven.behavState,XevenLabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
cm.Title = 'KNN Classifier Confusion Matrix';
% Pull data out of confusion matrix
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
KNN_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
disp(['K-nearest neighbor model prediction accuracy: ' num2str(KNN_accuracy) '%']); disp(' ')
% save model and figure
save([saveLoc '\IOS_KNN_SleepScoringModel.mat\'],'KNN_MDL')
savefig(KNN_confMat,[dirpath 'IOS_KNN_ConfusionMatrix']);
close(KNN_confMat)

%% Naive Bayes classifier
disp('Training naive Bayes Classifier...'); disp(' ')
NB_MDL = fitcnb(Xodd,Yodd,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
[XevenLabels,~] = predict(NB_MDL,Xeven);
% apply a logical patch on the REM events
REMindex = strcmp(XevenLabels,'REM Sleep');
reshapedREMindex = reshape(REMindex,dataLength,numFiles);
patchedREMindex = [];
% patch missing REM indeces due to theta band falling off
for c = 1:size(reshapedREMindex,2)
    remArray = reshapedREMindex(:,c);
    patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
    patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
end
% change labels for each event
for d = 1:length(XevenLabels)
    if patchedREMindex(d,1) == 1
        XevenLabels{d,1} = 'REM Sleep';
    end
end
% confusion matrix
NB_confMat = figure;
cm = confusionchart(Yeven.behavState,XevenLabels);
cm.ColumnSummary = 'column-normalized';
cm.RowSummary = 'row-normalized';
cm.Title = 'Naive Bayes Classifier Confusion Matrix';
% Pull data out of confusion matrix
confVals = cm.NormalizedValues;
totalScores = sum(confVals(:));
NB_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
disp(['Naive Bayes model prediction accuracy: ' num2str(NB_accuracy) '%']); disp(' ')
% save model and figure
save([saveLoc '\IOS_NB_SleepScoringModel.mat\'],'NB_MDL')
savefig(NB_confMat,[dirpath 'IOS_NB_ConfusionMatrix']);
close(NB_confMat)

%% Restatement of accuracies for viewing
disp(['SVM model prediction accuracy: ' num2str(SVM_accuracy) '%']); disp(' ')
disp(['Ensemble model prediction accuracy: ' num2str(EC_accuracy) '%']); disp(' ')
disp(['Decision Tree model prediction accuracy: ' num2str(DT_accuracy) '%']); disp(' ')
disp(['Random Forest model prediction accuracy: ' num2str(RF_accuracy) '%']); disp(' ')
disp(['K-nearest neighbor model prediction accuracy: ' num2str(KNN_accuracy) '%']); disp(' ')
disp(['Naive Bayes model prediction accuracy: ' num2str(NB_accuracy) '%']); disp(' ')

end
