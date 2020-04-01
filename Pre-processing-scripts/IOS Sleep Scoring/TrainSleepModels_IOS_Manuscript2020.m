function [] = TrainSleepModels_IOS_Manuscript2020(animalIDs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Train several machine learning techniques on manually scored sleep data, and evaluate each model's accuracy
%________________________________________________________________________________________________________________________

%% load in all the data to create a table of values
for aa = 1:length(animalIDs)
    clearvars -except animalIDs aa ConfusionData
    startingDirectory = cd;
    trainingDirectory = [animalIDs{1,aa} '\Training Data\'];
    cd(trainingDirectory)
    % character list of all training files
    trainingDataFileStruct = dir('*_TrainingData.mat');
    trainingDataFiles = {trainingDataFileStruct.name}';
    trainingDataFileIDs = char(trainingDataFiles);
    % Load each updated training set and concatenate the data into table
    for bb = 1:size(trainingDataFileIDs,1)
        trainingTableFileID = trainingDataFileIDs(bb,:);
        if bb == 1
            load(trainingTableFileID)
            dataLength = size(trainingTable,1);
            joinedTableOdd = trainingTable;
        elseif bb == 2
            load(trainingTableFileID)
            joinedTableEven = trainingTable;
        elseif rem(bb,2) == 1
            load(trainingTableFileID)
            joinedTableOdd = vertcat(joinedTableOdd,trainingTable);
        elseif rem(bb,2) == 0
            load(trainingTableFileID)
            joinedTableEven = vertcat(joinedTableEven,trainingTable);
        end
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
    % save model in desired location
    dirpath = [startingDirectory '\' animalIDs{1,aa} '\Figures\Sleep Models\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    save([dirpath animalIDs{1,aa} '_IOS_SVM_SleepScoringModel.mat'],'SVM_MDL')
    % determine k-fold loss of the model
    disp('Cross-validating (3-fold) the support vector machine classifier...'); disp(' ')
    CV_SVM_MDL = crossval(SVM_MDL,'kfold',3);
    loss = kfoldLoss(CV_SVM_MDL);
    disp(['k-fold loss classification error: ' num2str(loss*100) '%']); disp(' ')
    % save cross validation in desired location
    save([dirpath animalIDs{1,aa} '_IOS_SVM_ModelCrossValidation.mat'],'CV_SVM_MDL')
    % use the model to generate a set of scores for the even set of data
    [XevenLabels,~] = predict(SVM_MDL,Xeven);
    % apply a logical patch on the REM events
    REMindex = strcmp(XevenLabels,'REM Sleep');
    numFiles = length(XevenLabels)/dataLength;
    reshapedREMindex = reshape(REMindex,dataLength,numFiles);
    patchedREMindex = [];
    % patch missing REM indeces due to theta band falling off
    for cc = 1:size(reshapedREMindex,2)
        remArray = reshapedREMindex(:,cc);
        patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
        patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
    end
    % change labels for each event
    for dd = 1:length(XevenLabels)
        if patchedREMindex(dd,1) == 1
            XevenLabels{dd,1} = 'REM Sleep';
        end
    end
    ConfusionData.SVM.Ylabels{aa,1} = Yeven.behavState;
    ConfusionData.SVM.Xlabels{aa,1} = XevenLabels;
    % confusion matrix
    confMat = figure;
    cm = confusionchart(Yeven.behavState,XevenLabels);
    cm.ColumnSummary = 'column-normalized';
    cm.RowSummary = 'row-normalized';
    cm.Title = 'Support Vector Machine Classifier Confusion Matrix';
    % pull data out of confusion matrix
    confVals = cm.NormalizedValues;
    totalScores = sum(confVals(:));
    SVM_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
    disp(['SVM model prediction accuracy: ' num2str(SVM_accuracy) '%']); disp(' ')
    savefig(confMat,[dirpath animalIDs{1,aa} '_IOS_SVM_ConfusionMatrix']);
    close(confMat)
    
    %% Ensemble classification - AdaBoostM2, Subspace, Bag, LPBoost,RUSBoost, TotalBoost
    disp('Training Ensemble Classifier...'); disp(' ')
    t = templateTree('Reproducible',true);
    EC_MDL = fitcensemble(Xodd,Yodd,'OptimizeHyperparameters','auto','Learners',t,'HyperparameterOptimizationOptions',...
        struct('AcquisitionFunctionName','expected-improvement-plus'),'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_EC_SleepScoringModel.mat'],'EC_MDL')
    % use the model to generate a set of scores for the even set of data
    [XevenLabels,~] = predict(EC_MDL,Xeven);
    % apply a logical patch on the REM events
    REMindex = strcmp(XevenLabels,'REM Sleep');
    reshapedREMindex = reshape(REMindex,dataLength,numFiles);
    patchedREMindex = [];
    % patch missing REM indeces due to theta band falling off
    for ee = 1:size(reshapedREMindex,2)
        remArray = reshapedREMindex(:,ee);
        patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
        patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
    end
    % change labels for each event
    for ff = 1:length(XevenLabels)
        if patchedREMindex(ff,1) == 1
            XevenLabels{ff,1} = 'REM Sleep';
        end
    end
    ConfusionData.EC.Ylabels{aa,1} = Yeven.behavState;
    ConfusionData.EC.Xlabels{aa,1} = XevenLabels;
    % confusion matrix
    EC_confMat = figure;
    cm = confusionchart(Yeven.behavState,XevenLabels);
    cm.ColumnSummary = 'column-normalized';
    cm.RowSummary = 'row-normalized';
    cm.Title = 'Ensemble Classifier Confusion Matrix';
    % pull data out of confusion matrix
    confVals = cm.NormalizedValues;
    totalScores = sum(confVals(:));
    EC_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
    disp(['Ensemble model prediction accuracy: ' num2str(EC_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(EC_confMat,[dirpath animalIDs{1,aa} '_IOS_EC_ConfusionMatrix']);
    close(EC_confMat)
    
    %% Decision Tree classification
    disp('Training Decision Tree Classifier...'); disp(' ')
    DT_MDL = fitctree(Xodd,Yodd,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_DT_SleepScoringModel.mat'],'DT_MDL')
    % use the model to generate a set of scores for the even set of data
    [XevenLabels,~] = predict(DT_MDL,Xeven);
    % apply a logical patch on the REM events
    REMindex = strcmp(XevenLabels,'REM Sleep');
    reshapedREMindex = reshape(REMindex,dataLength,numFiles);
    patchedREMindex = [];
    % patch missing REM indeces due to theta band falling off
    for gg = 1:size(reshapedREMindex,2)
        remArray = reshapedREMindex(:,gg);
        patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
        patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
    end
    % change labels for each event
    for hh = 1:length(XevenLabels)
        if patchedREMindex(hh,1) == 1
            XevenLabels{hh,1} = 'REM Sleep';
        end
    end
    ConfusionData.DT.Ylabels{aa,1} = Yeven.behavState;
    ConfusionData.DT.Xlabels{aa,1} = XevenLabels;
    % confusion matrix
    DT_confMat = figure;
    cm = confusionchart(Yeven.behavState,XevenLabels);
    cm.ColumnSummary = 'column-normalized';
    cm.RowSummary = 'row-normalized';
    cm.Title = 'Decision Tree Classifier Confusion Matrix';
    % pull data out of confusion matrix
    confVals = cm.NormalizedValues;
    totalScores = sum(confVals(:));
    DT_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
    disp(['Decision Tree model prediction accuracy: ' num2str(DT_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(DT_confMat,[dirpath animalIDs{1,aa} '_IOS_DT_ConfusionMatrix']);
    close(DT_confMat)
    
    %% Random forest
    disp('Training Random Forest Classifier...'); disp(' ')
    numTrees = 128;
    RF_MDL = TreeBagger(numTrees,Xodd,Yodd,'Method','Classification','Surrogate','all','ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_RF_SleepScoringModel.mat'],'RF_MDL')
    % use the model to generate a set of scores for the even set of data
    [XevenLabels,~] = predict(RF_MDL,Xeven);
    % apply a logical patch on the REM events
    REMindex = strcmp(XevenLabels,'REM Sleep');
    reshapedREMindex = reshape(REMindex,dataLength,numFiles);
    patchedREMindex = [];
    % patch missing REM indeces due to theta band falling off
    for ii = 1:size(reshapedREMindex,2)
        remArray = reshapedREMindex(:,ii);
        patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
        patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
    end
    % change labels for each event
    for jj = 1:length(XevenLabels)
        if patchedREMindex(jj,1) == 1
            XevenLabels{jj,1} = 'REM Sleep';
        end
    end
    ConfusionData.RF.Ylabels{aa,1} = Yeven.behavState;
    ConfusionData.RF.Xlabels{aa,1} = XevenLabels;
    % confusion matrix
    RF_confMat = figure;
    cm = confusionchart(Yeven.behavState,XevenLabels);
    cm.ColumnSummary = 'column-normalized';
    cm.RowSummary = 'row-normalized';
    cm.Title = 'Random Forest Classifier Confusion Matrix';
    % pull data out of confusion matrix
    confVals = cm.NormalizedValues;
    totalScores = sum(confVals(:));
    RF_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
    disp(['Random Forest model prediction accuracy: ' num2str(RF_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(RF_confMat,[dirpath animalIDs{1,aa} '_IOS_RF_ConfusionMatrix']);
    close(RF_confMat)
    
    %% k-nearest neighbor classifier
    disp('Training k-nearest neighbor Classifier...'); disp(' ')
    t = templateKNN('NumNeighbors',5,'Standardize',1);
    KNN_MDL = fitcecoc(Xodd,Yodd,'Learners',t);
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_KNN_SleepScoringModel.mat'],'KNN_MDL')
    % use the model to generate a set of scores for the even set of data
    [XevenLabels,~] = predict(KNN_MDL,Xeven);
    % apply a logical patch on the REM events
    REMindex = strcmp(XevenLabels,'REM Sleep');
    reshapedREMindex = reshape(REMindex,dataLength,numFiles);
    patchedREMindex = [];
    % patch missing REM indeces due to theta band falling off
    for kk = 1:size(reshapedREMindex,2)
        remArray = reshapedREMindex(:,kk);
        patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
        patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
    end
    % change labels for each event
    for ll = 1:length(XevenLabels)
        if patchedREMindex(ll,1) == 1
            XevenLabels{ll,1} = 'REM Sleep';
        end
    end
    ConfusionData.KNN.Ylabels{aa,1} = Yeven.behavState;
    ConfusionData.KNN.Xlabels{aa,1} = XevenLabels;
    % confusion matrix
    KNN_confMat = figure;
    cm = confusionchart(Yeven.behavState,XevenLabels);
    cm.ColumnSummary = 'column-normalized';
    cm.RowSummary = 'row-normalized';
    cm.Title = 'KNN Classifier Confusion Matrix';
    % pull data out of confusion matrix
    confVals = cm.NormalizedValues;
    totalScores = sum(confVals(:));
    KNN_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
    disp(['K-nearest neighbor model prediction accuracy: ' num2str(KNN_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(KNN_confMat,[dirpath animalIDs{1,aa} '_IOS_KNN_ConfusionMatrix']);
    close(KNN_confMat)
    
    %% Naive Bayes classifier
    disp('Training naive Bayes Classifier...'); disp(' ')
    NB_MDL = fitcnb(Xodd,Yodd,'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
    % save model in desired location
    save([dirpath animalIDs{1,aa} '_IOS_NB_SleepScoringModel.mat'],'NB_MDL')
    [XevenLabels,~] = predict(NB_MDL,Xeven);
    % apply a logical patch on the REM events
    REMindex = strcmp(XevenLabels,'REM Sleep');
    reshapedREMindex = reshape(REMindex,dataLength,numFiles);
    patchedREMindex = [];
    % patch missing REM indeces due to theta band falling off
    for mm = 1:size(reshapedREMindex,2)
        remArray = reshapedREMindex(:,mm);
        patchedREMarray = LinkBinaryEvents_IOS_Manuscript2020(remArray',[5,0]);
        patchedREMindex = vertcat(patchedREMindex,patchedREMarray'); %#ok<*AGROW>
    end
    % change labels for each event
    for nn = 1:length(XevenLabels)
        if patchedREMindex(nn,1) == 1
            XevenLabels{nn,1} = 'REM Sleep';
        end
    end
    ConfusionData.NB.Ylabels{aa,1} = Yeven.behavState;
    ConfusionData.NB.Xlabels{aa,1} = XevenLabels;
    % confusion matrix
    NB_confMat = figure;
    cm = confusionchart(Yeven.behavState,XevenLabels);
    cm.ColumnSummary = 'column-normalized';
    cm.RowSummary = 'row-normalized';
    cm.Title = 'Naive Bayes Classifier Confusion Matrix';
    % pull data out of confusion matrix
    confVals = cm.NormalizedValues;
    totalScores = sum(confVals(:));
    NB_accuracy = (sum(confVals([1,5,9])/totalScores))*100;
    disp(['Naive Bayes model prediction accuracy: ' num2str(NB_accuracy) '%']); disp(' ')
    % save model and figure
    savefig(NB_confMat,[dirpath animalIDs{1,aa} '_IOS_NB_ConfusionMatrix']);
    close(NB_confMat)
    cd(startingDirectory)
end
% save confusion matrix results
saveLoc = [startingDirectory '\Support Files\'];
save([saveLoc 'ConfusionData.mat'],'ConfusionData')

end
