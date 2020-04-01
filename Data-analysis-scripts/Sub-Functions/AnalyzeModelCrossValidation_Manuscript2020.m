function [AnalysisResults] = AnalyzeModelCrossValidation_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Determine the probability of different resting durations including a sleeping event
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};

% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    modelLocation = [rootFolder '\' animalID '\Figures\Sleep Models\'];
    cd(modelLocation)
    % load the random forest model for evaluating the cross-validation
    modelName = [animalID '_IOS_EC_SleepScoringModel.mat'];
    load(modelName,'-mat')
    iterations = 100;
    X = EC_MDL.X;
    Y = EC_MDL.Y;
    % run a 3-fold cross-validation on the model 100 times
    for aa = 1:iterations
        crossVal_MDL = crossval(EC_MDL,'kfold',3);
        AnalysisResults.(animalID).ModelCrossValidation.MDL_Loss(aa,1) = kfoldLoss(crossVal_MDL);
    end
    % re-create the model 100 times with shuffled data and run a single 3-fold cross-validation on the model
    for bb = 1:iterations
        shuffYIdx = randperm(numel(Y));
        shuffY = Y(shuffYIdx);
        t = templateTree('Reproducible',true);
        shuffEC_MDL = fitcensemble(X,shuffY,'OptimizeHyperparameters','auto','Learners',t,'HyperparameterOptimizationOptions',...
            struct('AcquisitionFunctionName','expected-improvement-plus'),'ClassNames',{'Not Sleep','NREM Sleep','REM Sleep'});
        shuffCV_MDL = crossval(shuffEC_MDL,'kfold',3);
        AnalysisResults.(animalID).ModelCrossValidation.shuffMDL_Loss(bb,1) = kfoldLoss(shuffCV_MDL);
        close all
    end
    % save figure if desired
    if strcmp(saveFigs,'y') == true
        distributionFig = figure;
        sgtitle([animalID ' Ensemble classifier cross-validation'])
        subplot(1,2,1)
        histogram(AnalysisResults.(animalID).ModelCrossValidation.MDL_Loss,'Normalization','probability','FaceColor','k')
        title('3-fold cross-validation of real data')
        ylabel('Probability')
        xlabel('kfoldLoss')
        axis square
        set(gca,'box','off')
        subplot(1,2,2)
        histogram(AnalysisResults.(animalID).ModelCrossValidation.shuffMDL_Loss,'Normalization','probability','FaceColor','k')
        title('3-fold cross-validation of shuffled data')
        ylabel('Probability')
        xlabel('kfoldLoss')
        axis square
        set(gca,'box','off')
        % Save the figure to directory.
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Figures/Cross Validation/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(distributionFig,[dirpath animalID '_CrossValidationDistribution']);
        close(distributionFig)
    end
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
