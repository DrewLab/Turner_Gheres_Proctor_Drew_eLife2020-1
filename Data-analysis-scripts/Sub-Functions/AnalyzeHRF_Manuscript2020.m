function [AnalysisResults] = AnalyzeHRF_Manuscript2020(animalID,saveFigs,rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Analyze the spectral coherence between bilateral hemodynamic and neural signals.
%________________________________________________________________________________________________________________________

%% function parameters
animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
HRF_hemDataTypes = {'adjLH','adjRH'};
HRF_neuralBands = {'gammaBandPower','muaPower'};
behaviors = {'Contra','Whisk','Rest','NREM','REM'};

%% only run analysis for valid animal IDs
if any(strcmp(animalIDs,animalID))
    dataLocation = [rootFolder '/' animalID '/Bilateral Imaging/'];
    cd(dataLocation)
    for a = 1:length(HRF_hemDataTypes)
        hemDataType = HRF_hemDataTypes{1,a};
        for b = 1:length(HRF_neuralBands)
            neuralBand = HRF_neuralBands{1,b};
            [AnalysisResults] = CalculateHRFDeconvolution_Manuscript2020(animalID,neuralBand,hemDataType,'Contra',saveFigs,AnalysisResults);
            for c = 1:length(behaviors)
                behavior = behaviors{1,c};
                [AnalysisResults] = EvaluateCBVPredictionAccuracy_Manuscript2020(animalID,neuralBand,hemDataType,behavior,AnalysisResults);
            end
        end
    end
    % save data
    cd(rootFolder)
    save('AnalysisResults.mat','AnalysisResults')
end

end
