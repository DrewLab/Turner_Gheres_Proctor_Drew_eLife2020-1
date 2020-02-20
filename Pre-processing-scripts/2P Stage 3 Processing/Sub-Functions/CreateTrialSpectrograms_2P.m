function CreateTrialSpectrograms_2P(mergedDataFiles,neuralDataTypes)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Analyzes the raw neural data from each MergedData.mat file and calculates two different spectrograms.
%________________________________________________________________________________________________________________________
%
%   Inputs: List of MergedData.mat files.
%
%   Outputs: Saves a SpecData.mat file of the same name as MergedData containing the analysis.
%
%   Last Revised: February 21st, 2019
%________________________________________________________________________________________________________________________

for a = 1:size(mergedDataFiles, 1)
    specDataFileID = [mergedDataFiles(a,1:end - 15) '_SpecData.mat'];
    if ~exist(specDataFileID, 'file') == true
        mergedDataFileID = mergedDataFiles(a, :);
        load(mergedDataFileID);
        duration = MergedData.notes.trialDuration_Sec;
        anFs = MergedData.notes.anFs;
        expectedLength = duration*anFs;
        [animalID,hem,~,fileID,imageID,vesselID] = GetFileInfo2_2P(mergedDataFileID);
        for z = 1:length(neuralDataTypes)
            neuralDT = neuralDataTypes{1,z};
            rawNeuro = detrend(MergedData.data.(neuralDT)(1:expectedLength),'constant');
            
            w0 = 60/(anFs/2);
            bw = w0/35;
            [num,den] = iirnotch(w0, bw);
            rawNeuro2 = filtfilt(num, den, rawNeuro);
            
            % Spectrogram parameters
            params1.tapers = [1 1];
            params1.Fs = anFs;
            params1.fpass = [1 100];
            movingwin1 = [1 1/10];
            
            params5.tapers = [5 9];
            params5.Fs = anFs;
            params5.fpass = [1 100];
            movingwin5 = [5 1/5];
            
            disp(['Creating spectrogram for file number ' num2str(a) ' of ' num2str(size(mergedDataFiles, 1)) '...']); disp(' ')
            
            [S1, T1, F1] = mtspecgramc_2P(rawNeuro2, movingwin1, params1);
            [S5, T5, F5] = mtspecgramc_2P(rawNeuro2, movingwin5, params5);
            
            SpecData.(neuralDT).fiveSec.S = S5';
            SpecData.(neuralDT).fiveSec.T = T5;
            SpecData.(neuralDT).fiveSec.F = F5;
            SpecData.(neuralDT).fiveSec.params = params5;
            SpecData.(neuralDT).fiveSec.movingwin = movingwin5;
            
            SpecData.(neuralDT).oneSec.S = S1';
            SpecData.(neuralDT).oneSec.T = T1;
            SpecData.(neuralDT).oneSec.F = F1;
            SpecData.(neuralDT).oneSec.params = params1;
            SpecData.(neuralDT).oneSec.movingwin = movingwin1;
            
            save([animalID '_' hem '_' fileID '_' imageID '_' vesselID '_SpecData.mat'], 'SpecData');
        end
    end
    disp('File exists. continuing...'); disp(' ')
end

end
