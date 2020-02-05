rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
for a = 1:size(procDataFileIDs,1)
    load(procDataFileIDs(a,:))
    load(rawDataFileIDs(a,:))
    disp(num2str(a))
    ProcData.notes.dsFs = 30;   % downsampled Fs
    analogExpectedLength = ProcData.notes.trialDuration_sec*ProcData.notes.analogSamplingRate;
    neuralDataTypes = {'cortical_LH','cortical_RH','hippocampus'};
    for c = 1:length(neuralDataTypes)
        neuralDataType = neuralDataTypes{1,c};
        % MUA Band [300 - 3000]
        [muaPower,~] = ProcessNeuro_IOS_Manuscript2020(RawData,analogExpectedLength,'MUA',neuralDataType);
        ProcData.data.(neuralDataType).muaPower = muaPower;
        % Gamma Band [40 - 100]
        [gammaBandPower,~] = ProcessNeuro_IOS_Manuscript2020(RawData,analogExpectedLength,'Gam',neuralDataType);
        ProcData.data.(neuralDataType).gammaBandPower = gammaBandPower;
        % Beta [13 - 30 Hz]
        [betaBandPower,~] = ProcessNeuro_IOS_Manuscript2020(RawData,analogExpectedLength,'Beta',neuralDataType);
        ProcData.data.(neuralDataType).betaBandPower = betaBandPower;
        % Alpha [8 - 12 Hz]
        [alphaBandPower,~] = ProcessNeuro_IOS_Manuscript2020(RawData,analogExpectedLength,'Alpha',neuralDataType);
        ProcData.data.(neuralDataType).alphaBandPower = alphaBandPower;
        % Theta [4 - 8 Hz]
        [thetaBandPower,~] = ProcessNeuro_IOS_Manuscript2020(RawData,analogExpectedLength,'Theta',neuralDataType);
        ProcData.data.(neuralDataType).thetaBandPower = thetaBandPower;
        % Delta [1 - 4 Hz]
        [deltaBandPower,~] = ProcessNeuro_IOS_Manuscript2020(RawData,analogExpectedLength,'Delta',neuralDataType);
        ProcData.data.(neuralDataType).deltaBandPower = deltaBandPower;
    end
    save(procDataFileIDs(a,:),'ProcData')
end
