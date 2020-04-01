manualFileStruct = dir('*_ManualBaselineFileList.mat');
manualFile = {manualFileStruct.name}';
manualFileID = char(manualFile);
load(manualFileID)

procDataFileStruct = dir('*_ProcData.mat'); 
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);

for a = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(a,:);
    if strcmp(procDataFileID,ManualDecisions.fileIDs{a,1}) == true
        load(procDataFileID)      
        ProcData.manualBaselineInfo.fileDecision = ManualDecisions.validFiles{a,1};
        ProcData.manualBaselineInfo.startTime = ManualDecisions.startTimes{a,1};
        ProcData.manualBaselineInfo.endTime = ManualDecisions.endTimes{a,1};
        save(procDataFileID,'ProcData')
        clear ProcData
    else
        keyboard
    end
end