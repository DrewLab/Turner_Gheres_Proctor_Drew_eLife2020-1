function [decData] = DecimateRestData_Manuscript2020(data,fileIDs,durations,eventTimes,ManualDecisions)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Remove resting events from the various fields that aren't in the manual selection
%________________________________________________________________________________________________________________________

trialDuration_sec = 900;   % sec
offset = 0.5;   % sec
x = 1;
for a = 1:length(data)
    fileID = fileIDs{a,1};
    startTime = eventTimes(a,1);
    endTime = startTime + durations(a,1);
    manualStartTime = [];
    manualEndTime = [];
    for b = 1:length(ManualDecisions.fileIDs)
        [~,~,manualFileID] = GetFileInfo_IOS_Manuscript2020(ManualDecisions.fileIDs{b,1});
        if strcmp(fileID,manualFileID) == true
            manualStartTime = ManualDecisions.startTimes{b,1};
            manualEndTime = ManualDecisions.endTimes{b,1};
            break
        end
    end
    % check that the event falls within appropriate bounds
    if startTime >= manualStartTime && endTime <= manualEndTime
        if startTime >= offset && endTime <= (trialDuration_sec - offset)
            decData{x,1} = data{a,1}; %#ok<*AGROW>
            x = x + 1;
        end
    end
end

end

