function [singleTrialFig] = GenerateSingleFigures_2P(mergedDataFiles,RestingBaselines,saveFigs)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
%   Purpose: Create a summary figure for a single n minute two photon trial
%________________________________________________________________________________________________________________________

for a = 1:size(mergedDataFiles,1)
    mergedDataFile = mergedDataFiles(a,:);
    load(mergedDataFile)
    disp(['Creating single trial summary figure ' num2str(a) ' of ' num2str(size(mergedDataFiles,1)) '...']); disp(' ')
    [animalID,hem,fileDate,fileID,imageID,vesselID] = GetFileInfo2_2P(mergedDataFile);
    strDay = ConvertDate_2P(fileDate);
    
    %% Filter the whisker angle and identify the solenoid timing and location.
    % Setup butterworth filter coefficients for a 10 Hz lowpass based on the sampling rate
    [B, A] = butter(3,10/(MergedData.notes.dsFs/2),'low');
    % Whisker angle
    filteredWhiskerAngle = filtfilt(B,A,MergedData.data.whiskerAngle);
    binWhiskers = MergedData.data.binWhiskerAngle;
    % Force sensor
    filtForceSensor = filtfilt(B,A,MergedData.data.forceSensorL);
    binForce = MergedData.data.binForceSensorM;
    % EMG
    EMG = MergedData.data.EMG;
    % Solenoids
    LPadSol = MergedData.data.solenoids.LPadSol;
    RPadSol = MergedData.data.solenoids.RPadSol;
    AudSol = MergedData.data.solenoids.AudSol;
    
    %% Vessel diameter - normalize and then lowpass filer
    % Setup butterworth filter coefficients for a 1 Hz lowpass based on the sampling rate
    [D,C] = butter(3,1/(MergedData.notes.p2Fs/2),'low');
    vesselDiameter = MergedData.data.vesselDiameter;
    normVesselDiameter = (vesselDiameter - RestingBaselines.(vesselID).(strDay).vesselDiameter.baseLine)./(RestingBaselines.(vesselID).(strDay).vesselDiameter.baseLine);
    filtVesselDiameter = filtfilt(D,C,normVesselDiameter)*100;
    
    %% Normalized neural spectrogram
    specDataFile = [animalID '_' hem '_' fileID '_' imageID '_' vesselID '_SpecData.mat'];
    load(specDataFile,'-mat');
    cortNormS = SpecData.corticalNeural.fiveSec.normS.*100;
    hipNormS = SpecData.hippocampalNeural.fiveSec.normS.*100;
    T = SpecData.corticalNeural.fiveSec.T;
    F = SpecData.corticalNeural.fiveSec.F;
    
    %% Yvals for behavior Indices
    whisking_YVals = 1.10*max(filtVesselDiameter)*ones(size(binWhiskers));
    force_YVals = 1.20*max(filtVesselDiameter)*ones(size(binForce));
    LPad_Yvals = 1.30*max(filtVesselDiameter)*ones(size(LPadSol));
    RPad_Yvals = 1.30*max(filtVesselDiameter)*ones(size(RPadSol));
    Aud_Yvals = 1.30*max(filtVesselDiameter)*ones(size(AudSol));
    whiskInds = binWhiskers.*whisking_YVals;
    forceInds = binForce.*force_YVals;
    for x = 1:length(whiskInds)
        if whiskInds(1,x) == 0
            whiskInds(1,x) = NaN;
        end
        
        if forceInds(1,x) == 0
            forceInds(1,x) = NaN;
        end
    end
    
    %% Figure
    singleTrialFig = figure;
    fileID2 = strrep(fileID,'_',' ');
    sgtitle([animalID ' Two-photon behavioral characterization and vessel ' vesselID ' diameter changes for ' fileID2 ' image ' imageID])
    % XCorr
    ax1 = subplot(7,1,1);
    plot(MergedData.notes.MScan.shiftLags,MergedData.notes.MScan.shiftXCorr,'k')
    title('Time delay cross-correlation shift')
    ylabel('XCorr (A.U.)')
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % Force sensor
    ax2 = subplot(7,1,2);
    plot((1:length(filtForceSensor))/MergedData.notes.dsFs,filtForceSensor,'color',colors_IOS('sapphire'))
    title('Force sensor')
    ylabel('Force (V)')
    xlim([0 MergedData.notes.trialDuration_Sec])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % EMG
    ax3 = subplot(7,1,3);
    plot((1:length(EMG))/MergedData.notes.dsFs,EMG,'color', colors_IOS('deep carrot orange'))
    title('EMG')
    ylabel('EMG (Volts^2)')
    xlim([0 MergedData.notes.trialDuration_Sec])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % Whisker angle
    ax4 = subplot(7,1,4);
    plot((1:length(filteredWhiskerAngle))/MergedData.notes.dsFs,-filteredWhiskerAngle,'color','k')
    title('Whisker angle')
    ylabel('Whisker angle (deg)')
    xlim([0 MergedData.notes.trialDuration_Sec])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % Vessel diameter
    ax5 = subplot(7,1,5);
    plot((1:length(filtVesselDiameter))/MergedData.notes.p2Fs, detrend(filtVesselDiameter, 'constant'), 'color', colors_IOS('dark candy apple red'))
    hold on;
    s1 = scatter((1:length(binForce))/MergedData.notes.dsFs,forceInds,'.','MarkerEdgeColor',colors_IOS('sapphire'));
    s2 = scatter((1:length(binWhiskers))/MergedData.notes.dsFs,whiskInds,'.','MarkerEdgeColor',colors_IOS('electric purple'));
    s3 = scatter(LPadSol,LPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','c');
    s4 = scatter(RPadSol,RPad_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','m');
    s5 = scatter(AudSol,Aud_Yvals,'v','MarkerEdgeColor','k','MarkerFaceColor','g');
    title('Vessel diameter')
    ylabel('\DeltaD/D (%)')
    legend([s1,s2,s3,s4,s5],'Movement','Whisking','LPadSol','RPadSol','AudSol')
    xlim([0 MergedData.notes.trialDuration_Sec])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % Cortical LFP
    ax6 = subplot(7,1,6);
    semilog_imagesc_IOS(T,F,cortNormS,'y')
    axis xy
    c5 = colorbar;
    ylabel(c5,'\DeltaP/P (%)')
    caxis([-100 100])
    title('Cortical LFP')
    ylabel('Frequency (Hz)')
    xlim([0 MergedData.notes.trialDuration_Sec])
    set(gca,'Xticklabel',[])
    set(gca,'box','off')
    axis tight
    % Hippocampal LFP
    ax7 = subplot(7,1,7);
    semilog_imagesc_IOS(T,F,hipNormS,'y')
    axis xy
    c6 = colorbar;
    ylabel(c6,'\DeltaP/P (%)')
    caxis([-100 100])
    title('Hippocampal LFP')
    xlabel('Time (sec)')
    ylabel('Frequency (Hz)')
    xlim([0 MergedData.notes.trialDuration_Sec])
    set(gca,'box','off')
    axis tight
    % Axes properties
    linkaxes([ax2,ax3,ax4,ax5,ax6,ax7],'x')
    ax2Pos = get(ax2,'position');
    ax6Pos = get(ax6,'position');
    ax7Pos = get(ax7,'position');
    ax6Pos(3:4) = ax2Pos(3:4);
    ax7Pos(3:4) = ax2Pos(3:4);
    set(ax6,'position',ax6Pos);
    set(ax7,'position',ax7Pos);

    %% Save the file to directory.
    if strcmp(saveFigs,'y') == true
        [pathstr,~,~] = fileparts(cd);
        dirpath = [pathstr '/Combined Imaging/Figures/Single Trial Summary/'];
        if ~exist(dirpath,'dir')
            mkdir(dirpath);
        end
        savefig(singleTrialFig,[dirpath animalID '_' hem '_' fileID '_' imageID '_' vesselID '_SingleTrialFig']);
        close(singleTrialFig)
    end
end

end