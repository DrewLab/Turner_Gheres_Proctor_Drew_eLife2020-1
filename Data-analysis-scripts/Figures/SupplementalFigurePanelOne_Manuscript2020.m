function [] = SupplementalFigurePanelOne_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: 
%________________________________________________________________________________________________________________________

dataDir = [rootFolder '\Summary Figures and Structures\Cross Correlation ROI\'];
cd(dataDir)
% character list of RawData files
rawDataFileStruct = dir('*_RawData.mat');
rawDataFiles = {rawDataFileStruct.name}';
rawDataFileIDs = char(rawDataFiles);
% character list of ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
% character list of WindowCam files
windowCamFileStruct = dir('*_WindowCam.bin');
windowCamFiles = {windowCamFileStruct.name}';
windowCamFileIDs = char(windowCamFiles);
% animal ID
[animalID,fileDate,~] = GetFileInfo_IOS_Manuscript2020(rawDataFileIDs(1,:));
strDay = ConvertDate_IOS_Manuscript2020(fileDate);
hem = 'Barrels';
if isfield(AnalysisResults.(animalID),'CrossCorrExample') == false
    % take the cross correlation between the neural data and each pixel
    for aa = 1:size(rawDataFileIDs,1)
        rawDataFileID = rawDataFileIDs(aa,:);
        load(rawDataFileID,'-mat')
        procDataFileID = procDataFileIDs(aa,:);
        load(procDataFileID,'-mat')
        windowCamFileID = windowCamFileIDs(aa,:);
        % read camera frames from binary file
        [frames] = ReadDalsaBinary_IOS_Manuscript2020(animalID,windowCamFileID);
        if aa == 1
            % open figure for ROI drawing
            windowFig = figure;
            imagesc(frames{1})
            title([animalID ' CBV camera'])
            xlabel('Image size (pixels)')
            ylabel('Image size (pixels)')
            colormap gray
            colorbar
            axis image
            caxis([0,2^RawData.notes.CBVCamBitDepth])
            % draw ROI for the mask over the entire windows
            isok = false;
            while isok == false
                disp('Drag a rectangular ROI over the window '); disp(' ')
                [~,rect] = imcrop;
                hold on;
                ROIoutline = rectangle('Position',rect,'EdgeColor','r');
                checkMask = input('Is the ROI okay? (y/n): ','s'); disp(' ')
                if strcmp(checkMask,'y') == true
                    isok = true;
                    ROIs.([hem '_' strDay]).rect = rect;
                end
                delete(ROIoutline);
            end
            close(windowFig)
        end
        % extract the pixel values from the window ROIs
        imageMask = nan(size(frames{1}));
        rectMask = ROIs.([hem '_' strDay]).rect;
        rectMask = round(rectMask);
        imageMask(rectMask(2):(rectMask(2) + rectMask(4)),rectMask(1):(rectMask(1) + rectMask(3))) = 1;
        for c = 1:length(frames) - 1
            frame = frames{1,c};
            frameHold = double(frame).*imageMask;
            imageStack.(hem)(:,c) = frameHold(~isnan(frameHold));
        end
        % extract and process gamma band power
        [z,p,k] = butter(4,1/(ProcData.notes.dsFs/2),'low');
        [sos,g] = zp2sos(z,p,k);
        LH_gammaBandPower = detrend(filtfilt(sos,g,ProcData.data.cortical_LH.gammaBandPower - ProcData.data.cortical_LH.gammaBandPower(1)) + ProcData.data.cortical_LH.gammaBandPower(1),'constant');
        RH_gammaBandPower = detrend(filtfilt(sos,g,ProcData.data.cortical_RH.gammaBandPower - ProcData.data.cortical_RH.gammaBandPower(1)) + ProcData.data.cortical_RH.gammaBandPower(1),'constant');
        LH_muaPower = detrend(filtfilt(sos,g,ProcData.data.cortical_LH.muaPower - ProcData.data.cortical_LH.muaPower(1)) + ProcData.data.cortical_LH.muaPower(1),'constant');
        RH_muaPower = detrend(filtfilt(sos,g,ProcData.data.cortical_RH.muaPower - ProcData.data.cortical_RH.muaPower(1)) + ProcData.data.cortical_RH.muaPower(1),'constant');
        % cross correlation
        lagTime = 5;   % seconds
        maxLag = lagTime*RawData.notes.CBVCamSamplingRate;
        singleHem = ProcData.notes.hemisphere;
        if strcmp(singleHem,'LH') == true
            gammaBandArray = LH_gammaBandPower;
            muaArray = LH_muaPower;
        elseif strcmp(singleHem,'RH') == true
            gammaBandArray = RH_gammaBandPower;
            muaArray = RH_muaPower;
        end
        % extract pixel values from each numel index in matrix image
        for e = 1:size(imageStack.(hem),1)
            disp(['Analyzing matrix numel ' num2str(e) ' of ' num2str(size(imageStack.(hem),1))]); disp(' ')
            pixelArray = imageStack.(hem)(e,:);
            pixelArray = detrend(filtfilt(sos,g,pixelArray - pixelArray(1)) + pixelArray(1),'constant');
            % gamma band power cross-correlation
            [gammaXCorrVals,lags] = xcorr(pixelArray,gammaBandArray,maxLag,'coeff');
            zeroPoint = find(lags == 0);
            validGammaVals = gammaXCorrVals(zeroPoint:zeroPoint + 45);
            gammaMaxCorr = min(validGammaVals);
            if isnan(gammaMaxCorr) == true
                gammaCorrMatrix.(hem)(1,e) = 0;
            else
                gammaCorrMatrix.(hem)(1,e) = gammaMaxCorr;
            end
            % multi unit activity cross-correlation
            [muaXCorrVals,~] = xcorr(pixelArray,muaArray,maxLag,'coeff');
            muaValidVals = muaXCorrVals(zeroPoint:zeroPoint + 45);
            muaMaxCorr = min(muaValidVals);
            if isnan(muaMaxCorr) == true
                muaCorrMatrix.(hem)(1,e) = 0;
            else
                muaCorrMatrix.(hem)(1,e) = muaMaxCorr;
            end
        end
        % proper size of the ROI based on camera/lens magnification
        circRadius = 30;   % pixels ~= 1 mm diameter
        % place circle along the most correlation region of each hemisphere
        rectMask = ROIs.([hem '_' strDay]).rect;
        rectMask = round(rectMask);
        imgWidth = rectMask(3) + 1;
        imgHeight = rectMask(4) + 1;
        gammaCorrImg(:,:,aa) = reshape(gammaCorrMatrix.(hem),imgHeight,imgWidth); %#ok<*AGROW>
        muaCorrImg(:,:,aa) = reshape(muaCorrMatrix.(hem),imgHeight,imgWidth); %#ok<*AGROW>
    end
    % generate image
    isok = false;
    while isok == false
        windowFig = figure;
        imagesc(AnalysisResults.(animalID).CrossCorrExample.gammaCorrImgAvg)
        title([animalID ' peak pixel correlations'])
        xlabel('Image size (pixels)')
        ylabel('Image size (pixels)')
        colormap parula
        colorbar
        axis image
        disp('Move the ROI over the most correlated region'); disp(' ')
        circ = drawcircle('Center',[0,0],'Radius',circRadius,'Color','r');
        checkCircle = input('Is the ROI okay? (y/n): ','s'); disp(' ')
        circPosition = round(circ.Center);
        if strcmp(checkCircle,'y') == true
            isok = true;
            rectBottomLeftCorner = [rectMask(1),rectMask(2) + rectMask(4)];
            rectTopLeftCorner = [rectMask(1),rectMask(2)];
            circPositionEdit = [rectBottomLeftCorner(1) + circPosition(1),rectTopLeftCorner(2) + circPosition(2)];
            ROIs.([hem '_' strDay]).circPosition = circPositionEdit;
            ROIs.([hem '_' strDay]).circRadius = circRadius;
        end
        delete(windowFig);
    end
    % save results
    AnalysisResults.(animalID).CrossCorrExample.gammaCorrImgAvg = mean(gammaCorrImg,3);
    AnalysisResults.(animalID).CrossCorrExample.muaCorrImgAvg = mean(muaCorrImg,3);
    AnalysisResults.(animalID).CrossCorrExample.frame = frames{1};
    AnalysisResults.(animalID).CrossCorrExample.ROIs = ROIs;
    save([rootFolder '\AnalysisResults.mat\'],'AnalysisResults')
end
% image mask to set background pixels to zero
testFig = figure;
imagesc(AnalysisResults.(animalID).CrossCorrExample.gammaCorrImgAvg.*-1)
title('Gamma-band power vs. pixel reflectance')
xlabel('Image size (pixels)')
ylabel('Image size (pixels)')
caxis([.25,0.5])
axis image
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
imgMask = roipoly;
close(testFig)
%% Supplemental figure panel one
summaryFigure = figure;
sgtitle('Supplemental Figure Panel 1 - Turner Manuscript 2020')
%% [A] gamma cross correlation image
ax1 = subplot(2,2,1);
gammaImg = AnalysisResults.(animalID).CrossCorrExample.gammaCorrImgAvg.*-1;
gammaImg(imgMask == 0) = 0;
imagesc(gammaImg)
title({'[A] Gamma-band power vs. pixel reflectance','Cross-correlation'})
xlabel('Image width (pixels)')
ylabel('Image height (pixels)')
colormap(ax1,parula)
c1 = colorbar;
ylabel(c1,{'Corr. coefficient';'Gamma-band [30-100 Hz] vs. reflectance'},'rotation',-90,'VerticalAlignment','bottom')
caxis([.25,0.5])
axis image
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
%% [B] MUA cross correlation image
ax2 = subplot(2,2,3);
muaImg = AnalysisResults.(animalID).CrossCorrExample.muaCorrImgAvg*-1;
muaImg(imgMask == 0) = 0;
imagesc(muaImg)
title({'[B] MUA power vs. pixel reflectance','Cross-correlation'})
xlabel('Image width (pixels)')
ylabel('Image height (pixels)')
colormap(ax2,parula)
c2 = colorbar;
ylabel(c2,{'Corr. coefficient';'MUA [300-3000 Hz] vs. reflectance'},'rotation',-90,'VerticalAlignment','bottom')
caxis([0.1,0.35])
axis image
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
%% [C] original image with circular ROI
ax3 = subplot(2,2,[2,4]);
imagesc(AnalysisResults.(animalID).CrossCorrExample.frame)
hold on;
drawcircle('Center',AnalysisResults.(animalID).CrossCorrExample.ROIs.(['Barrels_' strDay]).circPosition,'Radius',AnalysisResults.(animalID).CrossCorrExample.ROIs.(['Barrels_' strDay]).circRadius,'Color','r');
title({'[C] 1 mm OD ROI placement','on example window'})
xlabel('Image width (pixels)')
ylabel('Image height (pixels)')
set(gca,'YDir','reverse')
colormap(ax3,gray)
c3 = colorbar;
ylabel(c3,'Pixel Intensity (12-bit)','rotation',-90,'VerticalAlignment','bottom')
caxis([0,2^12])
axis image
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
cd(rootFolder)
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Supplemental Figure Panel 1']);
set(summaryFigure,'PaperPositionMode','auto');
print('-painters','-dpdf','-bestfit',[dirpath 'Supplemental Figure Panel 1'])

end
