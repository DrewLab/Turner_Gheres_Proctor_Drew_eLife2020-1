function [AnalysisResults] = FigS1_Manuscript2020(rootFolder,saveFigs,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Generate figure panel S1 for Turner_Kederasetti_Gheres_Proctor_Costanzo_Drew_Manuscript2020
%________________________________________________________________________________________________________________________

%% set-up and process data
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
    AnalysisResults.(animalID).CrossCorrExample.imgMask = roipoly;
    close(testFig)
    save([rootFolder '\AnalysisResults.mat\'],'AnalysisResults')
end
cd(rootFolder)
%% set-up and process data
dataDir = [rootFolder '\Summary Figures and Structures\Pixel Drift Correction\'];
cd(dataDir)
% character list of ProcData files
procDataFileStruct = dir('*_ProcData.mat');
procDataFiles = {procDataFileStruct.name}';
procDataFileIDs = char(procDataFiles);
catLH_CBVdata = [];
catRH_CBVdata = [];
catCement_cementData = [];
% load the processed CBV/cement data from each file and concat it into one array
for aa = 1:size(procDataFileIDs,1)
    procDataFileID = procDataFileIDs(aa,:);
    load(procDataFileID,'-mat')
    samplingRate = ProcData.notes.CBVCamSamplingRate;
    LH_CBVdata = ProcData.data.CBV.LH;
    RH_CBVdata = ProcData.data.CBV.RH;
    Cement_cementData = ProcData.data.CBV.Cement;
    catLH_CBVdata = horzcat(catLH_CBVdata,LH_CBVdata); %#ok<*AGROW>
    catRH_CBVdata = horzcat(catRH_CBVdata,RH_CBVdata);
    catCement_cementData = horzcat(catCement_cementData,Cement_cementData);
end
% establish whether a slow exponential trend exists for the data
[B,A] = butter(3,0.01/(samplingRate/2),'low');
filtCatCement_cementData = filtfilt(B,A,catCement_cementData);
x = ((1:length(filtCatCement_cementData))/samplingRate)';
% create a weight vector for the trend
Cement_weightVec = ones(1,length(x));
Cement_secondHalfMean = mean(filtCatCement_cementData(floor(length(filtCatCement_cementData/2)):end));
for t = 1:length(Cement_weightVec)
    if filtCatCement_cementData(t) > Cement_secondHalfMean
        Cement_weightVec(t) = 10;
    end
end
% compare weighted models
Cement_modelFit = fit(x,filtCatCement_cementData','exp2','Weight',Cement_weightVec);
Cement_modelFit_Y = Cement_modelFit(x);
Cement_modelFit_norm = (Cement_modelFit_Y - min(Cement_modelFit_Y))./min(Cement_modelFit_Y);
Cement_modelFit_flip = 1 - Cement_modelFit_norm;
% apply exponential correction to original data
LH_adjCatC_CBVdata = catLH_CBVdata.*Cement_modelFit_flip';
RH_adjCatC_CBVdata = catRH_CBVdata.*Cement_modelFit_flip';
cd(rootFolder)
%% Fig. S1
summaryFigure = figure('Name','FigS1 (a-c,f-k)'); %#ok<*NASGU>
sgtitle('Figure Panel S1 (a-c,f-k) Turner Manuscript 2020')
%% [S1a] original image with circular ROI
ax1 = subplot(3,3,1);
imagesc(AnalysisResults.(animalID).CrossCorrExample.frame)
hold on;
drawcircle('Center',AnalysisResults.(animalID).CrossCorrExample.ROIs.(['Barrels_' strDay]).circPosition,'Radius',AnalysisResults.(animalID).CrossCorrExample.ROIs.(['Barrels_' strDay]).circRadius,'Color','r');
title({'[S1a] 1 mm OD ROI placement','on example window'})
xlabel('Image width (pixels)')
ylabel('Image height (pixels)')
set(gca,'YDir','reverse')
colormap(ax1,gray)
c1 = colorbar;
ylabel(c1,'Pixel intensity (12-bit)','rotation',-90,'VerticalAlignment','bottom')
caxis([0,2^12])
axis image
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
%% [S1b] gamma cross correlation image
ax2 = subplot(3,3,2);
gammaImg = AnalysisResults.(animalID).CrossCorrExample.gammaCorrImgAvg.*-1;
gammaImg(AnalysisResults.(animalID).CrossCorrExample.imgMask == 0) = 0;
imagesc(gammaImg)
title({'[S1b] Gamma-band vs. \DeltaR','Cross correlation'})
xlabel('Image width (pixels)')
ylabel('Image height (pixels)')
colormap(ax2,parula)
c2 = colorbar;
ylabel(c2,{'Corr. coefficient';'Gamma-band vs. \DeltaR'},'rotation',-90,'VerticalAlignment','bottom')
caxis([.25,0.5])
axis image
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
%% [S1c] MUA cross correlation image
ax3 = subplot(3,3,3);
muaImg = AnalysisResults.(animalID).CrossCorrExample.muaCorrImgAvg*-1;
muaImg(AnalysisResults.(animalID).CrossCorrExample.imgMask == 0) = 0;
imagesc(muaImg)
title({'[S1c] MUA vs. \DeltaR','Cross-correlation'})
xlabel('Image width (pixels)')
ylabel('Image height (pixels)')
colormap(ax3,parula)
c3 = colorbar;
ylabel(c3,{'Corr. coefficient';'MUA vs. \DeltaR'},'rotation',-90,'VerticalAlignment','bottom')
caxis([0.1,0.35])
axis image
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
%% [S1f] Data and exponential fit for cement ROI
ax4 = subplot(3,3,4);
p4a = plot(x,filtCatCement_cementData,'color',colors_Manuscript2020('deep carrot orange'),'LineWidth',1);
hold on
p4b = plot(x,Cement_modelFit_Y,'color',colors_Manuscript2020('electric purple'),'LineWidth',1);
title('[S1f] Cement ROI pixel drift')
xlabel('Time (s)')
ylabel('Pixel intensity (12-bit)')
legend([p4a,p4b],'cement ROI drift','exp2 fit')
axis tight
axis square
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [S1g] original left hemisphere
ax5 = subplot(3,3,5);
plot(x,catLH_CBVdata,'color','r','LineWidth',0.5)
title({'[S1g] Left hemisphere','original data'})
xlabel('Time (s)')
ylabel('Pixel intensity (12-bit)')
axis tight
axis square
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [S1h] original right hemisphere
ax6 = subplot(3,3,6);
plot(x,catRH_CBVdata,'color','c','LineWidth',0.5)
title({'[S1h] Right hemisphere','original data'})
xlabel('Time (s)')
ylabel('Pixel intensity (12-bit)')
axis tight
axis square
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [S1i] correction profile
ax7 = subplot(3,3,7);
plot(x,Cement_modelFit_flip,'color',colors_Manuscript2020('electric purple'),'LineWidth',1)
title('[S1i] Correction profile')
xlabel('Time (s)')
ylabel('Correction profile (%)')
axis tight
axis square
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [S1j] left hemisphere correction
ax8 = subplot(3,3,8);
plot(x,catLH_CBVdata,'color','r','LineWidth',0.5)
hold on
p8 = plot(x,LH_adjCatC_CBVdata,'color',colors_Manuscript2020('electric purple'),'LineWidth',0.5);
title({'[S1j] Left hemisphere','original vs. corrected data'})
xlabel('Time (s)')
ylabel('Pixel intensity (12-bit)')
legend(p8,'corrected')
axis tight
axis square
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% [S1k] right hemisphere correction
ax9 = subplot(3,3,9);
plot(x,catRH_CBVdata,'color','c','LineWidth',0.5)
hold on
plot(x,RH_adjCatC_CBVdata,'color',colors_Manuscript2020('electric purple'),'LineWidth',0.5);
title({'[S1k] Right hemisphere','original vs. corrected data'})
xlabel('Time (s)')
ylabel('Pixel intensity (12-bit)')
axis tight
axis square
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% save figure(s)
if strcmp(saveFigs,'y') == true
    dirpath = [rootFolder '\Summary Figures and Structures\'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(summaryFigure,[dirpath 'FigS1']);
    set(summaryFigure,'PaperPositionMode','auto');
    print('-painters','-dpdf','-bestfit',[dirpath 'FigS1'])
end

end
