function [] = PixelDriftExample_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%________________________________________________________________________________________________________________________
%
% Purpose: Analyze the pixel drift and show how it is corrected
%________________________________________________________________________________________________________________________

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
%% comparison showing original LH data and the corrected data
suppFigure = figure;
sgtitle('Correction of CBV camera pixel drift')
% exponential fit
subplot(2,3,1)
p1 = plot(x,filtCatCement_cementData,'k');
hold on
p2 = plot(x,Cement_modelFit_Y,'color',colors_Manuscript2020('electric purple'),'LineWidth',1);
title('Cement ROI pixel drift')
xlabel('Time (sec)')
ylabel('12-bit pixel val')
legend([p1,p2],'cement ROI drift','exp2 fit')
axis tight
axis square
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
% original left hemisphere
subplot(2,3,2)
plot(x,catLH_CBVdata,'k')
title({'Left hemisphere', 'original data'})
xlabel('Time (sec)')
ylabel('12-bit pixel val')
axis tight
axis square
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
% original right hemisphere
subplot(2,3,3)
plot(x,catRH_CBVdata,'k')
title({'Right hemisphere', 'original data'})
xlabel('Time (sec)')
ylabel('12-bit pixel val')
axis tight
axis square
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
% correction profile
subplot(2,3,4)
plot(x,Cement_modelFit_flip,'color',colors_Manuscript2020('electric purple'),'LineWidth',3)
title('Correction profile')
xlabel('Time (sec)')
ylabel('Normalized val')
axis tight
axis square
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
% left hemisphere correction
subplot(2,3,5)
plot(x,catLH_CBVdata,'k')
hold on
p3 = plot(x,LH_adjCatC_CBVdata,'color',colors_Manuscript2020('electric purple'));
title({'Left hemisphere', 'original vs. corrected data'})
xlabel('Time (sec)')
ylabel('12-bit pixel val')
legend([p3],'corrected')
axis tight
axis square
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
% right hemisphere correction
subplot(2,3,6)
plot(x,catRH_CBVdata,'k')
hold on
p4 = plot(x,RH_adjCatC_CBVdata,'color',colors_Manuscript2020('electric purple'));
title({'Right hemisphere', 'original vs. corrected data'})
xlabel('Time (sec)')
ylabel('12-bit pixel val')
legend([p4],'corrected')
axis tight
axis square
set(gca,'Ticklength',[0,0])
set(gca,'box','off')
% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(suppFigure,[dirpath 'Supplemental Figure - Pixel Drift Correction']);

end
