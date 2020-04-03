function [] = AvgLaserDopplerFlow_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
%   Purpose: Calculate the average laser doppler flow during different behavioral states
%________________________________________________________________________________________________________________________

animalIDs = {'T109','T110','T111','T119','T120','T121','T122','T123'};
behavFields = {'Whisk','Rest','NREM','REM'};
colorA = [(51/256),(160/256),(44/256)];   % rest color
colorB = [(192/256),(0/256),(256/256)];   % NREM color
colorC = [(255/256),(140/256),(0/256)];   % REM color
colorD = [(31/256),(120/256),(180/256)];  % whisk color

%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for b = 1:length(behavFields)
        behavField = behavFields{1,b};
        data.(behavField).LDFlow.flowMeans(a,1) = mean(AnalysisResults.(animalID).LDFlow.(behavField));
    end
end
% 
for c = 1:length(behavFields)
    behavField = behavFields{1,c};
    data.(behavField).LDFlow.behavMean = mean(data.(behavField).LDFlow.flowMeans);
    data.(behavField).LDFlow.behavStD = std(data.(behavField).LDFlow.flowMeans,0,1);
end

%% summary figure(s)
summaryFigure = figure;
xInds = ones(1,length(animalIDs));
%% CBV HbT
s1 = scatter(xInds*1,data.Whisk.LDFlow.flowMeans,'MarkerEdgeColor','k','MarkerFaceColor',colorD,'jitter','on','jitterAmount',0.25);
hold on
e1 = errorbar(1,data.Whisk.LDFlow.behavMean,data.Whisk.LDFlow.behavStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e1.Color = 'black';
s2 = scatter(xInds*2,data.Rest.LDFlow.flowMeans,'MarkerEdgeColor','k','MarkerFaceColor',colorA,'jitter','on','jitterAmount',0.25);
e2 = errorbar(2,data.Rest.LDFlow.behavMean,data.Rest.LDFlow.behavStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e2.Color = 'black';
s3 = scatter(xInds*3,data.NREM.LDFlow.flowMeans,'MarkerEdgeColor','k','MarkerFaceColor',colorB,'jitter','on','jitterAmount',0.25);
e3 = errorbar(3,data.NREM.LDFlow.behavMean,data.NREM.LDFlow.behavStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e3.Color = 'black';
s4 = scatter(xInds*4,data.REM.LDFlow.flowMeans,'MarkerEdgeColor','k','MarkerFaceColor',colorC,'jitter','on','jitterAmount',0.25);
e4 = errorbar(4,data.REM.LDFlow.behavMean,data.REM.LDFlow.behavStD,'d','MarkerEdgeColor','k','MarkerFaceColor','k');
e4.Color = 'black';
title('Mean Laser Doppler Flow')
ylabel('Flow Increase (%)')
legend([s1,s2,s3,s4],'Whisking','Awake Rest','NREM','REM','Location','NorthWest')
set(gca,'xtick',[])
set(gca,'xticklabel',[])
axis square
xlim([0,length(behavFields) + 1])
set(gca,'box','off')

% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath,'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Summary Figure - Laser Doppler Flow']);

end
