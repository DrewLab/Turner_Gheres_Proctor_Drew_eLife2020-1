function [] = SupplementalFigurePanelFive_Manuscript2020(rootFolder,AnalysisResults)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% Purpose:
%________________________________________________________________________________________________________________________

animalIDs = {'T99','T101','T102','T103','T105','T108','T109','T110','T111','T119','T120','T121','T122','T123'};
whiskDataTypes = {'ShortWhisks','IntermediateWhisks','LongWhisks'};
%% cd through each animal's directory and extract the appropriate analysis results
for a = 1:length(animalIDs)
    animalID = animalIDs{1,a};
    for c = 1:length(whiskDataTypes)
        whiskDataType = whiskDataTypes{1,c};
        % LH cortical
        data.(whiskDataType).adjLH.HbT(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).CBV_HbT.HbT;
        data.(whiskDataType).adjLH.CBV(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).CBV.CBV;
        data.(whiskDataType).adjLH.cortMUA(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).MUA.corticalData;
        data.(whiskDataType).adjLH.cortS(:,:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.corticalS;
        data.(whiskDataType).adjLH.cortT(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.T;
        data.(whiskDataType).adjLH.cortF(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.F;
        % RH cortical
        data.(whiskDataType).adjRH.HbT(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).CBV_HbT.HbT;
        data.(whiskDataType).adjRH.CBV(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).CBV.CBV;
        data.(whiskDataType).adjRH.cortMUA(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).MUA.corticalData;
        data.(whiskDataType).adjRH.cortS(:,:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.corticalS;
        data.(whiskDataType).adjRH.cortT(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.T;
        data.(whiskDataType).adjRH.cortF(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjRH.(whiskDataType).LFP.F;
        % hippocampal
        data.(whiskDataType).Hip.hipMUA(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).MUA.hippocampalData;
        data.(whiskDataType).Hip.hipS(:,:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.hippocampalS;
        data.(whiskDataType).Hip.hipT(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.T;
        data.(whiskDataType).Hip.hipF(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).LFP.F;
        % time vector
        data.(whiskDataType).timeVector(:,a) = AnalysisResults.(animalID).EvokedAvgs.Whisk.adjLH.(whiskDataType).timeVector;
    end
end
% concatenate the data from the contra and ipsi data
for e = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,e};
    data.(whiskDataType).HbT = cat(2,data.(whiskDataType).adjLH.HbT,data.(whiskDataType).adjRH.HbT);
    data.(whiskDataType).CBV = cat(2,data.(whiskDataType).adjLH.CBV,data.(whiskDataType).adjRH.CBV);
    data.(whiskDataType).cortMUA = cat(2,data.(whiskDataType).adjLH.cortMUA,data.(whiskDataType).adjRH.cortMUA);
    data.(whiskDataType).cortS = cat(3,data.(whiskDataType).adjLH.cortS,data.(whiskDataType).adjRH.cortS);
    data.(whiskDataType).cortT = cat(2,data.(whiskDataType).adjLH.cortT,data.(whiskDataType).adjRH.cortT);
    data.(whiskDataType).cortF = cat(2,data.(whiskDataType).adjLH.cortF,data.(whiskDataType).adjRH.cortF);
end
% concatenate the data from the contra and ipsi data
for e = 1:length(whiskDataTypes)
    whiskDataType = whiskDataTypes{1,e};
    data.(whiskDataType).meanHbT = mean(data.(whiskDataType).HbT,2);
    data.(whiskDataType).stdHbT = std(data.(whiskDataType).HbT,0,2);
    data.(whiskDataType).meanCBV = mean(data.(whiskDataType).CBV,2);
    data.(whiskDataType).stdCBV = std(data.(whiskDataType).CBV,0,2);
    data.(whiskDataType).meanCortMUA = mean(data.(whiskDataType).cortMUA,2);
    data.(whiskDataType).stdCortMUA = std(data.(whiskDataType).cortMUA,0,2);
    data.(whiskDataType).meanCortS = mean(data.(whiskDataType).cortS,3).*100;
    data.(whiskDataType).meanCortT = mean(data.(whiskDataType).cortT,2);
    data.(whiskDataType).meanCortF = mean(data.(whiskDataType).cortF,2);
    data.(whiskDataType).meanHipMUA = mean(data.(whiskDataType).Hip.hipMUA,2);
    data.(whiskDataType).stdHipMUA = std(data.(whiskDataType).Hip.hipMUA,0,2);
    data.(whiskDataType).meanHipS = mean(data.(whiskDataType).Hip.hipS,3).*100;
    data.(whiskDataType).meanHipT = mean(data.(whiskDataType).Hip.hipT,2);
    data.(whiskDataType).meanHipF = mean(data.(whiskDataType).Hip.hipF,2);
    data.(whiskDataType).meanTimeVector = mean(data.(whiskDataType).timeVector(:,a),2);
end
%% summary figure(s)
summaryFigure = figure;
sgtitle('Turner Manuscript 2020 - Supplemental Figure Panel Five')
%% [A] ShortWhisks whisks cortical MUA
ax1 = subplot(6,3,1);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA + data.ShortWhisks.stdCortMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCortMUA - data.ShortWhisks.stdCortMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[A] Short whisk cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')   
axis square
set(gca,'box','off')
ax1.TickLength = [0.03,0.03];
%% [B] IntermediateWhisks whisks cortical MUA
ax2 = subplot(6,3,2);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA + data.IntermediateWhisks.stdCortMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCortMUA - data.IntermediateWhisks.stdCortMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[B] Intermed whisk cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')   
axis square
set(gca,'box','off')
ax2.TickLength = [0.03,0.03];
%% [C] LongWhisks whisks cortical MUA
ax3 = subplot(6,3,3);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA + data.LongWhisks.stdCortMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCortMUA - data.LongWhisks.stdCortMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[C] Long whisk cortical MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')   
axis square
set(gca,'box','off')
ax3.TickLength = [0.03,0.03];
%% [D] ShortWhisks whisks cortical LFP
ax4 = subplot(6,3,4);
imagesc(data.ShortWhisks.meanCortT,data.ShortWhisks.meanCortF,data.ShortWhisks.meanCortS)
title('[D] Short whisk cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')   
c4 = colorbar;
ylabel(c4,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax4.TickLength = [0.03,0.03];
%% [E] IntermediateWhisks whisks cortical LFP
ax5 = subplot(6,3,5);
imagesc(data.IntermediateWhisks.meanCortT,data.IntermediateWhisks.meanCortF,data.IntermediateWhisks.meanCortS)
title('[E] Intermed whisk cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')   
c5 = colorbar;
ylabel(c5,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax5.TickLength = [0.03,0.03];
%% [F] LongWhisks whisks cortical LFP
ax6 = subplot(6,3,6);
imagesc(data.LongWhisks.meanCortT,data.LongWhisks.meanCortF,data.LongWhisks.meanCortS)
title('[F] Long whisk cortical LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')   
c6 = colorbar;
ylabel(c6,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax6.TickLength = [0.03,0.03];
%% [G] Short whisks hippocampal MUA
ax7 = subplot(6,3,7);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA + data.ShortWhisks.stdHipMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHipMUA - data.ShortWhisks.stdHipMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[G] Short whisk hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')   
axis square
set(gca,'box','off')
ax7.TickLength = [0.03,0.03];
%% [H] Intermediate whisks hippocampal MUA
ax8 = subplot(6,3,8);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA + data.IntermediateWhisks.stdHipMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHipMUA - data.IntermediateWhisks.stdHipMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[H] Intermed whisk hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')   
axis square
set(gca,'box','off')
ax8.TickLength = [0.03,0.03];
%% [I] Long whisks hippocampal MUA
ax9 = subplot(6,3,9);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA + data.LongWhisks.stdHipMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHipMUA - data.LongWhisks.stdHipMUA,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[I] Long whisk hippocampal MUA')
ylabel('\DeltaP/P (%)')
xlabel('Peri-whisk time (s)')   
axis square
set(gca,'box','off')
ax9.TickLength = [0.03,0.03];
%% [J] Short whisks hippocampal LFP
ax10 = subplot(6,3,10);
imagesc(data.ShortWhisks.meanHipT,data.ShortWhisks.meanHipF,data.ShortWhisks.meanHipS)
title('[J] Short whisk hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')   
c10 = colorbar;
ylabel(c10,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax10.TickLength = [0.03,0.03];
%% [K] Intermediate whisks hippocampal LFP
ax11 = subplot(6,3,11);
imagesc(data.IntermediateWhisks.meanHipT,data.IntermediateWhisks.meanHipF,data.IntermediateWhisks.meanHipS)
title('[K] Intermed whisk hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')   
c11 = colorbar;
ylabel(c11,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0 0])
axis square
axis xy
set(gca,'box','off')
ax11.TickLength = [0.03,0.03];
%% [L] Long whisks hippocampal LFP
ax12 = subplot(6,3,12);
imagesc(data.LongWhisks.meanHipT,data.LongWhisks.meanHipF,data.LongWhisks.meanHipS)
title('[L] Long whisk hippocampal LFP')
ylabel('Freq (Hz)')
xlabel('Peri-whisk time (s)')   
c12 = colorbar;
ylabel(c12,'\DeltaP/P (%)','rotation',-90,'VerticalAlignment','bottom')
caxis([-25,25])
set(gca,'Ticklength',[0,0])
axis square
axis xy
set(gca,'box','off')
ax12.TickLength = [0.03,0.03];
%% [M] Short whisks HbT
ax13 = subplot(6,3,13);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHbT,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHbT + data.ShortWhisks.stdHbT,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanHbT - data.ShortWhisks.stdHbT,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[M] Short whisk \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')  
axis square
set(gca,'box','off')
ax13.TickLength = [0.03,0.03];
%% [N] Intermediate whisks HbT
ax14 = subplot(6,3,14);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHbT,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHbT + data.IntermediateWhisks.stdHbT,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanHbT - data.IntermediateWhisks.stdHbT,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[N] Intermed whisk \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')  
axis square
set(gca,'box','off')
ax14.TickLength = [0.03,0.03];
%% [O] Long whisks HbT
ax15 = subplot(6,3,15);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHbT,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHbT + data.LongWhisks.stdHbT,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanHbT - data.LongWhisks.stdHbT,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[O] Long whisk \DeltaHbT (\muM)')
ylabel('\DeltaHbT (\muM)')
xlabel('Peri-whisk time (s)')  
axis square
set(gca,'box','off')
ax15.TickLength = [0.03,0.03];
%% [P] Short whisks refl
ax16 = subplot(6,3,16);
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCBV,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCBV + data.ShortWhisks.stdCBV,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.ShortWhisks.meanTimeVector,data.ShortWhisks.meanCBV - data.ShortWhisks.stdCBV,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[P] Short whisk reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peri-whisk time (s)')  
axis square
set(gca,'box','off')
ax16.TickLength = [0.03,0.03];
%% [Q] Intermediate whisks refl
ax17 = subplot(6,3,17);
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCBV,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCBV + data.IntermediateWhisks.stdCBV,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.IntermediateWhisks.meanTimeVector,data.IntermediateWhisks.meanCBV - data.IntermediateWhisks.stdCBV,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[Q] Intermed whisk reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peri-whisk time (s)')  
axis square
set(gca,'box','off')
ax17.TickLength = [0.03,0.03];
%% [R] Long whisks refl
ax18 = subplot(6,3,18);
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCBV,'color',colors_Manuscript2020('rich black'),'LineWidth',1);
hold on
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCBV + data.LongWhisks.stdCBV,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
plot(data.LongWhisks.meanTimeVector,data.LongWhisks.meanCBV - data.LongWhisks.stdCBV,'color',colors_Manuscript2020('battleship grey'),'LineWidth',0.5)
title('[R] Long whisk reflectance')
ylabel('\DeltaR/R (%)')
xlabel('Peri-whisk time (s)')  
axis square
set(gca,'box','off')
ax18.TickLength = [0.03,0.03];
%% axes positions
linkaxes([ax1,ax2,ax3,ax7,ax8,ax9],'xy')
linkaxes([ax4,ax5,ax6,ax10,ax11,ax12],'xy')
linkaxes([ax13,ax14,ax15],'xy')
linkaxes([ax16,ax17,ax18],'xy')
ax1Pos = get(ax1,'position');
ax2Pos = get(ax2,'position');
ax3Pos = get(ax3,'position');
ax4Pos = get(ax4,'position');
ax5Pos = get(ax5,'position');
ax6Pos = get(ax6,'position');
ax10Pos = get(ax10,'position');
ax11Pos = get(ax11,'position');
ax12Pos = get(ax12,'position');
ax4Pos(3:4) = ax1Pos(3:4);
ax5Pos(3:4) = ax2Pos(3:4);
ax6Pos(3:4) = ax3Pos(3:4);
ax10Pos(3:4) = ax1Pos(3:4);
ax11Pos(3:4) = ax2Pos(3:4);
ax12Pos(3:4) = ax3Pos(3:4);
set(ax4,'position',ax4Pos);
set(ax5,'position',ax5Pos);
set(ax6,'position',ax6Pos);
set(ax10,'position',ax10Pos);
set(ax11,'position',ax11Pos);
set(ax12,'position',ax12Pos);
%% save figure(s)
dirpath = [rootFolder '\Summary Figures and Structures\'];
if ~exist(dirpath, 'dir')
    mkdir(dirpath);
end
savefig(summaryFigure,[dirpath 'Supplemental Figure Panel 5']);
print('-painters','-dpdf','-fillpage',[dirpath 'Supplemental Figure Panel 5'])

end
