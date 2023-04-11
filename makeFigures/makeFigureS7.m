%% pRF relations: high-broadband VS low-broadband

% 20220620 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));

figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%%% 2D histogram %%%
ecog_APRFF_10k2b_visualizePRFrelations_lowbb % pRF results with lowbb

%-- modify axis range
rfsizelim = [0 5];
rfsizeax  = 1;
hF.Children.Children(rfsizeax).YLabel.Units = 'normalized';
nBins = hF.Children.Children(rfsizeax).Children(end).NumBins;
hF.Children.Children(rfsizeax).Children(end).XBinLimits = rfsizelim;
hF.Children.Children(rfsizeax).Children(end).YBinLimits = rfsizelim;
hF.Children.Children(rfsizeax).Children(end).NumBins = nBins;
hF.Children.Children(rfsizeax).XLim = rfsizelim;
hF.Children.Children(rfsizeax).YLim = rfsizelim;
    hF.Children.Children(rfsizeax).XTick = rfsizelim(1):rfsizelim(2);
    hF.Children.Children(rfsizeax).YTick = rfsizelim(1):rfsizelim(2);

savefigauto(hF,fullfile(plotsavedir,[figureID]),'-vector');
