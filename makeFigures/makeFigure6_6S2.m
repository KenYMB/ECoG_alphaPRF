%% pRF relations: high-broadband VS low-broadband

% 20220620 Yuasa
% 20240503 Yuasa - update

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));

ecog_APRFF_10k2b_visualizePRFrelations_lowbb % pRF results with lowbb

%%% Size Comparison %%%
figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

savefigauto(hF(2),fullfile(plotsavedir,[figureID]),'-vector');

%%% 2D histogram %%%
figureID = figureIDs{2};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%-- modify axis range
rfsizelim = [0 5];
rfsizeax  = 1;
hF(1).Children.Children(rfsizeax).YLabel.Units = 'normalized';
nBins = hF(1).Children.Children(rfsizeax).Children(end).NumBins;
hF(1).Children.Children(rfsizeax).Children(end).XBinLimits = rfsizelim;
hF(1).Children.Children(rfsizeax).Children(end).YBinLimits = rfsizelim;
hF(1).Children.Children(rfsizeax).Children(end).NumBins = nBins;
hF(1).Children.Children(rfsizeax).XLim = rfsizelim;
hF(1).Children.Children(rfsizeax).YLim = rfsizelim;
    hF(1).Children.Children(rfsizeax).XTick = rfsizelim(1):rfsizelim(2);
    hF(1).Children.Children(rfsizeax).YTick = rfsizelim(1):rfsizelim(2);

savefigauto(hF(1),fullfile(plotsavedir,[figureID 'B']),'-vector');

%%% pRF locations %%%
ecog_APRFF_10k2c_visualizePRFlocations_lowbb % pRRF locations with lowbb
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'A']),'-vector');