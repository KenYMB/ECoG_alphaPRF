%% pRF relations: broadband VS alpha

% 20220223 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));
figureIDs = strcat(figureIDs,["A","B"]);

ecog_APRFF_10g_visualizePRFrelations;

%%% 2D histogram %%%
figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(1),fullfile(plotsavedir,[figureID]),'-vector');

%%% Size scaling %%%
figureID = figureIDs{2};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(2),fullfile(plotsavedir,[figureID]),'-vector');
