%% pRF relations

% 20241121 Yuasa

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

%%%% Variance explaind VS pRF size
ecog_APRFF_10g1c_visualizePRFrelations_xval
savefigauto(hF(1),fullfile(plotsavedir,[figureID]),'-vector');
