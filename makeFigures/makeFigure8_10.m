%% makeFigure8,10

% 20220223 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure

ecog_APRFF_10g_visualizePRFrelations;

plotsavedir    = fullfile(plotsavePthP, 'Figure8');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(1),fullfile(plotsavedir,'Figure8'),'-vector');

plotsavedir    = fullfile(plotsavePthP, 'Figure10');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(2),fullfile(plotsavedir,'Figure10'),'-vector');
