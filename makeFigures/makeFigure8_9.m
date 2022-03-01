%% makeFigure8,9

% 20220223 Yuasa

%% Initialize
close all;

SetDefaultAnalysisPath;
plotsavepth    = fullfile(figPth, 'Publication');
issaveplot     = false;

%% Figure

ecog_APRF_10g_visualizePRFrelations

plotsavedir    = fullfile(plotsavepth, 'Figure8');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(1),fullfile(plotsavedir,'Figure8'));

plotsavedir    = fullfile(plotsavepth, 'Figure9');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(2),fullfile(plotsavedir,'Figure9'));
