%% makeFigure3

% 20220223 Yuasa

%% Initialize
close all;

SetDefaultAnalysisPath;
plotsavepth    = fullfile(figPth, 'Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavepth, 'Figure3');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

ecog_APRF_10c2_outputSpectrumModel;

savefigauto(hF(1),fullfile(plotsavedir,'Figure3a'));
savefigauto(hF(2),fullfile(plotsavedir,'Figure3b'));
savefigauto(hF(3),fullfile(plotsavedir,'Figure3c'));
savefigauto(hF(4),fullfile(plotsavedir,'Figure3d'));
