%% makeFigure2

% 20220223 Yuasa

%% Initialize
close all;

SetDefaultAnalysisPath;
plotsavepth    = fullfile(figPth, 'Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavepth, 'Figure2');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

ecog_APRF_10c_outputSpectrumPanel;

savefigauto(hF,fullfile(plotsavedir,'Figure2a'));
