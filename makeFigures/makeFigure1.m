%% makeFigure1

% 20220223 Yuasa

%% Initialize
close all;

SetDefaultAnalysisPath;
plotsavepth    = fullfile(figPth, 'Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavepth, 'Figure1');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

ecog_APRF_10b_outputSpectrumBar;

savefigauto(hF,fullfile(plotsavedir,'Figure1cd'));
