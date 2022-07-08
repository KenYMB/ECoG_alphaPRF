%% makeFigure1

% 20220223 Yuasa

%% Initialize
close all;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavePthP, 'Figure1');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

ecog_APRFF_10c2_outputSpectrumModel;

savefigauto(hF(1),fullfile(plotsavedir,'Figure1a'),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,'Figure1b'),'-vector');
savefigauto(hF(3),fullfile(plotsavedir,'Figure1c'),'-vector');
savefigauto(hF(4),fullfile(plotsavedir,'Figure1d'),'-vector');
