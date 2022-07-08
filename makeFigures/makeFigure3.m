%% makeFigure3

% 20220223 Yuasa

%% Initialize
close all;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavePthP, 'Figure3');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%-- Spectra
ecog_APRFF_10c_outputSpectrumPanel;

savefigauto(hF,fullfile(plotsavedir,'Figure3a'),'-vector');

%-- Surface
selsbj      = 2;
plotelecrad = 2;
ecog_APRFF_10a_plotSurface;
view(hF(1).Children(end),[-90 0]);
clear selsbj plotelecrad

savefigauto(hF,fullfile(plotsavedir,'Figure3b'),'Renderer','auto');

