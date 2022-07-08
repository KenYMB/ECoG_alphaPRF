%% makeFigureS5

% 20220620 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavePthP, 'FigureS5');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%%% Figure S5A %%%
ecog_APRFF_10k4_inPRFspectra

savefigauto(hF,fullfile(plotsavedir,'FigureS5a'),'-vector');


%%% Figure S5B %%%
ecog_APRFF_10k3_inPRFresponses               % show ECoG power in/out pRFs

savefigauto(hF,fullfile(plotsavedir,'FigureS5b'),'-vector');
