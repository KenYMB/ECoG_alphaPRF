%% makeFigure9

% 20220620 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavePthP, 'Figure9');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%%% Figure 9A %%%
ecog_APRFF_10k1a_checkalpha                  % FaC vs Fa

savefigauto(hF,fullfile(plotsavedir,'Figure9a'),'-vector');


%%% Figure 9B %%%
ecog_APRFF_10k1b_visualizePRFrelations_noC   % pRF results with Fa

savefigauto(hF,fullfile(plotsavedir,'Figure9b'),'-vector');

