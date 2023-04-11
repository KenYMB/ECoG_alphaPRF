%% pRF relations: model-based alpha VS model-free alpha

% 20220620 Yuasa

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

%%% Violin plot %%%
ecog_APRFF_10k1e_checkalpha_bothgaind          % FaC vs Fa

savefigauto(hF(1),fullfile(plotsavedir,[figureID]),'-vector');
