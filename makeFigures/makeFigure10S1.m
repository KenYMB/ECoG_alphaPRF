%% Detailed spatial profiles of exogenous attention

% 20230215 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;
plotavg  = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));

figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%%% Exogenous Attantion %%%
ecog_APRFF_10n1_ExogenousAttention;
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'a']),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'b']),'-vector');
