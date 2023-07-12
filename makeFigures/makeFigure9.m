%% alpha computation model

% 20220223 Yuasa

%% Initialize
close all;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));

figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%%% Model components %%%
ecog_APRFF_10c2_outputSpectrumModel;

savefigauto(hF(1),fullfile(plotsavedir,[figureID 'a']),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'b']),'-vector');
savefigauto(hF(3),fullfile(plotsavedir,[figureID 'c']),'-vector');
savefigauto(hF(4),fullfile(plotsavedir,[figureID 'd']),'-vector');
