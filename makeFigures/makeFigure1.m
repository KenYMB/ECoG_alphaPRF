%% representative PSD & electrode location

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

%%% Spectra %%%
ecog_APRFF_10b_outputSpectrumBar;

savefigauto(hF,fullfile(plotsavedir,[figureID 'cd']),'-vector');

%%% Surface %%%
selsbj      = 2;
plotelecrad = 2;
ecog_APRFF_10a_plotSurface;
view(hF(1).Children(end),[-45 0]);
clear selsbj plotelecrad

savefigauto(hF,fullfile(plotsavedir,[figureID 'b']),'Renderer','auto');

