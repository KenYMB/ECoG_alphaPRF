%% representative PSD corrections & electrode locations

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
% ecog_APRFF_10c_outputSpectrumPanel;
ecog_APRFF_10c_outputSpectrumPanelError;

savefigauto(hF,fullfile(plotsavedir,[figureID 'a']),'-vector');

%%% Surface %%%
selsbj      = 2;
plotelecrad = 2;
ecog_APRFF_10a_plotSurface;
view(hF(1).Children(end),[-90 0]);
clear selsbj plotelecrad

savefigauto(hF,fullfile(plotsavedir,[figureID 'b']),'Renderer','auto');

