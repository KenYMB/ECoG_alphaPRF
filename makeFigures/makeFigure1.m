%% representative PSD & electrode location

% 20220223 Yuasa
% 20231101 Yuasa - update
% 20241119 Yuasa - update

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

%%% Evoked %%%
% ecog_APRFF_10b2_outputVoltageBar;         % Averaged trials
ecog_APRFF_10b3_outputVoltageBarSingle;     % Single trial

savefigauto(hF,fullfile(plotsavedir,[figureID 'cd']),'-vector');

%%% Spectra %%%
ecog_APRFF_10b_outputSpectrumBar;

savefigauto(hF,fullfile(plotsavedir,[figureID 'ef']),'-vector');

%%% Surface %%%
selsbj      = 2;
plotelecrad = 2;
ecog_APRFF_10a_plotSurface;
view(hF(1).Children(end),[-45 0]);
clear selsbj plotelecrad

savefigauto(hF,fullfile(plotsavedir,[figureID 'b']),'Renderer','auto');

