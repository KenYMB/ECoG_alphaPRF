%% Spatial profiles with exogenous attention - individual

% 20230215 Yuasa

%% Initialize
close all;
clear modeldataID prfID useChans;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;
plotavg        = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));

figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%%% Exogenous Attantion %%%
iswideWin = false;
ecog_APRFF_10h2_coherenceDist_AllFreq;

figure(hF(1)); ylim([min(ylim) 0.24]);
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'a']),'-vector');
figure(hF(2)); ylim([min(ylim) 0.24]);
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'b']),'-vector');
