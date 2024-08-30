%% PSDs and freq-band averaged spectral powers

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

%%% PSD %%%
ecog_APRFF_10k4_inPRFspectra

savefigauto(hF,fullfile(plotsavedir,[figureID 'a']),'-vector');

%%% Averaged powers %%%
ecog_APRFF_10k3_inPRFresponses               % show ECoG power in/out pRFs

savefigauto(hF,fullfile(plotsavedir,[figureID 'b']),'-vector');
