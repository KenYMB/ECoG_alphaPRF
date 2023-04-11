%% coherence from HD grids

% 20220225 Yuasa

%% Initialize
close all;
clear modeldataID prfID useChans;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));

figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%%% Individual location %%%
iswideWin = false;
ecog_APRFF_10h_coherenceDistBootPair_both;       % save figure with different arrange

savefigauto(hF(1),fullfile(plotsavedir,[figureID 'a']),'-vector');
% savefigauto(hF(2),fullfile(plotsavedir,[figureID 'a-sub']),'-vector');
savefigauto(hF(3),fullfile(plotsavedir,[figureID 'b']),'-vector');
% savefigauto(hF(4),fullfile(plotsavedir,[figureID 'b-sub']),'-vector');
