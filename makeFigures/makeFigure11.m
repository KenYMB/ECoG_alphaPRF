%% makeFigure11

% 20220225 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavePthP, 'Figure11');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%-- Individual location
ecog_APRFF_10h_coherenceDistBootPair_both;       % save figure with different arrange

savefigauto(hF(1),fullfile(plotsavedir,'Figure11a'),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,'Figure11a-sub'),'-vector');
savefigauto(hF(3),fullfile(plotsavedir,'Figure11b'),'-vector');
savefigauto(hF(4),fullfile(plotsavedir,'Figure11b-sub'),'-vector');
