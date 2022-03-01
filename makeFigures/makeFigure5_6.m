%% makeFigure5,6

% 20220223 Yuasa

%% Initialize
close all;

SetDefaultAnalysisPath;
plotsavepth    = fullfile(figPth, 'Publication');
issaveplot     = false;

%% Figure

ecog_APRF_10f_channelSelection;

%%% Figure 5
plotsavedir    = fullfile(plotsavepth, 'Figure5');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(1),fullfile(plotsavedir,'Figure5'));

%%% Figure 6
plotsavedir    = fullfile(plotsavepth, 'Figure6');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(2),fullfile(plotsavedir,'Figure6a'));
savefigauto(hF(3),fullfile(plotsavedir,'Figure6b'));
savefigauto(hF(4),fullfile(plotsavedir,'Figure6c'));
savefigauto(hF(5),fullfile(plotsavedir,'Figure6d'));
