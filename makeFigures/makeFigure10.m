%% makeFigure10

% 20220225 Yuasa

%% Initialize
close all;

SetDefaultAnalysisPath;
plotsavepth    = fullfile(figPth, 'Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavepth, 'Figure10');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%-- Individual location
ecog_APRF_10h5_coherenceDistBootPair_both       % save figure with different arrange

savefigauto(hF(1),fullfile(plotsavedir,'Figure10a'));
savefigauto(hF(2),fullfile(plotsavedir,'Figure10b'));
savefigauto(hF(3),fullfile(plotsavedir,'Figure10c'));
savefigauto(hF(4),fullfile(plotsavedir,'Figure10d'));
