%% makeFigure7

% 20220223 Yuasa

%% Initialize
close all;

SetDefaultAnalysisPath;
plotsavepth    = fullfile(figPth, 'Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavepth, 'Figure7');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%-- Individual location
ecog_APRF_10g2_visualizePRFlocations       % save figure with different arrange

savefigauto(hF{1},fullfile(plotsavedir,'Figure7a'));
savefigauto(hF{2},fullfile(plotsavedir,'Figure7c'));

%-- Normalized location
ecog_APRF_10g3_visualizePRFbootlocation    % normalized bootstrapping location (normalized by ecc)

savefigauto(hF,fullfile(plotsavedir,'Figure7bd'));
