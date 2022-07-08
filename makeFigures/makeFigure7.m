%% makeFigure7

% 20220223 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavePthP, 'Figure7');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%-- Individual location
ecog_APRFF_10g2_visualizePRFlocations;       % save figure with different arrange

savefigauto(hF(1),fullfile(plotsavedir,'Figure7a'),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,'Figure7c'),'-vector');

%-- Normalized location
ecog_APRFF_10g3_visualizePRFbootlocation;    % normalized bootstrapping location (normalized by ecc)

savefigauto(hF(1),fullfile(plotsavedir,'Figure7bd'),'-vector');
    exportgraphics(hF(1).Children.Children(2),[fullfile(plotsavedir,'Figure7b') '.eps'],'BackgroundColor','none');
    exportgraphics(hF(1).Children.Children(1),[fullfile(plotsavedir,'Figure7d') '.eps'],'BackgroundColor','none');
    
    