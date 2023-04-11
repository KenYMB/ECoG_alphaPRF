%% pRF locations

% 20220223 Yuasa

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

%%% Individual location %%%
ecog_APRFF_10g2_visualizePRFlocations;       % save figure with different arrange

savefigauto(hF(1),fullfile(plotsavedir,[figureID 'a']),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'c']),'-vector');

%%% Normalized location %%%
ecog_APRFF_10g3_visualizePRFbootlocation;    % normalized bootstrapping location (normalized by ecc)

savefigauto(hF(1),fullfile(plotsavedir,[figureID 'bd']),'-vector');
    exportgraphics(hF(1).Children.Children(2),[fullfile(plotsavedir,[figureID 'b']) '.eps'],'BackgroundColor','none');
    exportgraphics(hF(1).Children.Children(1),[fullfile(plotsavedir,[figureID 'd']) '.eps'],'BackgroundColor','none');
    
    