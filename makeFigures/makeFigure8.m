%% Spatial profiles with exogenous attention

% 20230215 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;
plotavg  = true;  plotovl  = true;  plotallrois  = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));

figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%%% Alpha pRF %%%
ecog_APRFF_10n2_AlphaPRFprofile;
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'a']),'-vector');

%%% Exogenous Attantion %%%
ecog_APRFF_10n1_ExogenousAttention;
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'b']),'-vector');
