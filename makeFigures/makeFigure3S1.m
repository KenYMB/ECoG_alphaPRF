%% OG vs DoG model

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

%%% Variance explained %%%
ecog_APRFF_10j_checkmodel;

%-- modify axis range
axmin = -5;
hF.Children.Children(end-1).XLim(1) = axmin;
hF.Children.Children(end-1).YLim(1) = axmin;
hF.Children.Children(end-1).Children(1).YData(1) = -5;
hF.Children.Children(end-1).Children(2).XData(1) = -5;
hF.Children.Children(end-1).Children(3).YData(1) = -5;
hF.Children.Children(end-1).Children(4).XData(1) = -5;
hF.Children.Children(end-1).Children(5).YData(1) = -5;
hF.Children.Children(end-1).Children(5).XData(1) = -5;
savefigauto(hF,fullfile(plotsavedir,[figureID]),'-vector');
