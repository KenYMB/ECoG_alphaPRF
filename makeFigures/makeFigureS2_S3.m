%% makeFigure10

% 20220620 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
%-- Visual areas
ecog_APRFF_10l1_visualizeAreaGrid;

plotsavedir    = fullfile(plotsavePthP, 'FigureS2');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(1),fullfile(plotsavedir,'FigureS2b'),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,'FigureS2b-Border'),'-vector');
    exportgraphics(hF(1).Children(end),[fullfile(plotsavedir,'FigureS2b') '.png']);
    exportgraphics(hF(2).Children(end),[fullfile(plotsavedir,'FigureS2b-Border') '.png'],'BackgroundColor','none');

plotsavedir    = fullfile(plotsavePthP, 'FigureS3');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(3),fullfile(plotsavedir,'FigureS3b'),'-vector');
savefigauto(hF(4),fullfile(plotsavedir,'FigureS3b-Border'),'-vector');
    exportgraphics(hF(3).Children(end),[fullfile(plotsavedir,'FigureS3b') '.png']);
    exportgraphics(hF(4).Children(end),[fullfile(plotsavedir,'FigureS3b-Border') '.png'],'BackgroundColor','none');


%-- pRF maps
ecog_APRFF_10l2_visualizePRFsGrid;

plotsavedir    = fullfile(plotsavePthP, 'FigureS2');
savefigauto(hF(1),fullfile(plotsavedir,'FigureS2c_bb'),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,'FigureS2d_bb'),'-vector');
savefigauto(hF(3),fullfile(plotsavedir,'FigureS2e_bb'),'-vector');
savefigauto(hF(4),fullfile(plotsavedir,'FigureS2c_a'),'-vector');
savefigauto(hF(5),fullfile(plotsavedir,'FigureS2d_a'),'-vector');
savefigauto(hF(6),fullfile(plotsavedir,'FigureS2e_a'),'-vector');
    exportgraphics(hF(1).Children(end),[fullfile(plotsavedir,'FigureS2c_bb') '.png']);
    exportgraphics(hF(2).Children(end),[fullfile(plotsavedir,'FigureS2d_bb') '.png']);
    exportgraphics(hF(3).Children(end),[fullfile(plotsavedir,'FigureS2e_bb') '.png']);
    exportgraphics(hF(4).Children(end),[fullfile(plotsavedir,'FigureS2c_a') '.png']);
    exportgraphics(hF(5).Children(end),[fullfile(plotsavedir,'FigureS2d_a') '.png']);
    exportgraphics(hF(6).Children(end),[fullfile(plotsavedir,'FigureS2e_a') '.png']);

plotsavedir    = fullfile(plotsavePthP, 'FigureS3');
savefigauto(hF(7),fullfile(plotsavedir,'FigureS3c_bb'),'-vector');
savefigauto(hF(8),fullfile(plotsavedir,'FigureS3d_bb'),'-vector');
savefigauto(hF(9),fullfile(plotsavedir,'FigureS3e_bb'),'-vector');
savefigauto(hF(10),fullfile(plotsavedir,'FigureS3c_a'),'-vector');
savefigauto(hF(11),fullfile(plotsavedir,'FigureS3d_a'),'-vector');
savefigauto(hF(12),fullfile(plotsavedir,'FigureS3e_a'),'-vector');
    exportgraphics(hF(7).Children(end),[fullfile(plotsavedir,'FigureS3c_bb') '.png']);
    exportgraphics(hF(8).Children(end),[fullfile(plotsavedir,'FigureS3d_bb') '.png']);
    exportgraphics(hF(9).Children(end),[fullfile(plotsavedir,'FigureS3e_bb') '.png']);
    exportgraphics(hF(10).Children(end),[fullfile(plotsavedir,'FigureS3c_a') '.png']);
    exportgraphics(hF(11).Children(end),[fullfile(plotsavedir,'FigureS3d_a') '.png']);
    exportgraphics(hF(12).Children(end),[fullfile(plotsavedir,'FigureS3e_a') '.png']);


%-- Surface
SetDefault('subjectList_fname','subjectlist.tsv');
selsbj      = find(ismember(SetSubjectsList(subjectList_fname),...
               SetSubjectsList(subjectList_fname, 'hasHDgrid', 'yes')));
plotelecrad = 1;
ecog_APRFF_10a_plotSurface;
view(hF(1).Children(end),[50 0]);
view(hF(2).Children(end),[-50 0]);
clear selsbj plotelecrad

plotsavedir    = fullfile(plotsavePthP, 'FigureS2');
savefigauto(hF(1),fullfile(plotsavedir,'FigureS2a'),'Renderer','auto');

plotsavedir    = fullfile(plotsavePthP, 'FigureS3');
savefigauto(hF(2),fullfile(plotsavedir,'FigureS3a'),'Renderer','auto');

