%% pRF results on HD grids

% 20220620 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));
figureIDs = strcat(figureIDs,["A","B"]);

%%% Visual areas %%%
ecog_APRFF_10l1_visualizeAreaGrid;

figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'b']),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'b-Border']),'-vector');
    exportgraphics(hF(1).Children(end),[fullfile(plotsavedir,[figureID 'b']) '.png']);
    exportgraphics(hF(2).Children(end),[fullfile(plotsavedir,[figureID 'b-Border']) '.png'],'BackgroundColor','none');

figureID = figureIDs{2};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(3),fullfile(plotsavedir,[figureID 'b']),'-vector');
savefigauto(hF(4),fullfile(plotsavedir,[figureID 'b-Border']),'-vector');
    exportgraphics(hF(3).Children(end),[fullfile(plotsavedir,[figureID 'b']) '.png']);
    exportgraphics(hF(4).Children(end),[fullfile(plotsavedir,[figureID 'b-Border']) '.png'],'BackgroundColor','none');


%%% pRF maps %%%
ecog_APRFF_10l2_visualizePRFsGrid;

figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'c_bb']),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'd_bb']),'-vector');
savefigauto(hF(3),fullfile(plotsavedir,[figureID 'e_bb']),'-vector');
savefigauto(hF(4),fullfile(plotsavedir,[figureID 'c_a']),'-vector');
savefigauto(hF(5),fullfile(plotsavedir,[figureID 'd_a']),'-vector');
savefigauto(hF(6),fullfile(plotsavedir,[figureID 'e_a']),'-vector');
    exportgraphics(hF(1).Children(end),[fullfile(plotsavedir,[figureID 'c_bb']) '.png']);
    exportgraphics(hF(2).Children(end),[fullfile(plotsavedir,[figureID 'd_bb']) '.png']);
    exportgraphics(hF(3).Children(end),[fullfile(plotsavedir,[figureID 'e_bb']) '.png']);
    exportgraphics(hF(4).Children(end),[fullfile(plotsavedir,[figureID 'c_a']) '.png']);
    exportgraphics(hF(5).Children(end),[fullfile(plotsavedir,[figureID 'd_a']) '.png']);
    exportgraphics(hF(6).Children(end),[fullfile(plotsavedir,[figureID 'e_a']) '.png']);

figureID = figureIDs{2};
plotsavedir    = fullfile(plotsavePthP, figureID);
savefigauto(hF(7),fullfile(plotsavedir,[figureID 'c_bb']),'-vector');
savefigauto(hF(8),fullfile(plotsavedir,[figureID 'd_bb']),'-vector');
savefigauto(hF(9),fullfile(plotsavedir,[figureID 'e_bb']),'-vector');
savefigauto(hF(10),fullfile(plotsavedir,[figureID 'c_a']),'-vector');
savefigauto(hF(11),fullfile(plotsavedir,[figureID 'd_a']),'-vector');
savefigauto(hF(12),fullfile(plotsavedir,[figureID 'e_a']),'-vector');
    exportgraphics(hF(7).Children(end),[fullfile(plotsavedir,[figureID 'c_bb']) '.png']);
    exportgraphics(hF(8).Children(end),[fullfile(plotsavedir,[figureID 'd_bb']) '.png']);
    exportgraphics(hF(9).Children(end),[fullfile(plotsavedir,[figureID 'e_bb']) '.png']);
    exportgraphics(hF(10).Children(end),[fullfile(plotsavedir,[figureID 'c_a']) '.png']);
    exportgraphics(hF(11).Children(end),[fullfile(plotsavedir,[figureID 'd_a']) '.png']);
    exportgraphics(hF(12).Children(end),[fullfile(plotsavedir,[figureID 'e_a']) '.png']);


%%% Surface %%%
SetDefault('subjectList_fname','subjectlist.tsv');
selsbj      = find(ismember(SetSubjectsList(subjectList_fname),...
               SetSubjectsList(subjectList_fname, 'hasHDgrid', 'yes')));
plotelecrad = 1;
ecog_APRFF_10a_plotSurface;
view(hF(1).Children(end),[50 0]);
view(hF(2).Children(end),[-50 0]);
clear selsbj plotelecrad

figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'a']),'Renderer','auto');

figureID = figureIDs{2};
plotsavedir    = fullfile(plotsavePthP, figureID);
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'a']),'Renderer','auto');

