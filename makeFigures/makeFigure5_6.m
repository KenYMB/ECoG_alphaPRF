%% makeFigure5,6

% 20220223 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure

ecog_APRFF_10f_channelSelection;

%%% Figure 5
plotsavedir    = fullfile(plotsavePthP, 'Figure5');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(1),fullfile(plotsavedir,'Figure5'),'-vector');

%%% Figure 6
plotsavedir    = fullfile(plotsavePthP, 'Figure6');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
axis([hF([2,4]).Children],'ij');
  btmpad = 0.03;
  for iax = [2,4]
    savepos = hF(iax).Children.Position;
    axshift = savepos(2) - btmpad;
    hF(iax).Children.XAxisLocation = 'origin';
    hF(iax).Children.Position = savepos - [0 axshift 0 0];
    hF(iax).Children.Title.Position = hF(iax).Children.Title.Position - ...
        [0 axshift./diff(hF(iax).Children.Position([2,4])).*diff(ylim(hF(iax).Children)) 0];
  end
savefigauto(hF(2),fullfile(plotsavedir,'Figure6a'),'-vector');
savefigauto(hF(3),fullfile(plotsavedir,'Figure6b'),'-vector');
savefigauto(hF(4),fullfile(plotsavedir,'Figure6c'),'-vector');
savefigauto(hF(5),fullfile(plotsavedir,'Figure6d'),'-vector');
