%% threshold for channels selection

% 20220223 Yuasa
% 20240527 Yuasa - update

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));

%%%% All patients
ecog_APRFF_10f_channelSelection;

%%% variance explained of alpha pRFs
figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(1),fullfile(plotsavedir,[figureID]),'-vector');

%%% bootstrap distribution
figureID = figureIDs{2};
plotsavedir    = fullfile(plotsavePthP, figureID);
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
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'a']),'-vector');
savefigauto(hF(3),fullfile(plotsavedir,[figureID 'b']),'-vector');
savefigauto(hF(4),fullfile(plotsavedir,[figureID 'c']),'-vector');
savefigauto(hF(5),fullfile(plotsavedir,[figureID 'd']),'-vector');


%%%% Individual patients
ecog_APRFF_10f3_channelSelection_ind;

%%% variance explained of alpha pRFs
figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'b']),'-vector');
