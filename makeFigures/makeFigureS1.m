%% brain surfaces

% 20220223 Yuasa

%% Initialize
close all;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
figureIDs = strcat('Figure',strsplit(strrep(mfilename,'makeFigure',''),'_'));

figureID = figureIDs{1};
plotsavedir    = fullfile(plotsavePthP, figureID);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%%% Surface plots %%%
ecog_APRFF_10a_plotSurface;

for ii = 1:length(hF)
    sbjid   = strtok(hF(ii).Name,'_');
    figname = sprintf('%s_%s',figureID,sbjid);
    
    %-- fig
    view(hF(ii).Children(end),[0 0]);
    saveas(hF(ii),fullfile(plotsavedir,figname),'fig');

    %-- right
    view(hF(ii).Children(end),[30 0]);
    savefigauto(hF(ii),fullfile(plotsavedir,[figname '-right']),{'png','eps'},'Renderer','auto');

    %-- left
    view(hF(ii).Children(end),[-30 0]);
    savefigauto(hF(ii),fullfile(plotsavedir,[figname '-left']),{'png','eps'},'Renderer','auto');

end