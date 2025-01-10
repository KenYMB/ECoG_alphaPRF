%% effect of ERP-regression removal

% 20241122 Yuasa

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

%%%% Example electrodes
ecog_APRFF_10k5b_checkRegress_RefCh
%-- Adjust axes
hF(1).Children.Children(end).YLim = [-150 1300];
hF(2).Children.Children(end).YLim = [-150 1300];
set(hF(1).Children.Children(end).Children(1:9),'YData',hF(1).Children.Children(end).YLim);
set(hF(2).Children.Children(end).Children(1:9),'YData',hF(2).Children.Children(end).YLim);
hF(1).Children.Children(end-1).YLim = [-120 1050];
hF(2).Children.Children(end-1).YLim = [-120 1050];
set(hF(1).Children.Children(end-1).Children(1:9),'YData',hF(1).Children.Children(end-1).YLim);
set(hF(2).Children.Children(end-1).Children(1:9),'YData',hF(2).Children.Children(end-1).YLim);
hF(1).Children.Children(end-2).YLim = [-70 600];
hF(2).Children.Children(end-2).YLim = [-70 600];
set(hF(1).Children.Children(end-2).Children(1:9),'YData',hF(1).Children.Children(end-2).YLim);
set(hF(2).Children.Children(end-2).Children(1:9),'YData',hF(2).Children.Children(end-2).YLim);
hF(3).Children.Children(end).YLim = [-1.3 0.5];
hF(4).Children.Children(end).YLim = [-1.3 0.5];
set(hF(3).Children.Children(end).Children(1:9),'YData',hF(3).Children.Children(end).YLim);
set(hF(4).Children.Children(end).Children(1:9),'YData',hF(4).Children.Children(end).YLim);
hF(3).Children.Children(end-1).YLim = [-1.1 0.42];
hF(4).Children.Children(end-1).YLim = [-1.1 0.42];
set(hF(3).Children.Children(end-1).Children(1:9),'YData',hF(3).Children.Children(end-1).YLim);
set(hF(4).Children.Children(end-1).Children(1:9),'YData',hF(4).Children.Children(end-1).YLim);
hF(3).Children.Children(end-2).YLim = [-0.8 0.31];
hF(4).Children.Children(end-2).YLim = [-0.8 0.31];
set(hF(3).Children.Children(end-2).Children(1:9),'YData',hF(3).Children.Children(end-2).YLim);
set(hF(4).Children.Children(end-2).Children(1:9),'YData',hF(4).Children.Children(end-2).YLim);
%-- Update axis labels
for ifig = 1:4
for iax = reshape(hF(ifig).Children.Children(end+(-2:0)),1,[])
    yticklabels(iax,iax.YTick/100+1);
    if ifig<3
        set(iax,'YTickLabel',num2str(get(iax,'YTick')'./100+1));
    else
        if min(ylim(iax))>-1.2
            set(iax,'YTick',log10([1/8 1/4 1/2 1 2 4 8]));
            set(iax,'YTickLabel',{'1/8','1/4','1/2','1','2','4','8'});
        else
            set(iax,'YTick',log10([1/27 1/9 1/3 1 3 9 27]));
            set(iax,'YTickLabel',{'1/27','1/9','1/3','1','3','9','27'});
        end
    end
end
end
%-- Save fig
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'a_broadband_ts1']),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'a_broadband_ts2']),'-vector');
savefigauto(hF(5),fullfile(plotsavedir,[figureID 'a_broadband_loc']),'-vector');
savefigauto(hF(3),fullfile(plotsavedir,[figureID 'a_alpha_ts1']),'-vector');
savefigauto(hF(4),fullfile(plotsavedir,[figureID 'a_alpha_ts2']),'-vector');
savefigauto(hF(6),fullfile(plotsavedir,[figureID 'a_alpha_loc']),'-vector');


%%%% Summary
ecog_APRFF_10k5a_checkRegress
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'b']),'-vector');
