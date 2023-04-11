%% representative pRF fitting

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

%%% Spectra %%%
ecog_APRFF_10e_representative_pRFspectrum;

savefigauto(hF(1),fullfile(plotsavedir,[figureID 'a_broadband']),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'a_alpha']),'-vector');

%%% Timecourse %%%
ecog_APRFF_10d_representative_pRF;

figure(hF(1));    ylim([min(ylim) 850]);
yticks(0:200:max(ylim));    yticklabels(yticks/100+1);
savefigauto(hF(1),fullfile(plotsavedir,[figureID 'b_broadband']),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,[figureID 'b_alpha']),'-vector');

%%% pRF location %%%
ecog_APRFF_10d2_representative_pRF_boot;

savefigauto(hF(1),fullfile(plotsavedir,[figureID 'c']),'-vector');
