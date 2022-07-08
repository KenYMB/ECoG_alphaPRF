%% makeFigure4

% 20220223 Yuasa

%% Initialize
close all;
clear modeldataID prfID;

run_checkPath;
plotsavePthP   = SetDefaultAnalysisPath('FIG','Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavePthP, 'Figure4');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%-- Spectra
ecog_APRFF_10e_representative_pRFspectrum;

savefigauto(hF(1),fullfile(plotsavedir,'Figure4a_broadband'),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,'Figure4a_alpha'),'-vector');

%-- Timecourse
ecog_APRFF_10d_representative_pRF;

figure(hF(1));    ylim([min(ylim) 850]);
savefigauto(hF(1),fullfile(plotsavedir,'Figure4b_broadband'),'-vector');
savefigauto(hF(2),fullfile(plotsavedir,'Figure4b_alpha'),'-vector');

%-- pRF location
ecog_APRFF_10d2_representative_pRF_boot;

savefigauto(hF(1),fullfile(plotsavedir,'Figure4c'),'-vector');
