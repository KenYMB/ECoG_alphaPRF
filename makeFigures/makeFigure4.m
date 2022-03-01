%% makeFigure4

% 20220223 Yuasa

%% Initialize
close all;

SetDefaultAnalysisPath;
plotsavepth    = fullfile(figPth, 'Publication');
issaveplot     = false;

%% Figure
plotsavedir    = fullfile(plotsavepth, 'Figure4');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%-- Spectra
ecog_APRF_10e_representative_pRFspectrum;

savefigauto(hF(1),fullfile(plotsavedir,'Figure4a_broadband'));
savefigauto(hF(2),fullfile(plotsavedir,'Figure4a_alpha'));

%-- Timecourse
ecog_APRF_10d_representative_pRF;

savefigauto(hF{1},fullfile(plotsavedir,'Figure4b_broadband'));
savefigauto(hF{2},fullfile(plotsavedir,'Figure4b_alpha'));

%-- pRF location
ecog_APRF_10d2_representative_pRF_boot;

savefigauto(hF{1},fullfile(plotsavedir,'Figure4c'));
