%% s1_analyzePRF
% Preprocessing ECoG signals
% Fitting pRF model on ECoG timecourse

% 20220223 Yuasa

%% Initialize
run_tbUse;
clearvars -except analysisRootPath
if exist('gcp','file') && isempty(gcp('nocreate')), parpool;  end

%% Set test dataset
%-- Input path
dataPth     = fullfile(bidsRootPath,'derivatives',filesep);
%-- Output path 
savePth     = fullfile(analysisRootPath,'Data',filesep);
figPth      = fullfile(analysisRootPath, 'Figures',filesep);
%-- Subject list
subjectList_fname = 'subjectlistTEST.tsv';

%% Compute power spectral density
ecog_APRF_01a_preprocessing

%% pRF analysis
ecog_APRF_01b_analyzePRF
ecog_APRF_01c_visualizePRF

%% Finish session
if exist('gcp','file'), delete(gcp('nocreate'));  end
