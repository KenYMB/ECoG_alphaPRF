%% s1_analyzePRF
% Preprocessing ECoG signals
% Fitting pRF model on ECoG timecourse

% 20220223 Yuasa

%% Initialize
run_checkPath;
clearvars;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool;  end

%% Set dataset
% %-- Input path
% bidsPth     = fullfile(bidsRootPath,'derivatives',filesep);
% %-- Output path 
% datPth      = fullfile(analysisRootPath,'Data',filesep);
% figPth      = fullfile(analysisRootPath, 'Figures',filesep);

%-- Subject list
subjectList_fname = 'subjectlist.tsv';

%% Compute power spectral density
ecog_APRFF_01a_preprocessing

%% pRF analysis
ecog_APRFF_01b_analyzePRF
ecog_APRFF_01c_visualizePRF

%% Finish session
if exist('gcp','file'), delete(gcp('nocreate'));  end
