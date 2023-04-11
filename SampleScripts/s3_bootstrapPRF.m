%% s3_bootstrapPRF
% Apply bootstrap on pRF analysis

% 20220518 Yuasa
% 20230404 Yuasa - simplify

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

%% pRF analysis
% ecog_APRFF_01d_analyzePRFboot;
ecog_APRFF_10d0_analyzePRFboot_REP;

%% Finish session
if exist('gcp','file'), delete(gcp('nocreate'));  end
