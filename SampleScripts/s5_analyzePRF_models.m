%% s5_analyzePRF_models
% Fitting supplementary pRF models on ECoG timecourse

% 20220223 Yuasa
% 20241223 Yuasa - updates for new figures

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

%% no smoothing
prfmodelList = {'linear'};
gaussianList = {'gs'};
smoothingMode = 'none';
smoothingN = 0;
ecog_APRFF_01b_analyzePRF;

%% pRF analysis for other properties
clear smoothingMode smoothingN
prfmodelList = {'linear'};
gaussianList = {'gs'};
ecog_APRFF_01bb2_analyzePRF_alphachange
ecog_APRFF_01bb2_analyzePRF_lowbroadband
ecog_APRFF_01bb2_analyzePRF_allownegativegain

%% pRF analysis with other models
clear smoothingMode smoothingN
% prfmodelList = {'linear','css'};
% gaussianList = {'dog','og','gs'};
prfmodelList = {'linear'};
gaussianList = {'og'};
ecog_APRFF_01bb_analyzePRF_allprf

%% pRF analysis without ERP removal
clear smoothingMode smoothingN
ecog_APRFF_01a_preprocessing_noregress
ecog_APRFF_01b_analyzePRF_noregress

%% Finish session
if exist('gcp','file'), delete(gcp('nocreate'));  end
