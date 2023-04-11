%% s6_PRFrelations_models
% Compute threshold of pRFs estimated based on supplementary models

% 20220518 Yuasa

%% Initialize
run_checkPath;
clearvars;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool;  end

%% Set dataset
% %-- Input path
% prfPth         = 'pRFmodel';
% prfstatPth     = 'pRFanalysis';
% %-- Output path
% plotsavePth    = 'pRFrelations';

%-- Subject list
subjectList_fname = 'subjectlist.tsv';

%% Compute threshold with other models
clear allow*
noalphathresh = false;
noACorrect    = false;
prfmodel= 'linear';  gaussianmode = 'og';
ecog_APRFF_02a_distribution_xR2_halves          % estimate threshold
ecog_APRFF_02b_computePRFparamsFPM_bootall      % compute bootstrap

% prfmodel= 'linear';  gaussianmode = 'dog';
% ecog_APRFF_02a_distribution_xR2_halves          % estimate threshold
% ecog_APRFF_02b_computePRFparamsFPM_bootall      % compute bootstrap
% 
% prfmodel= 'css';  gaussianmode = 'gs';
% ecog_APRFF_02a_distribution_xR2_halves          % estimate threshold
% ecog_APRFF_02b_computePRFparamsFPM_bootall      % compute bootstrap
% 
% prfmodel= 'css';  gaussianmode = 'og';
% ecog_APRFF_02a_distribution_xR2_halves          % estimate threshold
% ecog_APRFF_02b_computePRFparamsFPM_bootall      % compute bootstrap
% 
% prfmodel= 'css';  gaussianmode = 'dog';
% ecog_APRFF_02a_distribution_xR2_halves          % estimate threshold
% ecog_APRFF_02b_computePRFparamsFPM_bootall      % compute bootstrap

%% Compute threshold for other properties
%-- without model for alpha computation
clear allow*
noalphathresh = false;
noACorrect    = true;
prfmodel= 'linear';  gaussianmode = 'gs';
ecog_APRFF_02a_distribution_xR2_halves          % estimate threshold
ecog_APRFF_02b_computePRFparamsFPM_bootall      % compute bootstrap

%% Compute prf relations w/o alpha threshold
clear allow*
noalphathresh = true;
noACorrect    = false;
prfmodel= 'linear';  gaussianmode = 'gs';
ecog_APRFF_02b_computePRFparamsFPM_bootall      % compute bootstrap

%% Finish session
if exist('gcp','file'), delete(gcp('nocreate'));  end
