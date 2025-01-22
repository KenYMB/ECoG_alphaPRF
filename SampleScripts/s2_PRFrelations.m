%% s2_PRFrelations
% Compute threshold of pRFs
% Visualize relations of broadband and alpha pRFs

% 20220518 Yuasa
% 20241223 Yuasa - updates for new figures

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

%% Compute threshold
ecog_APRFF_02a_distribution_xR2_halves        % estimate threshold
ecog_APRFF_02b_computePRFparamsFPM_bootall    % compute bootstrap

%% Compute selected channels
ecog_APRFF_02e_selectedchannels
ecog_APRFF_02e2_selectedchannels_ind

%% visualize relations
% ecog_APRFF_02c_visualizePRFrelations
% ecog_APRFF_02d_visualizePRFlocations_add

%% Finish session
if exist('gcp','file'), delete(gcp('nocreate'));  end
