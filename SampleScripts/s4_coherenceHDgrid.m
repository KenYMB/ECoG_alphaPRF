%% s4_coherenceHDgrid
% Compute inter-electrodes coherence in HD grids
% Visualize the coherence across distance

% 20220518 Yuasa

%% Initialize
run_checkPath;
clearvars;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool;  end

%% Set dataset
%-- Subject list
subjectList_fname = 'subjectlist.tsv';

%% Compute coherence
ecog_APRFF_03a_connectivityanalysis_widewin;

%% Compute error bar and fitting curve
ecog_APRFF_03d_CoherenceAcrossDist
ecog_APRFF_03e_bootstrapCoherence;
ecog_APRFF_03g_fitCoherence;

%% Visualize
% ecog_APRFF_03cA_connectivityanalysisTS;    % plot coherence across distance for alpha
% ecog_APRFF_03cB_connectivityanalysisTS;    % plot coherence across distance for broadband
% 
% ecog_APRFF_03f_plotbootstrapCoherence;

%% Finish session
if exist('gcp','file'), delete(gcp('nocreate'));  end
