%% ECoG Alpha pRF
% ecog_APRFF_01aa_preprocessing_lag
%   load ECoG data from BIDS files

% 20200126 reference: tde_run.m
% 20210811 Yuasa: allow lag in regression
% 20220413 major update
% %% without ERP %%

%% Prepare parallel computation
isstartpar = false;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool([1 40]); isstartpar = true;  end

%% Define dataset
%-- Set path
checkPath;

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

%% load raw data
opts = [];
opts.compute        = false;
opts.doplots        = false;
opts.outputDir      = 'Preprocessed';
opts.plotsavedir    = 'dataselection';
opts.fileid         = 'data_visualelecs';

[data] = ecog_prf_regressData(subjectList,opts);

%% Regress out ERP 
opts = [];
opts.baseline_time  = [-0.2 0];
opts.issave         = true;
opts.doplots        = true;
opts.outputDir      = 'Preprocessed';
opts.plotsavedir    = 'regression';
opts.plot.XLim      = [-0.2 0.81];
opts.plot.fontSize  = 16;

otps.allowlag       = true;
opts.maxlag         = 0.01;
opts.reference_time = [0 0.5];
opts.fileid         = 'data_regresslag';
[data] = ecog_prf_regressData(data, opts);

%% compute spectrum
opts = [];
opts.target_time    = [0 0.5];
opts.issave         = true;
opts.doplots        = true;
opts.doplots_alpha  = true;
opts.outputDir      = 'Spectrum';
opts.plotsavedir    = 'spectrum';
opts.plot.XLim      = [1 150];
opts.plot.fontSize  = 16;

opts.fileid         = 'freq_spectra_regresslag';
[freq] = ecog_prf_spectra(data, opts);

%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
