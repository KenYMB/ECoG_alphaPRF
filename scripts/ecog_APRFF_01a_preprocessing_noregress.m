%% ECoG Alpha pRF
% ecog_APRF_01a_preprocessing
%   load ECoG data from BIDS files

% 20200126 reference: tde_run.m
% 20220413 major update
% %% without ERP %%

%% Prepare parallel computation
isstartpar = false;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool([1 40]); isstartpar = true;  end

%% Define dataset
%-- Set path
checkPath;

%-- Working space
SetDefaultAnalysisPath('SHOW');

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

%-- modifiable parameters
SetDefault('skipexist',true);

%% Load BIDS data
%-- automatically selected visual channels and GB channels
compute     = [];
inputDir    = 'ECoGCAR';
outputDir   = 'Raw';
epochTime   = [-1.0 1.5];
fsample     = 512;

[data] = ecog_prf_getData(compute, subjectList, [], [], [], epochTime, fsample, inputDir, outputDir);

%% select epochs and channels
opts = [];
opts.baseline_time  = [-0.2 0];
opts.stim_on_time   = [0 0.5];
opts.epoch_time     = [-0.2 0.8];
opts.doplots        = false;
opts.outputDir      = 'Preprocessed';
opts.issave         = true;
opts.skipexist      = skipexist;

[data] = ecog_prf_selectData(data,[],opts);

%% compute spectrum
opts = [];
opts.target_time    = [0 0.5];
opts.doplots        = true;
opts.fileid         = 'freq_spectra_noregress';
opts.outputDir      = 'Spectrum';
opts.plotsavedir    = 'spectrum_noregress';
opts.plot.XLim      = [1 180];
opts.plot.fontSize  = 16;
opts.issave         = true;
opts.skipexist      = skipexist;

[freq] = ecog_prf_spectra(data, opts);

%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
