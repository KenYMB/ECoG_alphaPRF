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
opts.doplots        = true;
opts.outputDir      = 'Preprocessed';
opts.plotsavedir    = 'dataselection';
% opts.plot.XLim      = [-0.2 0.81];
opts.plot.fontSize  = 16;
opts.issave         = true;
opts.skipexist      = skipexist;

[data] = ecog_prf_selectData(data,[],opts);

%% Regress out ERP 
opts = [];
opts.baseline_time  = [-0.2 0];
opts.epoch_time     = [-0.0 0.8];
opts.doplots        = true;
opts.outputDir      = 'Preprocessed';
opts.plotsavedir    = 'regression';
% opts.plot.XLim      = [-0.2 0.81];
opts.plot.fontSize  = 16;
opts.issave         = true;
opts.skipexist      = skipexist;

[data] = ecog_prf_regressData(data, opts);

%% compute spectrum
opts = [];
% opts.target_time    = [-0.15 0.85];     % 500ms hann window is 0.65 at t=0, 1.0 at t=100ms
opts.target_time    = [-0.2 0.8];       % 500ms hann window is 0.90 at t=0, 1.0 at t=50ms
opts.f_wintyp       = @hann;
opts.f_winlng       = 0.5;      % 0.5s (=Fs/2)
opts.f_winovl       = 0.75;     % 75% overlap
opts.doplots        = true;
opts.outputDir      = 'Spectrum';
opts.plotsavedir    = 'spectrum';
opts.plot.XLim      = [1 180];
opts.plot.fontSize  = 16;
opts.issave         = true;
opts.skipexist      = skipexist;

[freq] = ecog_prf_spectra(data, opts);

%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
