%% ECoG Alpha pRF
% ecog_APRF_01a_preprocessing
%   load ECoG data from BIDS files

% 20200126 reference: tde_run.m

%% Prepare parallel computation
isstartpar = false;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool; isstartpar = true;  end

%% Define paths and dataset
%-- Input & Output path
checkPath;
SetDefaultAnalysisPath;

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%% Load BIDS data
%-- automatically selected visual channels and GB channels
compute     = [];
inputDir    = fullfile(dataPth, 'ECoGCAR');
outputDir   = fullfile(savePth,'Raw');
sessionList = [];
epochTime   = [-1.0 1.5];
fsample     = 512;

[data] = ecog_prf_getData(compute, inputDir, outputDir, subjectList, sessionList, epochTime, fsample);

%% select epochs and channels
opts = [];
opts.baseline_time  = [-0.2 0];
opts.stim_on        = [0 0.5];
opts.issave         = true;
opts.doplots        = true;
opts.outputDir      = fullfile(savePth, 'Preprocessed');
opts.plotsavedir    = fullfile(figPth, 'dataselection');
opts.plot.XLim      = [-0.2 0.81];
opts.plot.fontSize  = 16;

[data] = ecog_prf_selectData(data,[],opts);

%% Regress out ERP 
opts = [];
opts.baseline_time  = [-0.2 0];
opts.issave         = true;
opts.doplots        = true;
opts.outputDir      = fullfile(savePth, 'Preprocessed');
opts.plotsavedir    = fullfile(figPth, 'regression');
opts.plot.XLim      = [-0.2 0.81];
opts.plot.fontSize  = 16;

[data] = ecog_prf_regressData(data, opts);

%% compute spectrum
opts = [];
opts.target_time    = [0 0.5];
opts.issave         = true;
opts.doplots        = true;
opts.doplots_alpha  = true;
opts.outputDir      = fullfile(savePth, 'Spectrum');
opts.plotsavedir    = fullfile(figPth, 'spectrum');
opts.plot.XLim      = [1 150];
opts.plot.fontSize  = 16;

[freq] = ecog_prf_spectra(data, opts);

%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
