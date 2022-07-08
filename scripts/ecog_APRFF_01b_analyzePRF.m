%% ECoG Alpha pRF
% ecog_APRFF_01b_analyzePRF
%   load ECoG data from BIDS files
% 
% Prev, ecog_APRFF_01a_preprocessing

% 20200126 reference: tde_run.m
% 20210913 Yuasa: apply beta options based on subjects
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

%-- modifiable parameters
SetDefault('skipexist',true);
SetDefault('allowlag',false);
SetDefault('average','runs');
SetDefault('smoothingMode','decimate');
SetDefault('smoothingN',3);
SetDefault('prfmodel','linear');
SetDefault('gaussianmode', 'gs');
SetDefault('isbetawide', ismember(alphaFitTypes(subjectList,'name'),'betawide'));

%% Load spctrum

opts = [];
opts.outputDir      = 'Spectrum';
opts.compute        = false;
opts.doplots        = false;
opts.allowlag       = allowlag;

[freq] = ecog_prf_spectra(subjectList, opts);

%% fit alpha

opts = [];
opts.outputDir      = 'Spectrum';
opts.average        = average;
opts.gammafit       = false;
opts.estimateIAF    = true;
opts.allownegfit    = true;
opts.issave         = true;
opts.skipexist      = skipexist;
opts.allowbetafit   = isbetawide;
opts.allowwidefit   = isbetawide;
spcrm_params = ecog_prf_fitalpha(freq, opts);

%% construct time series data for model fit
opts = [];
opts.outputDir      = 'pRFmodel';
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;
opts.issave         = true;
opts.skipexist      = skipexist;

opts.targetBAND     ='bbS';
modeldata_bb = ecog_prf_constructTimeSeries(spcrm_params, opts);

opts.targetBAND = cell(size(spcrm_params));
opts.targetBAND(~isbetawide)  ={'FaCLb'};
opts.targetBAND(isbetawide)   ={'FaCLbBW'};
modeldata_a = ecog_prf_constructTimeSeries(spcrm_params, opts);

%% analyzePRF
opts = [];
opts.outputDir      = 'pRFmodel';
opts.prfmodel       = prfmodel;
opts.gaussianmode   = gaussianmode;
opts.issave         = true;
opts.skipexist      = skipexist;

prf_params_bb = ecog_prf_analyzePRF(modeldata_bb, opts);
prf_params_a  = ecog_prf_analyzePRF(modeldata_a, opts);

%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
