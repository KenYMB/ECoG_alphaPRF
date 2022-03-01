%% ECoG Alpha pRF
% ecog_APRF_01b_analyzePRF
%   pRF analysis

% 20200126 reference: tde_run.m
% 20210913 Yuasa: apply beta options based on subjects

%% Prepare parallel computation
isstartpar = false;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool; isstartpar = true;  end

%% Define paths and dataset
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');
isbetawide = ismember(alphaFitTypes(subjectList,'name'),'betawide');
%-- modifiable parameters
SetDefault('allowlag',false);

%% Load spctrum

opts = [];
opts.compute        = false;
opts.doplots        = false;
opts.allowlag       = allowlag;

[freq] = ecog_prf_spectra(subjectList, opts);

%% fit alpha

opts = [];
opts.target_time    = [0 0.5];
opts.issave         = true;
opts.average        = 'runs';
opts.gammafit       = false;
opts.estimateIAF    = true;
opts.allownegfit    = true;

spcrm_params = cell(size(freq));
if any(~isbetawide)
opts.allowbetafit   = false;
opts.allowwidefit   = false;
tmp = ecog_prf_fitalpha(freq(~isbetawide), opts);
spcrm_params(~isbetawide) = tmp;
end
if any(isbetawide)
opts.allowbetafit   = true;
opts.allowwidefit   = true;
tmp = ecog_prf_fitalpha(freq(isbetawide), opts);
spcrm_params(isbetawide) = tmp;
end

%% construct time series data for model fit
opts = [];
opts.smoothingMode  ='decimate';
opts.smoothingN     = 3;
opts.issave         = true;

opts.targetBAND     ='bbS';
modeldata_bb = ecog_prf_constructTimeSeries(spcrm_params, opts);

opts.targetBAND = cell(size(spcrm_params));
opts.targetBAND(~isbetawide)  ={'FaCLb'};
opts.targetBAND(isbetawide)   ={'FaCLbBW'};
modeldata_a = ecog_prf_constructTimeSeries(spcrm_params, opts);

%% analyzePRF
opts = [];
opts.prfmodel       = 'linear';
opts.gaussianmode   = 'gs';
opts.issave         = true;

prf_params_bb = ecog_prf_analyzePRF(modeldata_bb, opts);
prf_params_a  = ecog_prf_analyzePRF(modeldata_a, opts);

%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
