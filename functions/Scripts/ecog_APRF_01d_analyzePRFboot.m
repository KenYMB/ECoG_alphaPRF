%% ECoG Alpha pRF (simple broadband computation)
% ecog_APRF_01d_analyzePRFboot
%   load ECoG spectra & compute broadband and alpha power in timecourse
%   compute broadband power without linear regression
%   with bootstrapping

% 20210916 - Yuasa

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
opts.fileid         ='freq_spectra-params-boot';
opts.issave         = true;
opts.skipexist      = true;
opts.bootstrap      = 1000;
opts.target_time    = [0 0.5];
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

%% time series data for model fit
opts = [];
opts.fileid         ='freq_spectra-timeseries-boot';
opts.issave         = true;
opts.skipexist      = true;
opts.smoothingMode  ='decimate';
opts.smoothingN     = 3;

opts.targetBAND     ='bbS';
modeldata_bb = ecog_prf_constructTimeSeries(spcrm_params, opts);

opts.targetBAND = cell(size(spcrm_params));
opts.targetBAND(~isbetawide)  ={'FaCLb'};
opts.targetBAND(isbetawide)   ={'FaCLbBW'};
modeldata_a = ecog_prf_constructTimeSeries(spcrm_params, opts);

%% analyzePRF
opts = [];
opts.fileid         ='prf-boot';
opts.issave         = true;
opts.skipexist      = true;
opts.prfmodel       = 'linear';
opts.gaussianmode   = 'gs';

prf_params_bb = ecog_prf_analyzePRF(modeldata_bb, opts);
prf_params_a = ecog_prf_analyzePRF(modeldata_a, opts);


%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
