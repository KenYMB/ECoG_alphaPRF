%% ECoG Alpha pRF
% ecog_APRF_01a_preprocessing
%   load ECoG data from BIDS files

% 20210907 Yuasa: run for allbeta options
% %% without ERP %%

%% prefix
% close all; clear all;
if isempty(gcp('nocreate')),  parpool([1 40]); end
% startupToolboxToolbox;

%% Define paths and dataset
checkPath;

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%-- modifiable parameters
SetDefault('allowlag',false);
SetDefault('average','runs');
SetDefault('smoothingMode','decimate');
SetDefault('smoothingN',3);
SetDefault('prfmodelList',{'linear','css'});
SetDefault('gaussianlList', {'dog','og','gs'});

SetDefault('allowbetafit',false);
SetDefault('allowwidefit',false);
broadbandType = 'bbS';
alphaType     = 'FaCLb';
if allowbetafit, alphaType =[alphaType 'B']; end
if allowwidefit, alphaType =[alphaType 'W']; end

bootstrap      = 1000;

try
%% load time series
opts = [];
opts.outputDir      = 'pRFmodel';
opts.fileid         ='freq_spectra-timeseries-boot';
opts.compute        = false;
opts.allowlag       = allowlag;
opts.average        = average;
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;

% opts.targetBAND     = broadbandType;
% modeldata_bb = ecog_prf_constructTimeSeries(subjectList, opts);
opts.targetBAND     = alphaType;
modeldata_a   = ecog_prf_constructTimeSeries(subjectList, opts);

catch
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
opts.fileid         ='freq_spectra-params-boot';
opts.issave         = true;
opts.skipexist      = true;
opts.bootstrap      = bootstrap;
opts.average        = average;
opts.gammafit       = false;
opts.estimateIAF    = true;
opts.allownegfit    = true;
opts.allowbetafit   = allowbetafit;
opts.allowwidefit   = allowwidefit;
spcrm_params   = ecog_prf_fitalpha(freq, opts);

%% construct time series data for model fit
opts = [];
opts.outputDir      = 'pRFmodel';
opts.fileid         ='freq_spectra-timeseries-boot';
opts.issave         = true;
opts.skipexist      = true;
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;

% opts.targetBAND     = broadbandType;
% modeldata_bb  = ecog_prf_constructTimeSeries(spcrm_params, opts);
opts.targetBAND     = alphaType;
modeldata_a   = ecog_prf_constructTimeSeries(spcrm_params, opts);
end

%% analyzePRF
opts = [];
opts.outputDir      = 'pRFmodel';
opts.fileid         ='prf-boot';
opts.issave         = true;
opts.skipexist      = true;

for prfmodel=prfmodelList
    for gaussianmode=gaussianlList

opts.prfmodel       = prfmodel{:};
opts.gaussianmode   = gaussianmode{:};

% prf_params_bb   = ecog_prf_analyzePRF(modeldata_bb, opts);
prf_params_a    = ecog_prf_analyzePRF(modeldata_a, opts);

    end
end