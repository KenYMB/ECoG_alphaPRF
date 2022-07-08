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
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

%-- modifiable parameters
SetDefault('skipexist',true);
SetDefault('allowlag',false);
SetDefault('average','runs');
SetDefault('smoothingMode','decimate');
SetDefault('smoothingN',3);
SetDefault('prfmodelList',{'linear'});
SetDefault('gaussianList', {'gs'});

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
opts.issave         = true;
opts.skipexist      = skipexist;
opts.average        = average;
opts.gammafit       = false;
opts.estimateIAF    = true;
opts.allownegfit    = true;
opts.allowbetafit   = false;
opts.allowwidefit   = false;
spcrm_params   = ecog_prf_fitalpha(freq, opts);

opts.allowbetafit   = true;
opts.allowwidefit   = false;
spcrm_paramsB  = ecog_prf_fitalpha(freq, opts);

opts.allowbetafit   = true;
opts.allowwidefit   = true;
spcrm_paramsBW = ecog_prf_fitalpha(freq, opts);

%% construct time series data for model fit
opts = [];
opts.outputDir      = 'pRFmodel';
opts.issave         = true;
opts.skipexist      = skipexist;
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;

opts.targetBAND     ='bbS';
modeldata_bb  = ecog_prf_constructTimeSeries(spcrm_params, opts);

opts.targetBAND     ='FaCLb';
modeldata_a   = ecog_prf_constructTimeSeries(spcrm_params, opts);

opts.targetBAND     ='FaCLbB';
modeldata_aB  = ecog_prf_constructTimeSeries(spcrm_paramsB, opts);

opts.targetBAND     ='FaCLbBW';
modeldata_aBW = ecog_prf_constructTimeSeries(spcrm_paramsBW, opts);

%% analyzePRF
opts = [];
opts.outputDir      = 'pRFmodel';
opts.issave         = true;
opts.skipexist      = skipexist;

for prfmodel=prfmodelList
    for gaussianmode=gaussianList

opts.prfmodel       = prfmodel{:};
opts.gaussianmode   = gaussianmode{:};

prf_params_bb   = ecog_prf_analyzePRF(modeldata_bb, opts);
prf_params_a    = ecog_prf_analyzePRF(modeldata_a, opts);
prf_params_aB   = ecog_prf_analyzePRF(modeldata_aB, opts);
prf_params_aBW  = ecog_prf_analyzePRF(modeldata_aBW, opts);

    end
end
