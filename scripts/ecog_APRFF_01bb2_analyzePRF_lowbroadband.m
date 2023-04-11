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
SetDefault('allowlag',false);
SetDefault('average','runs');
SetDefault('smoothingMode','decimate');
SetDefault('smoothingN',3);
SetDefault('prfmodelList',{'linear'});
SetDefault('gaussianList', {'gs'});

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
opts.skipexist      = true;
opts.average        = average;
opts.gammafit       = false;
opts.estimateIAF    = true;
opts.allownegfit    = true;
opts.allowbetafit   = false;
opts.allowwidefit   = false;
spcrm_params   = ecog_prf_fitalpha(freq, opts);

%% construct time series data for model fit
opts = [];
opts.issave         = true;
opts.skipexist      = true;
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;

opts.targetBAND     ='bbL';
modeldata_bb   = ecog_prf_constructTimeSeries(spcrm_params, opts);

%% analyzePRF
opts = [];
opts.issave         = true;
opts.skipexist      = true;

for prfmodel=prfmodelList
    for gaussianmode=gaussianList

opts.prfmodel       = prfmodel{:};
opts.gaussianmode   = gaussianmode{:};

modeldata_bb    = ecog_prf_analyzePRF(modeldata_bb, opts);

    end
end
