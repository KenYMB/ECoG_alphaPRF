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
SetDefault('isbetawide', ismember(alphaFitTypes(subjectList,'name'),'betawide'));

SetDefault('skip2prf',false);
SetDefault('computenoflip',false);

if ~skip2prf
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
opts.skipexist      = skipexist;
opts.average        = average;
opts.gammafit       = false;
opts.estimateIAF    = true;
opts.allownegfit    = true;
opts.allowbetafit   = false;
opts.allowwidefit   = false;
spcrm_params   = ecog_prf_fitalpha(freq, opts);
else
spcrm_params = subjectList;
end

%% construct time series data for model fit
opts = [];
opts.issave         = true;
opts.skipexist      = skipexist;
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;
    opts.compute    = ~skip2prf;
    opts.issave     = ~skip2prf;

opts.targetBAND     ='bbS';
modeldata_bb   = ecog_prf_constructTimeSeries(spcrm_params, opts);

opts.targetBAND = cell(size(spcrm_params));
opts.targetBAND(~isbetawide)  ={'FaCLb'};
opts.targetBAND(isbetawide)   ={'FaCLbBW'};
modeldata_a    = ecog_prf_constructTimeSeries(spcrm_params, opts);

opts.targetBAND     ='FaLb';
modeldata_aU   = ecog_prf_constructTimeSeries(spcrm_params, opts);

%-- no-flip for alpha
if computenoflip
    modeldata_a_N   = modeldata_a;
    modeldata_aU_N  = modeldata_aU;
    for isbj = 1:length(modeldata_a_N)
        modeldata_a_N{isbj}.targetBAND  = strrep(modeldata_a_N{isbj}.targetBAND,'b','Nb');
        modeldata_aU_N{isbj}.targetBAND = strrep(modeldata_aU_N{isbj}.targetBAND,'b','Nb');
        modeldata_a_N{isbj}.datats   = cellfun(@(D) -D, modeldata_a_N{isbj}.datats,...
            'UniformOutput',false);
        modeldata_aU_N{isbj}.datats  = cellfun(@(D) -D, modeldata_aU_N{isbj}.datats,...
            'UniformOutput',false);
    end
end

%% analyzePRF
opts = [];
opts.issave         = true;
opts.skipexist      = skipexist;
opts.fileid         = 'prfbidirgain';

opts.noneggain      = false;

for prfmodel=prfmodelList
    for gaussianmode=gaussianList

opts.prfmodel       = prfmodel{:};
opts.gaussianmode   = gaussianmode{:};

prf_params_bb   = ecog_prf_analyzePRF(modeldata_bb, opts);
prf_params_a    = ecog_prf_analyzePRF(modeldata_a, opts);
prf_params_aU   = ecog_prf_analyzePRF(modeldata_aU, opts);

%-- no-flip for alpha
if computenoflip
    prf_params_a_N    = ecog_prf_analyzePRF(modeldata_a_N, opts);
    prf_params_aU_N   = ecog_prf_analyzePRF(modeldata_aU_N, opts);
end

    end
end
