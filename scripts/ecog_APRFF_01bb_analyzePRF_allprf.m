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
SetDefault('prfmodelList',{'linear','css'});
SetDefault('gaussianList', {'dog','og','gs'});
SetDefault('isbetawide', ismember(alphaFitTypes(subjectList,'name'),'betawide'));

%% Load time series data for model fit

opts = [];
opts.compute        = false;
opts.allowlag       = allowlag;
opts.average        = average;
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;

opts.targetBAND     ='bbS';
modeldata_bb = ecog_prf_constructTimeSeries(subjectList, opts);

opts.targetBAND = cell(size(subjectList));
opts.targetBAND(~isbetawide)  ={'FaCLb'};
opts.targetBAND(isbetawide)   ={'FaCLbBW'};
modeldata_a = ecog_prf_constructTimeSeries(subjectList, opts);

%% analyzePRF
% prfmodelList  = {'linear','fixexpt','css'};
% gaussianlList = {'gs','dog','og'};

opts = [];
opts.issave         = true;
opts.skipexist      = skipexist;

for prfmodel=prfmodelList
    for gaussianmode=gaussianList

opts.prfmodel       = prfmodel{:};
opts.gaussianmode   = gaussianmode{:};

prf_params_bb   = ecog_prf_analyzePRF(modeldata_bb, opts);
prf_params_a    = ecog_prf_analyzePRF(modeldata_a, opts);

    end
end
