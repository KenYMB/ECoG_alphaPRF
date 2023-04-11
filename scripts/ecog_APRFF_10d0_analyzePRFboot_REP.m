%% ECoG Alpha pRF (simple broadband computation)
% ecog_APRF_01gg_analyzePRFboot
%   load ECoG spectra & compute broadband and alpha power in timecourse
%   compute broadband power without linear regression
%   with bootstrapping

% 20210916 - Yuasa
% 20230404 - Yuasa: reduce channels
% 
% %% without ERP %%

%% Prepare parallel computation
isstartpar = false;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool([1 40]); isstartpar = true;  end

%% Define paths and dataset
checkPath;

repelec = table({'Oc17','Oc18','GB102','GB103'}',{'p02','p02','p10','p10'}','VariableNames',{'name','subject_name'});

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
% SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, 'all');
selsbj = find(ismember(subjectList,repelec.subject_name));
subjectList = subjectList(selsbj);

%-- modifiable parameters
SetDefault('allowlag',false);
SetDefault('average','runs');
SetDefault('smoothingMode','decimate');
SetDefault('smoothingN',3);
SetDefault('prfmodel','linear');
SetDefault('gaussianmode', 'gs');
SetDefault('isbetawide', ismember(alphaFitTypes(subjectList,'name'),'betawide'));
% numbootstrap = 1000;
numbootstrap = 200;

%% Load spctrum

opts = [];
opts.outputDir      = 'Spectrum';
opts.compute        = false;
opts.doplots        = false;
opts.allowlag       = allowlag;

[freq] = ecog_prf_spectra(subjectList, opts);

%-- select representative data
for isbj = 1:length(freq)
    repchs = ismember(freq{isbj}.channels.name,repelec.name) &...
             ismember(freq{isbj}.channels.subject_name,repelec.subject_name);
    freq{isbj}.channels(~repchs,:)      = [];
    freq{isbj}.spectra(~repchs,:,:)     = [];
    freq{isbj}.spectra_off(~repchs,:,:) = [];
end

%% fit alpha
opts = [];
opts.fileid         = 'freq_spectra-params-boot';
opts.outputDir      = 'Spectrum';
opts.bootstrap      = numbootstrap;
opts.average        = average;
opts.gammafit       = false;
opts.estimateIAF    = true;
opts.allownegfit    = true;
opts.issave         = true;
opts.skipexist      = true;
opts.allowbetafit   = isbetawide;
opts.allowwidefit   = isbetawide;
spcrm_params = ecog_prf_fitalpha(freq, opts);

%% time series data for model fit
opts = [];
opts.fileid         = 'freq_spectra-timeseries-bootrep';
opts.outputDir      = 'pRFmodel';
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;
opts.issave         = true;
opts.skipexist      = true;

opts.targetBAND     ='bbS';
modeldata_bb = ecog_prf_constructTimeSeries(spcrm_params, opts);

opts.targetBAND = cell(size(spcrm_params));
opts.targetBAND(~isbetawide)  ={'FaCLb'};
opts.targetBAND(isbetawide)   ={'FaCLbBW'};
modeldata_a = ecog_prf_constructTimeSeries(spcrm_params, opts);

%% analyzePRF
opts = [];
opts.fileid         = 'prf-bootrep';
opts.outputDir      = 'pRFmodel';
opts.prfmodel       = prfmodel;
opts.gaussianmode   = gaussianmode;
opts.issave         = true;
opts.skipexist      = true;

prf_params_bb = ecog_prf_analyzePRF(modeldata_bb, opts);
prf_params_a = ecog_prf_analyzePRF(modeldata_a, opts);

%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
