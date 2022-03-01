%% ECoG Alpha pRF
% ecog_APRF_01c_visualizePRF
%   show pRF results

% 20200221 Yuasa
% 20220223 Yuasa: Updates

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

%% Load pRF results
clear alphaType broadbandType

average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';
selectchs      = 'allchs';
    allowmixbeta   = true;
va_area = 'wangarea';

ecog_APRF_INITa_loaddata;

%% plot pRF results
opts = [];
opts.plotsavedir = fullfile(figPth, 'pRF');
opts.closefig    = 'yes';
ecog_prf_plotPRFs(modeldata_bb,prf_params_bb,modeldata_a,prf_params_a,opts);


