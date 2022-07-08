%% ECoG Alpha pRF
% ecog_APRFF_01c_visualizePRF
%   show pRF results

% 20200221 Yuasa
% 20220223 Yuasa: Updates

%% Define dataset
%-- Set path
checkPath;

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

%% Load pRF results
clear alphaType broadbandType

%-- modifiable parameters
SetDefault('allowlag',false);

SetDefault('average','runs');
SetDefault('smoothingMode','decimate');
SetDefault('smoothingN',3);
SetDefault('prfmodel','linear');
SetDefault('gaussianmode', 'gs');

selectchs      = 'allchs';
    allowmixbeta   = true;
va_area = 'wangarea';

ecog_APRFF_INITa_loaddata;

%% plot pRF results
opts = [];
opts.outputDir   = 'pRFmodel';
opts.plotsavedir = 'pRF';
opts.closefig    = 'yes';
ecog_prf_plotPRFs(modeldata_bb,prf_params_bb,modeldata_a,prf_params_a,opts);


