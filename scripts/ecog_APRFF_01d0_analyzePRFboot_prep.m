%% ECoG Alpha pRF
% ecog_APRF_01a_preprocessing
%   load ECoG data from BIDS files

% 20210907 Yuasa: run for allbeta options
% %% without ERP %%

% See also: ecog_APRFF_01d0_analyzePRFboot_merge

%% Prepare parallel computation
isstartpar = false;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool([1 40]); isstartpar = true;  end

%% Define paths and dataset
checkPath;
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

%-- modifiable parameters
SetDefault('allowlag',false);
SetDefault('average','runs');
SetDefault('allowbetafit',contains(alphaFitTypes(subjectList,'name'),'beta'));
SetDefault('allowwidefit',contains(alphaFitTypes(subjectList,'name'),'wide'));
SetDefault('bootstrapindex',nan);

%% Load spctrum
opts = [];
opts.outputDir      = 'Spectrum';
opts.compute        = false;
opts.doplots        = false;
opts.allowlag       = allowlag;

[freq] = ecog_prf_spectra(subjectList, opts);

%% fit alpha
if isnan(bootstrapindex)
    fileidpostfix = '';
    numbootstrap = 1000;
else
    fileidpostfix = sprintf('-%d',bootstrapindex);
    nchans          = sum(cellfun(@(C) height(C.channels),freq));
    if     nchans<6,  numbootstrap = 1000;
    elseif nchans<10,  numbootstrap = 500;
    elseif nchans<20,  numbootstrap = 200;
    elseif nchans<40,  numbootstrap = 100;
    elseif nchans<80,  numbootstrap = 50;
    elseif nchans<160, numbootstrap = 20;
    elseif nchans<320, numbootstrap = 10;
    else,              numbootstrap = 5;
    end
end

opts = [];
opts.fileid         =['freq_spectra-params-boot' fileidpostfix];
opts.outputDir      = 'Spectrum';
opts.bootstrap      = numbootstrap;
opts.average        = average;
opts.gammafit       = false;
opts.estimateIAF    = true;
opts.allownegfit    = true;
opts.issave         = true;
opts.skipexist      = true;
opts.allowbetafit   = allowbetafit;
opts.allowwidefit   = allowwidefit;
spcrm_params = ecog_prf_fitalpha(freq, opts);

%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
