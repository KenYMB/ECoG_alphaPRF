%% s3_bootstrapPRF
% Apply bootstrap on pRF analysis

% 20220518 Yuasa

%% Initialize
run_checkPath;
clearvars;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool;  end

%% Set dataset
% %-- Input path
% bidsPth     = fullfile(bidsRootPath,'derivatives',filesep);
% %-- Output path 
% datPth      = fullfile(analysisRootPath,'Data',filesep);
% figPth      = fullfile(analysisRootPath, 'Figures',filesep);

%-- Subject list
subjectList_fname = 'subjectlist.tsv';

%-- Set numbers of segmentation of bootstrapping
% nbootloop = [5 10 20 20 10 5 10 50 50];
nbootloop = [5 10 10 10 10 5 5 50 50];
%{
%-- Modify nbootloop in the next cell refering numbers of channels
allowlag    = false;
subjectList = SetSubjectsList(subjectList_fname, 'all');

opts = [];
opts.outputDir      = 'Spectrum';
opts.compute        = false;
opts.doplots        = false;
opts.allowlag       = allowlag;
[freq] = ecog_prf_spectra(subjectList, opts);

nchans = cellfun(@(D) height(D.channels),freq)';
fprintf('# of channels are');
fprintf(' %d',nchans);
fprintf('\n');

nbootloop = getelement([1 2 5 10 20 50 100 200] ,sum(double(nchans'>=[6 10 20 40 80 160 320]),2)+1)';
%}
%% Compute power spectral density
for selsbj = 1:length(nbootloop)
    parfor bootstrapindex = 1:nbootloop(selsbj)
        ecog_APRFF_01d0_analyzePRFboot_prep;
    end
end

%% pRF analysis
ecog_APRFF_01d0_analyzePRFboot_merge;
ecog_APRFF_01d_analyzePRFboot;

%% Finish session
if exist('gcp','file'), delete(gcp('nocreate'));  end
