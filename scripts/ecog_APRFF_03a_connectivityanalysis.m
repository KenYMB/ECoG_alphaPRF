% compute connectivities across electrodes for subjects with HDgrid

% 20201111 Yuasa
% 20210515 Yuasa - modify 

%%
% close all; clear all;
if isempty(gcp('nocreate')),  parpool([4 10]); end
% startupToolboxToolbox;

%% Define paths and dataset
checkPath;
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
if isnumeric(selsbj)
subjectList = SetSubjectsList(subjectList_fname, selsbj);
else
subjectList = SetSubjectsList(subjectList_fname, 'hasHDgrid','yes');
end

%-- modifiable parameters
SetDefault('skipexist',true);

%% Load signal
opts = [];
opts.compute        = false;
opts.issave         = false;
opts.doplots        = false;
[data] = ecog_prf_regressData(subjectList, opts);

%% Compute cross-spectrum
opts = [];
opts.outputDir      = 'xSpectrum';
opts.target_time    = [0 0.5];
opts.channels       = 'GB*';
opts.inclauto       = true;
opts.issave         = true;
opts.skipexist      = skipexist;
opts.targetBAND     = 'broadband';
opts.stimNames      = 'HORIZONTAL*';
[xfreq_bb1] = ecog_prf_crossspectra(data, opts);
opts.stimNames      = 'VERTICAL*';
[xfreq_bb2] = ecog_prf_crossspectra(data, opts);
opts.stimNames      = 'DIAGONAL*';
[xfreq_bb3] = ecog_prf_crossspectra(data, opts);
opts.stimNames      = 'BLANK';
[xfreq_bb0] = ecog_prf_crossspectra(data, opts);

opts.targetBAND     = 'alpha';
opts.stimNames      = 'HORIZONTAL*';
[xfreq_a1] = ecog_prf_crossspectra(data, opts);
opts.stimNames      = 'VERTICAL*';
[xfreq_a2] = ecog_prf_crossspectra(data, opts);
opts.stimNames      = 'DIAGONAL*';
[xfreq_a3] = ecog_prf_crossspectra(data, opts);
opts.stimNames      = 'BLANK';
[xfreq_a0] = ecog_prf_crossspectra(data, opts);

%% Compute Coherence
opts = [];
opts.outputDir      = 'xSpectrum';
opts.average    = 'runs';
opts.method     = 'mscoh';
opts.issave         = true;
opts.skipexist      = skipexist;
opts.isavgfreq  = true;
[coh_bb1]   = ecog_prf_connectivity(xfreq_bb1, opts);
[coh_bb2]   = ecog_prf_connectivity(xfreq_bb2, opts);
[coh_bb3]   = ecog_prf_connectivity(xfreq_bb3, opts);
[coh_bb0]   = ecog_prf_connectivity(xfreq_bb0, opts);
opts.isavgfreq  = false;
[coh_a1]    = ecog_prf_connectivity(xfreq_a1, opts);
[coh_a2]    = ecog_prf_connectivity(xfreq_a2, opts);
[coh_a3]    = ecog_prf_connectivity(xfreq_a3, opts);
[coh_a0]    = ecog_prf_connectivity(xfreq_a0, opts);

opts.average    = 'stimuli';
opts.isavgfreq  = true;
[coh_bb1]   = ecog_prf_connectivity(xfreq_bb1, opts);
[coh_bb2]   = ecog_prf_connectivity(xfreq_bb2, opts);
[coh_bb3]   = ecog_prf_connectivity(xfreq_bb3, opts);
[coh_bb0]   = ecog_prf_connectivity(xfreq_bb0, opts);
opts.isavgfreq  = false;
[coh_a1]    = ecog_prf_connectivity(xfreq_a1, opts);
[coh_a2]    = ecog_prf_connectivity(xfreq_a2, opts);
[coh_a3]    = ecog_prf_connectivity(xfreq_a3, opts);
[coh_a0]    = ecog_prf_connectivity(xfreq_a0, opts);
