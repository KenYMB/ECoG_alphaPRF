% compute connectivities across electrodes for subjects with HDgrid

% 20201111 Yuasa
% 20210515 Yuasa - modify 

%% Define paths and dataset
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
plotsavePth    = fullfile(figPth, 'pRFrelations');
if ~exist(plotsavePth,'dir'), mkdir(plotsavePth); end
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%% Load signal
opts = [];
opts.compute        = false;
opts.issave         = false;
opts.doplots        = false;
[data] = ecog_prf_regressData(subjectList, opts);

%% Compute cross-spectrum
opts = [];
opts.target_time    = [0 0.5];
opts.channels       = 'GB*';
opts.inclauto       = true;
opts.issave         = true;
opts.outputDir      = fullfile(savePth, 'xSpectrum');
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
opts.issave     = true;
opts.average    = 'runs';
opts.method     = 'mscoh';
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
