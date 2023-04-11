% compute connectivities across electrodes for subjects with HDgrid

% 20201111 Yuasa
% 20210515 Yuasa - modify 
% 20230307 Yuasa - for low-broadband
% 20230320 Yuasa - for wide-window

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
corrmethod = 'mscoh';                   % Compute ms-coherence
SetDefault('skipexist',true);

%% Load signal
opts = [];
opts.compute        = false;
opts.issave         = false;
opts.doplots        = false;
[data] = ecog_prf_regressData(subjectList, opts);

%% Compute high-freq-resolution Coherence
opts = [];
opts.outputDir      = 'xSpectrum';
% opts.target_time    = [-0.15 0.85];     % 500ms hann window is 0.65 at t=0, 1.0 at t=100ms
opts.target_time    = [-0.2 0.8];       % 500ms hann window is 0.90 at t=0, 1.0 at t=50ms
opts.f_winlng       = 0.5;    % 0.5s = Fs/2 
opts.f_winovl       = 3/4;    % 75%  = Fs/2*(3/4)
opts.channels       = 'GB*';
opts.inclauto       = true;
opts.issave         = true;
opts.skipexist      = skipexist;
opts.targetBAND     = 'all';
opts.calcmode       = 'coh';
opts.stimNames      = 'HORIZONTAL*';
[coh_all1] = ecog_prf_crossspectra(data, opts);
opts.stimNames      = 'VERTICAL*';
[coh_all2] = ecog_prf_crossspectra(data, opts);
opts.stimNames      = 'DIAGONAL*';
[coh_all3] = ecog_prf_crossspectra(data, opts);
opts.stimNames      = 'BLANK';
[coh_all0] = ecog_prf_crossspectra(data, opts);

opts.stimNames      = 'SHUFFLE';
[coh_sfl] = ecog_prf_crossspectra(data, opts);

%% Reshape data
opts = [];
opts.outputDir      = 'xSpectrum';
opts.average      = 'runs';
opts.method         = corrmethod;
opts.issave         = true;
opts.skipexist      = skipexist;
opts.isavgfreq      = true;
opts.passthrough    = true;
opts.targetBAND   = 'broadband';
opts.isavgfreq  = true;
[coh_bb1]   = ecog_prf_connectivity(coh_all1, opts);
[coh_bb2]   = ecog_prf_connectivity(coh_all2, opts);
[coh_bb3]   = ecog_prf_connectivity(coh_all3, opts);
[coh_bb0]   = ecog_prf_connectivity(coh_all0, opts);
opts.targetBAND   = 'alpha';
opts.isavgfreq  = false;
[coh_a1]    = ecog_prf_connectivity(coh_all1, opts);
[coh_a2]    = ecog_prf_connectivity(coh_all2, opts);
[coh_a3]    = ecog_prf_connectivity(coh_all3, opts);
[coh_a0]    = ecog_prf_connectivity(coh_all0, opts);
