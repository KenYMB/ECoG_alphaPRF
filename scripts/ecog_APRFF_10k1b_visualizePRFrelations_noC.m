% visualize pRF correlations computed in 06c (equivalant 03c2 + 03d2)
%   Update xR2 with full time-series (no decimation)
%   Folk from ecog_APRF_06d_visualizePRFparamsFPM_bootall_fullts.m
%         and ecog_APRF_06e_visualizePRFsFPM_bootall_fullts.m

% 20210107 Yuasa - Update from 06d,06e
% 20210809 Yuasa - modify for paper
% 20210916 Yuasa - for APRFF

%% prefix
% close all; clearvars;
% startupToolboxToolbox;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'pRFrelations-representative');
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
else
plotsavePth    = 'pRFrelations-representative';
prfPth         = 'pRFmodel';
prfstatPth     = 'pRFanalysis';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'noACorr');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%-- Plotting Setting
FntSiz  = 18;
LFntSiz = 22;

%% Load dataset
clear alphaType broadbandType

%%% load analyzePRF & recompute full-ts R2
average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';
selectchs      = 'wangprobchs';
    allowlag       = false;
    allowbeta      = false;
    allowwide      = false;
    noACorrect     = true;

va_area = 'wangarea';
usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

%% Load correlations
% nboot     = 5000;
boottype  = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling
robustfit = true;

[prfs,prfs_corr,prfs_corrx,rois] = ...
    ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit,1,0);

%%% prepare visualize
% nsubjects = length(prf_params_bb);
% nroi  = length(rois);
% selroi = (nroi-1):nroi;
% nvisarea  = length(prf_va_bb);
% areaList  = unique(prf_all_bb.channels.(va_area));

res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');
hF = gobjects(0);

%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations
% close all;

% fitting = '-ellipse';
fitting = '-robustellipse';     % exclude outliers

%% plot 2D histogram (replication around 90 deg for angle)
%% %%%%%%%%%%%%%%%%%%%%
%%% arrange in one figure
if issaveplot

relt    = 'bb_a';
flds    = {'ang','ecc','rfsize'};
selroi  = {'V1-V3','Dorsolateral'};
opts = [];
opts.mothod      = 'hist';
opts.ang_limpad  = 90;       % degree for angle
opts.pix2deg     = cfactor;
opts.FontSize    = FntSiz;
opts.showfitted  = strtok(fitting,'-');
opts.fit_alpha   = alpha;
h = ecog_prf_plotPRFrelations(prfs,prfs_corr,relt,flds,selroi,rois,opts);
hF = [hF h];
    
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_2Dhist%s_%s',threshold_bb,threshold_a,eclimit,fitting,'all');
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');   end

end

%% %%%%%%%%%%%%%%%%%%%%
%%% cross-correlations
% close all;

showref = true;
domean  = true;

%% plot fitted line in ROIs
if issaveplot

plotmethod = 'line';

relt    = 'ecc_rfsize';
flds    = {'bb','a'};
selroi  = {'V1-V3','Dorsolateral'};
opts = [];
opts.mothod        = plotmethod;
opts.overplot      = true;
opts.pix2deg       = cfactor;
opts.FontSize      = LFntSiz;
opts.showfitted    = 'fitlm';
opts.fit_alpha     = alpha;
opts.fitlm_domean  = domean;
opts.fitlm_showref = showref;

h = ecog_prf_plotPRFrelations(prfs,prfs_corrx,relt,flds,selroi,rois,opts);
hF = [hF h];

if showref,     isrefname = '_withref';
else,           isrefname = '';
end
    
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_line-roi_%s%s-largeFont',threshold_bb,threshold_a,eclimit,relt,isrefname);
if domean,  figname = sprintf('%s-mean',figname);   end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));    end

end


%% %%%%%%%%%%%%%%%%%%%%
%% Compare models
%% %%%%%%%%%%%%%%%%%%%%
if ~issaveplot
%% Load results
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
prf_all_a_noAC  = prf_all_a;
%% Load prfs
%-- Set parameter for main data
clear alphaType broadbandType
    allowmixbeta   = true;
    noACorrect     = false;
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

%-- Load main data
% nboot         = 5000;
boottype      = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling
robustfit     = true;
noalphathresh = true;
noeccthresh   = true;          % if true, avoid to apply eccentricity restriction

if noalphathresh,       threshold_a = nan;  end   % just apply broadband threshold
[prfsO,prfsO_corr,prfsO_corrx,rois,nroi] = ...
    ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit,1,1);

%% Update prfs to plot same electrodes

prfsO.gain_a  = cell(size(prfsO.chanidx,1),1);
prfsO.gain_bb = cell(size(prfsO.chanidx,1),1);

prfsU = prfsO;
prfsU_corr = prfsO_corr;
for iroi=1:nroi
    
    %-- get bootstrap channel index
    if iscell(prfsO.chanidx)
      ichx = [prfsO.chanidx{iroi,:}];
    else
      chidx = [];
      chidx(:,1) = nearlyeq(prf_all_bb.ecc,prfsO.ecc_bb{iroi}');
      chidx(:,2) = nearlyeq(prf_all_bb.ang,prfsO.ang_bb{iroi}');
      chidx(:,3) = nearlyeq(prf_all_bb.rfsize,prfsO.rfsize_bb{iroi}');
      chidx(:,4) = nearlyeq(prf_all_a.ecc,prfsO.ecc_a{iroi}');
      chidx(:,5) = nearlyeq(prf_all_a.ang,prfsO.ang_a{iroi}');
      chidx(:,6) = nearlyeq(prf_all_a.rfsize,prfsO.rfsize_a{iroi}');
      ichx = mode(chidx,2);
    end
    
    %-- collect data
    prfsU.R2_a{iroi,1}       = prf_all_a_noAC.xval(ichx)';
    prfsU.ecc_a{iroi,1}      = prf_all_a_noAC.ecc(ichx)';
    prfsU.ang_a{iroi,1}      = prf_all_a_noAC.ang(ichx)';
    prfsU.rfsize_a{iroi,1}   = prf_all_a_noAC.rfsize(ichx)';

    %-- collect gain
    prfsO.gain_a{iroi,1}     = prf_all_a.gain(ichx)';
    prfsO.gain_bb{iroi,1}    = prf_all_bb.gain(ichx)';
    prfsU.gain_a{iroi,1}     = prf_all_a_noAC.gain(ichx)';
    prfsU.gain_bb{iroi,1}    = prfsO.gain_bb{iroi,1};

end

%% Merge prfs
tarroi = 'V1-V3';
iroi   = find(ismember(rois,tarroi),1);

prfsC = [];
paramnames = fieldnames(prfsO);
for iprm = reshape(paramnames,1,[])
    if size(prfsO.(iprm{:}),1) == nroi
        prfsC.(iprm{:}) = vertcat(prfsO.(iprm{:})(iroi,:),prfsU.(iprm{:})(iroi,:));
    else
        prfsC.(iprm{:}) = prfsO.(iprm{:});
    end
end

prfsC_corr = [];
% paramnames = fieldnames(prfsO_corr);
% for iprm = reshape(paramnames,1,[])
%     if size(prfsO_corr.(iprm{:}),1) == nroi
%         prfsC_corr.(iprm{:}) = vertcat(prfsO_corr.(iprm{:})(iroi,:),prfsU_corr.(iprm{:})(iroi,:));
%     else
%         prfsC_corr.(iprm{:}) = prfsO_corr.(iprm{:});
%     end
% end

%-- Exclude electrodes with too large eccentricity
if ~noeccthresh
% badch = (prfsC.ecc_a{1}>eclimit)|(prfsC.ecc_a{2}>eclimit);
for iroi = 1:2
badch = (prfsC.ecc_a{iroi}>eclimit);
for iprm = reshape(paramnames,1,[])
    if iscell(prfsC.(iprm{:})) && size(prfsC.(iprm{:}){iroi},2) == length(badch)
        prfsC.(iprm{:}){iroi}(:,badch) = [];
    end
end
end
end

%-- Set ROI names
roisC = {sprintf('Model-based\nAlpha Suppression');sprintf('Power Change\nin Alpha')};

%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations
% close all;

% fitting = '-ellipse';
fitting = '-robustellipse';     % exclude outliers

%% plot 2D histogram (replication around 90 deg for angle)
%% %%%%%%%%%%%%%%%%%%%%
%%% arrange in one figure

relt    = 'bb_a';
flds    = {'ang','ecc','rfsize'};
opts = [];
opts.mothod      = 'hist';
opts.ang_limpad  = 90;       % degree for angle
opts.pix2deg     = cfactor;
opts.FontSize    = FntSiz;
opts.showfitted  = strtok(fitting,'-');
opts.fit_alpha   = alpha;
h = ecog_prf_plotPRFrelations(prfsC,prfsC_corr,relt,flds,roisC,roisC,opts);
hF = [hF h];
    
figname = sprintf('prf-modelcomparison_2Dhist%s_%s',fitting,'all');
if noeccthresh, figname = sprintf('%s-noecc',figname);  end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');   end


%% plot gain in 2D histogram
%%% arrange in one figure

relt    = 'bb_a';
flds    = {'gain'};
opts = [];
opts.mothod      = 'hist';
opts.ang_limpad  = 90;       % degree for angle
opts.pix2deg     = cfactor;
opts.FontSize    = FntSiz;
opts.showfitted  = strtok(fitting,'-');
opts.fit_alpha   = alpha;
opts.showdiag = false;
h = ecog_prf_plotPRFrelations(prfsC,prfsC_corr,relt,flds,roisC,roisC,opts);
hF = [hF h];
    
figname = sprintf('prf-modelcomparison_gain2Dhist%s_%s',fitting,'all');
if noeccthresh, figname = sprintf('%s-noecc',figname);  end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');   end

end
