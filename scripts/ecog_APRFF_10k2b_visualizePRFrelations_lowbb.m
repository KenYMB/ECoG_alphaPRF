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
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'lowbb');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%-- Plotting Setting
FntSiz  = 18;
LFntSiz = 22;

%% Load dataset
%%% load analyzePRF & recompute full-ts R2
average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';
selectchs      = 'wangprobchs';
    allowlag       = false;
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;
    noACorrect     = false;

va_area = 'wangarea';
usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

%%% load low broadband
[modeldata_bbl, prf_params_bbl] = ecog_prf_loadprfs(subjectList,'bbL',prfPth,modeldataID,prfID,average,smoothingMode,smoothingN,prfmodel,gaussianmode,selectchs,selectch_exFEF,selectch_thresh,usefulltsR2,usefulltsxR2,skipsummarizeROIs);

%% Load correlations

% nboot     = 5000;
boottype  = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling
robustfit = true;

[prfs,prfs_corr,prfs_corrx,rois,nroi,nboot] = ...
    ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit,1,0);


%%% prepare visualize
% nsubjects = length(prf_params_bb);
% nroi  = length(rois);
selroi = (nroi-1):nroi;
% nvisarea  = length(prf_va_bb);
% areaList  = unique(prf_all_bb.channels.(va_area));

res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');
hF = gobjects(0);

%% Get bootstrap information
if ~issaveplot
precisefit = false;
else
precisefit = true;
end

if iscell(prfs.chanidx)         % is false, prfs.chanidx doesn't work as intended
  rearropt = {};
else
  ecog_APRFF_INITb_mergedata;
end
[model_all_bbl, prf_all_bbl]    = ecog_prf_mergeprfs(modeldata_bbl,prf_params_bbl,va_area,usexvalparams,rearropt);

if robustfit,  fitopt = {'RobustOpts','on'};
else,          fitopt = {'RobustOpts','off'};
end

prfs_corrx.ecc_rfsize_bbl_fit  = cell(nroi,1);
prfs_corrx.ecc_rfsize_bbl_corr = cell(nroi,1);
prfs_corrx.ecc_rfsize_bbl_R2   = cell(nroi,1);
for iroi=1:nroi
    
    %-- get bootstrap channel index
    if iscell(prfs.chanidx)
      ichx = [prfs.chanidx{iroi,:}];
    else
      chidx = [];
      chidx(:,1) = nearlyeq(prf_all_bb.ecc,prfs.ecc_bb{iroi}');
      chidx(:,2) = nearlyeq(prf_all_bb.ang,prfs.ang_bb{iroi}');
      chidx(:,3) = nearlyeq(prf_all_bb.rfsize,prfs.rfsize_bb{iroi}');
      chidx(:,4) = nearlyeq(prf_all_a.ecc,prfs.ecc_a{iroi}');
      chidx(:,5) = nearlyeq(prf_all_a.ang,prfs.ang_a{iroi}');
      chidx(:,6) = nearlyeq(prf_all_a.rfsize,prfs.rfsize_a{iroi}');
      ichx = mode(chidx,2);
    end
    
    %-- collect data
    prfs.R2_bbl{iroi,1}       = prf_all_bbl.xval(ichx)';
    prfs.ecc_bbl{iroi,1}      = prf_all_bbl.ecc(ichx)';
    prfs.ang_bbl{iroi,1}      = prf_all_bbl.ang(ichx)';
    prfs.rfsize_bbl{iroi,1}   = prf_all_bbl.rfsize(ichx)';
    
    %-- line fit
    if precisefit
        iprfs_num    = prfs.num(iroi,:);
        iprfs_ecc    = prfs.ecc_bbl{iroi};
        iprfs_rfsize = prfs.rfsize_bbl{iroi};
        iprfs_fit  = {};    iprfs_corr  = {};    iprfs_R2  = {};
        parfor iboot = 1:nboot
        dnum = iprfs_num(iboot);
        doffset = sum(iprfs_num(1:(iboot-1)));
          if dnum > 3
            iidx = (1:dnum)+doffset;
            fitparms = fitlm(iprfs_ecc(iidx),iprfs_rfsize(iidx),fitopt{:});
            iprfs_fit{iboot}   = fitparms.Coefficients.Estimate;
            iprfs_corr{iboot}  = sqrt(fitparms.Rsquared.Ordinary);
            iprfs_R2{iboot}    = fitparms.Rsquared.Adjusted;
          end
        end
        prfs_corrx.ecc_rfsize_bbl_fit{iroi}    = horzcat(iprfs_fit{:});
        prfs_corrx.ecc_rfsize_bbl_corr{iroi}   = horzcat(iprfs_corr{:});
        prfs_corrx.ecc_rfsize_bbl_R2{iroi}     = horzcat(iprfs_R2{:});
    else
        fitparms = fitlm(prfs.ecc_bbl{iroi},prfs.rfsize_bbl{iroi},fitopt{:});
        prfs_corrx.ecc_rfsize_bbl_fit{iroi}   = fitparms.Coefficients.Estimate;
        prfs_corrx.ecc_rfsize_bbl_corr{iroi}  = sqrt(fitparms.Rsquared.Ordinary);
        prfs_corrx.ecc_rfsize_bbl_R2{iroi}    = fitparms.Rsquared.Adjusted;
    end
end

%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations
% close all;

fitting = '-ellipse';

%% plot 2D histogram (replication around 90 deg for angle)
%% %%%%%%%%%%%%%%%%%%%%
%%% arrange in one figure
if issaveplot
    relt    = {'bb_bbl','a_bbl'};
    selroi  = {'V1-V3','Dorsolateral'};
else
    relt    = {'bb_bbl'};
    selroi  = 'V1-V3';
end
flds    = {'ang','ecc','rfsize'};
mdls    = {'bb','a','bbl'};
ModelNames = {'Broadband (70–180 Hz)','Alpha','Broadband (3–26 Hz)'};
opts = [];
opts.mothod      = 'hist';
opts.ang_limpad  = 90;       % degree for angle
opts.pix2deg     = cfactor;
opts.FontSize    = FntSiz;
opts.showfitted  = strtok(fitting,'-');
opts.fit_alpha   = alpha;
opts.condNames   = ModelNames;
h = ecog_prf_plotPRFrelations(prfs,prfs_corr,relt,flds,mdls,selroi,rois,opts);
hF = [hF h];
    
for iplt = 1:length(h)
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_2Dhist%s_%s_%s',threshold_bb,threshold_a,eclimit,fitting,'all',relt{iplt});
if issaveplot,   savefigauto(h(iplt),fullfile(plotsavedir, figname),'-vector');   end
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
mdls    = {'bb','a','bbl'};
ModelNames = {'Broadband (70–180 Hz)','Alpha','Broadband (3–26 Hz)'};
selroi  = {'V1','V2','V3','V1-V3','Dorsolateral'};
opts = [];
opts.mothod        = plotmethod;
opts.overplot      = true;
opts.pix2deg       = cfactor;
opts.FontSize      = LFntSiz;
opts.showfitted    = 'fitlm';
opts.fit_alpha     = alpha;
opts.fitlm_domean  = domean;
opts.fitlm_showref = showref;
opts.condNames   = ModelNames;

h = ecog_prf_plotPRFrelations(prfs,prfs_corrx,relt,mdls,selroi,rois,opts);
hF = [hF h];

if showref,     isrefname = '_withref';
else,           isrefname = '';
end
    
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_line-roi_%s%s-largeFont',threshold_bb,threshold_a,eclimit,relt,isrefname);
if domean,  figname = sprintf('%s-mean',figname);   end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));    end

end
