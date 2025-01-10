% visualize pRF cross-correlation to xval
%   Folk from 10g

% 20241120 Yuasa - compute correlation between xval and rfsize

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
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth));
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%-- Plotting Setting
FntSiz  = 18;
LFntSiz = 22;

%%
%%% load analyzePRF & recompute full-ts R2
clear alphaType broadbandType
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

va_area = 'wangarea';
usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

%% Load correlations

% nboot     = 5000;
boottype  = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling
robustfit = true;

[prfs,prfs_corr,prfs_corrx,rois,nroi] = ...
    ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit);


%%% prepare visualize
% nsubjects = length(prf_params_bb);
% nroi  = length(rois);
% selroi = (nroi-1):nroi;
% nvisarea  = length(prf_va_bb);
% areaList  = unique(prf_all_bb.channels.(va_area));

res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

%-- Obtain channel index
if ~iscell(prfs.chanidx)
  prfs.chanidx = {};
  for iroi = 1:nroi
      %-- get channel index
      chidx = [];
      chidx(:,1) = nearlyeq(prf_all_bb.ecc,prfs.ecc_bb{iroi}');
      chidx(:,2) = nearlyeq(prf_all_bb.ang,prfs.ang_bb{iroi}');
      chidx(:,3) = nearlyeq(prf_all_bb.rfsize,prfs.rfsize_bb{iroi}');
      chidx(:,4) = nearlyeq(prf_all_a.ecc,prfs.ecc_a{iroi}');
      chidx(:,5) = nearlyeq(prf_all_a.ang,prfs.ang_a{iroi}');
      chidx(:,6) = nearlyeq(prf_all_a.rfsize,prfs.rfsize_a{iroi}');
      ichx = mode(chidx,2)';
      prfs.chanidx{iroi,1} = ichx;
  end
end

%-- ROIs of interest
if issaveplot
    roioi = [1,2,3,nroi-1,nroi];        % V1,V2,V3,low,high
else
    roioi = [nroi-1 nroi];              % low,high
end

%% Set xval params
if ~isfield(prfs,'xval_a')
    for iroi=1:size(prfs.chanidx,1)
        prfs.xval_bb{iroi,1} = prf_all_bb.xval([prfs.chanidx{iroi,:}])';
        prfs.xval_a{iroi,1}  = prf_all_a.xval([prfs.chanidx{iroi,:}])';
    end
end

%-- Compute cross-correlation between xval and size (<3min)
if ~isfield(prfs_corrx,'xval_rfsize_a_corr')
    nboot = size(prfs.chanidx,2);
    fitopt = {'RobustOpts','on'};

    flds  = {'xval_rfsize'};
 
     for sbfld = flds
       prfs_corrx.([sbfld{:} '_bb' '_fit'])    = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_bb' '_corr'])   = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_bb' '_R2'])     = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_a' '_fit'])     = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_a' '_corr'])    = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_a' '_R2'])      = cell(nroi,1);
     end

    fprintf(1,'Computing');
    warnst = warning('off');
    try
    for iroi = roioi
        for iboot = 1:nboot
            dnum = prfs.num(iroi,iboot);
            doffset = sum(prfs.num(iroi,1:(iboot-1)));
            if dnum > 3
                for sbfld = flds
                    sbflds = strsplit(sbfld{:},'_');
                    dat1 = prfs.([sbflds{1} '_bb']){iroi}((1:dnum)+doffset);
                    dat2 = prfs.([sbflds{2} '_bb']){iroi}((1:dnum)+doffset);

                    fitparms = fitlm(dat1,dat2,fitopt{:});
                    Rval = corr(dat1',dat2');
                    prfs_corrx.([sbfld{:} '_bb' '_fit']){iroi}  = [prfs_corrx.([sbfld{:} '_bb' '_fit']){iroi} fitparms.Coefficients.Estimate];
                    prfs_corrx.([sbfld{:} '_bb' '_R2']){iroi}   = [prfs_corrx.([sbfld{:} '_bb' '_R2']){iroi} fitparms.Rsquared.Ordinary];
                    prfs_corrx.([sbfld{:} '_bb' '_corr']){iroi} = [prfs_corrx.([sbfld{:} '_bb' '_corr']){iroi} Rval];
                    
                    dat1 = prfs.([sbflds{1} '_a']){iroi}((1:dnum)+doffset);
                    dat2 = prfs.([sbflds{2} '_a']){iroi}((1:dnum)+doffset);

                    fitparms = fitlm(dat1,dat2,fitopt{:});
                    Rval = corr(dat1',dat2');
                    prfs_corrx.([sbfld{:} '_a' '_fit']){iroi}  = [prfs_corrx.([sbfld{:} '_a' '_fit']){iroi} fitparms.Coefficients.Estimate];
                    prfs_corrx.([sbfld{:} '_a' '_R2']){iroi}   = [prfs_corrx.([sbfld{:} '_a' '_R2']){iroi} fitparms.Rsquared.Ordinary];
                    prfs_corrx.([sbfld{:} '_a' '_corr']){iroi} = [prfs_corrx.([sbfld{:} '_a' '_corr']){iroi} Rval];
                end
            end
        end
        fprintf(1,'.');
    end

    warning(warnst);
    catch ME
    warning(warnst);
    rethrow(ME);
    end
    fprintf(1,'\n');
end

%% Load low-broadband & Set xval params

bblList = {'V1','V2','V3','V1-V3','low'};
[modeldata_bbl, prf_params_bbl] = ecog_prf_loadprfs(subjectList,'bbL',prfPth,modeldataID,prfID,average,smoothingMode,smoothingN,prfmodel,gaussianmode,selectchs,selectch_exFEF,selectch_thresh,usefulltsR2,usefulltsxR2,skipsummarizeROIs);
[model_all_bbl, prf_all_bbl]    = ecog_prf_mergeprfs(modeldata_bbl,prf_params_bbl,va_area,usexvalparams,rearropt);

if ~isfield(prfs,'xval_bbl')
    for iroi=1:size(prfs.chanidx,1)
        prfs.R2_bbl{iroi,1}     = prf_all_bbl.R2([prfs.chanidx{iroi,:}])';
        prfs.ecc_bbl{iroi,1}    = prf_all_bbl.ecc([prfs.chanidx{iroi,:}])';
        prfs.ang_bbl{iroi,1}    = prf_all_bbl.ang([prfs.chanidx{iroi,:}])';
        prfs.rfsize_bbl{iroi,1} = prf_all_bbl.rfsize([prfs.chanidx{iroi,:}])';
        prfs.xval_bbl{iroi,1}   = prf_all_bbl.xval([prfs.chanidx{iroi,:}])';
        if ~ismember(rois{iroi},bblList)     % exclude pRF results from high visual areas, where low-broadband does not exist
            prfs.R2_bbl{iroi,1}(:)     = nan;
            prfs.ecc_bbl{iroi,1}(:)    = nan;
            prfs.ang_bbl{iroi,1}(:)    = nan;
            prfs.rfsize_bbl{iroi,1}(:) = nan;
            prfs.xval_bbl{iroi,1}(:)   = nan;
        end
    end
end

%-- Compute cross-correlation between xval and size (<1.5min)
if issaveplot,     mixedbb4err = false;
else,              mixedbb4err = true;
end
if ~isfield(prfs_corrx,'xval_rfsize_bbl_corr')
    nboot = size(prfs.chanidx,2);
    fitopt = {'RobustOpts','on'};

    flds  = {'xval_rfsize'};
 
     for sbfld = flds
       prfs_corrx.([sbfld{:} '_bbl' '_fit'])    = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_bbl' '_corr'])   = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_bbl' '_R2'])     = cell(nroi,1);
     end

    fprintf(1,'Computing');
    warnst = warning('off');
    try
    for iroi = roioi
        if ismember(rois{iroi},bblList)
          for iboot = 1:nboot
            dnum = prfs.num(iroi,iboot);
            doffset = sum(prfs.num(iroi,1:(iboot-1)));
            if dnum > 3
                for sbfld = flds
                    sbflds = strsplit(sbfld{:},'_');
                    dat1 = prfs.([sbflds{1} '_bbl']){iroi}((1:dnum)+doffset);
                    dat2 = prfs.([sbflds{2} '_bbl']){iroi}((1:dnum)+doffset);
                    if mixedbb4err
                        rplel = find(rand(1,dnum)>0.5);
                        dat1(rplel) = prfs.([sbflds{1} '_bb']){iroi}(rplel+doffset);
                        dat2(rplel) = prfs.([sbflds{2} '_bb']){iroi}(rplel+doffset);
                    end

                    fitparms = fitlm(dat1,dat2,fitopt{:});
                    Rval = corr(dat1',dat2');
                        
                    prfs_corrx.([sbfld{:} '_bbl' '_fit']){iroi}  = [prfs_corrx.([sbfld{:} '_bbl' '_fit']){iroi} fitparms.Coefficients.Estimate];
                    prfs_corrx.([sbfld{:} '_bbl' '_R2']){iroi}   = [prfs_corrx.([sbfld{:} '_bbl' '_R2']){iroi} fitparms.Rsquared.Ordinary];
                    prfs_corrx.([sbfld{:} '_bbl' '_corr']){iroi} = [prfs_corrx.([sbfld{:} '_bbl' '_corr']){iroi} Rval];
                end
            end
          end
          if mixedbb4err      
              prfs_corrx.([sbfld{:} '_bb' '_fit']){iroi}   = [prfs_corrx.([sbfld{:} '_bbl' '_fit']){iroi}];
              prfs_corrx.([sbfld{:} '_bb' '_R2']){iroi}    = [prfs_corrx.([sbfld{:} '_bbl' '_R2']){iroi}];
              prfs_corrx.([sbfld{:} '_bb' '_corr']){iroi}  = [prfs_corrx.([sbfld{:} '_bbl' '_corr']){iroi}];
              prfs_corrx.([sbfld{:} '_bbl' '_fit']){iroi}  = [nan;nan];
              prfs_corrx.([sbfld{:} '_bbl' '_R2']){iroi}   = [nan];
              prfs_corrx.([sbfld{:} '_bbl' '_corr']){iroi} = [nan];
          end
          fprintf(1,'.');
        else   
            prfs_corrx.([sbfld{:} '_bbl' '_fit']){iroi}  = [nan;nan];
            prfs_corrx.([sbfld{:} '_bbl' '_R2']){iroi}   = [nan];
            prfs_corrx.([sbfld{:} '_bbl' '_corr']){iroi} = [nan];
        end
    end

    warning(warnst);
    catch ME
    warning(warnst);
    rethrow(ME);
    end
    fprintf(1,'\n');
end

%% Comstruct Maximum-Probability prfs
prfs_mpm = prfs;
prfs_mpm.chanidx = {};
for iroi = 1:nroi
    %-- concatenate channel index
    ichx = [prfs.chanidx{iroi,:}];
    prfs_mpm.chanidx{iroi,1} = ichx;
end

%-- set wangprob
wangprob = channels(:,startsWith(channels.Properties.VariableNames,'wangprob')&endsWith(channels.Properties.VariableNames,rois));
%-- low
catrois = {'V1';'V2';'V3'};
catidx  = endsWith(wangprob.Properties.VariableNames,catrois);
wangprob.wangprob_low = sum(wangprob{:,catidx},2);
%-- high
catrois = {'V3a';'V3b';'LO1';'LO2';'TO';'IPS'};
catidx  = endsWith(wangprob.Properties.VariableNames,catrois);
wangprob.wangprob_high = sum(wangprob{:,catidx},2);

%-- select maximum-probability roi
flds = string(fieldnames(prfs_mpm))';
plroiset = {[1,2,3],[6:(nroi-2)],[nroi-1 nroi]};        % V1–V3;dorsolateral;low-high
for plroi = plroiset
    plroi = plroi{:};
for iroi = plroi
    [ichx,chidx] = unique(prfs_mpm.chanidx{iroi});
    chidx(wangprob{ichx,iroi} < max(wangprob{ichx,plroi},[],2)) = [];
    for ifld = flds
    if iscell(prfs_mpm.(ifld)) && size(prfs_mpm.(ifld),1)==nroi && ~all(isnan(prfs_mpm.(ifld){iroi,1}))
    prfs_mpm.(ifld){iroi,1} = prfs_mpm.(ifld){iroi,1}(chidx);
    end
    end
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize
%% %%%%%%%%%%%%%%%%%%%%%%%%

alpha = 0.32;   % 0.32 for 1sd, 0.05 for 2sd, 0.003 for 3sd
% plcol = get(groot,'defaultAxesColorOrder');
hF = gobjects(0);

    plroiset = {[nroi-1 nroi]};        % low,high
    plotmethod = 'scatter';
    %%% MPM plot
    prfs_pl  = prfs_mpm;
    MkrSiz     = 30;

showref = true;
domean  = true;

%% %%%%%%%%%%%%%%%%%%%%
%%% train-R2 vs Size

if issaveplot
relt    = 'R2_rfsize';
flds    = {'bb','a'};
opts = [];
opts.mothod        = plotmethod;
opts.overplot      = true;
opts.pix2deg       = cfactor;
opts.FontSize      = LFntSiz;
if showref
opts.showfitted    = 'fitlm';
else
opts.showfitted    = 'none';
end
opts.fit_alpha     = alpha;
opts.showdiag = false;
opts.fitlm_domean  = domean;
opts.fitlm_showref = showref;
opts.MarkerSize    = MkrSiz;

for plroi = plroiset
plroi = plroi{:};
    
selroi  = rois(plroi);
h = ecog_prf_plotPRFrelations(prfs_pl,prfs_corrx,relt,flds,selroi,rois,opts);
hF = [hF h];

if showref,     isrefname = '_withref';
else,           isrefname = '';
end
    
figureName = sprintf('prf-%02d%%-%02d%%-ecc%02d_line-roi_%s%s-largeFont',threshold_bb,threshold_a,eclimit,relt,isrefname);
if domean,  figureName = sprintf('%s-mean',figureName);   end
if plroi(1)~=1,   figureName = sprintf('%s-alt',figureName);   end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName));    end
end
end

%% %%%%%%%%%%%%%%%%%%%%
%%% xvalidated-R2 vs Size

if issaveplot
relt    = 'xval_rfsize';
flds    = {'bb','a'};
opts = [];
opts.mothod        = plotmethod;
opts.overplot      = true;
opts.pix2deg       = cfactor;
opts.FontSize      = LFntSiz;
if showref
opts.showfitted    = 'fitlm';
% opts.showfitted    = 'tls';
else
opts.showfitted    = 'none';
end
opts.fit_alpha     = alpha;
opts.showdiag = false;
opts.fitlm_domean  = domean;
opts.fitlm_showref = showref;
opts.MarkerSize    = MkrSiz;

for plroi = plroiset
plroi = plroi{:};
    
selroi  = rois(plroi);
h = ecog_prf_plotPRFrelations(prfs_pl,prfs_corrx,relt,flds,selroi,rois,opts);
hF = [hF h];

if showref,     isrefname = '_withref';
else,           isrefname = '';
end
    
figureName = sprintf('prf-%02d%%-%02d%%-ecc%02d_line-roi_%s%s-largeFont',threshold_bb,threshold_a,eclimit,relt,isrefname);
if domean,  figureName = sprintf('%s-mean',figureName);   end
if plroi(1)~=1,   figureName = sprintf('%s-alt',figureName);   end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName));    end
end
end

%% %%%%%%%%%%%%%%%%%%%%
%%% xvalidated-R2 vs Size with low-bb

if issaveplot    
    plroiset = {[nroi-1]};             % low
else
    plroiset = {[nroi-1 nroi]};        % low,high
end
    
relt    = 'xval_rfsize';
flds    = {'bb','a','bbl'};
fldNames = {'Broadband (70–180 Hz)','Alpha','Broadband (3–26 Hz)'};
opts = [];
opts.mothod        = plotmethod;
opts.overplot      = true;
opts.pix2deg       = cfactor;
opts.FontSize      = LFntSiz;
if showref
opts.showfitted    = 'fitlm';
% opts.showfitted    = 'tls';
else
opts.showfitted    = 'none';
end
opts.fit_alpha     = alpha;
opts.showdiag = false;
opts.fitlm_domean  = domean;
opts.fitlm_showref = showref;
opts.MarkerSize    = MkrSiz;
opts.condNames     = fldNames;
opts.Color         = plcol([1,2,5],:);

for plroi = plroiset
plroi = plroi{:};
    
selroi  = rois(plroi);
h = ecog_prf_plotPRFrelations(prfs_pl,prfs_corrx,relt,flds,selroi,rois,opts);
hF = [hF h];

if showref,     isrefname = '_withref';
else,           isrefname = '';
end
    
figureName = sprintf('prf-%02d%%-%02d%%-ecc%02d_line-roi_%s%s-bbl-largeFont',threshold_bb,threshold_a,eclimit,relt,isrefname);
if domean,  figureName = sprintf('%s-mean',figureName);   end
if plroi(1)~=1,   figureName = sprintf('%s-alt',figureName);   end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName));    end
end

