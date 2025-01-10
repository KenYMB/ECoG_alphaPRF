% ouput R-vales of pRF correlations computed in 06c (equivalant 03c2 + 03d2)
%   Folk from 10g

% 20240513 Yuasa - Output r-values for pRF sets

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

alpha = 0.32;   % 0.32 for 1sd, 0.05 for 2sd, 0.003 for 3sd

plcol = get(groot,'defaultAxesColorOrder');

hF = gobjects(0);

%% Construct Maximum-Probability prfs (each roi include channels which has max probabirity to be assigned to them)
%-- get channel index
prfs_mpm = prfs;
prfs_mpm.chanidx = {};
for iroi = 1:nroi
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
      ichx = mode(chidx,2)';
    end
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
plroiset = {[1,2,3],[nroi-1 nroi]};        % V1–V3;low-high
for plroi = plroiset
    plroi = plroi{:};
for iroi = plroi
    [ichx,chidx] = unique(prfs_mpm.chanidx{iroi});
    chidx(wangprob{ichx,iroi} < max(wangprob{ichx,plroi},[],2)) = [];
    for ifld = flds
    if iscell(prfs_mpm.(ifld)) && size(prfs_mpm.(ifld),1)==nroi
    prfs_mpm.(ifld){iroi,1} = prfs_mpm.(ifld){iroi,1}(chidx);
    end
    end
end
end



%% ang
%-- use all bootstrapped data
for iroi = 12:13  % V1-3, dorso
datX = prfs.ang_bb{iroi};
datY = prfs.ang_a{iroi};
    datX = datX - 360.*( (datX - datY) > 180 );
    datY = datY - 360.*( (datY - datX) > 180 );
fprintf(1,'%s   \tR = %.3f\n',rois{iroi},corr(datX',datY'));
end


%-- compute in iterations
boot_corr = [];
for iroi = 12:13  % V1-3, dorso
datX = prfs.ang_bb{iroi};
datY = prfs.ang_a{iroi};
    datX = datX - 360.*( (datX - datY) > 180 );
    datY = datY - 360.*( (datY - datX) > 180 );
prf_nums = [0,cumsum(prfs.num(iroi,:))]; % get index numbers of each iteration
for itr = 1:size(prfs.num,2)
    selitr = (prf_nums(itr)+1):prf_nums(itr+1);
    boot_corr{iroi}(itr) = corr(datX(selitr)',datY(selitr)'); 
end
fprintf(1,'%s   \tR = %.3f (p=%.4f)\n',rois{iroi},mean(boot_corr{iroi}),...
          2*mean((abs(boot_corr{iroi}) - mean(boot_corr{iroi}))>= abs(corr(datX',datY'))));
end

%% ecc
isrobustfit = false;
%-- use all bootstrapped data
for iroi = 12:13  % V1-3, dorso
datX = prfs.ecc_bb{iroi};
datY = prfs.ecc_a{iroi};
    if isrobustfit
        robustalpha = 0.0455;   % 2sigma
        datDist = sqrt((datX - median(datX)).^2 + (datY - median(datY)).^2);
        datDist(datDist<=0) = max(min(datDist).*1e-3,eps);
        pd   = fitdist(datDist','Lognormal');
        ngch = datDist > pd.icdf(1-robustalpha);
        datX = datX(~ngch);
        datY = datY(~ngch);
    end
fprintf(1,'%s   \tR = %.3f\n',rois{iroi},corr(datX',datY'));
end

%-- compute in iterations
boot_corr = [];
for iroi = 12:13  % V1-3, dorso
datX = prfs.ecc_bb{iroi};
datY = prfs.ecc_a{iroi};
    if isrobustfit
        robustalpha = 0.0455;   % 2sigma
        datDist = sqrt((datX - median(datX)).^2 + (datY - median(datY)).^2);
        datDist(datDist<=0) = max(min(datDist).*1e-3,eps);
        pd   = fitdist(datDist','Lognormal');
        ngch = datDist > pd.icdf(1-robustalpha);
        datX = datX(~ngch);
        datY = datY(~ngch);
    end
prf_nums = [0,cumsum(prfs.num(iroi,:))]; % get index numbers of each iteration
for itr = 1:size(prfs.num,2)
    selitr = (prf_nums(itr)+1):prf_nums(itr+1);
    boot_corr{iroi}(itr) = corr(datX(selitr)',datY(selitr)'); 
end
fprintf(1,'%s   \tR = %.3f (p=%.4f)\n',rois{iroi},mean(boot_corr{iroi}),...
          2*mean((abs(boot_corr{iroi}) - mean(boot_corr{iroi}))>= abs(corr(datX',datY'))));
end

%% size
isrobustfit = false;
% isrobustfit = true;
%-- use all bootstrapped data
for iroi = 12:13  % V1-3, dorso
datX = prfs.rfsize_bb{iroi};
datY = prfs.rfsize_a{iroi};
    if isrobustfit
        robustalpha = 0.0455;   % 2sigma
        datDist = sqrt((datX - median(datX)).^2 + (datY - median(datY)).^2);
        datDist(datDist<=0) = max(min(datDist).*1e-3,eps);
        pd   = fitdist(datDist','Lognormal');
        ngch = datDist > pd.icdf(1-robustalpha);
        datX = datX(~ngch);
        datY = datY(~ngch);
    end
fprintf(1,'%s   \tR = %.3f\n',rois{iroi},corr(datX',datY'));
end

%-- compute in iterations
boot_corr = [];
for iroi = 12:13  % V1-3, dorso
datX = prfs.rfsize_bb{iroi};
datY = prfs.rfsize_a{iroi};
    if isrobustfit
        robustalpha = 0.0455;   % 2sigma
        datDist = sqrt((datX - median(datX)).^2 + (datY - median(datY)).^2);
        datDist(datDist<=0) = max(min(datDist).*1e-3,eps);
        pd   = fitdist(datDist','Lognormal');
        ngch = datDist > pd.icdf(1-robustalpha);
        datX = datX(~ngch);
        datY = datY(~ngch);
    end
prf_nums = [0,cumsum(prfs.num(iroi,:))]; % get index numbers of each iteration
for itr = 1:size(prfs.num,2)
    selitr = (prf_nums(itr)+1):prf_nums(itr+1);
    boot_corr{iroi}(itr) = corr(datX(selitr)',datY(selitr)'); 
end
fprintf(1,'%s   \tR = %.3f (p=%.4f)\n',rois{iroi},mean(boot_corr{iroi}),...
          2*mean((abs(boot_corr{iroi}) - mean(boot_corr{iroi}))>= abs(corr(datX',datY'))));
end

%% ecc vs size
isrobustfit = false;
% isrobustfit = true;
%-- use all bootstrapped data
for iroi = 12:13  % V1-3, dorso
%-- bb
datX = prfs.ecc_bb{iroi};
datY = prfs.rfsize_bb{iroi};
    if isrobustfit
        robustalpha = 0.0027;   % 3sigma
        pd   = fitdist(datX','Gamma');
        ngch = datX > pd.icdf(1-robustalpha);
        pd   = fitdist(datY','Gamma');
        ngch = ngch | datY > pd.icdf(1-robustalpha);
        datX = datX(~ngch);
        datY = datY(~ngch);
    end
fprintf(1,'%s   \tR = %.3f\n',rois{iroi},corr(datX',datY'));
%-- a
datX = prfs.ecc_a{iroi};
datY = prfs.rfsize_a{iroi};
    if isrobustfit
        robustalpha = 0.0027;   % 3sigma
        pd   = fitdist(datX','Gamma');
        ngch = datX > pd.icdf(1-robustalpha);
        pd   = fitdist(datY','Gamma');
        ngch = ngch | datY > pd.icdf(1-robustalpha);
        datX = datX(~ngch);
        datY = datY(~ngch);
    end
fprintf(1,'%s   \tR = %.3f\n',rois{iroi},corr(datX',datY'));
end

