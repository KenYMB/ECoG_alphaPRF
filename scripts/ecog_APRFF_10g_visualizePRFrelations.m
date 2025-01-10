% visualize pRF correlations computed in 06c (equivalant 03c2 + 03d2)
%   Update xR2 with full time-series (no decimation)
%   Folk from ecog_APRF_06d_visualizePRFparamsFPM_bootall_fullts.m
%         and ecog_APRF_06e_visualizePRFsFPM_bootall_fullts.m

% 20210107 Yuasa - Update from 06d,06e
% 20210809 Yuasa - modify for paper
% 20210916 Yuasa - Update for new environment

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

%% Comstruct Maximum-Probability prfs
prfs_mpm = prfs;
prfs_mpm.chanidx = {};
for iroi = 1:nroi
    %-- get channel index
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

%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations
% close all;

fitting = '-ellipse';
% fitting = '-robustellipse';

%% plot 2D histogram (replication around 90 deg for angle)
%% %%%%%%%%%%%%%%%%%%%%
%%% arrange in one figure

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

figureName = sprintf('prf-%02d%%-%02d%%-ecc%02d_2Dhist%s_%s',threshold_bb,threshold_a,eclimit,fitting,'all');
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName),'-vector');    end

%% %%%%%%%%%%%%%%%%%%%%
%%% cross-correlations

showref = true;
domean  = true;

%% plot fitted line in ROIs

if issaveplot
    plroiset = {[1 2 3 nroi],...       % V1,V2,V3,high
                [2 3 nroi-1 nroi]};    % V2,V3,low,high
    plotmethod = 'line';
    prfs_pl  = prfs;
    MkrSiz     = [];
else
    plroiset = {[nroi-1 nroi]};        % low,high
    plotmethod = 'scatter';
    %%% MPM plot
    prfs_pl  = prfs_mpm;
    MkrSiz     = 30;
%     %%% Probabilistic plot
%     prfs_pl  = prfs;
%     MkrSiz     = 10;
end

relt    = 'ecc_rfsize';
flds    = {'bb','a'};
opts = [];
opts.mothod        = plotmethod;
opts.overplot      = true;
opts.pix2deg       = cfactor;
opts.FontSize      = LFntSiz;
opts.showfitted    = 'fitlm';
opts.fit_alpha     = alpha;
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

%%
%% ColorBar for 2Dhist
if issaveplot
    
meshres = 300;
maxrad  = meshres/2;

hF(end+1) = figure;

[X,Y] = meshgrid(0:meshres,0:meshres/10);
Z = X;
surf(X,Y,Z,'EdgeColor','none');
axis equal;
colormap([[0,0,0];hot].^0.3);
set(gcf, 'InvertHardcopy', 'off', 'Color', [1 1 1], 'Position', [200 100 50 250]);
view([90 0.1]); box on;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);

figureName = 'colorbar';
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figureName),{'png','eps'},'-vector');  end

end

