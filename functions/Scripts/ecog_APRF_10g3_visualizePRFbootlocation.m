% visualize pRF correlations computed in 06c (equivalant 03c2 + 03d2)
%   Update xR2 with full time-series (no decimation)

% 20210107 Yuasa - Update from 06d,06e
% 20210806 Yuasa - show bootstrapping alpha pRFs normalized for broadband pRFs
%                  Use different computation for normalize

%% Prepare parallel computation
isstartpar = false;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool; isstartpar = true;  end

%% Define paths and dataset
% close all; clear all;
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
SetDefault('issaveplot',true);
if issaveplot
    plotsavePth    = fullfile(figPth, 'pRFrelations-representative');
end
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%% load analyzePRF & recompute full-ts R2
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
% va_area = 'wangarea';

usefulltsR2   = false;
usefulltsxR2  = false;

    va_area = 'wangarea';
    modeldataID  = 'freq_spectra-timeseries-boot';
    prfID        = 'prf-boot';
    usefulltsxR2 = false;
    ecog_APRF_INITa_loaddata;
    ecog_APRF_INITb_mergedata;

    model_all_bb_boot = model_all_bb;
    model_all_a_boot  = model_all_a;
    prf_all_bb_boot   = prf_all_bb;
    prf_all_a_boot    = prf_all_a;
    
    ecog_APRF_INITc_postfix;
    
modeldataID = [];
prfID       = [];

ecog_APRF_INITa_loaddata;
ecog_APRF_INITb_mergedata;

%% Load correlations
ecog_APRF_INITd_threshold

% nboot     = 5000;
boottype  = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling

filename = sprintf('all_prfparams-%sall%s-%s%s_thresh%02d',boottype,R2mode,selectchs,postfix,threshold_bb);
if exist('threshold_a','var')&&~isnan(threshold_a)
    filename = sprintf('%s-%02d',filename,threshold_a); end
if exist('eclimit','var')&&~isnan(eclimit)
    filename = sprintf('%s-ecc%02d',filename,eclimit);     end
filepath = fullfile(savePth,'pRFanalysis',filename);

if exist([filepath '.mat'],'file')
 fprintf('loading %s...',filename);
 load(filepath);
 fprintf('\n');
 nroi = length(rois);
 nboot = size(prfs.num,2);
 SetDefault('threshold_a',nan); SetDefault('eclimit',nan);
else
 error('Failed to find %s',filename);
end

%-- update ROI labels
rois(ismember(rois,'low'))={'V1-V3'};
rois(ismember(rois,'high'))={'Dorsolateral'};
rois(ismember(rois,'dorsolateral'))={'Dorsolateral'};

%-- reject inaccurate ROIs
rejthr = nboot .* 2;
for ifld = fieldnames(prfs)'
    for iroi = 1:nroi
        if ~ismember(ifld{:},{'num','chanidx'}) && length(prfs.(ifld{:}){iroi})<rejthr
        prfs.(ifld{:}){iroi} = nan;
        end
    end
end

%%% prepare visualize
nsubjects = length(prf_params_bb);
nroi  = length(rois);
selroi = (nroi-1):nroi;
% nvisarea  = length(prf_va_bb);
% areaList  = unique(prf_all_bb.channels.(va_area));

res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');

if issaveplot
    plotsavedir    = fullfile(plotsavePth, postfix(2:end),'boot');
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations
% close all;
fontSiz = 18;

%--- Normalize parameter
normres = 100;
normby   = 'ecc';   % 'size','ecc','each'
% % normcntr = [0,0];    % center of reference pRF
% normcntr = [1,0];    % center of reference pRF

%-- for figure name
normpostfix = sprintf('-BY%s',normby);
% if ~any(normcntr), normpostfix = sprintf('%s-%s',normpostfix,'shift'); end

%% normalize randomized pRF (random ROI assign & random sampling = full sample of bootstrap)
switch normby
    case 'size',    normeccby   = prfs.rfsize_bb;
                    normsizby   = prfs.rfsize_bb;
    case 'ecc',     normeccby   = prfs.ecc_bb;
                    normsizby   = prfs.ecc_bb;
    case 'both',    normeccby   = prfs.ecc_bb;
                    normsizby   = prfs.rfsize_bb;
end
prfs_norm   = [];
prfs_norm.rfsize = cellfun(@(x,y) x./y.*normres,prfs.rfsize_a,normsizby,'UniformOutput',false);
prfs_norm.ecc    = cellfun(@(x,y) x./y.*normres,prfs.ecc_a,normeccby,'UniformOutput',false);
prfs_norm.ang    = cellfun(@minus ,prfs.ang_a,prfs.ang_bb,'UniformOutput',false);
[prfs_norm.R, prfs_norm.C] = cellfun(@(r,t) pol2cart(t./180.*pi,r),...
                                prfs_norm.ecc,prfs_norm.ang,'UniformOutput',false);
                            
    prfs_norm0   = [];
    prfs_norm0.rfsize = cellfun(@(x,y) x./y.*normres,prfs.rfsize_bb,normsizby,'UniformOutput',false);
    prfs_norm0.ecc    = cellfun(@(x,y) x./y.*normres,prfs.ecc_bb,normeccby,'UniformOutput',false);
    prfs_norm0.ang    = cellfun(@minus ,prfs.ang_bb,prfs.ang_bb,'UniformOutput',false);
    [prfs_norm0.R, prfs_norm0.C] = cellfun(@(r,t) pol2cart(t./180.*pi,r),...
                                    prfs_norm0.ecc,prfs_norm0.ang,'UniformOutput',false);

%%% plot parameters
prfparams = [];
pickupboot = 1:100;

plscale   = 1./normres;
plres     = normres.*2;
plsigma   = [1 1];
plotlng   = 10;

%%% plot prep
axplt = linspace(-plotlng./plscale,plotlng./plscale,plres);
[xi, yi]=meshgrid(axplt,axplt);
          
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot PRF locations
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot bootstrapping PRF locations with another nomalization
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- average pRFs across electrodes in each iteration
plroi  = selroi;
nplroi = length(plroi);
    %%
    %% plot bootstrapping norm pRFs
    %%
    nonlin    = false;
    normamp   = true;
    shiftcntr = true;
    isrot     = true;
    
    %% constract bootstrapping norm pRFs (same height)    
    %-- construct bootstrap pRF images
    if nonlin,  meanfun = @(x,y) geomean(x,y,'omitnan');
    else,       meanfun = @(x,y) mean(x,y,'omitnan');
    end
    imPRFb  = zeros([size(xi),nplroi,nboot]);
    imPRF0  = zeros([size(xi),nplroi,nboot]);
    refcntr = zeros([2,nplroi,nboot]);
    for ii = 1:nplroi
       iroi = selroi(ii);
       fprintf('Bootstrapping normalized pRFs to broadband in %s...',rois{iroi});
       parfor iboot = 1:nboot
          bootidx = [1:prfs.num(iroi,iboot)]+sum(prfs.num(iroi,1:(iboot-1)));
          
          %-- normalize pRF
          PRFi = permute([prfs_norm.R{iroi}(bootidx); ...
                          prfs_norm.C{iroi}(bootidx); ...
                          prfs_norm.rfsize{iroi}(bootidx)],[3,1,2]);
          PRFexp = (bsxfun(@minus,xi,PRFi(1,1,:)).^2 + bsxfun(@minus,yi,PRFi(1,2,:)).^2)./(2.*PRFi(1,3,:).^2);
          PRFamp = 1;
          imPRFi = PRFamp  .* exp(-PRFexp);
          %%-- average
          imPRFb(:,:,ii,iboot) = meanfun(imPRFi,3);
          
          %-- reference pRF
          PRFi = permute([prfs_norm0.R{iroi}(bootidx); ...
                          prfs_norm0.C{iroi}(bootidx); ...
                          prfs_norm0.rfsize{iroi}(bootidx)],[3,1,2]);
          PRFexp = (bsxfun(@minus,xi,PRFi(1,1,:)).^2 + bsxfun(@minus,yi,PRFi(1,2,:)).^2)./(2.*PRFi(1,3,:).^2);
          PRFamp = 1;
          imPRFi = PRFamp  .* exp(-PRFexp);
          %%-- average
          imPRF0(:,:,ii,iboot) = meanfun(imPRFi,3);
          refcntr(:,ii,iboot)  = meanfun(PRFi(:,1:2,:),3);
       end
       fprintf('\n');
    end
    refcntr = meanfun(refcntr,3);
    
    %% plot bootstrapping norm pRFs (same height)
    dispratio = 2;
    
    hF = figure('Position',[150 100 250*nplroi 255],'MenuBar','none');
    tiledlayout(1,nplroi,'Padding','compact','TileSpacing','compact');
    ht = [];
    subcol = mean([plcol(2,:);ones(1,3).*0.85],1);
    refcol = mean([plcol(1,:);ones(1,3).*0.85],1);
    imPRFcur = imPRFb;
    imPRFref = imPRF0;
    if normamp
        imPRFcur = imPRFcur ./ max(imPRFcur,[],[1,2],'omitnan');
        imPRFref = imPRFref ./ max(imPRFref,[],[1,2],'omitnan');
    end
    
    %-- plot some sample of bootstrap pRFs
    for ii = 1:nplroi
      iroi = selroi(ii);
      ht(ii) = nexttile;
      plot([-1 1]*plotlng,[0 0],'k:',[0 0],[-1 1]*plotlng,'k:');
      hold on;
       for iboot = pickupboot % 1:nboot
           hcnt= contour(axplt.*plscale,axplt.*plscale,imPRFcur(:,:,ii,iboot),...
               exp(-plsigma.^2./2), 'LineColor', subcol,'LineWidth',0.7,'LineStyle',':');
       end
       for iboot = pickupboot % 1:nboot
           hcnt= contour(axplt.*plscale,axplt.*plscale,imPRFref(:,:,ii,iboot),...
               exp(-plsigma.^2./2), 'LineColor', refcol,'LineWidth',0.7,'LineStyle',':');
       end
       axis xy equal
       if shiftcntr
           set(gca,'XLim',[-1 1]*plotlng./dispratio+refcntr(1,ii).*plscale,...
               'YLim',[-1 1]*plotlng./dispratio+refcntr(2,ii).*plscale,...
               'FontSize', 14);
       else
           set(gca,'XLim',[-1 1]*plotlng./dispratio,'YLim',[-1 1]*plotlng./dispratio,...
               'FontSize', 14);
       end
       title(rois{iroi});
       
       %%-- delete axis
          hA = gca;
          hA.XTick = [];   hA.YTick = [];   hA.ZTick = [];
          hA.XAxis.Visible = 'off';
          hA.YAxis.Visible = 'off';
          hA.ZAxis.Visible = 'off';
    end
    
    %-- plot averaged bootstrap pRFs
    imPRFave = meanfun(imPRFcur,4);
    imPRFave0 = meanfun(imPRFref,4);
    if normamp
        imPRFave  = imPRFave ./ max(imPRFave,[],[1,2],'omitnan');
        imPRFave0 = imPRFave0 ./ max(imPRFave0,[],[1,2],'omitnan');
    end
    for ii = 1:nplroi
      hcnt= contour(ht(ii),axplt.*plscale,axplt.*plscale,imPRFave0(:,:,ii),...
                  exp(-plsigma.^2./2), 'LineColor', plcol(1,:),'LineWidth',1.5);
      hcnt= contour(ht(ii),axplt.*plscale,axplt.*plscale,imPRFave(:,:,ii),...
                  exp(-plsigma.^2./2), 'LineColor', plcol(2,:),'LineWidth',1.5);
    end
    
    %-- plot center
    xscale = .2 ./ dispratio;
    for ii = 1:nplroi
        plot(ht(ii),refcntr(1,ii).*plscale + [-1 1]*xscale,refcntr(2,ii).*plscale + [-1 1]*xscale,'k-','LineWidth',1);
        plot(ht(ii),refcntr(1,ii).*plscale + [-1 1]*xscale,refcntr(2,ii).*plscale + [1 -1]*xscale,'k-','LineWidth',1);
    end
    
    %-- axis rotation
    if isrot
        for ii = 1:nplroi
            view(ht(ii),[-90 90]);
        end
    end
    
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_location-polarnormboot-merge%s_size-%s',threshold_bb,threshold_a,eclimit,normpostfix,'All');
if nonlin,    figname = sprintf('%s_geomean',figname);    end
if normamp,   figname = sprintf('%s_normamp',figname);    end
if shiftcntr, figname = sprintf('%s_shiftcntr',figname);  end
if isrot,     figname = sprintf('%s_rot',figname);        end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));  end

%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
