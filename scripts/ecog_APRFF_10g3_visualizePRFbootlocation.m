% visualize normalized pRF locations based on computation in 06c (equivalant 03c2 + 03d2)
%   Update xR2 with full time-series (no decimation)
%   Folk from ecog_APRF_06d_visualizePRFparamsFPM_bootall_fullts.m
%         and ecog_APRF_06e_visualizePRFsFPM_bootall_fullts.m

% 20210107 Yuasa - Update from 06d,06e
% 20210806 Yuasa - show bootstrapping alpha pRFs normalized for broadband pRFs
%                  Use different computation for normalize

%% prefix
% close all; clear all;
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
FntSiz  = 14;

%% Prepare parallel computation
isstartpar = false;
if exist('gcp','file') && isempty(gcp('nocreate')), parpool([1 40]); isstartpar = true;  end

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
% va_area = 'wangarea';

usefulltsR2   = false;
usefulltsxR2  = false;

    va_area = 'wangarea';
    modeldataID  = 'freq_spectra-timeseries-boot';
    prfID        = 'prf-boot';
    usefulltsxR2 = false;
%     ecog_APRFF_INITa_loaddata;
%     ecog_APRFF_INITb_mergedata;
% 
%     model_all_bb_boot = model_all_bb;
%     model_all_a_boot  = model_all_a;
%     prf_all_bb_boot   = prf_all_bb;
%     prf_all_a_boot    = prf_all_a;
    
    ecog_APRFF_INITc_postfix;
    
% modeldataID = [];
% prfID       = [];
% 
% ecog_APRFF_INITa_loaddata;
% ecog_APRFF_INITb_mergedata;

%% Load correlations
ecog_APRFF_INITd_threshold

% nboot     = 5000;
boottype  = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling
robustfit = true;

[prfs,prfs_corr,prfs_corrx,rois,~,nboot] = ...
    ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit);

%%% prepare visualize
% nsubjects = length(prf_params_bb);
nroi  = length(rois);
selroi = (nroi-1):nroi;
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
    
    hF(end+1) = figure('Position',[150 100 250*nplroi 255],'MenuBar','none');
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
               'FontSize', FntSiz);
       else
           set(gca,'XLim',[-1 1]*plotlng./dispratio,'YLim',[-1 1]*plotlng./dispratio,...
               'FontSize', FntSiz);
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
    
figureName = sprintf('prf-%02d%%-%02d%%-ecc%02d_location-polarnormboot-merge%s_size-%s',threshold_bb,threshold_a,eclimit,normpostfix,'All');
if nonlin,    figureName = sprintf('%s_geomean',figureName);    end
if normamp,   figureName = sprintf('%s_normamp',figureName);    end
if shiftcntr, figureName = sprintf('%s_shiftcntr',figureName);  end
if isrot,     figureName = sprintf('%s_rot',figureName);        end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName));    end

%%
%% plot bootstrapping norm pRFs <reverse>
%%
if issaveplot

nonlin  = false;
normamp = true;
    
    %% reversal norm
    switch normby
        case 'size',    normeccby   = prfs.rfsize_a;
                        normsizby   = prfs.rfsize_a;
        case 'ecc',     normeccby   = prfs.ecc_a;
                        normsizby   = prfs.ecc_a;
        case 'both',    normeccby   = prfs.ecc_a;
                        normsizby   = prfs.rfsize_a;
    end
    prfs_normR   = [];
    prfs_normR.rfsize = cellfun(@(x,y) x./y.*normres,prfs.rfsize_bb,normsizby,'UniformOutput',false);
    prfs_normR.ecc    = cellfun(@(x,y) x./y.*normres,prfs.ecc_bb,normeccby,'UniformOutput',false);
    prfs_normR.ang    = cellfun(@minus ,prfs.ang_bb,prfs.ang_a,'UniformOutput',false);
    [prfs_normR.R, prfs_normR.C] = cellfun(@(r,t) pol2cart(t./180.*pi,r),...
                                    prfs_normR.ecc,prfs_normR.ang,'UniformOutput',false);

        prfs_normR0   = [];
        prfs_normR0.rfsize = cellfun(@(x,y) x./y.*normres,prfs.rfsize_a,normsizby,'UniformOutput',false);
        prfs_normR0.ecc    = cellfun(@(x,y) x./y.*normres,prfs.ecc_a,normeccby,'UniformOutput',false);
        prfs_normR0.ang    = cellfun(@minus ,prfs.ang_a,prfs.ang_a,'UniformOutput',false);
        [prfs_normR0.R, prfs_normR0.C] = cellfun(@(r,t) pol2cart(t./180.*pi,r),...
                                        prfs_normR0.ecc,prfs_normR0.ang,'UniformOutput',false);

    %% constract bootstrapping norm pRFs (same height)    
    %-- construct bootstrap pRF images
    if nonlin,  meanfun = @(x,y) geomean(x,y,'omitnan');
    else,       meanfun = @(x,y) mean(x,y,'omitnan');
    end
    imPRFb = zeros([size(xi),nplroi,nboot]);
    imPRF0 = zeros([size(xi),nplroi,nboot]);
    for ii = 1:nplroi
       iroi = selroi(ii);
       fprintf('Bootstrapping normalized pRFs to alpha in %s...',rois{iroi});
       parfor iboot = 1:nboot
          bootidx = [1:prfs.num(iroi,iboot)]+sum(prfs.num(iroi,1:(iboot-1)));
          
          %-- normalize pRF
          PRFi = permute([prfs_normR.R{iroi}(bootidx); ...
                          prfs_normR.C{iroi}(bootidx); ...
                          prfs_normR.rfsize{iroi}(bootidx)],[3,1,2]);
          PRFexp = (bsxfun(@minus,xi,PRFi(1,1,:)).^2 + bsxfun(@minus,yi,PRFi(1,2,:)).^2)./(2.*PRFi(1,3,:).^2);
          PRFamp = 1;
          imPRFi = PRFamp  .* exp(-PRFexp);
          %%-- average
          imPRFb(:,:,ii,iboot) = meanfun(imPRFi,3);
          
          %-- reference pRF
          PRFi = permute([prfs_normR0.R{iroi}(bootidx); ...
                          prfs_normR0.C{iroi}(bootidx); ...
                          prfs_normR0.rfsize{iroi}(bootidx)],[3,1,2]);
          PRFexp = (bsxfun(@minus,xi,PRFi(1,1,:)).^2 + bsxfun(@minus,yi,PRFi(1,2,:)).^2)./(2.*PRFi(1,3,:).^2);
          PRFamp = 1;
          imPRFi = PRFamp  .* exp(-PRFexp);
          %%-- average
          imPRF0(:,:,ii,iboot) = meanfun(imPRFi,3);
       end
       fprintf('\n');
    end
    
    %% plot bootstrapping norm pRFs (same height)
    dispratio = 5;
    
    hF(end+1) = figure('Position',[150 100 250*nplroi 255],'MenuBar','none');
    tiledlayout(1,nplroi,'Padding','compact','TileSpacing','compact');
    ht = [];
    subcol = mean([plcol(1,:);ones(1,3).*0.85],1);
    refcol = mean([plcol(2,:);ones(1,3).*0.85],1);
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
       for iboot = pickupboot % nboot
           hcnt= contour(axplt.*plscale,axplt.*plscale,imPRFcur(:,:,ii,iboot),...
               exp(-plsigma.^2./2), 'LineColor', subcol,'LineWidth',0.7,'LineStyle',':');
       end
       for iboot = pickupboot % nboot
           hcnt= contour(axplt.*plscale,axplt.*plscale,imPRFref(:,:,ii,iboot),...
               exp(-plsigma.^2./2), 'LineColor', refcol,'LineWidth',0.7,'LineStyle',':');
       end
       axis xy equal
       if shiftcntr
           set(gca,'XLim',[-1 1]*plotlng./dispratio+refcntr(1,ii).*plscale,...
               'YLim',[-1 1]*plotlng./dispratio+refcntr(2,ii).*plscale,...
               'FontSize', FntSiz);
       else
           set(gca,'XLim',[-1 1]*plotlng./dispratio,'YLim',[-1 1]*plotlng./dispratio,...
               'FontSize', FntSiz);
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
                  exp(-plsigma.^2./2), 'LineColor', plcol(2,:),'LineWidth',1.5);
      hcnt= contour(ht(ii),axplt.*plscale,axplt.*plscale,imPRFave(:,:,ii),...
                  exp(-plsigma.^2./2), 'LineColor', plcol(1,:),'LineWidth',1.5);
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
    
figureName = sprintf('prf-%02d%%-%02d%%-ecc%02d_location-polarnormbootRev-merge%s_size-%s',threshold_bb,threshold_a,eclimit,normpostfix,'All');
if nonlin,    figureName = sprintf('%s_geomean',figureName);    end
if normamp,   figureName = sprintf('%s_normamp',figureName);    end
if shiftcntr, figureName = sprintf('%s_shiftcntr',figureName);  end
if isrot,     figureName = sprintf('%s_rot',figureName);        end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName));    end

end

%% Close parallel pool
if isstartpar, delete(gcp('nocreate'));  end
