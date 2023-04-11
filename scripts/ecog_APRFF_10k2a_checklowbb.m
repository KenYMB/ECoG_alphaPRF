% Compare Gaussian Models, Channel Selection
%   Folk from ecog_APRF_07b_checkmodel_fullts
 
% 20210510 Yuasa
 
%% prefix
% close all; clearvars;
% startupToolboxToolbox;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'pRFselections-representative');
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
else
plotsavePth    = 'pRFselections-representative';
prfPth         = 'pRFmodel';
prfstatPth     = 'pRFanalysis';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth), 'modelselection', 'lowbb');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');
 
%-- Plotting Setting
FntSiz    = 20;

%% Load dataset
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
    noACorrect     = false;
va_area = 'wangarea';

usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

%%% load low broadband
[modeldata_bbl, prf_params_bbl] = ecog_prf_loadprfs(subjectList,'bbL',prfPth,modeldataID,prfID,average,smoothingMode,smoothingN,prfmodel,gaussianmode,selectchs,selectch_exFEF,selectch_thresh,usefulltsR2,usefulltsxR2,skipsummarizeROIs);
[model_all_bbl, prf_all_bbl]    = ecog_prf_mergeprfs(modeldata_bbl,prf_params_bbl,va_area,usexvalparams,rearropt);

%% load analyzePRF @DOG,OG,LFS

modeldata  = struct();
prf_params = struct();
model_all  = struct();
prf_all    = struct();
threshold  = struct();

%-- merge loaded data
ii = 1;
modeldata(ii).a   = modeldata_a;
modeldata(ii).bb  = modeldata_bb;
prf_params(ii).a  = prf_params_a;
prf_params(ii).bb = prf_params_bb;
 
model_all(ii).a   = model_all_a;
model_all(ii).bb  = model_all_bb;
prf_all(ii).a     = prf_all_a;
prf_all(ii).bb    = prf_all_bb;

threshold(ii).a   = threshold_a;
threshold(ii).bb  = threshold_bb;

ii = 2;
modeldata(ii).a   = modeldata_bbl;
modeldata(ii).bb  = modeldata_bbl;
prf_params(ii).a  = prf_params_bbl;
prf_params(ii).bb = prf_params_bbl;
 
model_all(ii).a   = model_all_bbl;
model_all(ii).bb  = model_all_bbl;
prf_all(ii).a     = prf_all_bbl;
prf_all(ii).bb    = prf_all_bbl;

threshold(ii).a   = threshold_a;
threshold(ii).bb  = threshold_bb;

%% %%%%%%%%%%%%%%
%% Visualize
%% %%%%%%%%%%%%%%
%%% prepare visualize
res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;
 
alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');
hF = gobjects(0);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FaCLb(BW) vs FaLb
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelsets = [1 2];      % FaCLb(BW) vs FaLb
targetBAND = ('lowBroadband');

% eclimit = 50;
idat = 1;
okch1bb = ~(prf_all(modelsets(idat)).bb.xval<=threshold(modelsets(idat)).bb | prf_all(modelsets(idat)).bb.ecc >= eclimit);  % use tilde for nan
okch1a  = ~(prf_all(modelsets(idat)).a.xval<=threshold(modelsets(idat)).a | prf_all(modelsets(idat)).a.ecc >= eclimit);  % use tilde for nan
idat = 2;
okch2bb = ~(prf_all(modelsets(idat)).bb.xval<=threshold(modelsets(idat)).bb | prf_all(modelsets(idat)).bb.ecc >= eclimit);  % use tilde for nan
okch2a  = ~(prf_all(modelsets(idat)).a.xval<=threshold(modelsets(idat)).a | prf_all(modelsets(idat)).a.ecc >= eclimit);  % use tilde for nan

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Test Performance
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
%% Difference of Gaussian Models (scatter)
roundunit = 20;
hF(end+1) = figure('Menubar','none','Position',[200 200 850 420]);
subplot_er(1,2,1);
    scatter(prf_all(modelsets(2)).bb.xval,prf_all(modelsets(1)).bb.xval);
    ll = max(-100,floor(nanmin([prf_all(modelsets(2)).bb.xval(:);prf_all(modelsets(1)).bb.xval(:)])/roundunit)*roundunit);
    ul = min(100,ceil(nanmax([prf_all(modelsets(2)).bb.xval(:);prf_all(modelsets(1)).bb.xval(:)])/roundunit)*roundunit);
    axis([ll ul ll ul],'square');
    hold on;
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    set(gca,'FontSize',FntSiz);
    xlabel('Broadband (3–26 Hz)'); ylabel('Broadband (70–180 Hz)');
subplot_er(1,2,2);
    scatter(prf_all(modelsets(2)).a.xval,prf_all(modelsets(1)).a.xval);
    ll = max(-100,floor(nanmin([prf_all(modelsets(2)).a.xval(:);prf_all(modelsets(1)).a.xval(:)])/roundunit)*roundunit);
    ul = min(100,ceil(nanmax([prf_all(modelsets(2)).a.xval(:);prf_all(modelsets(1)).a.xval(:)])/roundunit)*roundunit);
    axis([ll ul ll ul],'square');
    hold on;
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    set(gca,'FontSize',FntSiz);
    xlabel('Broadband (3–26 Hz)'); ylabel('Alpha');
 
figname =sprintf('xR2-%s%s_%s_all',targetBAND,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Test Performance w/ threshold in each ROI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;

%% Difference of Gaussian Model in training accuracy (scatter) w/ rough ROI
plcol_wang = [ones(3,1)*plcol(1,:);...
              ones(6,1)*plcol(4,:);];
plshp_wang = [repmat('o',1,3) repmat('^',1,6)];
          
roundunit = 20;
% hF(end+1) = figure('Menubar','none','Position',[200 200 860 420],'defaultAxesColorOrder',plcol_wang);
hF(end+1) = figure('Menubar','none','Position',[200 200 900 440],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
roilist = {'V1','V2','V3','V3a','V3b','LO1','LO2','TO','IPS'}';
nexttile;
    ll = 100; ul = -100;
    ii=1;
    for iroi = roilist'
    okch = ((okch1bb)|(okch2bb)) & ismember(prf_all(modelsets(2)).bb.channels.wangarea,iroi);
    scatter(prf_all(modelsets(2)).bb.xval(okch),prf_all(modelsets(1)).bb.xval(okch),plshp_wang(ii),'LineWidth',1.6);
    if sum(okch)
    ll = min(ll,floor(nanmin([prf_all(modelsets(2)).bb.xval(okch);prf_all(modelsets(1)).bb.xval(okch)])/roundunit)*roundunit);
    ul = max(ul,ceil(nanmax([prf_all(modelsets(2)).bb.xval(okch);prf_all(modelsets(1)).bb.xval(okch)])/roundunit)*roundunit);
    end
    hold on;
    ii=ii+1;
    end
    hs = flipud(findobj(get(gca,'Children'),'Type','Scatter'));
    ll = max(-100,ll);
    ul = min(100,ul);
    axis([ll ul ll ul],'square');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    plot([getelement(xlim,1) threshold(modelsets(2)).bb],[1 1].*threshold(modelsets(1)).bb,'k-.',...
         [1 1].*threshold(modelsets(2)).bb,[getelement(ylim,1) threshold(modelsets(1)).bb],'k-.');  % threshold
    hl=legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southwest');
%     hl.Layout.Tile = 'east';
    
    set(gca,'FontSize',FntSiz);
    xlabel('Broadband (3–26 Hz)'); ylabel('Broadband (70–180 Hz)');
    xticks(yticks); xtickangle(0);
    
nexttile;
    ll = 100; ul = -100;
    ii=1;
    for iroi = roilist'
    okch = ((okch1a)|(okch2a)) & ismember(prf_all(modelsets(2)).a.channels.wangarea,iroi);
    scatter(prf_all(modelsets(2)).a.xval(okch),prf_all(modelsets(1)).a.xval(okch),plshp_wang(ii),'LineWidth',1.6);
    if sum(okch)
    ll = min(ll,floor(nanmin([prf_all(modelsets(2)).a.xval(okch);prf_all(modelsets(1)).a.xval(okch)])/roundunit)*roundunit);
    ul = max(ul,ceil(nanmax([prf_all(modelsets(2)).a.xval(okch);prf_all(modelsets(1)).a.xval(okch)])/roundunit)*roundunit);
    end
    hold on;
    ii=ii+1;
    end
    hs = flipud(findobj(get(gca,'Children'),'Type','Scatter'));
    ll = max(-100,ll);
    ul = min(100,ul);
    axis([ll ul ll ul],'square');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    plot([getelement(xlim,1) threshold(modelsets(2)).a],[1 1].*threshold(modelsets(1)).a,'k-.',...
         [1 1].*threshold(modelsets(2)).a,[getelement(ylim,1) threshold(modelsets(1)).a],'k-.');  % threshold
%     hl=legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
% %     hl.Layout.Tile = 'east';
    
    set(gca,'FontSize',FntSiz);
    xlabel('Broadband (3–26 Hz)'); ylabel('Alpha');
    xticks(yticks); xtickangle(0);
 
figname =sprintf('xR2-%s%s_%s_ROI',targetBAND,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Test Performance w/ threshold in V1–V3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;

ispltall = true;
%% Difference of Gaussian Model in training accuracy (scatter) w/ rough ROI
plcol_wang = [ones(3,1)*plcol(1,:);...
              ones(6,1)*plcol(4,:);];
plshp_wang = [repmat('o',1,3) repmat('^',1,6)];
if ispltall
    plcol_gray = [plcol_wang .* 0  +  ones(1,3).*0.7 .* 2]./2;
end
          
roundunit = 20;
% hF(end+1) = figure('Menubar','none','Position',[200 200 860 420],'defaultAxesColorOrder',plcol_wang);
hF(end+1) = figure('Menubar','none','Position',[200 200 900 440],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
roilist = {'V1','V2','V3'}';
nexttile;
    ll = 100; ul = -100;
    ii=1;
    for iroi = roilist'
    roich = ismember(prf_all(modelsets(2)).bb.channels.wangarea,iroi);
    okch = ((okch1bb)|(okch2bb)) & roich;
    scatter(prf_all(modelsets(2)).bb.xval(okch),prf_all(modelsets(1)).bb.xval(okch),plshp_wang(ii),'LineWidth',1.6);
    if sum(okch)
    ll = min(ll,floor(nanmin([prf_all(modelsets(2)).bb.xval(okch);prf_all(modelsets(1)).bb.xval(okch)])/roundunit)*roundunit);
    ul = max(ul,ceil(nanmax([prf_all(modelsets(2)).bb.xval(okch);prf_all(modelsets(1)).bb.xval(okch)])/roundunit)*roundunit);
    end
    hold on;
    ii=ii+1;
    end
    hs = flipud(findobj(get(gca,'Children'),'Type','Scatter'));
    if ispltall
        ii=1;
        for iroi = roilist'
        roich = ismember(prf_all(modelsets(2)).bb.channels.wangarea,iroi);
        okch = ~((okch1bb)|(okch2bb)) & roich;
        scatter(prf_all(modelsets(2)).bb.xval(okch),prf_all(modelsets(1)).bb.xval(okch),plshp_wang(ii),'LineWidth',1.6,'MarkerEdgeColor',plcol_gray(ii,:));
        if sum(okch)
        ll = min(ll,floor(nanmin([prf_all(modelsets(2)).bb.xval(okch);prf_all(modelsets(1)).bb.xval(okch)])/roundunit)*roundunit);
        ul = max(ul,ceil(nanmax([prf_all(modelsets(2)).bb.xval(okch);prf_all(modelsets(1)).bb.xval(okch)])/roundunit)*roundunit);
        end
        hold on;
        ii=ii+1;
        end
    end
    ll = max(-100,ll);
    ul = min(100,ul);
    axis([ll ul ll ul],'square');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    plot([getelement(xlim,1) threshold(modelsets(2)).bb],[1 1].*threshold(modelsets(1)).bb,'k-.',...
         [1 1].*threshold(modelsets(2)).bb,[getelement(ylim,1) threshold(modelsets(1)).bb],'k-.');  % threshold
    hl=legend(hs([1]),["V1-V3"],'Location','southeast');
%     hl.Layout.Tile = 'east';
    
    set(gca,'FontSize',FntSiz);
    xlabel('Broadband (3–26 Hz)'); ylabel('Broadband (70–180 Hz)');
    xticks(yticks); xtickangle(0);
    
nexttile;
    ll = 100; ul = -100;
    ii=1;
    for iroi = roilist'
    roich = ismember(prf_all(modelsets(2)).a.channels.wangarea,iroi);
    okch = ((okch1a)|(okch2a)) & roich;
    scatter(prf_all(modelsets(2)).a.xval(okch),prf_all(modelsets(1)).a.xval(okch),plshp_wang(ii),'LineWidth',1.6);
    if sum(okch)
    ll = min(ll,floor(nanmin([prf_all(modelsets(2)).a.xval(okch);prf_all(modelsets(1)).a.xval(okch)])/roundunit)*roundunit);
    ul = max(ul,ceil(nanmax([prf_all(modelsets(2)).a.xval(okch);prf_all(modelsets(1)).a.xval(okch)])/roundunit)*roundunit);
    end
    hold on;
    ii=ii+1;
    end
    hs = flipud(findobj(get(gca,'Children'),'Type','Scatter'));
    if ispltall
        ii=1;
        for iroi = roilist'
        roich = ismember(prf_all(modelsets(2)).a.channels.wangarea,iroi);
        okch = ~((okch1a)|(okch2a)) & roich;
        scatter(prf_all(modelsets(2)).a.xval(okch),prf_all(modelsets(1)).a.xval(okch),plshp_wang(ii),'LineWidth',1.6,'MarkerEdgeColor',plcol_gray(ii,:));
        if sum(okch)
        ll = min(ll,floor(nanmin([prf_all(modelsets(2)).a.xval(okch);prf_all(modelsets(1)).a.xval(okch)])/roundunit)*roundunit);
        ul = max(ul,ceil(nanmax([prf_all(modelsets(2)).a.xval(okch);prf_all(modelsets(1)).a.xval(okch)])/roundunit)*roundunit);
        end
        hold on;
        ii=ii+1;
        end
    end
    ll = max(-100,ll);
    ul = min(100,ul);
    axis([ll ul ll ul],'square');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    plot([getelement(xlim,1) threshold(modelsets(2)).a],[1 1].*threshold(modelsets(1)).a,'k-.',...
         [1 1].*threshold(modelsets(2)).a,[getelement(ylim,1) threshold(modelsets(1)).a],'k-.');  % threshold
%     hl=legend(hs([1]),["V1-V3"],'Location','southeast');
% %     hl.Layout.Tile = 'east';
    
    set(gca,'FontSize',FntSiz);
    xlabel('Broadband (3–26 Hz)'); ylabel('Alpha');
    xticks(yticks); xtickangle(0);
 
figname =sprintf('xR2-%s%s_%s_ROI-low',targetBAND,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

