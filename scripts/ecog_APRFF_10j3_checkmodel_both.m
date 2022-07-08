% Compare pRF models & Gaussian models, Channel Selection
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
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth), 'modelselection', 'pRFmodel');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');
 
%-- Plotting Setting
FntSiz    = 20;

%% general parameter
clear alphaType broadbandType

average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
% prfmodel       ='linear';
% gaussianmode   ='gs';
selectchs      = 'wangprobchs';
    allowlag       = false;
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;
va_area = 'wangarea';

usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

%% load analyzePRF @DOG,OG,LFS
gaussianModes = {'gs','og'};
prfmodelModels = {'linear','css'};
modelNames = ["Difference of Gaussians","CSS"];
%   gaussianModes = {'gs','dog'};
%   prfmodelModels = {'linear','css'};
%   modelNames = ["Linear - DoG","CSS - DoG(FULL)"];

modeldata  = struct();
prf_params = struct();
model_all  = struct();
prf_all    = struct();
threshold  = struct();
for ii=length(prfmodelModels):-1:1
gaussianmode = gaussianModes{ii};
prfmodel     = prfmodelModels{ii};
 
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold
 
%-- merge loaded data
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
end

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
%% Model Comparison
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelsets = [1 2];      % linear-gs vs css-dog

modelVS = [strjoin([prfmodelModels(modelsets);gaussianModes(modelsets)],{'-','VS','-'})];
ecog_APRFF_INITc_postfix;

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
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
    title('Broadband');
subplot_er(1,2,2);
    scatter(prf_all(modelsets(2)).a.xval,prf_all(modelsets(1)).a.xval);
    ll = max(-100,floor(nanmin([prf_all(modelsets(2)).a.xval(:);prf_all(modelsets(1)).a.xval(:)])/roundunit)*roundunit);
    ul = min(100,ceil(nanmax([prf_all(modelsets(2)).a.xval(:);prf_all(modelsets(1)).a.xval(:)])/roundunit)*roundunit);
    axis([ll ul ll ul],'square');
    hold on;
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    set(gca,'FontSize',FntSiz);
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
    title('Alpha');
 
figname =sprintf('xR2-%s%s_%s_all',modelVS,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Test Performance w/ threshold
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
%% Histogram of difference
% roundunit = 5;

% okch = (okch1bb&okch1a)|(okch2bb&okch2a);
%-- broadband
hF(end+1) = figure('Menubar','none','Position',[200 200 800 420]);
ht = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile;  okch = (okch1bb)|(okch2bb);
xval_diff = prf_all(modelsets(1)).bb.xval(okch) - prf_all(modelsets(2)).bb.xval(okch);
h=histogram(xval_diff);
    h.NumBins = 20;
    roundunit = h.BinWidth;
    llul = xlim;
    ll   = floor(min(llul)/roundunit)*roundunit;
    ul   = ceil(max(llul)/roundunit)*roundunit;
    xlim([ll ul]);
    h.BinWidth = roundunit;
set(gca,'FontSize',FntSiz);
hold on;
    hz=plot([0 0],ylim,'k-','LineWidth',1.2);   % axis
    hv=plot([1 1].*mean(xval_diff),ylim,'b--','LineWidth',1.6);  % threshold
title('Broadband');
%-- alpha
nexttile;  okch = (okch1a)|(okch2a);
xval_diff = prf_all(modelsets(1)).a.xval(okch) - prf_all(modelsets(2)).a.xval(okch);
h(2)=histogram(xval_diff);
    h(2).NumBins = 20;
    roundunit = h(2).BinWidth;
    llul = xlim;
    ll   = floor(min(llul)/roundunit)*roundunit;
    ul   = ceil(max(llul)/roundunit)*roundunit;
    xlim([ll ul]);
    h(2).BinWidth = roundunit;
set(gca,'FontSize',FntSiz);
hold on;
    hz(2)=plot([0 0],ylim,'k-','LineWidth',1.2);   % axis
    hv(2)=plot([1 1].*mean(xval_diff),ylim,'b--','LineWidth',1.6);  % threshold
title('Alpha');
tt=title(ht,sprintf('%s – %s',modelNames(modelsets(1)),modelNames(modelsets(2))),...
    'HorizontalAlignment','center','FontSize',FntSiz*1.1,'FontWeight','bold');

figname =sprintf('xR2-%s%s_%s_thresh-hist',modelVS,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

%-- shift bar position to centerize zero
h(1).BinLimits = h(1).BinLimits + h(1).BinWidth./2;
h(2).BinLimits = h(2).BinLimits + h(2).BinWidth./2;
    hz(1).YData = ht.Children(end).YLim;    hv(1).YData = ht.Children(end).YLim;
    hz(2).YData = ht.Children(end-1).YLim;    hv(2).YData = ht.Children(end-1).YLim;

figname =sprintf('xR2-%s%s_%s_thresh-histcenter',modelVS,R2mode,selectchs);
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
lll = -100; ull = 100;
nexttile;
    ll = ull; ul = lll;
    ii=1;
    for iroi = roilist'
    okch = ((okch1bb)|(okch2bb)) & ismember(prf_all(modelsets(2)).bb.channels.wangarea,iroi);
    scatter(prf_all(modelsets(2)).bb.xval(okch),prf_all(modelsets(1)).bb.xval(okch),plshp_wang(ii),'LineWidth',1.6);
    if sum(okch)
    dats = floor([prf_all(modelsets(2)).bb.xval(okch);prf_all(modelsets(1)).bb.xval(okch)]/roundunit)*roundunit;
    ll = nanmin([ll;dats(dats>=lll)]);
    dats = ceil([prf_all(modelsets(2)).bb.xval(okch);prf_all(modelsets(1)).bb.xval(okch)]/roundunit)*roundunit;
    ul = nanmax([ul;dats(dats<=ull)]);
    end
    hold on;
    ii=ii+1;
    end
    hs = flipud(findobj(get(gca,'Children'),'Type','Scatter'));
    axis([ll ul ll ul],'square');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    plot([getelement(xlim,1) threshold(modelsets(2)).bb],[1 1].*threshold(modelsets(1)).bb,'k-.',...
         [1 1].*threshold(modelsets(2)).bb,[getelement(ylim,1) threshold(modelsets(1)).bb],'k-.');  % threshold
%     legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
    
    set(gca,'FontSize',FntSiz);
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
    title('Broadband');
    xticks(yticks); xtickangle(0);
    
nexttile;
    ll = ull; ul = lll;
    ii=1;
    for iroi = roilist'
    okch = ((okch1a)|(okch2a)) & ismember(prf_all(modelsets(2)).a.channels.wangarea,iroi);
    scatter(prf_all(modelsets(2)).a.xval(okch),prf_all(modelsets(1)).a.xval(okch),plshp_wang(ii),'LineWidth',1.6);
    if sum(okch)
    dats = floor([prf_all(modelsets(2)).a.xval(okch);prf_all(modelsets(1)).a.xval(okch)]/roundunit)*roundunit;
    ll = nanmin([ll;dats(dats>=lll)]);
    dats = ceil([prf_all(modelsets(2)).a.xval(okch);prf_all(modelsets(1)).a.xval(okch)]/roundunit)*roundunit;
    ul = nanmax([ul;dats(dats<=ull)]);
    end
    hold on;
    ii=ii+1;
    end
    hs = flipud(findobj(get(gca,'Children'),'Type','Scatter'));
    axis([ll ul ll ul],'square');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    plot([getelement(xlim,1) threshold(modelsets(2)).a],[1 1].*threshold(modelsets(1)).a,'k-.',...
         [1 1].*threshold(modelsets(2)).a,[getelement(ylim,1) threshold(modelsets(1)).a],'k-.');  % threshold
    hl=legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
%     hl.Layout.Tile = 'east';
    
    set(gca,'FontSize',FntSiz);
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
    title('Alpha');
    xticks(yticks); xtickangle(0);
 
figname =sprintf('xR2-%s%s_%s_ROI',modelVS,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
