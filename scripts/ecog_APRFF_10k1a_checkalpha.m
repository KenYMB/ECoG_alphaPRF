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
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth), 'modelselection', 'noACorr');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');
 
%-- Plotting Setting
FntSiz    = 20;
SFntSiz   = 14;

%% general parameter
clear alphaType broadbandType

average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';
selectchs      = 'wangprobchs';
    allowlag       = false;
va_area = 'wangarea';

usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

%% load analyzePRF @DOG,OG,LFS
noACorrects = [false, true];
targetBANDs = {'FaCLb','FaLb'};
% modelNames  = ["Alpha Suppression","Power Change in Alpha"];
modelNames  = ["Model-based Alpha Suppression","Power Change in Alpha"];

modeldata  = struct();
prf_params = struct();
model_all  = struct();
prf_all    = struct();
threshold  = struct();
for ii=length(noACorrects):-1:1
noACorrect   = noACorrects(ii);
if noACorrect
    allowbeta      = false;
    allowwide      = false;
    allowmixbeta   = false;
else
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;
end
 
clear alphaType broadbandType
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
%% FaCLb(BW) vs FaLb
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelsets = [1 2];      % FaCLb(BW) vs FaLb
targetBAND = strjoin(targetBANDs,'VS');

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
if issaveplot
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
 
figname =sprintf('xR2-%s%s_%s_all',targetBAND,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Test Performance w/ threshold
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
%% Histogram of difference
if issaveplot
% roundunit = 5;

% okch = (okch1bb&okch1a)|(okch2bb&okch2a);
%-- broadband
hF(end+1) = figure('Menubar','none','Position',[200 200 800 420]);
ht = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile;  okch = (okch1a)|(okch2a);
xval_diff = prf_all(modelsets(1)).a.xval(okch) - prf_all(modelsets(2)).a.xval(okch);
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
title('Alpha');
%-- alpha
nexttile;  okch = (okch1a&okch1bb)|(okch2a&okch2bb);
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
title('Alpha with Broadband threshold');
tt=title(ht,sprintf('%s – %s',modelNames(modelsets(1)),modelNames(modelsets(2))),...
    'HorizontalAlignment','center','FontSize',FntSiz*1.1,'FontWeight','bold');

figname =sprintf('xR2-%s%s_%s_thresh-hist',targetBAND,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

%-- shift bar position to centerize zero
h(1).BinLimits = h(1).BinLimits + h(1).BinWidth./2;
h(2).BinLimits = h(2).BinLimits + h(2).BinWidth./2;
    hz(1).YData = ht.Children(end).YLim;    hv(1).YData = ht.Children(end).YLim;
    hz(2).YData = ht.Children(end-1).YLim;    hv(2).YData = ht.Children(end-1).YLim;

figname =sprintf('xR2-%s%s_%s_thresh-histcenter',targetBAND,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Test Performance w/ threshold in each ROI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;

%% Difference of Gaussian Model in training accuracy (scatter) w/ rough ROI
if issaveplot
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
%     legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
    
    set(gca,'FontSize',FntSiz);
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
    title('Alpha');
    xticks(yticks); xtickangle(0);
    
nexttile;
    ll = 100; ul = -100;
    ii=1;
    for iroi = roilist'
    okch = ((okch1a&okch1bb)|(okch2a&okch2bb)) & ismember(prf_all(modelsets(2)).a.channels.wangarea,iroi);
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
    hl=legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
%     hl.Layout.Tile = 'east';
    
    set(gca,'FontSize',FntSiz);
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
    title('Alpha with Broadband threshold');
    xticks(yticks); xtickangle(0);
 
figname =sprintf('xR2-%s%s_%s_ROI',targetBAND,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Test Performance w/ threshold in V1-V3
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;

ispltall = false;
%% Difference of Gaussian Model in training accuracy (scatter) w/ rough ROI
if issaveplot,      npanel = 2;
else,               npanel = 1;
end

plcol_wang = [ones(3,1)*plcol(1,:);...
              ones(6,1)*plcol(4,:);];
plshp_wang = [repmat('o',1,3) repmat('^',1,6)];
if ispltall
    plcol_gray = [plcol_wang .* 0  +  ones(1,3).*0.7 .* 2]./2;
end
          
roundunit = 20;
% hF(end+1) = figure('Menubar','none','Position',[200 200 860 420],'defaultAxesColorOrder',plcol_wang);
hF(end+1) = figure('Menubar','none','Position',[200 200 450*npanel 440],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(1,npanel,'Padding','compact','TileSpacing','compact');
roilist = {'V1','V2','V3'}';
nexttile;
    ll = 100; ul = -100;
    ii=1;
    for iroi = roilist'
    roich = ismember(prf_all(modelsets(2)).a.channels.wangarea,iroi);
    okch  = ((okch1a)|(okch2a)) & roich;
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
        okch  = ~((okch1a)|(okch2a)) & roich;
        scatter(prf_all(modelsets(2)).a.xval(okch),prf_all(modelsets(1)).a.xval(okch),plshp_wang(ii),'LineWidth',1.6,'MarkerEdgeColor',plcol_gray(ii,:));
        hold on;
        ii=ii+1;
        end
    end
    uistack(hs,'top');
    ll = max(-100,ll);
    ul = min(100,ul);
    axis([ll ul ll ul],'square');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    plot([getelement(xlim,1) threshold(modelsets(2)).a],[1 1].*threshold(modelsets(1)).a,'k-.',...
         [1 1].*threshold(modelsets(2)).a,[getelement(ylim,1) threshold(modelsets(1)).a],'k-.');  % threshold
%     legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
    
    set(gca,'FontSize',FntSiz);
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
    title('Alpha');
    xticks(yticks); xtickangle(0);
    
if npanel > 1
nexttile;
    ll = 100; ul = -100;
    ii=1;
    for iroi = roilist'
    roich = ismember(prf_all(modelsets(2)).a.channels.wangarea,iroi);
    okch = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich;
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
        okch = ~((okch1a&okch1bb)|(okch2a&okch2bb)) & roich;
        scatter(prf_all(modelsets(2)).a.xval(okch),prf_all(modelsets(1)).a.xval(okch),plshp_wang(ii),'LineWidth',1.6,'MarkerEdgeColor',plcol_gray(ii,:));
        hold on;
        ii=ii+1;
        end
    end
    uistack(hs,'top');
    ll = max(-100,ll);
    ul = min(100,ul);
    axis([ll ul ll ul],'square');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    plot([getelement(xlim,1) threshold(modelsets(2)).a],[1 1].*threshold(modelsets(1)).a,'k-.',...
         [1 1].*threshold(modelsets(2)).a,[getelement(ylim,1) threshold(modelsets(1)).a],'k-.');  % threshold
%     hl=legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
%     hl.Layout.Tile = 'east';
    
    set(gca,'FontSize',FntSiz);
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
    title('Alpha with Broadband threshold');
    xticks(yticks); xtickangle(0);
end
 
figname =sprintf('xR2-%s%s_%s_ROI-low',targetBAND,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

%%
% close all;

%% %%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%
%% pRF plots
%% %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
if issaveplot
goodch = (okch1a&okch1bb)&~(okch2a&okch2bb);
badch  = ~(okch1a&okch1bb)&(okch2a&okch2bb);
bothgoodch = (okch1a&okch1bb)&(okch2a&okch2bb);

%-- bad in correction
opts = [];
opts.plot.pix2deg = cfactor;
opts.plot.XLim    = [-1 1].*12;
opts.plot.YLim    = [-1 1].*12;
opts.plot.addChsToTitle     = 'yes';
opts.plot.addSbjToTitle     = 'yes';
opts.plot.addBensonToTitle  = 'no';
opts.plot.addWangToTitle    = 'no';
opts.plot.fontSize          = SFntSiz;
opts.plot.nSubPlots = [0 min(7,numel(badch))];

opts.plot.showaxis = 0;         % if show X & Y axis
ecog_plotGridPRF(badch, opts, prf_all(1).bb,prf_all(1).a,prf_all(2).a);

set(gcf,'MenuBar','none');

figname =sprintf('pRF-%s%s_%s_location-%s',targetBAND,R2mode,selectchs,targetBANDs{2});
set(gcf,'Name',figname);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

%-- good in correction
opts = [];
opts.plot.pix2deg = cfactor;
opts.plot.XLim    = [-1 1].*12;
opts.plot.YLim    = [-1 1].*12;
opts.plot.addChsToTitle     = 'yes';
opts.plot.addSbjToTitle     = 'yes';
opts.plot.addBensonToTitle  = 'no';
opts.plot.addWangToTitle    = 'no';
opts.plot.fontSize          = SFntSiz;
opts.plot.nSubPlots = [0 min(7,numel(goodch))];

opts.plot.showaxis = 0;         % if show X & Y axis
ecog_plotGridPRF(goodch, opts, prf_all(1).bb,prf_all(1).a,prf_all(2).a);

set(gcf,'MenuBar','none');

figname =sprintf('pRF-%s%s_%s_location-%s',targetBAND,R2mode,selectchs,targetBANDs{1});
set(gcf,'Name',figname);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

%-- both good
opts = [];
opts.plot.pix2deg = cfactor;
opts.plot.XLim    = [-1 1].*12;
opts.plot.YLim    = [-1 1].*12;
opts.plot.addChsToTitle     = 'yes';
opts.plot.addSbjToTitle     = 'yes';
opts.plot.addBensonToTitle  = 'no';
opts.plot.addWangToTitle    = 'no';
opts.plot.fontSize          = SFntSiz;
opts.plot.nSubPlots = [0 min(7,numel(bothgoodch))];

opts.plot.showaxis = 0;         % if show X & Y axis
ecog_plotGridPRF(bothgoodch, opts, prf_all(1).bb,prf_all(1).a,prf_all(2).a);

set(gcf,'MenuBar','none');

figname =sprintf('pRF-%s%s_%s_location-%s',targetBAND,R2mode,selectchs,'Both');
set(gcf,'Name',figname);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

end
