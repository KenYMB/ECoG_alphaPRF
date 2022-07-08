% Plot coherence between electrodes in HD grid across distance
%   for alpha & broadband with bootstrap

% 20210601 Yuasa - compute time series
%                  for alpha
% 20220207 Yuasa - bootstrap
% 20220208 Yuasa - bootstrap all paird electrodes
%                  need to run ecog_APRFF_03d_bootstrapCoherence in advance

%% %%%%%%%%%%%%%%%%%%%%
%% test
%% %%%%%%%%%%%%%%%%%%%%
%% prefix
% close all; clear all;
% startupToolboxToolbox;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'Connectivity-representative');
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
SetDefault('xspctrmPth',    'xSpectrum');
else
plotsavePth    = 'Connectivity-representative';
prfPth         = 'pRFmodel';
prfstatPth     = 'pRFanalysis';
xspctrmPth     = 'xSpectrum';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth));
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
HDsubjectList = SetSubjectsList(subjectList_fname, 'hasHDgrid','yes');

%-- Plotting Setting
FntSiz = 20;

%%
disttype    = 'norm';            % 'square','diamond','norm'
% useChans = 'pRFchs';        % pRFchs, SELchs, ALLchs
useChans = 'SELchs';        % pRFchs, SELchs, ALLchs
% arounddist  = [1 2 3 6];
    arounddist  = 1:6;

nboot = 5000;

if issaveplot,  plsbjs = 1:(length(HDsubjectList)+1);
else,           plsbjs = 1:length(HDsubjectList);
end
    
hF = gobjects(0);
for selsbj = plsbjs
%-- Dataset specs
if selsbj > length(HDsubjectList),  subjectList = HDsubjectList;
else,                               subjectList = HDsubjectList(selsbj);
end

%% load time series data
clear alphaType broadbandType
decN = 3;
% decN = 1;

average        ='runs';
prfmodel       = 'linear';
gaussianmode   = 'gs';
smoothingMode  ='decimate';
smoothingN     = decN;
% selectchs      = 'wangprobchs';     % only use wangprobchs as seed but use all grid channels for averaged coherence
selectchs      = 'GB*';             % use all grid channels
    allowlag       = false;
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;
    
ecog_APRFF_INITa_loaddata;

%%% Subject Name
nsbj = length(subjectList);
if nsbj==1
   subject = subjectList{1};
else
   subject = 'all';
end

%%% Bootstrapping file
filename = sprintf('cohboot-%s-%s-%s.mat',subject,useChans,disttype);
filepath = fullfile(SetDefaultAnalysisPath('DAT',xspctrmPth),filename);
assert(exist(filepath,'file'),'Need to run ecog_APRFF_03d_bootstrapCoherence in advance.');
load(filepath);

%%% Fitting file
filename = sprintf('cohboot-fit-%s-%s-%s.mat',subject,useChans,disttype);
filepath = fullfile(SetDefaultAnalysisPath('DAT',xspctrmPth),filename);
assert(exist(filepath,'file'),'Need to run ecog_APRFF_03f_fitCoherence in advance.');
load(filepath);

%%
%%
%%
%% %%%%%%%%%%%%%
%% Visualize


%% parameter for plot
   alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd
   alphaD= 0.05;   % 0.05 for 2sd, 0.32 for 1sd
   plcol = get(groot,'defaultAxesColorOrder');
   %%%
   
   switch useChans
       case {'pRFchs','SELchs'},   chantitle = 'pRF channels';
       case {'ALLchs'},            chantitle = 'All channels';
       otherwise,                  chantitle = '';
   end
   
    %%% color for blank
    mkblcl = @(c) mean([c;ones(3,3)*0.6],1);
%     mkblcl = @(c) c;

    %%% plot type
    tarBAND   = 'both';
    inPRFmode = 'AB';

%% Coherence across distanece
%%% Plot Raw w/ fitting
ylintick  = 0:0.1:1.0;
% ylinraw   = [0.25 0.6];
ylinraw   = [0.25 1.0];
xs = linspace(0,max(arounddist),100);
% cohfit = @(x,xdata)x(:,1).*exp(x(:,2).*xdata)+x(:,3);

hF(end+1) = figure('MenuBar','none');
%-- fitting curve
yd = cohbootBBall;
ys = cohfit(fitparamsBBall,xs);
    yd1 = mean(yd,1,'omitnan');
    yd2 = prctile(yd,(alpha/2)*100,1);
    yd3 = prctile(yd,(1-alpha/2)*100,1);
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl = plot(xs,y1,'-','Color',plcol(1,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([xs, fliplr(xs)], [y2, fliplr(y3)],plcol(1,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl,'top');
  hd = errorbar(arounddist,yd1,yd1-yd2,yd3-yd1,'.','Color',plcol(1,:),'LineWidth',1.5,'MarkerSize',35);
%-- fitting curve
yd = cohbootAall;
ys = cohfit(fitparamsAall,xs);
    yd1 = mean(yd,1,'omitnan');
    yd2 = prctile(yd,(alpha/2)*100,1);
    yd3 = prctile(yd,(1-alpha/2)*100,1);
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(2) = plot(xs,y1,'-','Color',plcol(2,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([xs, fliplr(xs)], [y2, fliplr(y3)],plcol(2,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl,'top');
  hd(2) = errorbar(arounddist,yd1,yd1-yd2,yd3-yd1,'.','Color',plcol(2,:),'LineWidth',1.5,'MarkerSize',35);
set(gca,'FontSize',FntSiz);
set(hd,{'CapSize'},{18});
legend(hl,{'Broadband','Alpha'},'Location','northeast');
xlabel('Distance (mm)'); ylabel('Coherence');
xlim(minmax([0 arounddist]));
xticklabels(xticks*3);      % multiply 3 to convert mm
ylim(ylinraw);
yticks(ylintick);

   figureName = sprintf('CoherenceTS_%s-%s-%s-%s_Distance%s-ALL',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figureName = sprintf('%s-decimate%d',figureName,decN); end
   set(gcf,'Name',figureName);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figureName));  end

   
%%% Plot Raw w/ fitting w/o band
if issaveplot
ylintick  = 0:0.1:1.0;
% ylinraw   = [0.25 0.6];
ylinraw   = [0.25 1.0];
xs = linspace(0,max(arounddist),100);

hF(end+1) = figure('MenuBar','none');
ys = cohbootBBall;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl = errorbar(arounddist,y1,y1-y2,y3-y1,'.','Color',plcol(1,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
        %-- fitting curve
        params0 = [max(y1),-1,min(y1)];
        params = lsqcurvefit(cohfit,params0,arounddist,y1);
        plot(xs,cohfit(params,xs),'--','Color',plcol(1,:),'LineWidth',1.5);
ys = cohbootAall;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(2) = errorbar(arounddist,y1,y1-y2,y3-y1,'.','Color',plcol(2,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
        %-- fitting curve
        params0 = [max(y1),-1,min(y1)];
        params = lsqcurvefit(cohfit,params0,arounddist,y1);
        plot(xs,cohfit(params,xs),'--','Color',plcol(2,:),'LineWidth',1.5);
%     hb = plot(xlim,[1,1].*chcavgtsAdat,'k--','LineWidth',1.2);
%     uistack(hb,'bottom');
set(hl,{'CapSize'},{18});
set(gca,'FontSize',FntSiz);
legend(hl,{'Broadband','Alpha'},'Location','northeast');
xlabel('Distance (mm)'); ylabel('Coherence');
xlim(minmax([0 arounddist]));
xticklabels(xticks*3);      % multiply 3 to convert mm
ylim(ylinraw);
yticks(ylintick);

   figureName = sprintf('CoherenceTS_%s-%s-%s-%s_Distance2%s-ALL',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figureName = sprintf('%s-decimate%d',figureName,decN); end
   set(gcf,'Name',figureName);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figureName));  end
end


%% Difference of Coherence between pRF conditions
if issaveplot
%%% difference plot @ x=3
irnd = 1;       % index of arounddist
ylindiff  = [-1 1].*0.04./nsbj;

hF(end+1) = figure('Menubar','none');
set(gcf,'Position',get(gcf,'Position').*[1 1 0.75 1]);
ys = [cohbootBBprf(:,irnd) - cohbootBBbsl(:,irnd), ...
      cohbootBBout(:,irnd) - cohbootBBbsl(:,irnd)];
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alphaD/2)*100,1);
    y3 = prctile(ys,(1-alphaD/2)*100,1);
    xs = [-1 1]*.2 +2;
hb=errorbar(xs,y1,y1-y2,y3-y1, 'o');
ys = [cohbootAprf(:,irnd) - cohbootAbsl(:,irnd), ...
      cohbootAout(:,irnd) - cohbootAbsl(:,irnd)];
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alphaD/2)*100,1);
    y3 = prctile(ys,(1-alphaD/2)*100,1);
    xs = [-1 1]*.2 +1;
hold on;
hb(2)=errorbar(xs,y1,y1-y2,y3-y1, 'o');
ylim(ylindiff); xlim([0.5 2.5]);
xticks(sort(cat(2,hb.XData)));
xticklabels([{'in'},{'out'},{'in'},{'out'}]);

hold on; hl = plot(xlim,[0 0],'k:','LineWidth',1.2);  uistack(hl,'bottom');
set(hb,{'linew'},{2.6},{'MarkerSize'},{9},{'MarkerFaceColor'},{'auto'},{'CapSize'},{10});
set(gca,'FontSize',FntSiz);
title(sprintf('%s',chantitle));
ylabel(sprintf('Different Coherence from BLANK'));
AX = gca; AX.XAxis.FontSize=AX.XAxis.FontSize-2;
legend(hb,{'Boradband','Alpha'},'Location','northwest');
   
   figureName = sprintf('CoherenceTS_%s-%s-%s-%s_Dot%s-diff',subject,useChans,'both',disttype,'AB');
   if decN>1, figureName = sprintf('%s-decimate%d',figureName,decN); end
   set(gcf,'Name',figureName);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figureName));  end


%%% difference plot @ x=3 (est)
if issaveplot
xtim = 1;       % distance
ylindiff  = [-1 1].*0.04./nsbj;

hF(end+1) = figure('Menubar','none');
set(gcf,'Position',get(gcf,'Position').*[1 1 0.75 1]);
ys = [cohfit(fitparamsBBprf,xtim) - cohfit(fitparamsBBbsl,xtim), ...
      cohfit(fitparamsBBout,xtim) - cohfit(fitparamsBBbsl,xtim)];
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alphaD/2)*100,1);
    y3 = prctile(ys,(1-alphaD/2)*100,1);
    xs = [-1 1]*.2 +2;
hb=errorbar(xs,y1,y1-y2,y3-y1, 'o');
ys = [cohfit(fitparamsAprf,xtim) - cohfit(fitparamsAbsl,xtim), ...
      cohfit(fitparamsAout,xtim) - cohfit(fitparamsAbsl,xtim)];
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alphaD/2)*100,1);
    y3 = prctile(ys,(1-alphaD/2)*100,1);
    xs = [-1 1]*.2 +1;
hold on;
hb(2)=errorbar(xs,y1,y1-y2,y3-y1, 'o');
ylim(ylindiff); xlim([0.5 2.5]);
xticks(sort(cat(2,hb.XData)));
xticklabels([{'in'},{'out'},{'in'},{'out'}]);

hold on; hl = plot(xlim,[0 0],'k:','LineWidth',1.2);  uistack(hl,'bottom');
set(hb,{'linew'},{2.6},{'MarkerSize'},{9},{'MarkerFaceColor'},{'auto'},{'CapSize'},{10});
set(gca,'FontSize',FntSiz);
title(sprintf('%s',chantitle));
ylabel(sprintf('Different Coherence from BLANK'));
AX = gca; AX.XAxis.FontSize=AX.XAxis.FontSize-2;
legend(hb,{'Boradband','Alpha'},'Location','northwest');
   
   figureName = sprintf('CoherenceTS_%s-%s-%s-%s_Dot%s-diff',subject,useChans,'both',disttype,'AB');
   if decN>1, figureName = sprintf('%s-decimate%d',figureName,decN); end
   set(gcf,'Name',figureName);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figureName));  end
end

%%% difference plot @ x=0 (est)
if issaveplot
xtim = 0;       % distance
ylindiff  = [-1 1].*0.2./nsbj;

hF(end+1) = figure('Menubar','none');
set(gcf,'Position',get(gcf,'Position').*[1 1 0.75 1]);
ys = [cohfit(fitparamsBBprf,xtim) - cohfit(fitparamsBBbsl,xtim), ...
      cohfit(fitparamsBBout,xtim) - cohfit(fitparamsBBbsl,xtim)];
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alphaD/2)*100,1);
    y3 = prctile(ys,(1-alphaD/2)*100,1);
    xs = [-1 1]*.2 +2;
hb=errorbar(xs,y1,y1-y2,y3-y1, 'o');
ys = [cohfit(fitparamsAprf,xtim) - cohfit(fitparamsAbsl,xtim), ...
      cohfit(fitparamsAout,xtim) - cohfit(fitparamsAbsl,xtim)];
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alphaD/2)*100,1);
    y3 = prctile(ys,(1-alphaD/2)*100,1);
    xs = [-1 1]*.2 +1;
hold on;
hb(2)=errorbar(xs,y1,y1-y2,y3-y1, 'o');
ylim(ylindiff); xlim([0.5 2.5]);
xticks(sort(cat(2,hb.XData)));
xticklabels([{'in'},{'out'},{'in'},{'out'}]);

hold on; hl = plot(xlim,[0 0],'k:','LineWidth',1.2);  uistack(hl,'bottom');
set(hb,{'linew'},{2.6},{'MarkerSize'},{9},{'MarkerFaceColor'},{'auto'},{'CapSize'},{10});
set(gca,'FontSize',FntSiz);
title(sprintf('%s',chantitle));
ylabel(sprintf('Different Coherence from BLANK'));
AX = gca; AX.XAxis.FontSize=AX.XAxis.FontSize-2;
legend(hb,{'Boradband','Alpha'},'Location','northwest');
   
   figureName = sprintf('CoherenceTS_%s-%s-%s-%s_Dot%s-diff-est0',subject,useChans,'both',disttype,'AB');
   if decN>1, figureName = sprintf('%s-decimate%d',figureName,decN); end
   set(gcf,'Name',figureName);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figureName));  end
end
end

end
