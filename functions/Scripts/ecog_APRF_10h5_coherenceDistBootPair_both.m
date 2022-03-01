% compute coherence across electrodes for subjects with HDgrid

% 20210601 Yuasa - compute time series (test)
%                  for alpha
% 20220207 Yuasa - bootstrap
% 20220208 Yuasa - bootstrap all paird electrodes
%                  need to run ecog_APRFF_03d_bootstrapCoherence in advance

%% Define paths and dataset
% close all; clear all;
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
SetDefault('issaveplot',true);
if issaveplot
    plotsavedir    = fullfile(figPth, 'Coherence-boot');
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectLists = reshape(SetSubjectsList(subjectList_fname, 'hasHDgrid','yes'),1,[]);

%% %%%%%%%%%%%%%%%%%%%%
%% test
%% %%%%%%%%%%%%%%%%%%%%
disttype    = 'norm';            % 'square','diamond','norm'
% useChans = 'pRFchs';        % pRFchs, SELchs, ALLchs
useChans = 'SELchs';        % pRFchs, SELchs, ALLchs
% arounddist  = [1 2 3 6];
    arounddist  = 1:6;

nboot = 5000;

hF = [];
for subjectList = cellfun(@(C) {{C}},subjectLists) %[{subjectLists},cellfun(@(C) {{C}},subjectLists)]
subjectList = subjectList{:};
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
    
ecog_APRF_INITa_loaddata;

%%% Subject Name
if length(subjectList)==1
   subject = subjectList{1};
else
   subject = 'all';
end

%%% Bootstrapping file
filename = sprintf('cohboot-%s-%s-%s.mat',subject,useChans,disttype);
filepath = fullfile(analysisRootPath, 'Data', 'xSpectrum',filename);

%%% Coherence
load(filepath);
%%
%%
% close all
%%
%% %%%%%%%%%%%%%
%% Visualize


%% parameter for plot
   FntSiz = 20;
   alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd
   plcol = get(groot,'defaultAxesColorOrder');
   %%%
   
    %%% color for blank
    mkblcl = @(c) mean([c;ones(3,3)*0.6],1);
%     mkblcl = @(c) c;

%%
tarBAND   = 'both';
inPRFmode = 'AB';

%%% Plot Raw
ylinraw   = [0.25 0.6];
ylintick  = 0.2:0.1:0.7;

hF = [hF,figure('MenuBar','none')];
ys = cohbootBBprf;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl = plot(arounddist,y1,'.-','Color',plcol(1,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(1,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl,'top');
ys = cohbootAprf;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(2) = plot(arounddist,y1,'.-','Color',plcol(2,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(2,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(2),'top');
%     hb = plot(xlim,[1,1].*chcavgtsAdat,'k--','LineWidth',1.2);
%     uistack(hb,'bottom');
set(gca,'FontSize',FntSiz);
legend(hl,{'Broadband','Alpha'},'Location','northeast');
xlabel('Distance'); ylabel('Coherence');
xlim(minmax(arounddist));
ylim(ylinraw);
yticks(ylintick);

   figname = sprintf('CoherenceTS_%s-%s-%s-%s_Distance%s',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));  end

%%% difference plot
irnd = 1;       % index of arounddist
ylindiff  = [-1 1].*0.04;

hF = [hF,figure('MenuBar','none')];
set(gcf,'Position',get(gcf,'Position').*[1 1 0.75 1]);
hb=boxplot([cohbootAprf(:,irnd)-cohbootAbsl(:,irnd);cohbootAout(:,irnd)-cohbootAbsl(:,irnd)],...
        [repmat({'in pRF - BLANK'},nboot,1);...
         repmat({'out pRF - BLANK'},nboot,1)],...
         'Width',0.5);
% set(gcf,'Position',get(gcf,'Position').*[1 1 1.05 1]);
% hb=boxplot([cohbootAprf(:,irnd)-cohbootAbsl(:,irnd);cohbootAout(:,irnd)-cohbootAbsl(:,irnd);cohbootAprf(:,irnd)-cohbootAout(:,irnd)],...
%         [repmat({'in pRF - BLANK'},nboot,1);...
%          repmat({'out pRF - BLANK'},nboot,1);...
%          repmat({'in pRF - out pRF'},nboot,1)],...
%          'Width',0.5);
ylim(ylindiff);
hold on; hl = plot(xlim,[0 0],'k:','LineWidth',1.2);  uistack(hl,'bottom');
set(hb,{'linew'},{1.6});
set(gca,'FontSize',FntSiz);
title(sprintf('%s','pRF channels'));
ylabel(sprintf('Coherence - Peak Alpha Frequency'));
AX = gca; AX.XAxis.FontSize=AX.XAxis.FontSize-2;
   
   figname = sprintf('CoherenceTS_%s-%s-%s-%s_BoxDist%s-diff2',subject,useChans,'alpha',disttype,'A');
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));  end


end
