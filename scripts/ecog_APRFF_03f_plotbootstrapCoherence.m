% compute connectivities across electrodes for subjects with HDgrid
%   investigate how to analyze based on a patient across all electrodes

% 20210601 Yuasa - compute time series
%                  for alpha
% 20220207 Yuasa - bootstrap
% 20220208 Yuasa - bootstrap all paird electrodes
%                  need to run ecog_APRFF_03d_bootstrapCoherence in advance
% 20221027 Yuasa - enable to change some variables from outside of the script

%%
close all; % clearvars;
% if isempty(gcp('nocreate')),  parpool([1 40]); end
% startupToolboxToolbox;
%% Define paths and dataset
checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
plotsavePth    = 'Connectivity-boot';
xspctrmPth     = 'xSpectrum';
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'Distance');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
HDsubjectList = SetSubjectsList(subjectList_fname, 'hasHDgrid','yes');

%%
SetDefault('cohmethod','mscoh');      % 'mscoh', 'imcoh'
SetDefault('disttype','norm');        % 'square','diamond','norm'
SetDefault('useChans','SELchs');      % 'pRFchs', 'SELchs', 'ALLchs'
    
for selsbj = 1:(length(HDsubjectList)+1)
%-- Dataset specs
if selsbj > length(HDsubjectList),  subjectList = HDsubjectList;
else,                               subjectList = HDsubjectList(selsbj);
end

%% load coherence
%%% Subject Name
nsbj = length(subjectList);
if nsbj==1
   subject = subjectList{1};
else
   subject = 'all';
end

%%% Bootstrapping file
filename = sprintf('%sboot_%s-%s-%s.mat',cohmethod,subject,useChans,disttype);
filepath = fullfile(SetDefaultAnalysisPath('DAT',xspctrmPth),filename);

%%% Coherence
load(filepath);
%%
%%
close all
%%
%% %%%%%%%%%%%%%
%% Visualize


%% parameter for plot
   FntSiz = 20;
   alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd
   plcol = get(groot,'defaultAxesColorOrder');
   decN = 1;
   %%%
   
   switch useChans
       case {'pRFchs','SELchs'},   chantitle = 'pRF channels';
       case {'ALLchs'},            chantitle = 'All channels';
       otherwise,                  chantitle = '';
   end
   
    %%% color for blank
    mkblcl = @(c) mean([c;ones(3,3)*0.6],1);
%     mkblcl = @(c) c;

%%
tarBAND   = 'both';
inPRFmode = 'AB';

%%% Plot Raw
ylinraw   = [0.25 0.6];
ylintick  = 0.2:0.1:0.7;

%%-- IN-pRF
figure('MenuBar','none');
ys = cohbootBBprf;  ih = 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl = plot(arounddist,y1,'.-','Color',plcol(1,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(1,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
ys = cohbootAprf;  ih = ih + 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(ih) = plot(arounddist,y1,'.-','Color',plcol(2,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(2,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
%     hb = plot(xlim,[1,1].*chcavgtsAdat,'k--','LineWidth',1.2);
%     uistack(hb,'bottom');
set(gca,'FontSize',FntSiz);
legend(hl,{'Broadband','Alpha'},'Location','northeast');
xlabel('Distance'); ylabel('Coherence');
xlim(minmax(arounddist));
ylim(ylinraw);
yticks(ylintick);

   figname = sprintf('CoherenceTS_%s-%s-%s-%s_Distance%s-inPRF',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figname),{'png','eps'});  end

%%-- BLANK
figure('MenuBar','none');
ys = cohbootBBbsl;  ih = 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl = plot(arounddist,y1,'.-','Color',plcol(1,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(1,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
ys = cohbootAbsl;  ih = ih + 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(ih) = plot(arounddist,y1,'.-','Color',plcol(2,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(2,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
%     hb = plot(xlim,[1,1].*chcavgtsAdat,'k--','LineWidth',1.2);
%     uistack(hb,'bottom');
set(gca,'FontSize',FntSiz);
legend(hl,{'Broadband','Alpha'},'Location','northeast');
xlabel('Distance'); ylabel('Coherence');
xlim(minmax(arounddist));
ylim(ylinraw);
yticks(ylintick);

   figname = sprintf('CoherenceTS_%s-%s-%s-%s_Distance%s-BLANK',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figname),{'png','eps'});  end   

%%-- ALL
figure('MenuBar','none');
ys = cohbootBBall;  ih = 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl = plot(arounddist,y1,'.-','Color',plcol(1,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(1,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
ys = cohbootAall;  ih = ih + 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(ih) = plot(arounddist,y1,'.-','Color',plcol(2,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(2,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
%     hb = plot(xlim,[1,1].*chcavgtsAdat,'k--','LineWidth',1.2);
%     uistack(hb,'bottom');
set(gca,'FontSize',FntSiz);
legend(hl,{'Broadband','Alpha'},'Location','northeast');
xlabel('Distance'); ylabel('Coherence');
xlim(minmax(arounddist));
ylim(ylinraw);
yticks(ylintick);

   figname = sprintf('CoherenceTS_%s-%s-%s-%s_Distance%s-ALL',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figname),{'png','eps'});  end   
   
   
%%-- Three
lighten = @(color,p) color.*p + ones(size(color)).*(1-p);
darken  = @(color,p) color.*p + zeros(size(color)).*(1-p);

figure('MenuBar','none');
ys = cohbootBBbsl;  ih = 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl = plot(arounddist,y1,'.-','Color',plcol(1,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(1,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
ys = cohbootBBout;  ih = ih + 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(ih) = plot(arounddist,y1,'.-','Color',darken(plcol(1,:),0.5),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],darken(plcol(1,:),0.5),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
ys = cohbootBBprf;  ih = ih + 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(ih) = plot(arounddist,y1,'.-','Color',lighten(plcol(1,:),0.5),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],lighten(plcol(1,:),0.5),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
  
ys = cohbootAbsl;  ih = ih + 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(ih) = plot(arounddist,y1,'.-','Color',plcol(2,:),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(2,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
ys = cohbootAout;  ih = ih + 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(ih) = plot(arounddist,y1,'.-','Color',darken(plcol(2,:),0.5),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],darken(plcol(2,:),0.5),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
ys = cohbootAprf;  ih = ih + 1;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(ih) = plot(arounddist,y1,'.-','Color',lighten(plcol(2,:),0.5),'LineWidth',1.5,'MarkerSize',35);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],lighten(plcol(2,:),0.5),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(ih),'top');
  
%     hb = plot(xlim,[1,1].*chcavgtsAdat,'k--','LineWidth',1.2);
%     uistack(hb,'bottom');
set(gca,'FontSize',FntSiz);
legend(hl,{'Broadband: BLANK','Broadband: out-pRF','Broadband: in-pRF',...
           'Alpha: BLANK','Alpha: out-pRF','Alpha: in-pRF'},'Location','northeast');
xlabel('Distance'); ylabel('Coherence');
xlim(minmax(arounddist));
ylim(ylinraw);
yticks(ylintick);

   figname = sprintf('CoherenceTS_%s-%s-%s-%s_Distance%s-3conditions',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figname),{'png','eps'});  end   
   
%% difference plot
irnd = 1;       % index of arounddist
ylindiff  = [-1 1].*0.04;

%%-- alpha
figure('Menubar','none');
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
title(sprintf('%s',chantitle));
ylabel(sprintf('Different Coherence'));
AX = gca; AX.XAxis.FontSize=AX.XAxis.FontSize-2;
   
   figname = sprintf('CoherenceTS_%s-%s-%s-%s_BoxDist%s-diff2',subject,useChans,'alpha',disttype,'A');
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figname),{'png','eps'});  end
   
%%-- broadband
figure('Menubar','none');
set(gcf,'Position',get(gcf,'Position').*[1 1 0.75 1]);
hb=boxplot([cohbootBBprf(:,irnd)-cohbootBBbsl(:,irnd);cohbootBBout(:,irnd)-cohbootBBbsl(:,irnd)],...
        [repmat({'in pRF - BLANK'},nboot,1);...
         repmat({'out pRF - BLANK'},nboot,1)],...
         'Width',0.5);
ylim(ylindiff);
hold on; hl = plot(xlim,[0 0],'k:','LineWidth',1.2);  uistack(hl,'bottom');
set(hb,{'linew'},{1.6});
set(gca,'FontSize',FntSiz);
title(sprintf('%s',chantitle));
ylabel(sprintf('Different Coherence'));
AX = gca; AX.XAxis.FontSize=AX.XAxis.FontSize-2;
   
   figname = sprintf('CoherenceTS_%s-%s-%s-%s_BoxDist%s-diff2',subject,useChans,'broadband',disttype,'BB');
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figname),{'png','eps'});  end

%% violin plot
% irnd = 1;       % index of arounddist
% ylindiff  = [-1 1].*0.04;
% 
% %%-- alpha
% figure('Menubar','none');
% set(gcf,'Position',get(gcf,'Position').*[1 1 0.75 1]);
% pltdat = cat(3,cat(2,cohbootBBprf(:,irnd)-cohbootBBbsl(:,irnd),cohbootBBout(:,irnd)-cohbootBBbsl(:,irnd)),...
%                cat(2,cohbootAprf(:,irnd)-cohbootAbsl(:,irnd),cohbootAout(:,irnd)-cohbootAbsl(:,irnd)));
% hb=violinplotsplit(pltdat,...
%         [{'in pRF - BLANK'},{'out pRF - BLANK'}],...
%          'ViolinAlpha',0.6 ,'Width',0.4, 'ShowData',false, 'ShowMean',true);
% ylim(ylindiff);
% set(gca,'FontSize',FntSiz);
% title(sprintf('%s',chantitle));
% ylabel(sprintf('Different Coherence'));
% AX = gca; AX.XAxis.FontSize=AX.XAxis.FontSize-2;
% hold on; hl = plot(xlim,[0 0],'k:','LineWidth',1.2);  uistack(hl,'bottom');
%    
%    figname = sprintf('CoherenceTS_%s-%s-%s-%s_Violin%s-diff2',subject,useChans,'both',disttype,'AB');
%    if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
%    set(gcf,'Name',figname);
%    if issaveplot,  savefigauto(gcf, fullfile(plotsavedir,figname),{'png','eps'});  end

end
