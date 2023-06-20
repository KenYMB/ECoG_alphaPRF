%% ECoG Alpha pRF (simple broadband computation)
% plot bootstrap results
% for figure 4

% 20201215
% 20210224 - update for finalize figure 4
% 20210508 - use another y axis
% %% without ERP %%

%%
% close all; clear all;
% startupToolboxToolbox;
run_checkPath;

%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'pRF-representative');
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
else
plotsavePth    = 'pRF-representative';
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
LFntSiz = 23;

%% Representative channels
if issaveplot
repelec = {'name',{'Oc18','GB103','Oc17','GB102'},...
            'subject_name',{'p02','p10','p02','p10'}};
else
% repelec = {'name',{'Oc17'},'subject_name',{'p02'}};
repelec = {'name',{'GB103'},'subject_name',{'p10'}};
end

selsbj      = ismember(subjectList,repelec{4});
subjectList = subjectList(selsbj);

selset      = ismember(repelec{4},subjectList);
repelec{2}  = repelec{2}(selset);
repelec{4}  = repelec{4}(selset);

%% load analyzePRF (w/ bootstrap)
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

%%% load original pRF data
modeldataID = [];
prfID       = [];
usefulltsxR2 = false;
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;

%% set threshold
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;
thresh_ecc = eclimit;
% thresh_ecc = nan;

elec_ok = ~(prf_all_bb.xval<=threshold_bb | prf_all_a.xval<=threshold_a ...
         | prf_all_bb.ecc >= thresh_ecc | prf_all_a.ecc >= thresh_ecc); 

%% visualize parameter
res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

plcol = get(groot,'defaultAxesColorOrder');

hF = gobjects(0);

%% plot averaged pRF time-series

for tarBANDs = {'broadband', 'alpha-noflip'}
tarBAND = tarBANDs{:};

switch tarBAND
    case {'broadband'}
        mdlall      = model_all_bb;
        prfall      = prf_all_bb;
        cntcol      = plcol(1,:);
        tname  = 'Broadband';
    case {'alpha','alpha-noflip'}
        mdlall      = model_all_a;
        prfall      = prf_all_a;
        cntcol      = plcol(2,:);
        tname  = 'Alpha';
end

el = findtable(prfall.channels,repelec{:});
for ii=1:numel(repelec{2})
figureTitle = sprintf('pRFplot-%s-%s_ts',prf_all_bb.channels.subject_name(el(ii)),prf_all_bb.channels.name{el(ii)});

opt = [];
opt.plot.fontSize   = FntSiz;
opt.plot.addWangToTitle   = 'no';
opt.plot.addBensonToTitle = 'no';
opt.plot.addSbjToTitle    = 'yes';
opt.plot.addR2ToTitle     = 'yes';
opt.plot.model_options    = {'Color',cntcol};

if contains(tarBAND,'-noflip')
    opt.plot.reverseYaxis = 'no';
end

% boundary = [0 28 40 56 84 96 112 140 152 168 196 208 224]+.5;     % for all event change
boundary = [0 40 56 96 112 152 168 208 224]+.5;         % for BLANK

%-- plot w/o polynomial projection
bndfactor = 1/224*75;
opt.plot.boundaries = boundary.*bndfactor;
opt.plot.options = {'XTick', boundary.*bndfactor,'XTickLabel',{}};
opt.plot.XLim    = [floor(boundary(1)) ceil(boundary(end))].*bndfactor;
    
opt.skipprojection = 'yes';
p = ecog_plotGridPRFts(mdlall.datats,mdlall.stimulus,prfall,el(ii), opt);
hF = [hF p{:}];
  set(gcf,'MenuBar','none','Position',get(gcf,'Position').*[1,1,2.1,1.6]);
  R2val = get(gca,'Title'); R2val = R2val.String(end-6:end-2);
  title(sprintf('\\fontsize{22}\\color[rgb]{%g %g %g}%s (%s)',cntcol,tname,R2val));
  adjustplot(gcf);
  
switch tarBAND
    case {'broadband'}
        set(gca,'YTick',0:mean(get(gca,'YTick')):max(get(gca,'YTick'))*1.5);
        set(gca,'YTickLabel',num2str(get(gca,'YTick')'./100+1));
    case {'alpha','alpha-noflip'}
%         set(gca,'YTick',log10(round(10.^get(gca,'YTick'),1,'significant')));
        if min(yticks)>-1.2
            set(gca,'YTick',log10([1/8 1/4 1/2 1 2 4 8]));
            set(gca,'YTickLabel',{'1/8','1/4','1/2','1','2','4','8'});
%             set(gcf,'Position',get(gcf,'Position').*[1,1,1.025,1]);
        else
%             set(gca,'YTick',log10([0.03 0.1 0.3 1 3 ]);
%             set(gca,'YTickLabel',num2str(10.^get(gca,'YTick')'));
            set(gca,'YTick',log10([1/27 1/9 1/3 1 3 9 27]));
            set(gca,'YTickLabel',{'1/27','1/9','1/3','1','3','9','27'});
        end
end
    
    %-- add gray-box for blank
    pb=plotbox(gca,mat2cell(boundary(2:end).*bndfactor,1,2*ones(1,4)),ylim,'gray','FaceAlpha',0.5);
    uistack(pb,'bottom')
    
      figureName = sprintf('%s_%s-%s',figureTitle,'orig-noproj-showblank-ratio',tarBAND);
      if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
      
        %-- show boundaries
        boundaries = ([28 84 140 196]+0.5).*bndfactor;
        hl = straightline(boundaries,'v','k--');
        uistack(hl,'bottom')

          figureName = sprintf('%s_%s-%s',figureTitle,'orig-noproj-showbound-ratio',tarBAND);
          if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
    
end
end


%% plot averaged pRF time-series (plot broadband and alpha in one panel) with large font
if issaveplot

el = findtable(model_all_bb.channels,repelec{:});

alignzero  = true;
ratioaxis  = true;
maskblank  = true;
showR2 = true;
trnsC = 0.7;

ndat = 224;
    t_dec = decimate(1:ndat,smoothingN);
    boundaries = ([0 28 40 56 84 96 112 140 152 168 196 208 224]+0.5);
    blankbnd   = boundaries([3,4;6,7;9,10;12,13]);
    medblank   = mean(blankbnd,2);

for ii=1:numel(repelec{2})
%% plot
figureTitle = sprintf('pRFplot-%s-%s_ts',prf_all_bb.channels.subject_name(el(ii)),prf_all_bb.channels.name{el(ii)});

%%% compute time series of model & data
opt = [];
opt.skipprojection = 'yes';
opt.catchan        = 'yes';
opt.whichelec      = el(ii);
opt.flipgain       = 'no';
[modelts_bb, datts_bb] = reconPRFdog(prf_all_bb,model_all_bb,opt);
opt.flipgain       = 'yes';
[modelts_a, datts_a]   = reconPRFdog(prf_all_a,model_all_a,opt);

%%% plot time series of model & data
datts = cat(1, datts_bb',datts_a');
mdlts = cat(1, modelts_bb',modelts_a');

hF(end+1) = figure('Menubar','none'); hF(end).Position(3)=720; hF(end).Position(4)=550;
yyaxis left
if alignzero,   plot(zeros([1 ndat]),':','LineWidth',3.0,'Color','k'); hold on;    end
h = plot(t_dec,datts(1,:),':o','LineWidth',2.6,'Color',plcol(1,:));  hold on;
h.Color = [h.Color trnsC];
h1 = plot(t_dec,mdlts(1,:),'-','LineWidth',2.6);
if ~alignzero,  plot(zeros([1 ndat]),':','LineWidth',2.0,'Color',h1.Color);    end
yy1 = ylim;
yyaxis right
h = plot(t_dec,datts(2,:),':o','LineWidth',2.6,'Color',plcol(2,:)); hold on;
h.Color = [h.Color trnsC];
h2 = plot(t_dec,mdlts(2,:),'-','LineWidth',2.6);
if ~alignzero,  plot(zeros([1 ndat]),':','LineWidth',2.0,'Color',h2.Color);    end
yy2 = ylim;

hA = get(gca,'YAxis');
if alignzero
    hA(1).Limits = [-1 1] .* round(max(datts(1,:))*1.2);
    hA(2).Limits = [-1 1] .* round(max(-datts(2,:))*1.1,2);
else
    hA(1).Limits = yy1+[-0 1.5].*0.05.*diff(yy1);
    hA(2).Limits = yy2+[-0 1.5].*0.05.*diff(yy2);
end
if ratioaxis
    hA(1).TickValues(hA(1).TickValues<0) = [];
    hA(1).TickLabel = arrayfun(@(x) sprintf('%g',x./100+1),hA(1).TickValues,'UniformOutput',false);
    if max(abs(hA(2).TickValues)) > 1
%         hA(2).TickValues = log10([0.03 0.1 0.3 1 3 10 30]);
%         hA(2).TickLabel = arrayfun(@(x) sprintf('%g',10.^x),hA(2).TickValues,'UniformOutput',false);
        hA(2).TickValues = log10([1/27 1/9 1/3 1 3 9 27]);
        hA(2).TickLabel = {'1/27','1/9','1/3','1','3','9','27'};
    else
        hA(2).TickValues = log10([1/8 1/4 1/2 1 2 4 8]);
        hA(2).TickLabel = {'1/8','1/4','1/2','1','2','4','8'};
    end
    set(gca,'Position',get(gca,'Position')+[-.02 0 0 0])
end

%-- parameters
% lng = legend([h1 h2],...
%             {sprintf('Broadband (%.1f%%)',prf_all_bb.xval(el(ii))),...
%              sprintf('Alpha (%.1f%%)',prf_all_a.xval(el(ii)))},...
%              'AutoUpdate','off');
% lng.NumColumns = 4;
set(gca,'FontSize',LFntSiz);
title(sprintf('%s - %s',prf_all_bb.channels.subject_name(el(ii)),prf_all_bb.channels.name{el(ii)}));
xlim([1 ndat]);
yyaxis left
ylabel('Power Ratio (Broadband)');
yyaxis right
A=ylabel('Power Ratio (Alpha)','Rotation',270,'VerticalAlignment','bottom');
set(gca,'Position',get(gca,'Position')+[-1 0 1 0].*3e-5.*hF(end).Position(3))

%-- Boundary
yyaxis left;
hl = straightline(boundaries,'v','k--');
set(hl,'LineWidth',1);
hlist = get(get(hl(1),'Parent'),'Children');
hlist = hlist([(length(hl)+1):end,1:length(hl)]);
yyaxis right;

%-- BLANK
if maskblank
  yyaxis left;
  xticks([]);
  pb=plotbox(gca,mat2cell(blankbnd,ones(1,4),2),ylim,'gray','FaceAlpha',0.5);
  pb(1).Parent.Children = cat(1,pb(1).Parent.Children((length(pb)+1):end),pb');
  yyaxis right;
else
  xticks(medblank);
  xticklabels(repmat({'BLANK'},4,1));
end

%-- R2
if showR2
  yyaxis right;
  text(max(xlim)-1600/hF(end).Position(3)*LFntSiz,max(ylim)*0.88,sprintf('R^2 = %.1f%%',prf_all_bb.xval(el(ii))),...
      'FontSize',LFntSiz,'FontWeight','bold','Color',plcol(1,:),'BackgroundColor','w','EdgeColor','k','LineWidth',1.4);
  text(max(xlim)-1600/hF(end).Position(3)*LFntSiz,max(ylim)*-0.88,sprintf('R^2 = %.1f%%',prf_all_a.xval(el(ii))),...
      'FontSize',LFntSiz,'FontWeight','bold','Color',plcol(2,:),'BackgroundColor','w','EdgeColor','k','LineWidth',1.4);
end
    if ratioaxis, postfix = '-ratio';   else, postfix = '';     end
    figureName = sprintf('%s_%s%s-%s',figureTitle,'orig-noproj-showblank-wide',postfix,'Both');
    set(gcf,'Name',figureName);
    if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end

end

end
