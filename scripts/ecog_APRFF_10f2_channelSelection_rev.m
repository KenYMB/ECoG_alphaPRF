% Compare decimation & full-ts, Channel Selection

% 20210408 - update to demonstrate reverse selection

%%
% close all; clear all;
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
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'channelselection');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%-- Plotting Setting
FntSiz  = 18;

%% load analyzePRF & recompute full-ts R2 & set data for bootstrap with half-trials
clear alphaType broadbandType

average        ='runs';
prfmodel       ='linear';
smoothingMode  ='decimate';
smoothingN     = 3;
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

%% %%%%%%%%%%%%%%
%% Visualize
%% %%%%%%%%%%%%%%
%%% prepare visualize
res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');
plcol(1:2,:) = [49 127 160; 183 20 51]./255;

hF = gobjects(0);

%% Channel selection with FPM

threshold = threshold_a;
nboot     = 5000;
%-- compute xR2_a in ROIs (wang_prob)
filename = sprintf('all_chansel-wangprob%s-%s%s_thresh%02d-rev',R2mode,selectchs,postfix,threshold);
filepath = fullfile(SetDefaultAnalysisPath('DAT',prfstatPth),filename);

if exist([filepath '.mat'],'file')
 fprintf('loading %s...',filename);
 load(filepath);
 fprintf('\n');
else
  error('Need to run ecog_APRFF_02e_selectedchannels in advance.');
end

%-- reject inaccurate ROIs
nroi = length(rois);
for iroi=1:nroi
    if sum(isnan(mn_keep_roi{iroi})) > nboot*0.85
        mn_keep_roi{iroi}(:) = nan;
        n_keep_roi{iroi}(:) = nan;
    else
        %-- set NaN for invalid condition for mean
        n_keep_roi{iroi}(n_keep_roi{iroi}==0) = nan;
    end
    if sum(isnan(mn_discard_roi{iroi})) > nboot*0.85
        mn_discard_roi{iroi}(:) = nan;
        n_discard_roi{iroi}(:) = nan;
    else
        %-- set NaN for invalid condition for mean
        n_discard_roi{iroi}(n_discard_roi{iroi}==0) = nan;
    end
end

%-- Set alpha
alpha  = 0.32;   % 0.05 for 2sd, 0.32 for 1sd
talpha = 0.05/nroi/2;

%-- update ROI labels
rois(ismember(rois,'low'))={'V1-V3'};
rois(ismember(rois,'high'))={'Dorsolateral'};
rois(ismember(rois,'dorsolateral'))={'Dorsolateral'};

%% Channel selection with FPM in low/high visual area

%-- reject inaccurate ROIs
shiftframe = 0.2;
shiftbar   = 0.02;

for iroi = 1:2
switch iroi
    case 1
        selroi = [1 2 3 12];
        roiset = 'low';
    case 2
        selroi = [6:11 13];
        roiset = 'high';
end
nsel = length(selroi);

%-- figure out
hF(end+1) = figure('defaultAxesColorOrder',plcol);
h3 = errorbar((1:nsel)-shiftbar,cellfun(@nanmean,mn_keep_roi(selroi)),...
    cellfun(@(C) nanmean(C)-prctile(C,alpha/2*100),mn_keep_roi(selroi)),...
    cellfun(@(C) prctile(C,(1-alpha/2)*100)-nanmean(C),mn_keep_roi(selroi)),...
    'o','LineStyle','none',...
    'LineWidth',2,'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
hold on, h4 = errorbar((1:nsel)+shiftbar,cellfun(@nanmean,mn_discard_roi(selroi)),...
    cellfun(@(C) nanmean(C)-prctile(C,alpha/2*100),mn_discard_roi(selroi)),...
    cellfun(@(C) prctile(C,(1-alpha/2)*100)-nanmean(C),mn_discard_roi(selroi)),...
    'o','LineStyle','none',...
    'LineWidth',2,'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
xlim([0.5,nsel+0.5]+shiftframe);  plot(xlim,[0,0],'k--');  plot([1 1].*(nsel-0.5)+shiftframe,ylim,'k:');
xticks(1:nsel); xticklabels(rois(selroi)); xtickangle(0);
set(gca,'FontSize',FntSiz);
ylabel('Variance Explained in Broadband (%)');
hl=legend({sprintf('Variance explained in Alpha > %0d%%',threshold_a),...
           sprintf('Variance explained in Alpha ≤ %0d%%',threshold_a)});
 
title('Cross-validated Variance Explained in Broadband');
ylim([-20 118]);
set(gcf,'Unit','pixels')
set(gcf,'Position',get(gcf,'Position').*[.5 .5 0 0]+[0 0 nsel*120 420]);

for ii=1:nsel
    mn_dif = nanmean(mn_keep_roi{selroi(ii)}) - nanmean(mn_discard_roi{selroi(ii)});
    if abs(mn_dif)<4,   mn_shift = 2*sign(mn_dif);
    else,               mn_shift = 0;  end
    text(ii+0.05-shiftbar, nanmean(mn_keep_roi{selroi(ii)})+mn_shift, sprintf('%5.1f',nanmean(n_keep_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(1,:));
    text(ii+0.05+shiftbar, nanmean(mn_discard_roi{selroi(ii)})-mn_shift, sprintf('%5.1f',nanmean(n_discard_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(2,:));
end

fprintf('%s\n',sprintf('%5s',rois{selroi}));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_keep_roi(selroi))));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_discard_roi(selroi))));
fprintf('%s\n',sprintf('%5d',cellfun(@(C) prctile(C,talpha*100),mn_keep_roi(selroi)) - cellfun(@(C) prctile(C,(1-talpha)*100),mn_discard_roi(selroi))>0));

figureName = sprintf('channelselection%s-rev-%02d%%-%s',R2mode,threshold,roiset);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end

%% Channel selection with FPM in low/high visual area in one panel

%-- reject inaccurate ROIs
shiftframe = 0.15;
shiftbar   = 0.02;
tnsel      = 11;

hF(end+1) = figure('defaultAxesColorOrder',plcol);
ht = tiledlayout(1,tnsel,'TileSpacing','compact','Padding','compact');

set(gcf,'Position',get(gcf,'Position').*[.3 .5 0 0]+[0 0 tnsel*105 420]);
for iroi = 1:2
switch iroi
    case 1
        selroi = [1 2 3 12];
        roiset = 'low';
    case 2
        selroi = [6:11 13];
        roiset = 'high';
end
nsel = length(selroi);
nexttile([1, nsel]);

%-- figure out
h3 = errorbar((1:nsel)-shiftbar,cellfun(@nanmean,mn_keep_roi(selroi)),...
    cellfun(@(C) nanmean(C)-prctile(C,alpha/2*100),mn_keep_roi(selroi)),...
    cellfun(@(C) prctile(C,(1-alpha/2)*100)-nanmean(C),mn_keep_roi(selroi)),...
    'o','LineStyle','none',...
    'LineWidth',2,'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
hold on, h4 = errorbar((1:nsel)+shiftbar,cellfun(@nanmean,mn_discard_roi(selroi)),...
    cellfun(@(C) nanmean(C)-prctile(C,alpha/2*100),mn_discard_roi(selroi)),...
    cellfun(@(C) prctile(C,(1-alpha/2)*100)-nanmean(C),mn_discard_roi(selroi)),...
    'o','LineStyle','none',...
    'LineWidth',2,'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
ylim([-20 118]);
xlim([0.5,nsel+0.5]+shiftframe);  plot(xlim,[0,0],'k--');
plot([1 1].*(nsel-0.5)+shiftframe,ylim,'k:','LineWidth',1.2);
xticks(1:nsel); xticklabels(rois(selroi)); xtickangle(0);
set(gca,'FontSize',FntSiz);
 
if iroi == 1
    ylabel('Variance Explained in Broadband (%)');
else
    yticklabels({});
end
if iroi == 2
    hl=legend({sprintf('Variance explained in Alpha > %0d%%',threshold_a),...
        sprintf('Variance explained in Alpha ≤ %0d%%',threshold_a)});
end
set(gcf,'Unit','pixels')

for ii=1:nsel
    mn_dif = nanmean(mn_keep_roi{selroi(ii)}) - nanmean(mn_discard_roi{selroi(ii)});
    if abs(mn_dif)<4,   mn_shift = 2*sign(mn_dif);
    else,               mn_shift = 0;  end
    text(ii+0.05-shiftbar, nanmean(mn_keep_roi{selroi(ii)})+mn_shift, sprintf('%5.1f',nanmean(n_keep_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(1,:));
    text(ii+0.05+shiftbar, nanmean(mn_discard_roi{selroi(ii)})-mn_shift, sprintf('%5.1f',nanmean(n_discard_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(2,:));
end

fprintf('%s\n',sprintf('%5s',rois{selroi}));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_keep_roi(selroi))));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_discard_roi(selroi))));
fprintf('%s\n',sprintf('%5d',cellfun(@(C) prctile(C,talpha*100),mn_keep_roi(selroi)) - cellfun(@(C) prctile(C,(1-talpha)*100),mn_discard_roi(selroi))>0));

end
htt = title(ht,'Cross-validated Variance Explained in Broadband','FontSize',FntSiz,'FontWeight','bold');
figureName = sprintf('channelselection%s-rev-%02d%%-%s',R2mode,threshold,'merge');
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end


%% Channel selection with FPM in low/high visual area in one panel -short-

%-- reject inaccurate ROIs
shiftframe = 0.15;
shiftbar   = 0.02;

hF(end+1) = figure('defaultAxesColorOrder',plcol);
ht = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');

selroi = [1 2 3 12 13];
nsel = length(selroi);
nexttile;

set(gcf,'Position',get(gcf,'Position').*[.3 .5 0 0]+[0 0 nsel*105 420]);

%-- figure out
h3 = errorbar((1:nsel)-shiftbar,cellfun(@nanmean,mn_keep_roi(selroi)),...
    cellfun(@(C) nanmean(C)-prctile(C,alpha/2*100),mn_keep_roi(selroi)),...
    cellfun(@(C) prctile(C,(1-alpha/2)*100)-nanmean(C),mn_keep_roi(selroi)),...
    'o','LineStyle','none',...
    'LineWidth',2,'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
hold on, h4 = errorbar((1:nsel)+shiftbar,cellfun(@nanmean,mn_discard_roi(selroi)),...
    cellfun(@(C) nanmean(C)-prctile(C,alpha/2*100),mn_discard_roi(selroi)),...
    cellfun(@(C) prctile(C,(1-alpha/2)*100)-nanmean(C),mn_discard_roi(selroi)),...
    'o','LineStyle','none',...
    'LineWidth',2,'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
ylim([-20 118]);
xlim([0.5,nsel+0.5]+shiftframe);  plot(xlim,[0,0],'k--');
plot([1 1].*(nsel-1.5)+shiftframe,ylim,'k:','LineWidth',1.5);
xticks(1:nsel); xticklabels(rois(selroi)); xtickangle(0);
set(gca,'FontSize',FntSiz);
 
    ylabel('Variance Explained in Broadband (%)');
    hl=legend({sprintf('Variance explained in Alpha > %0d%%',threshold_a),...
               sprintf('Variance explained in Alpha ≤ %0d%%',threshold_a)});
    
set(gcf,'Unit','pixels')

for ii=1:nsel
    mn_dif = nanmean(mn_keep_roi{selroi(ii)}) - nanmean(mn_discard_roi{selroi(ii)});
    if abs(mn_dif)<4,   mn_shift = 2*sign(mn_dif);
    else,               mn_shift = 0;  end
    text(ii+0.05-shiftbar, nanmean(mn_keep_roi{selroi(ii)})+mn_shift, sprintf('%5.1f',nanmean(n_keep_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(1,:));
    text(ii+0.05+shiftbar, nanmean(mn_discard_roi{selroi(ii)})-mn_shift, sprintf('%5.1f',nanmean(n_discard_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(2,:));
end

fprintf('%s\n',sprintf('%5s',rois{selroi}));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_keep_roi(selroi))));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_discard_roi(selroi))));
fprintf('%s\n',sprintf('%5d',cellfun(@(C) prctile(C,talpha*100),mn_keep_roi(selroi)) - cellfun(@(C) prctile(C,(1-talpha)*100),mn_discard_roi(selroi))>0));


htt = title(ht,'Cross-validated Variance Explained in Broadband','FontSize',FntSiz,'FontWeight','bold');
figureName = sprintf('channelselection%s-rev-%02d%%-%s',R2mode,threshold,'short');
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot xR2 relations in 2D histogram
sbfld = {'R2'};
plotwide = true;
% plotwide = false;

hF(end+1) = figure('Position',[1000/2 900/2 500 440],'MenuBar','none');

    %-- set data
    dat1 = prf_all_bb.xval;
    dat2 = prf_all_a.xval;
    thresh1 = threshold_bb;
    thresh2 = threshold_a;

    %-- plot histogram & set axis
    h1 = histogram2(dat1,dat2,'DisplayStyle','tile','ShowEmptyBins','on');

        if plotwide
         R2lim = [-150 100];
         namewide = '-wide';
        else
         R2lim = [-100 100];
         namewide = '';
        end
        xlim(R2lim);
        ylim(R2lim);
        Nbin = 25;
        
    h1.XBinLimits = xlim;
    h1.YBinLimits = ylim;
    h1.NumBins = [Nbin Nbin];
    h1.EdgeColor = 'none';
    view([0 90]); axis square;
    xlabel('Broadband cross-validated R^2 (%)'); ylabel('Alpha cross-validated R^2 (%)');
    title(sprintf('Distribution of Cross-Validated\nVariance Explained'));
    set(gca,'FontSize',FntSiz)
    
    %-- modify BinLimits separate bins at threshold
    mn_shift1 = mod((thresh1 - xlim*[1;0]),h1.BinWidth(1));
    mn_shift2 = mod((thresh2 - ylim*[1;0]),h1.BinWidth(2));
    
        h1.XBinLimits = xlim + mn_shift1 - [h1.BinWidth(1) 0];
        h1.YBinLimits = ylim + mn_shift2 - [h1.BinWidth(2) 0];
    
    hold on
    %-- show mask
    surf([xlim*[1;0] thresh1],[ylim*[1;0] thresh2],zeros(2),ones(2,2,3)*0.6,'FaceAlpha',0.5,'EdgeColor','none')
    surf([thresh1 xlim*[0;1]],[ylim*[1;0] thresh2],zeros(2),ones(2,2,3)*0.6,'FaceAlpha',0.5,'EdgeColor','none')
    surf([xlim*[1;0] thresh1],[thresh2 ylim*[0;1]],zeros(2),ones(2,2,3)*0.6,'FaceAlpha',0.5,'EdgeColor','none')
    
    %-- plot threshold
    plot([1 1]*thresh1,ylim,'-','LineWidth',1.2,'Color',[1 1 1]*0.6);
    plot(xlim,[1 1]*thresh2,'-','LineWidth',1.2,'Color',[1 1 1]*0.6);
    xticks(unique([yticks thresh1]));
    yticks(unique([yticks thresh2]));
    
    %-- plot zero-point
%     plot([0 0],ylim,'w-','LineWidth',1.0);
%     plot(xlim,[0 0],'w-','LineWidth',1.0);
    
    %-- set diagonal
    axlim = [min([xlim,ylim]), max([xlim,ylim])];
%     plot(axlim,axlim,'w--','LineWidth',0.8);

    %-- color set
    colormap([[0,0,0];hot].^0.3);
    hc=colorbar;
    hc.Ticks = hc.Limits;
    hc.TickLabels = {'0'; 'Max'};

set(gcf,'Name',sbfld{:});

figureName = sprintf('prf-%02d%%-%02d%%-ecc%02d_2Dhist_%s%s',threshold_bb,threshold_a,eclimit,sbfld{:},namewide);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end


%% plot xR2 relations in scatter plot
sbfld = {'R2'};
plotwide = true;
% plotwide = false;

hF(end+1) = figure('Position',[1000/2 900/2 500 440],'MenuBar','none','defaultAxesColorOrder',plcol);

    %-- set data
    dat1 = prf_all_bb.xval;
    dat2 = prf_all_a.xval;
    thresh1 = threshold_bb;
    thresh2 = threshold_a;
    R21 = prf_all_bb.xval;
    R22 = prf_all_a.xval;

    %-- plot histogram & set axis
    h1 = scatter(dat1(R21>thresh1 & R22>thresh2),dat2(R21>thresh1 & R22>thresh2));
    hold on
    h2 = scatter(dat1(~(R21>thresh1 & R22>thresh2)),dat2(~(R21>thresh1 & R22>thresh2)),'MarkerEdgeColor','k');

        if plotwide
         R2lim = [-150 100];
         namewide = '-wide';
        else
         R2lim = [-100 100];
         namewide = '';
        end
        xlim(R2lim);
        ylim(R2lim);
        
    view([0 90]); axis square;
    xlabel('Broadband cross-validated R^2 (%)'); ylabel('Alpha cross-validated R^2 (%)');
    title(sprintf('Distribution of Cross-Validated\nVariance Explained'));
    set(gca,'FontSize',FntSiz)
    
    %-- plot threshold
    plot([1 1]*thresh1,ylim,'-','LineWidth',1.2,'Color',[1 1 1]*0.6);
    plot(xlim,[1 1]*thresh2,'-','LineWidth',1.2,'Color',[1 1 1]*0.6);
    xticks(unique([yticks thresh1]));
    yticks(unique([yticks thresh2]));
    
    %-- plot zero-point
%     plot([0 0],ylim,'w-','LineWidth',1.0);
%     plot(xlim,[0 0],'w-','LineWidth',1.0);
    
    %-- set diagonal
    axlim = [min([xlim,ylim]), max([xlim,ylim])];
%     plot(axlim,axlim,'w--','LineWidth',0.8);

set(gcf,'Name',sbfld{:});
box on

figureName = sprintf('prf-%02d%%-%02d%%-ecc%02d_scatter_%s%s',threshold_bb,threshold_a,eclimit,sbfld{:},namewide);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end

