% Compare decimation & full-ts, Channel Selection

% 20201008 Yuasa
% 20220223 Yuasa: update bootstrap method
% 20240524 Yuasa: for individuals

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

plcol  = get(groot,'defaultAxesColorOrder');
plcol(1:2,:) = [49 127 160; 183 20 51]./255;

hF = gobjects(0);

%% Channel selection with FPM

threshold = threshold_bb;
nboot     = 5000;
%-- compute xR2_a in ROIs (wang_prob)
filename = sprintf('all_chansel-wangprob%s-%s%s_thresh%02d-ind',R2mode,selectchs,postfix,threshold);
filepath = fullfile(SetDefaultAnalysisPath('DAT',prfstatPth),filename);

if exist([filepath '.mat'],'file')
 fprintf('loading %s...',filename);
 load(filepath);
 fprintf('\n');
else
  error('Need to run ecog_APRFF_02e2_selectedchannels_ind in advance.');
end

nsbj = length(subjects);
nroi = length(rois);

%-- Set alpha
alpha  = 0.32;   % 0.05 for 2sd, 0.32 for 1sd
talpha = 0.05/nroi/2;

%-- update ROI labels
rois(ismember(rois,'low'))={'V1-V3'};
rois(ismember(rois,'high'))={'Dorsolateral'};
rois(ismember(rois,'dorsolateral'))={'Dorsolateral'};

%% %%%%%%%%%%%%%%%%
%% distribution
%% %%%%%%%%%%%%%%%%
%% Channel selection with FPM in low/high visual area in one panel -short-

%-- reject inaccurate ROIs
shiftframe = 0.30;
shiftbar   = 0.035;

selroi = [12 13];       figwidth = 100; showelecnum = ' n=%.1f';
nsel = length(subjects);

hF(end+1) = figure('defaultAxesColorOrder',plcol);
ht = tiledlayout(length(selroi),1,'TileSpacing','compact','Padding','compact');
set(gcf,'Position',get(gcf,'Position').*[.3 .5 0 0]+[0 0 nsel*figwidth 420*length(selroi)]);
for iroi=selroi
nexttile;

%-- figure out
h3 = errorbar((1:nsel)-shiftbar,cellfun(@nanmean,mn_keep_roi(:,iroi)),...
    cellfun(@(C) nanmean(C)-prctile(C,alpha/2*100),mn_keep_roi(:,iroi)),...
    cellfun(@(C) prctile(C,(1-alpha/2)*100)-nanmean(C),mn_keep_roi(:,iroi)),...
    'o','LineStyle','none',...
    'LineWidth',2,'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
hold on, h4 = errorbar((1:nsel)+shiftbar,cellfun(@nanmean,mn_discard_roi(:,iroi)),...
    cellfun(@(C) nanmean(C)-prctile(C,alpha/2*100),mn_discard_roi(:,iroi)),...
    cellfun(@(C) prctile(C,(1-alpha/2)*100)-nanmean(C),mn_discard_roi(:,iroi)),...
    'o','LineStyle','none',...
    'LineWidth',2,'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
ylim([-35 75]);
xlim([0.5,nsel+0.5]+shiftframe);  plot(xlim,[0,0],'k--');
xticks(1:nsel); xticklabels(subjects); xtickangle(0);
set(gca,'FontSize',FntSiz);
 
%     ylabel('Alpha cross-validated R^2 (%)');
    ylabel('Variance Explained in Alpha (%)');
    hl=legend({'Visually Selective Electrodes','Non-Visually Selective Electrodes'});
    
set(gcf,'Unit','pixels')

for ii=1:nsel
    mn_dif = nanmean(mn_keep_roi{ii,iroi}) - nanmean(mn_discard_roi{ii,iroi});
    if abs(mn_dif)<4,   mn_shift = 2*sign(mn_dif);
    else,               mn_shift = 0;  end
    text(ii+0.05-shiftbar, nanmean(mn_keep_roi{ii,iroi})+mn_shift, sprintf(showelecnum,nanmean(n_keep_roi{ii,iroi})),...
        'FontSize',FntSiz,'Color',plcol(1,:));
    text(ii+0.05+shiftbar, nanmean(mn_discard_roi{ii,iroi})-mn_shift, sprintf(showelecnum,nanmean(n_discard_roi{ii,iroi})),...
        'FontSize',FntSiz,'Color',plcol(2,:));
end

fprintf('%s\n',sprintf('%5s',rois{selroi}));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_keep_roi(selroi))));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_discard_roi(selroi))));
fprintf('%s\n',sprintf('%5d',cellfun(@(C) prctile(C,talpha*100),mn_keep_roi(selroi)) - cellfun(@(C) prctile(C,(1-talpha)*100),mn_discard_roi(selroi))>0));

title(gca,rois(iroi),'FontSize',FntSiz);
end

htt = title(ht,'Cross-validated Variance Explained in Alpha','FontSize',FntSiz,'FontWeight','bold');
figureName = sprintf('channelselection%s-%02d%%-%s',R2mode,threshold,'ind');
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
