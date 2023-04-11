% Compare decimation & full-ts, Channel Selection

% 20201008 Yuasa
% 20220223 Yuasa: update bootstrap method

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
filename = sprintf('all_chansel-wangprob%s-%s%s_thresh%02d',R2mode,selectchs,postfix,threshold);
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

%% %%%%%%%%%%%%%%%%
%% distribution
%% %%%%%%%%%%%%%%%%
%% Channel selection with FPM in low/high visual area

if issaveplot
%-- reject inaccurate ROIs
shiftframe = 0.15;
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
showelecnum = '%5.1f';
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
ylim([-20 60]);
xlim([0.5,nsel+0.5]+shiftframe);  plot(xlim,[0,0],'k--');
plot([1 1].*(nsel-0.5)+shiftframe,ylim,'k:','LineWidth',1.2);
xticks(1:nsel); xticklabels(rois(selroi)); xtickangle(0);
set(gca,'FontSize',FntSiz);
% ylabel('Alpha cross-validated R^2 (%)');
ylabel('Variance Explained in Alpha (%)');
% hl=legend({'Broadband cross-validated R^2 > 31%','Broadband cross-validated R^2 ≤ 31%'});
hl=legend({'Visually Selective Electrodes','Non-Visually Selective Electrodes'});
 
title('Cross-validated Variance Explained in Alpha');
set(gcf,'Unit','pixels')
set(gcf,'Position',get(gcf,'Position').*[.5 .5 0 0]+[0 0 nsel*120 420]);

for ii=1:nsel
    mn_dif = nanmean(mn_keep_roi{selroi(ii)}) - nanmean(mn_discard_roi{selroi(ii)});
    if abs(mn_dif)<4,   mn_shift = 2*sign(mn_dif);
    else,               mn_shift = 0;  end
    text(ii+0.05-shiftbar, nanmean(mn_keep_roi{selroi(ii)})+mn_shift, sprintf(showelecnum,nanmean(n_keep_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(1,:));
    text(ii+0.05+shiftbar, nanmean(mn_discard_roi{selroi(ii)})-mn_shift, sprintf(showelecnum,nanmean(n_discard_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(2,:));
end

fprintf('%s\n',sprintf('%5s',rois{selroi}));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_keep_roi(selroi))));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_discard_roi(selroi))));
fprintf('%s\n',sprintf('%5d',cellfun(@(C) prctile(C,talpha*100),mn_keep_roi(selroi)) - cellfun(@(C) prctile(C,(1-talpha)*100),mn_discard_roi(selroi))>0));

figureName = sprintf('channelselection%s-%02d%%-%s',R2mode,threshold,roiset);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end
end

%% Channel selection with FPM in low/high visual area in one panel

if issaveplot
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
showelecnum = '%5.1f';
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
ylim([-20 60]);
xlim([0.5,nsel+0.5]+shiftframe);  plot(xlim,[0,0],'k--');
plot([1 1].*(nsel-0.5)+shiftframe,ylim,'k:','LineWidth',1.2);
xticks(1:nsel); xticklabels(rois(selroi)); xtickangle(0);
set(gca,'FontSize',FntSiz);
 
if iroi == 1
%     ylabel('Alpha cross-validated R^2 (%)');
    ylabel('Variance Explained in Alpha (%)');
else
    yticklabels({});
end
if iroi == 2
%     hl=legend({'Broadband cross-validated R^2 > 31%','Broadband cross-validated R^2 ≤ 31%'});
    hl=legend({'Visually Selective Electrodes','Non-Visually Selective Electrodes'});
end
set(gcf,'Unit','pixels')

for ii=1:nsel
    mn_dif = nanmean(mn_keep_roi{selroi(ii)}) - nanmean(mn_discard_roi{selroi(ii)});
    if abs(mn_dif)<4,   mn_shift = 2*sign(mn_dif);
    else,               mn_shift = 0;  end
    text(ii+0.05-shiftbar, nanmean(mn_keep_roi{selroi(ii)})+mn_shift, sprintf(showelecnum,nanmean(n_keep_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(1,:));
    text(ii+0.05+shiftbar, nanmean(mn_discard_roi{selroi(ii)})-mn_shift, sprintf(showelecnum,nanmean(n_discard_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(2,:));
end

fprintf('%s\n',sprintf('%5s',rois{selroi}));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_keep_roi(selroi))));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_discard_roi(selroi))));
fprintf('%s\n',sprintf('%5d',cellfun(@(C) prctile(C,talpha*100),mn_keep_roi(selroi)) - cellfun(@(C) prctile(C,(1-talpha)*100),mn_discard_roi(selroi))>0));

end
htt = title(ht,'Cross-validated Variance Explained in Alpha','FontSize',FntSiz,'FontWeight','bold');
figureName = sprintf('channelselection%s-%02d%%-%s',R2mode,threshold,'merge');
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end


%% Channel selection with FPM in low/high visual area in one panel -short-

%-- reject inaccurate ROIs
shiftframe = 0.15;
shiftbar   = 0.02;

hF(end+1) = figure('defaultAxesColorOrder',plcol);
ht = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');

if issaveplot,      selroi = [1 2 3 12 13]; figwidth = 105; showelecnum = '%5.1f';
else,               selroi = [12 13];       figwidth = 200; showelecnum = ' n=%.0f';
end
nsel = length(selroi);
nexttile;

set(gcf,'Position',get(gcf,'Position').*[.3 .5 0 0]+[0 0 nsel*figwidth 420]);

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
ylim([-20 60]);
xlim([0.5,nsel+0.5]+shiftframe);  plot(xlim,[0,0],'k--');
plot([1 1].*(nsel-1.5)+shiftframe,ylim,'k:','LineWidth',1.5);
xticks(1:nsel); xticklabels(rois(selroi)); xtickangle(0);
set(gca,'FontSize',FntSiz);
 
%     ylabel('Alpha cross-validated R^2 (%)');
    ylabel('Variance Explained in Alpha (%)');
    hl=legend({'Visually Selective Electrodes','Non-Visually Selective Electrodes'});
    
set(gcf,'Unit','pixels')

for ii=1:nsel
    mn_dif = nanmean(mn_keep_roi{selroi(ii)}) - nanmean(mn_discard_roi{selroi(ii)});
    if abs(mn_dif)<4,   mn_shift = 2*sign(mn_dif);
    else,               mn_shift = 0;  end
    text(ii+0.05-shiftbar, nanmean(mn_keep_roi{selroi(ii)})+mn_shift, sprintf(showelecnum,nanmean(n_keep_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(1,:));
    text(ii+0.05+shiftbar, nanmean(mn_discard_roi{selroi(ii)})-mn_shift, sprintf(showelecnum,nanmean(n_discard_roi{selroi(ii)})),...
        'FontSize',FntSiz,'Color',plcol(2,:));
end

fprintf('%s\n',sprintf('%5s',rois{selroi}));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_keep_roi(selroi))));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_discard_roi(selroi))));
fprintf('%s\n',sprintf('%5d',cellfun(@(C) prctile(C,talpha*100),mn_keep_roi(selroi)) - cellfun(@(C) prctile(C,(1-talpha)*100),mn_discard_roi(selroi))>0));


htt = title(ht,'Cross-validated Variance Explained in Alpha','FontSize',FntSiz,'FontWeight','bold');
figureName = sprintf('channelselection%s-%02d%%-%s',R2mode,threshold,'short');
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end

%% %%%%%%%%%%%%%%%%
%% Violin
%% %%%%%%%%%%%%%%%%

% if issaveplot
% 
% %-- reject inaccurate ROIs
% 
% hF(end+1) = figure('defaultAxesColorOrder',plcol);
% ht = tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
% 
% selroi = [12 13];       figwidth = 200;
% 
% nsel = length(selroi);
% nexttile;
% 
% set(gcf,'Position',get(gcf,'Position').*[.3 .5 0 0]+[0 0 nsel*figwidth 420]);
% 
%     D = cat(3,cat(1,mn_keep_roi{selroi})',cat(1,mn_discard_roi{selroi})');
%     hV = violinplotsplit(D,rois(selroi));
%     ylim([-25 70]);
%     set(gca,'FontSize',FntSiz);
%     hl = plot(xlim,[0,0],'k--'); uistack(hl,'bottom');
%     legend([hV([1,3]).ViolinPlot],{'Visually Selective Electrodes','Non-Visually Selective Electrodes'});
%     ylabel('Variance Explained in Alpha (%)');
%     
% set(gcf,'Unit','pixels')
% 
% for ii=1:nsel
%     mn_dif = nanmean(mn_keep_roi{selroi(ii)}) - nanmean(mn_discard_roi{selroi(ii)});
%     if abs(mn_dif)<4,   mn_shift = 2*sign(mn_dif);
%     else,               mn_shift = 0;  end
%     text(ii+0.05, nanmean(mn_keep_roi{selroi(ii)})+mn_shift, sprintf('%5.1f',nanmean(n_keep_roi{selroi(ii)})),...
%         'FontSize',FntSiz,'Color',plcol(1,:));
%     text(ii-0.05, nanmean(mn_discard_roi{selroi(ii)})-mn_shift, sprintf('%5.1f',nanmean(n_discard_roi{selroi(ii)})),...
%         'FontSize',FntSiz,'Color',plcol(2,:),'HorizontalAlignment','right');
% end
% 
% figureName = sprintf('channelselection%s-%02d%%-%s',R2mode,threshold,'Violin');
% if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
% 
% end

%% %%%%%%%%%%%%%%%%
%% histogram
%% %%%%%%%%%%%%%%%%

%% permutation
filename = sprintf('all_cod-permhalves%s-%s%s',R2mode,selectchs,postfix);
filepath = fullfile(SetDefaultAnalysisPath('DAT',prfstatPth),filename);

if exist([filepath '.mat'],'file')
    clear Nresamp;
    load(filepath);
    nboot = size(bootchs,1);
    model_data_a   = cellfun(@(C) C(bootchs(:,1),:),model_all_a.datats,'UniformOutput',false);
    model_data_bb  = cellfun(@(C) C(bootchs(:,1),:),model_all_bb.datats,'UniformOutput',false);
    bootstimulus   = model_all_bb.stimulus(bootchs(:,1),:);
    params_perm_a  = prf_all_a.params(:,:,bootchs(:,2));
    params_perm_bb = prf_all_bb.params(:,:,bootchs(:,2));
else
    error('Need to run ecog_APRFF_02a_distribution_xR2_halves in advance.');
end

%-- reset prf_all
[~,prf_all_a]      = ecog_rearrangePRF(prf_params_a,va_area);
[~,prf_all_bb]     = ecog_rearrangePRF(prf_params_bb,va_area);

%%

Nresamp = 1:nboot;

isround = 1;
% Nresamp = 1:numel(cod_bb);
plotwide = true;
% plotwide = false;

%% histogram of permutation xR2 in broadband

barcol = [[1 1 1]*0.4; [1 1 1]*0.6];

% visch = find(~ismember(channels.(va_area),'none'));
visch = 1:height(channels);

cod     = cod_bb;
%-- plot histogram
hF(end+1) = figure('Position',[300 300 500 400],'MenuBar','none');
subplot_er(1,1,1);
alpha_bb = 0.05;
    thresh1 = prctile(cod(Nresamp),(1-alpha_bb)*100);
    if isround,     thresh1 = round(thresh1);     end
    barshift = mod(thresh1,10);
    medval  = median(cod(Nresamp));
%     histogram(cod_bb(Nresamp),'BinWidth',10,'FaceColor',[1 1 1]*0.6);
    hh1 = histogram(cod(Nresamp( cod(Nresamp)>=round(thresh1,0) )),'BinWidth',10,'FaceColor',barcol(1,:));
    hh1.BinLimits = hh1.BinLimits + barshift - [10 0];
    hold on;
    hh2 = histogram(cod(Nresamp( cod(Nresamp)<round(thresh1,0) )),'BinWidth',10,'FaceColor',barcol(2,:));
    hh2.BinLimits = hh2.BinLimits + barshift - [10 0];
    xlim([-150 100]);
    if plotwide
     xlim([-150 100]);
     namewide = '-wide';
    else
     xlim([-100 100]);
     namewide = '';
    end
    set(gca,'FontSize',18);
    plot(thresh1.*[1 1],ylim,'k--','LineWidth',1.2);
    plot(medval.*[1 1],ylim,'k-.','LineWidth',1.6);
    title(sprintf('Distribution of Shuffled\n Variance Explained in Broadband'));
    if thresh1>=70
        text(thresh1-5,diff(ylim)*0.8,sprintf('%.0f%%',(1-alpha_bb)*100),'FontSize',18,'HorizontalAlignment','right');
    else
        text(thresh1+5,diff(ylim)*0.8,sprintf('%.0f%%',(1-alpha_bb)*100),'FontSize',18);
    end
    xlabel('Broadband shuffled R^2 (%)');
    ylabel('Probability');
    xticks(unique([xticks thresh1]));
%     xticks(unique([xticks thresh1 round(medval,1)]));
    ytickbase = yticks; ytickbase = round(diff(ytickbase(1:2))/numel(Nresamp),2)*numel(Nresamp);
    yticks(0:ytickbase:(ytickbase*length(yticks)));
    yticklabels( num2str(yticks'/numel(Nresamp)) );
    
figureName = sprintf('prf-distribution-shuffle%s_permhalves%s_%.0f%%-N%d',namewide,R2mode,alpha_bb*100,length(Nresamp));
if isround, figureName = [figureName '-round'];   end
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end

fprintf('%s median is %.2f%%\n','Shuffled broadband',medval);

%% histogram of broadband xR2

barcol = [[1 1 1]*0.4; [1 1 1]*0.6];
% barcol = plcol;

% visch = find(~ismember(channels.(va_area),'none'));
visch = 1:height(channels);

cod     = cod_bb;
cod_raw = prf_all_bb.xval;
%-- plot histogram
hF(end+1) = figure('Position',[300 300 500 400],'MenuBar','none');
subplot_er(1,1,1);
alpha_bb = 0.05;
    thresh1  = prctile(cod(Nresamp),(1-alpha_bb)*100);
    if isround,     thresh1 = round(thresh1);     end
    barshift = mod(thresh1,10);
    medval  = median(cod_raw);
%     histogram(cod_raw,'BinWidth',10);
    hh1 = histogram(cod_raw( cod_raw>=round(thresh1,0) ),'BinWidth',10,'FaceColor',barcol(1,:));
    hh1.BinLimits = hh1.BinLimits + barshift - [10 0];
    hold on;
    hh2 = histogram(cod_raw( cod_raw<round(thresh1,0) ),'BinWidth',10,'FaceColor',barcol(2,:));
    hh2.BinLimits = hh2.BinLimits + barshift - [10 0];
    uistack(hh2,'down');
    if plotwide
     xlim([-150 100]);
     namewide = '-wide';
    else
     xlim([-100 100]);
     namewide = '';
    end
    set(gca,'FontSize',18);
    plot(thresh1.*[1 1],ylim,'k--','LineWidth',1.2);
    plot(medval.*[1 1],ylim,'k-.','LineWidth',1.6);
    title(sprintf('Distribution of Cross-Validated\n Variance Explained in Broadband'));
    xlabel('Broadband cross-validated R^2 (%)');
    ylabel('# of electrodes');
    xticks(unique([xticks thresh1]));
%     xticks(unique([xticks thresh1 round(medval,1)]));
    
figureName = sprintf('prf-distribution-raw%s_permhalves%s_%.0f%%-N%d',namewide,R2mode,alpha_bb*100,length(Nresamp));
if isround, figureName = [figureName '-round'];   end
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end

fprintf('%s median is %.2f%%\n','Broadband',medval);

%% histogram of permutation xR2 in alpha

barcol = [[1 1 1]*0.4; [1 1 1]*0.6];

% visch = find(~ismember(channels.(va_area),'none'));
visch = 1:height(channels);

cod     = cod_a;
%-- plot histogram
hF(end+1) = figure('Position',[300 300 500 400],'MenuBar','none');
subplot_er(1,1,1);
alpha_a = 0.05;
    thresh1 = prctile(cod(Nresamp),(1-alpha_a)*100);
    if isround,     thresh1 = round(thresh1);     end
    barshift = mod(thresh1,10);
    medval  = median(cod(Nresamp));
%     histogram(cod(Nresamp),'BinWidth',10,'FaceColor',[1 1 1]*0.6);
    hh1 = histogram(cod(Nresamp( cod(Nresamp)>=round(thresh1,0) )),'BinWidth',10,'FaceColor',barcol(1,:));
    hh1.BinLimits = hh1.BinLimits + barshift - [10 0];
    hold on;
    hh2 = histogram(cod(Nresamp( cod(Nresamp)<round(thresh1,0) )),'BinWidth',10,'FaceColor',barcol(2,:));
    hh2.BinLimits = hh2.BinLimits + barshift - [10 0];
    uistack(hh2,'down');
    if plotwide
     xlim([-150 100]);
     namewide = '-wide';
    else
     xlim([-100 100]);
     namewide = '';
    end
    set(gca,'FontSize',18);
    plot(thresh1.*[1 1],ylim,'k--','LineWidth',1.2);
    plot(medval.*[1 1],ylim,'k-.','LineWidth',1.6);
    title(sprintf('Distribution of Shuffled\n Variance Explained in Alpha'));
    if thresh1>=70
        text(thresh1-5,diff(ylim)*0.8,sprintf('%.0f%%',(1-alpha_a)*100),'FontSize',18,'HorizontalAlignment','right');
    else
        text(thresh1+5,diff(ylim)*0.8,sprintf('%.0f%%',(1-alpha_a)*100),'FontSize',18);
    end
    xlabel('Alpha shuffled R^2 (%)');
    ylabel('Probability');
    xticks(unique([xticks thresh1]));
%     xticks(unique([xticks thresh1 round(medval,1)]));
    ytickbase = yticks; ytickbase = round(diff(ytickbase(1:2))/numel(Nresamp),2)*numel(Nresamp);
    yticks(0:ytickbase:(ytickbase*length(yticks)));
    yticklabels( num2str(yticks'/numel(Nresamp)) );
    
figureName = sprintf('prf-distribution-shuffle-alpha%s_permhalves%s_%.0f%%-N%d',namewide,R2mode,alpha_a*100,length(Nresamp));
if isround, figureName = [figureName '-round'];   end
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end

fprintf('%s median is %.2f%%\n','Shuffled alpha',medval);

%% histogram of alpha xR2

barcol = [[1 1 1]*0.4; [1 1 1]*0.6];
% barcol = plcol;

% visch = find(~ismember(channels.(va_area),'none'));
visch = 1:height(channels);

cod     = cod_a;
cod_raw = prf_all_a.xval;
%-- plot histogram
hF(end+1) = figure('Position',[300 300 500 400],'MenuBar','none');
subplot_er(1,1,1);
alpha_a = 0.05;
    thresh1  = prctile(cod(Nresamp),(1-alpha_a)*100);
    if isround,     thresh1 = round(thresh1);     end
    barshift = mod(thresh1,10);
    medval  = median(cod_raw);
%     histogram(cod_raw,'BinWidth',10);
    hh1 = histogram(cod_raw( cod_raw>=round(thresh1,0) ),'BinWidth',10,'FaceColor',barcol(1,:));
    hh1.BinLimits = hh1.BinLimits + barshift - [10 0];
    hold on;
    hh2 = histogram(cod_raw( cod_raw<round(thresh1,0) ),'BinWidth',10,'FaceColor',barcol(2,:));
    hh2.BinLimits = hh2.BinLimits + barshift - [10 0];
    uistack(hh2,'down');
    if plotwide
     xlim([-150 100]);
     namewide = '-wide';
    else
     xlim([-100 100]);
     namewide = '';
    end
    set(gca,'FontSize',18);
    plot(thresh1.*[1 1],ylim,'k--','LineWidth',1.2);
    plot(medval.*[1 1],ylim,'k-.','LineWidth',1.6);
    title(sprintf('Distribution of Cross-Validated\n Variance Explained in Alpha'));
    xlabel('Alpha cross-validated R^2 (%)');
    ylabel('# of electrodes');
    xticks(unique([xticks thresh1]));
%     xticks(unique([xticks thresh1 round(medval,1)]));
    
figureName = sprintf('prf-distribution-raw-alpha%s_permhalves%s_%.0f%%-N%d',namewide,R2mode,alpha_a*100,length(Nresamp));
if isround, figureName = [figureName '-round'];   end
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end

fprintf('%s median is %.2f%%\n','Alpha',medval);
