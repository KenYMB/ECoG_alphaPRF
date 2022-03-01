% Compare decimation & full-ts, Channel Selection

% 20210408 - update to demonstrate reverse selection

%% Define paths and dataset
% close all; clear all;
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
SetDefault('issaveplot',true);
if issaveplot
    plotsavePth    = fullfile(figPth, 'pRFselections-representative');
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%% load analyzePRF & recompute full-ts R2 & set data for bootstrap with half-trials
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

clear alphaType broadbandType
ecog_APRF_INITa_loaddata;
ecog_APRF_INITb_mergedata;
ecog_APRF_INITc_postfix;
ecog_APRF_INITd_threshold;

%% %%%%%%%%%%%%%%
%% Visualize
%% %%%%%%%%%%%%%%
%%% prepare visualize
res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;
fontSiz   = 18;

alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');
plcol(1:2,:) = [49 127 160; 183 20 51]./255;

plotsavedir    = fullfile(plotsavePth, 'channelselection', postfix(2:end));
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

close all;

%% Channel selection with FPM

threshold = threshold_a;
nboot     = 5000;
%-- compute xR2_a in ROIs (wang_prob)
filename = sprintf('all_chansel-wangprob%s-%s%s_thresh%02d-rev',R2mode,selectchs,postfix,threshold);
filepath = fullfile(savePth,'pRFanalysis',filename);

if exist([filepath '.mat'],'file')
 fprintf('loading %s...',filename);
 load(filepath);
 fprintf('\n');
else
 fprintf('computing %s\n',filename);
 rois = categories(prf_all_bb.channels.(va_area));
 roiidx = ~ismember(rois,{'none','PHC','SPL1','FEF'});
 rois = [rois(roiidx);{'low';'high'}];      % add low/high visual area
 nroi = length(rois);
 mn_keep_roi = cell(1,nroi); mn_discard_roi = cell(1,nroi);
 n_keep_roi  = cell(1,nroi); n_discard_roi  = cell(1,nroi);
 for iboot=1:nboot       % ~3.4m for 2000 iteration
     %%-- assign roi label based on probability
     [~,prf_all_bb]   = ecog_rearrangePRF(prf_params_bb,va_area,'prob');
     %%-- randomly select channels
     boot_keep = randi(height(prf_all_bb.channels),1,height(prf_all_bb.channels));
     channels  = prf_all_bb.channels(boot_keep,:);
     xval      = prf_all_bb.xval(boot_keep,:);
     elec_ok   = prf_all_a.xval(boot_keep,:)>threshold;
     for iroi=1:nroi
         switch lower(rois{iroi})
             case {'low','v1-v3'}
                 catrois = {'V1';'V2';'V3'};
                 conbinedroi = true;
             case {'high','dorsolateral'}
                 catrois = {'V3a';'V3b';'LO1';'LO2';'TO';'IPS'};
                 conbinedroi = true;
             case {'high+'}
                 catrois = {'hV4','V3a';'V3b';'LO1';'LO2';'TO';'IPS'};
                 conbinedroi = true;
             otherwise
                 conbinedroi = false;
         end
         if conbinedroi
             elec_roi = false(height(channels),1);
             for jroi=1:length(catrois)
                 elec_roi = elec_roi | (channels.(va_area)==catrois{jroi});
             end
         else
             elec_roi = channels.(va_area)==rois{iroi};
         end
         elec_keep = elec_ok & elec_roi;
         elec_discard = ~elec_ok & elec_roi;
         if any(elec_keep)
             mn_keep_roi{iroi}(iboot) = mean(xval(elec_keep));
         else
             mn_keep_roi{iroi}(iboot) = nan;
         end
         n_keep_roi{iroi}(iboot)  = sum(elec_keep);         % electrode numbers (save even for invalid condition)
         if any(elec_discard)
             mn_discard_roi{iroi}(iboot) = mean(xval(elec_discard));
         else
             mn_discard_roi{iroi}(iboot) = nan;
         end
         n_discard_roi{iroi}(iboot)  = sum(elec_discard);   % electrode numbers (save even for invalid condition)
     end
 end
 
 saveauto(filepath,'rois','mn_keep_roi','mn_discard_roi','n_keep_roi','n_discard_roi','nboot','threshold');
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
figure('defaultAxesColorOrder',plcol),
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
set(gca,'FontSize',fontSiz);
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
        'FontSize',fontSiz,'Color',plcol(1,:));
    text(ii+0.05+shiftbar, nanmean(mn_discard_roi{selroi(ii)})-mn_shift, sprintf('%5.1f',nanmean(n_discard_roi{selroi(ii)})),...
        'FontSize',fontSiz,'Color',plcol(2,:));
end

fprintf('%s\n',sprintf('%5s',rois{selroi}));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_keep_roi(selroi))));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_discard_roi(selroi))));
fprintf('%s\n',sprintf('%5d',cellfun(@(C) prctile(C,talpha*100),mn_keep_roi(selroi)) - cellfun(@(C) prctile(C,(1-talpha)*100),mn_discard_roi(selroi))>0));

figname = sprintf('channelselection%s-rev-%02d%%-%s',R2mode,threshold,roiset);
savefigauto(gcf, fullfile(plotsavedir, figname));
end

%% Channel selection with FPM in low/high visual area in one panel

%-- reject inaccurate ROIs
shiftframe = 0.15;
shiftbar   = 0.02;
tnsel      = 11;

figure('defaultAxesColorOrder',plcol),
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
set(gca,'FontSize',fontSiz);
 
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
        'FontSize',fontSiz,'Color',plcol(1,:));
    text(ii+0.05+shiftbar, nanmean(mn_discard_roi{selroi(ii)})-mn_shift, sprintf('%5.1f',nanmean(n_discard_roi{selroi(ii)})),...
        'FontSize',fontSiz,'Color',plcol(2,:));
end

fprintf('%s\n',sprintf('%5s',rois{selroi}));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_keep_roi(selroi))));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_discard_roi(selroi))));
fprintf('%s\n',sprintf('%5d',cellfun(@(C) prctile(C,talpha*100),mn_keep_roi(selroi)) - cellfun(@(C) prctile(C,(1-talpha)*100),mn_discard_roi(selroi))>0));

end
htt = title(ht,'Cross-validated Variance Explained in Broadband','FontSize',fontSiz,'FontWeight','bold');
figname = fullfile(plotsavedir, sprintf('channelselection%s-rev-%02d%%-%s',R2mode,threshold,'merge'));
savefigauto(gcf,figname);


%% Channel selection with FPM in low/high visual area in one panel -short-

%-- reject inaccurate ROIs
shiftframe = 0.15;
shiftbar   = 0.02;

figure('defaultAxesColorOrder',plcol),
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
set(gca,'FontSize',fontSiz);
 
    ylabel('Variance Explained in Broadband (%)');
    hl=legend({sprintf('Variance explained in Alpha > %0d%%',threshold_a),...
               sprintf('Variance explained in Alpha ≤ %0d%%',threshold_a)});
    
set(gcf,'Unit','pixels')

for ii=1:nsel
    mn_dif = nanmean(mn_keep_roi{selroi(ii)}) - nanmean(mn_discard_roi{selroi(ii)});
    if abs(mn_dif)<4,   mn_shift = 2*sign(mn_dif);
    else,               mn_shift = 0;  end
    text(ii+0.05-shiftbar, nanmean(mn_keep_roi{selroi(ii)})+mn_shift, sprintf('%5.1f',nanmean(n_keep_roi{selroi(ii)})),...
        'FontSize',fontSiz,'Color',plcol(1,:));
    text(ii+0.05+shiftbar, nanmean(mn_discard_roi{selroi(ii)})-mn_shift, sprintf('%5.1f',nanmean(n_discard_roi{selroi(ii)})),...
        'FontSize',fontSiz,'Color',plcol(2,:));
end

fprintf('%s\n',sprintf('%5s',rois{selroi}));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_keep_roi(selroi))));
fprintf('%s\n',sprintf('%5.1f',cellfun(@nanmean,n_discard_roi(selroi))));
fprintf('%s\n',sprintf('%5d',cellfun(@(C) prctile(C,talpha*100),mn_keep_roi(selroi)) - cellfun(@(C) prctile(C,(1-talpha)*100),mn_discard_roi(selroi))>0));


htt = title(ht,'Cross-validated Variance Explained in Broadband','FontSize',fontSiz,'FontWeight','bold');
figname = fullfile(plotsavedir, sprintf('channelselection%s-rev-%02d%%-%s',R2mode,threshold,'short'));
savefigauto(gcf,figname);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot xR2 relations in 2D histogram
sbfld = {'R2'};
plotwide = true;
% plotwide = false;

figure('Position',[1000/2 900/2 500 440],'MenuBar','none');

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
    set(gca,'FontSize',fontSiz)
    
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

figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_2Dhist_%s%s',threshold_bb,threshold_a,eclimit,sbfld{:},namewide);
savefigauto(gcf, fullfile(plotsavedir, figname));


%% plot xR2 relations in scatter plot
sbfld = {'R2'};
plotwide = true;
% plotwide = false;

figure('Position',[1000/2 900/2 500 440],'MenuBar','none','defaultAxesColorOrder',plcol),

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
    set(gca,'FontSize',fontSiz)
    
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

figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_scatter_%s%s',threshold_bb,threshold_a,eclimit,sbfld{:},namewide);
savefigauto(gcf, fullfile(plotsavedir, figname));

%%
close all;
