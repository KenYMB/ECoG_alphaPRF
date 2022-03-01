% Compare decimation & full-ts, Channel Selection

% 20201008 Yuasa
% 20220223 Yuasa: update bootstrap method

%% Define paths and dataset
% close all; clear all;
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
SetDefault('issaveplot',true);
if issaveplot
    plotsavePth    = fullfile(figPth, 'pRFselections-representative');
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

plcol  = get(groot,'defaultAxesColorOrder');
plcol(1:2,:) = [49 127 160; 183 20 51]./255;

if issaveplot
    plotsavedir    = fullfile(plotsavePth, 'channelselection', postfix(2:end));
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

% close all;

%% Channel selection with FPM

threshold = threshold_bb;
nboot     = 5000;
%-- compute xR2_a in ROIs (wang_prob)
filename = sprintf('all_chansel-wangprob%s-%s%s_thresh%02d',R2mode,selectchs,postfix,threshold);
filepath = fullfile(savePth,'pRFanalysis',filename);

if exist([filepath '.mat'],'file')
 fprintf('loading %s...',filename);
 load(filepath);
 fprintf('\n');
else
 fprintf('computing %s\n',filename);
 rois = categories(prf_all_a.channels.(va_area));
 roiidx = ~ismember(rois,{'none','PHC','SPL1','FEF'});
 rois = [rois(roiidx);{'V1-V3';'dorsolateral'}];      % add low/high visual area
 nroi = length(rois);
 mn_keep_roi = cell(1,nroi); mn_discard_roi = cell(1,nroi);
 n_keep_roi  = cell(1,nroi); n_discard_roi  = cell(1,nroi);
 for iboot=1:nboot       % ~3.4m for 2000 iteration
     %%-- assign roi label based on probability
     [~,prf_all_a]   = ecog_rearrangePRF(prf_params_a,va_area,'prob');
     %%-- randomly select channels
     boot_keep = randi(height(prf_all_a.channels),1,height(prf_all_a.channels));
     channels  = prf_all_a.channels(boot_keep,:);
     xval      = prf_all_a.xval(boot_keep,:);
     elec_ok   = prf_all_bb.xval(boot_keep,:)>threshold;
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

%% %%%%%%%%%%%%%%%%
%% distribution
%% %%%%%%%%%%%%%%%%
%% Channel selection with FPM in low/high visual area in one panel -short-

%-- reject inaccurate ROIs
shiftframe = 0.15;
shiftbar   = 0.02;

hF = figure('defaultAxesColorOrder',plcol);
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
ylim([-20 60]);
xlim([0.5,nsel+0.5]+shiftframe);  plot(xlim,[0,0],'k--');
plot([1 1].*(nsel-1.5)+shiftframe,ylim,'k:','LineWidth',1.5);
xticks(1:nsel); xticklabels(rois(selroi)); xtickangle(0);
set(gca,'FontSize',fontSiz);
 
%     ylabel('Alpha cross-validated R^2 (%)');
    ylabel('Variance Explained in Alpha (%)');
    hl=legend({'Visually Selective Electrodes','Non-Visually Selective Electrodes'});
    
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


htt = title(ht,'Cross-validated Variance Explained in Alpha','FontSize',fontSiz,'FontWeight','bold');
figname = sprintf('channelselection%s-%02d%%-%s',R2mode,threshold,'short');
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));  end

%% %%%%%%%%%%%%%%%%
%% histogram
%% %%%%%%%%%%%%%%%%

%% permutation
filename = sprintf('all_cod-permhalves%s-%s%s',R2mode,selectchs,postfix);
filepath = fullfile(savePth,'pRFanalysis',filename);

if ~exist([filepath '.mat'],'file')
    %%% params
    nboot   = 5000;
    alpha    = .05;
    
    bootchs  = [randi(nchan,nboot,1) randi(nchan-1,nboot,1)];
        bootchs(diff(bootchs,[],2)==0,2) = nchan; % avoid correct pairs
    model_data_a   = cellfun(@(C) C(bootchs(:,1),:),model_all_a.datats,'UniformOutput',false);
    model_data_bb  = cellfun(@(C) C(bootchs(:,1),:),model_all_bb.datats,'UniformOutput',false);
    bootstimulus   = model_all_bb.stimulus(bootchs(:,1),:);
    params_perm_a  = prf_all_a.params(:,:,bootchs(:,2));
    params_perm_bb = prf_all_bb.params(:,:,bootchs(:,2));
    
    %%% scrambled R2
    tic;
    %-- 31.96s for all subjects in bb & a with 1 iteration (41.30s for 10 iterations => 1.04s/iter)
    [~, ~, cod_bb] = ecog_computePRFtimeseries(model_all_bb.stimulus,model_all_bb.datats,params_perm_bb,prf_all_bb.options);
    ptest_bb   = prf_all_bb.xval > prctile(cod_bb,(1-alpha)*100,3);
    fprintf('%.1f%% channels are significantly good for broadband\n',sum(ptest_bb)/length(ptest_bb)*100);
    [~, ~, cod_a] = ecog_computePRFtimeseries(model_all_a.stimulus,model_all_a.datats,params_perm_a,prf_all_a.options);
    ptest_a   = prf_all_a.xval > prctile(cod_a,(1-alpha)*100,3);
    fprintf('%.1f%% channels are significantly good for alpha\n',sum(ptest_a)/length(ptest_a)*100);
    toc;
    
    saveauto(filepath,'bootchs', 'alpha', 'cod_bb', 'cod_a', 'ptest_bb', 'ptest_a');
else
    clear Nresamp;
    load(filepath);
    nboot = size(bootchs,1);
    model_data_a   = cellfun(@(C) C(bootchs(:,1),:),model_all_a.datats,'UniformOutput',false);
    model_data_bb  = cellfun(@(C) C(bootchs(:,1),:),model_all_bb.datats,'UniformOutput',false);
    bootstimulus   = model_all_bb.stimulus(bootchs(:,1),:);
    params_perm_a  = prf_all_a.params(:,:,bootchs(:,2));
    params_perm_bb = prf_all_bb.params(:,:,bootchs(:,2));
end

%-- reset prf_all
[~,prf_all_a]      = ecog_rearrangePRF(prf_params_a,va_area);
[~,prf_all_bb]     = ecog_rearrangePRF(prf_params_bb,va_area);

%%
% close all;

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
hF(2) = figure('Position',[300 300 500 400],'MenuBar','none');
tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
nexttile;
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
    
figname = sprintf('prf-distribution-shuffle%s_permhalves%s_%.0f%%-N%d',namewide,R2mode,alpha_bb*100,length(Nresamp));
if isround, figname = [figname '-round'];   end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));  end


%% histogram of broadband xR2

barcol = [[1 1 1]*0.4; [1 1 1]*0.6];
% barcol = plcol;

% visch = find(~ismember(channels.(va_area),'none'));
visch = 1:height(channels);

cod     = cod_bb;
cod_raw = prf_all_bb.xval;
%-- plot histogram
hF(3) = figure('Position',[300 300 500 400],'MenuBar','none');
tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
nexttile;
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
    
figname = sprintf('prf-distribution-raw%s_permhalves%s_%.0f%%-N%d',namewide,R2mode,alpha_bb*100,length(Nresamp));
if isround, figname = [figname '-round'];   end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));  end

%% histogram of permutation xR2 in alpha

barcol = [[1 1 1]*0.4; [1 1 1]*0.6];

% visch = find(~ismember(channels.(va_area),'none'));
visch = 1:height(channels);

cod     = cod_a;
%-- plot histogram
hF(4) = figure('Position',[300 300 500 400],'MenuBar','none');
tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
nexttile;
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
    
figname = sprintf('prf-distribution-shuffle-alpha%s_permhalves%s_%.0f%%-N%d',namewide,R2mode,alpha_a*100,length(Nresamp));
if isround, figname = [figname '-round'];   end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));  end

%% histogram of alpha xR2

barcol = [[1 1 1]*0.4; [1 1 1]*0.6];
% barcol = plcol;

% visch = find(~ismember(channels.(va_area),'none'));
visch = 1:height(channels);

cod     = cod_a;
cod_raw = prf_all_a.xval;
%-- plot histogram
hF(5) = figure('Position',[300 300 500 400],'MenuBar','none');
tiledlayout(1,1,'TileSpacing','compact','Padding','compact');
nexttile;
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
    
figname = sprintf('prf-distribution-raw-alpha%s_permhalves%s_%.0f%%-N%d',namewide,R2mode,alpha_a*100,length(Nresamp));
if isround, figname = [figname '-round'];   end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));  end

