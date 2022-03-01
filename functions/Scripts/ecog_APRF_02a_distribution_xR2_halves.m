% Scrambled channel label to check xR2 distributions
%   Separate time series in halves (as in analyzePRF)
%   Update xR2 with full time-series (no decimation)

% 20210910 Yuasa - folk from ecog_APRF_06bb_distribution_xR2_halves_gen
% 20220201 Yuasa - reconstruct script for paper

%% Define paths and dataset
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
plotsavePth    = fullfile(figPth, 'pRFrelations');
if ~exist(plotsavePth,'dir'), mkdir(plotsavePth); end
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%% load analyzePRF & recompute full-ts R2 & set data for bootstrap with half-trials
% average        ='runs';
% smoothingMode  ='decimate';
% smoothingN     = 3;
% prfmodel       ='linear';
% gaussianmode   ='gs';
% selectchs      = 'wangprobchs';
%     allowlag       = false;
%     allowbeta      = true;
%     allowwide      = true;
%     allowmixbeta   = true;
va_area = 'wangarea';
    
    SetDefault('average','runs');
    SetDefault('smoothingMode','decimate');
    SetDefault('smoothingN',3);
    SetDefault('prfmodel','linear');
    SetDefault('gaussianmode','gs');
    SetDefault('selectchs','wangprobchs');
    SetDefault('allowlag',false);
    SetDefault('allowbeta',true);
    SetDefault('allowwide',true);
    SetDefault('allowmixbeta',true);
    SetDefault('justcompute',false);
    SetDefault('forcecompute',false);

usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = true;

clear alphaType broadbandType
ecog_APRF_INITa_loaddata;
ecog_APRF_INITb_mergedata;
ecog_APRF_INITc_postfix;

%% permutation
filename = sprintf('all_cod-permhalves%s-%s%s',R2mode,selectchs,postfix);
filepath = fullfile(savePth,'pRFanalysis',filename);
if ~exist(fileparts(filepath),'dir'), mkdir(fileparts(filepath)); end

if ~exist([filepath '.mat'],'file') || forcecompute
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
    [~, ~, cod_bb] = ecog_computePRFtimeseries(bootstimulus,model_data_bb,params_perm_bb,prf_all_bb.options);
    ptest_bb   = prf_all_bb.xval > prctile(cod_bb,(1-alpha)*100,'all');
    fprintf('%.1f%% channels are significantly good for broadband\n',sum(ptest_bb)/length(ptest_bb)*100);
    [~, ~, cod_a] = ecog_computePRFtimeseries(bootstimulus,model_data_a,params_perm_a,prf_all_a.options);
    ptest_a   = prf_all_a.xval > prctile(cod_a,(1-alpha)*100,'all');
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

%% %%%%%%%%%%%%%%%%%%%
%% Visualize
%% %%%%%%%%%%%%%%%%%%%
if ~justcompute
    
plotsavedir    = fullfile(plotsavePth, postfix(2:end),'threshold');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%% %%%%%%%%%%%%%%%%%
%% not set individual threshold for each electrode
close all;

Nresamp = 1:nboot;

for isround = 0:1
% isround = 1;

%% histogram

% visch = find(~ismember(channels.(va_area),'none'));
visch = 1:height(channels);

%-- plot histogram
figure('Position',[300 300 1000 400],'MenuBar','none');
subplot_er(1,2,1);
alpha_bb = 0.05;
    histogram(cod_bb(Nresamp),'BinWidth',10);
    xlim([-100 100]);
    set(gca,'FontSize',18);
    hold on;
    thresh1 = prctile(cod_bb(Nresamp),(1-alpha_bb)*100);
    if isround,     thresh1 = round(thresh1);     end
    plot(thresh1.*[1 1],ylim,'r--','LineWidth',1.2);
    title('broadband');
    if thresh1>=70
        text(thresh1-5,diff(ylim)*0.8,sprintf('%.2f%%',thresh1),'FontSize',18,'HorizontalAlignment','right');
    else
        text(thresh1+5,diff(ylim)*0.8,sprintf('%.2f%%',thresh1),'FontSize',18);
    end
    text((max(xlim)+thresh1)/2,diff(ylim)*0.6,sprintf('N=%d',sum(prf_all_bb.xval>thresh1)),'FontSize',18,'HorizontalAlignment','center');
subplot_er(1,2,2);
alpha_a = 0.05;
% alpha_a = 0.32;
% [f,x] = ecdf(cod_a(Nresamp)); alpha_a = 1 - f(nearest(x,0));
    histogram(cod_a(Nresamp),'BinWidth',10);
    xlim([-100 100]);
    set(gca,'FontSize',18);
    hold on;
    thresh2 = prctile(cod_a(Nresamp),(1-alpha_a)*100);
    if isround,     thresh2 = round(thresh2);     end
    plot(thresh2.*[1 1],ylim,'r--','LineWidth',1.2);
    title('alpha');
    if thresh2>=70
        text(thresh2-5,diff(ylim)*0.8,sprintf('%.2f%%',thresh2),'FontSize',18,'HorizontalAlignment','right');
    else
        text(thresh2+5,diff(ylim)*0.8,sprintf('%.2f%%',thresh2),'FontSize',18);
    end
    text((max(xlim)+thresh2)/2,diff(ylim)*0.6,sprintf('N=%d',sum(prf_all_a.xval>thresh2)),'FontSize',18,'HorizontalAlignment','center');
fprintf('N=%d electrodes satisfy both thresholds\n',...
    sum((prf_all_bb.xval>thresh1)&(prf_all_a.xval>thresh2)));
    
figname = fullfile(plotsavedir, sprintf('prf-distribution_permhalves%s-%s_%.0f%%-%.0f%%-N%d',R2mode,selectchs,alpha_bb*100,alpha_a*100,length(Nresamp)));
if isround, figname = [figname '-round'];   end
hgexport(gcf,figname,hgexport('factorystyle'),'Format','png');

%% show threshold percentage

alpha_bb = 0.05;
alpha_a  = 0.05;
% alpha_a  = 0.32;
% [f,x]  = ecdf(cod_a(Nresamp)); alpha_a = 1 - f(nearest(x,0));
eclimit = 50;
    fprintf('------------------------------\n');
    thresh1 = prctile(cod_bb(Nresamp),(1-alpha_bb)*100);
    if isround,     thresh1 = round(thresh1);     end
    fprintf('broadband: N = %d under R2>%.4f%%\n',sum(prf_all_bb.xval>thresh1),thresh1);
    thresh2 = prctile(cod_a(Nresamp),(1-alpha_a)*100);
    if isround,     thresh2 = round(thresh2);     end
    fprintf('alpha: N = %d under R2>%.4f%%\n',sum(prf_all_a.xval>thresh2),thresh2);
    fprintf('both: N=%d (N=%d under ecc<50pixel)\n',...
        sum((prf_all_bb.xval>thresh1)&(prf_all_a.xval>thresh2)),...
        sum((prf_all_bb.xval>thresh1)&(prf_all_a.xval>thresh2)...
             & prf_all_bb.ecc < eclimit & prf_all_a.ecc < eclimit));
    fprintf('------------------------------\n');
end
end
