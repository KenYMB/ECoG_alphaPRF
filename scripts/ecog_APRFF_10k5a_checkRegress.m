% Compare ERP regression
%   Folk from ecog_APRF_07b_checkmodel_fullts
 
% 20231129 Yuasa
 
%% prefix
% close all; clearvars;
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
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth), 'modelselection', 'regression');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');
 
%-- Plotting Setting
FntSiz    = 20;

%% general parameter
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
va_area = 'wangarea';

usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

%% load analyzePRF @DOG,OG,LFS
modelNames  = ["With ERP-Regression","Without ERP-Regression"];

modeldata  = struct();
prf_params = struct();
model_all  = struct();
prf_all    = struct();
threshold  = struct();
for ii=1:length(modelNames)
if contains(modelNames(ii),'Without')
    modeldataID = 'freq_spectra-timeseries_noregress';
    prfID       = 'prf_noregress';
else
    modeldataID = 'freq_spectra-timeseries';
    prfID       = 'prf';
end
 
clear alphaType broadbandType
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
if contains(modelNames(ii),'Without')
postfix = [postfix '_noregress'];
end
ecog_APRFF_INITd_threshold
 
%-- merge loaded data
modeldata(ii).a   = modeldata_a;
modeldata(ii).bb  = modeldata_bb;
prf_params(ii).a  = prf_params_a;
prf_params(ii).bb = prf_params_bb;
 
model_all(ii).a   = model_all_a;
model_all(ii).bb  = model_all_bb;
prf_all(ii).a     = prf_all_a;
prf_all(ii).bb    = prf_all_bb;

threshold(ii).a   = threshold_a;
threshold(ii).bb  = threshold_bb;


%-- flip gain for alpha (convert from ratio of suppression to power change)
% model_all(ii).a.datats = cellfun(@(C)-C,model_all(ii).a.datats,'UniformOutput',false);
prf_all(ii).a.gain     = -prf_all(ii).a.gain;
prf_all(ii).a.params   = prf_all(ii).a.params.*[1,1,1,-1,1,1,1];
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FaCLb vs FaCLb-withERP
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelsets = [2 1];      % w/o vs w/ Regress
targetBAND = ('noRegress');

% eclimit = 50;
idat = 1;
okch1bb = ~(prf_all(modelsets(idat)).bb.xval<=threshold(modelsets(idat)).bb | prf_all(modelsets(idat)).bb.ecc >= eclimit);  % use tilde for nan
okch1a  = ~(prf_all(modelsets(idat)).a.xval<=threshold(modelsets(idat)).a | prf_all(modelsets(idat)).a.ecc >= eclimit);  % use tilde for nan
idat = 2;
okch2bb = ~(prf_all(modelsets(idat)).bb.xval<=threshold(modelsets(idat)).bb | prf_all(modelsets(idat)).bb.ecc >= eclimit);  % use tilde for nan
okch2a  = ~(prf_all(modelsets(idat)).a.xval<=threshold(modelsets(idat)).a | prf_all(modelsets(idat)).a.ecc >= eclimit);  % use tilde for nan

okchORbb = okch1bb | okch2bb;
okchORa  = okch1a  | okch2a;

onlych1bb = okch1bb&~okch2bb;
onlych2bb = ~okch1bb&okch2bb;

% close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Responses (as N-1 fold increase or decrease: 100% == 0 == no change)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res       = [100 100];
stim = load('bar_apertures.mat','bar_apertures');
stim = imresize(stim.bar_apertures,res);

newparam = 'peakresp';
for ii=1:length(prf_all)
opt = [];
opt.flipgain      = 'no';
opt.skipprojection = 'yes';
opt.catchan        = 'yes';
modelts  = reconPRFdog(prf_all(ii).bb,stim,opt);
gainsign = sign(prf_all(ii).bb.gain);
prf_all(ii).bb.(newparam) = max(modelts'.*gainsign,[],2,'omitnan').*gainsign./100;  % % change -> ratio (N-1 fold change)
modelts  = reconPRFdog(prf_all(ii).a,stim,opt);
gainsign = sign(prf_all(ii).a.gain);
prf_all(ii).a.(newparam) = (10.^max(modelts'.*gainsign,[],2,'omitnan')-1).*gainsign;  % log ratio -> N-1 fold change
end

%% %%%%%%%%%%%%%%
%% Visualize
%% %%%%%%%%%%%%%%
%%% prepare visualize
% res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;
 
alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');
hF = gobjects(0);

ispltall = false;

%% pRF parameters
if issaveplot

if ispltall
    selch_bb = true(size(okch1bb));  selch_a = true(size(okch1a));
    chmode = 'all';
else
    selch_bb = okchORbb;   selch_a = okchORa;
    chmode = 'thrsh';
end

% selch_bb = onlych1bb;  selch_a = onlych1bb;
% selch_bb = onlych2bb;  selch_a = onlych2bb;

tarparam = 'xval';
hF(end+1) = figure; hT=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
set(gcf,'Position',get(gcf,'Position').*[0.8 1 1.6 1]);
nexttile; scatter(prf_all(modelsets(1)).bb.(tarparam)(selch_bb),prf_all(modelsets(2)).bb.(tarparam)(selch_bb)); hold on;
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom');
  hl = plot(xlim,[0,0],'k-',[0,0],ylim,'k-','LineWidth',1.5); uistack(hl,'bottom');
  xlim([max(min(xlim),-100),min(max(xlim),100)]);
  ylim([max(min(ylim),-100),min(max(ylim),100)]);
  title('Broadband'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
nexttile; scatter(prf_all(modelsets(1)).a.(tarparam)(selch_a),prf_all(modelsets(2)).a.(tarparam)(selch_a)); hold on;
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom');
  hl = plot(xlim,[0,0],'k-',[0,0],ylim,'k-','LineWidth',1.5); uistack(hl,'bottom');
  xlim([max(min(xlim),-100),min(max(xlim),100)]);
  ylim([max(min(ylim),-100),min(max(ylim),100)]);
  title('alpha'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
set(findall(gcf,'-property','FontSize'),'FontSize',FntSiz);
title(hT,'Variance Explained (%)','FontSize',FntSiz);
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-compare_ch-%s_%s',chmode,tarparam)));  end

tarparam = 'ang';
hF(end+1) = figure; hT=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
set(gcf,'Position',get(gcf,'Position').*[0.8 1 1.6 1]);
nexttile; scatter(prf_all(modelsets(1)).bb.(tarparam)(selch_bb),prf_all(modelsets(2)).bb.(tarparam)(selch_bb)); hold on;
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom');
  xlim([0,360]); xticks(0:90:360);
  ylim([0,360]); yticks(0:90:360);
  title('Broadband'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
nexttile; scatter(prf_all(modelsets(1)).a.(tarparam)(selch_a),prf_all(modelsets(2)).a.(tarparam)(selch_a)); hold on;
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom');
  xlim([0,360]); xticks(0:90:360);
  ylim([0,360]); yticks(0:90:360);
  title('Alpha'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
set(findall(gcf,'-property','FontSize'),'FontSize',FntSiz);
title(hT,'Angle (deg)','FontSize',FntSiz);
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-compare_ch-%s_%s',chmode,tarparam)));  end

tarparam = 'ecc';
hF(end+1) = figure; hT=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
set(gcf,'Position',get(gcf,'Position').*[0.8 1 1.6 1]);
nexttile; scatter(prf_all(modelsets(1)).bb.(tarparam)(selch_bb).*cfactor,prf_all(modelsets(2)).bb.(tarparam)(selch_bb).*cfactor); hold on;
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom');
  title('Broadband'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
nexttile; scatter(prf_all(modelsets(1)).a.(tarparam)(selch_a).*cfactor,prf_all(modelsets(2)).a.(tarparam)(selch_a).*cfactor); hold on;
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom');
  title('Alpha'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
set(findall(gcf,'-property','FontSize'),'FontSize',FntSiz);
title(hT,'Eccentricity (deg)','FontSize',FntSiz);
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-compare_ch-%s_%s',chmode,tarparam)));  end

tarparam = 'rfsize';
hF(end+1) = figure; hT=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
set(gcf,'Position',get(gcf,'Position').*[0.8 1 1.6 1]);
nexttile; scatter(prf_all(modelsets(1)).bb.(tarparam)(selch_bb).*cfactor,prf_all(modelsets(2)).bb.(tarparam)(selch_bb).*cfactor); hold on;
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom');
  title('Broadband'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
nexttile; scatter(prf_all(modelsets(1)).a.(tarparam)(selch_a).*cfactor,prf_all(modelsets(2)).a.(tarparam)(selch_a).*cfactor); hold on;
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom');
  title('Alpha'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
set(findall(gcf,'-property','FontSize'),'FontSize',FntSiz);
title(hT,'Size (deg)','FontSize',FntSiz);
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-compare_ch-%s_%s',chmode,tarparam)));  end

tarparam = 'peakresp';
hF(end+1) = figure; hT=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
set(gcf,'Position',get(gcf,'Position').*[0.8 1 1.6 1]);
nexttile; scatter(prf_all(modelsets(1)).bb.(tarparam)(selch_bb),prf_all(modelsets(2)).bb.(tarparam)(selch_bb)); hold on;
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom');
  title('Broadband'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
    xlim([max(min(xlim),-100),min(max(xlim),100)]);    ylim([max(min(ylim),-100),min(max(ylim),100)]);
    if length(xticks)>5,   if max(xlim)<=0  % negative axis
                               tickval = xticks; xticks(flip(tickval(length(tickval):-2:1)));
                           else
                               tickval = xticks; xticks(tickval(1:2:length(tickval)));
                           end
    end
    if length(yticks)>5,   if max(ylim)<=0  % negative axis
                               tickval = yticks; yticks(flip(tickval(length(tickval):-2:1)));
                           else
                               tickval = yticks; yticks(tickval(1:2:length(tickval)));
                           end
    end
    xticklabels(xticks+sign(xticks)+(xticks==0));    ticknames = xticklabels;
    ticknames(xticks<0)  = cellfun(@(d) sprintf('1/%d',-str2double(d)),ticknames(xticks<0),'UniformOutput',false);
    xticklabels(ticknames);
    yticklabels(yticks+sign(yticks)+(yticks==0));    ticknames = yticklabels;
    ticknames(yticks<0)  = cellfun(@(d) sprintf('1/%d',-str2double(d)),ticknames(yticks<0),'UniformOutput',false);
    yticklabels(ticknames);
nexttile; scatter(prf_all(modelsets(1)).a.(tarparam)(selch_a),prf_all(modelsets(2)).a.(tarparam)(selch_a)); hold on;
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom');
  title('Alpha'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
    xlim([max(min(xlim),-100),min(max(xlim),100)]);    ylim([max(min(ylim),-100),min(max(ylim),100)]);
    if length(xticks)>5,   if max(xlim)<=0  % negative axis
                               tickval = xticks; xticks(flip(tickval(length(tickval):-2:1)));
                           else
                               tickval = xticks; xticks(tickval(1:2:length(tickval)));
                           end
    end
    if length(yticks)>5,   if max(ylim)<=0  % negative axis
                               tickval = yticks; yticks(flip(tickval(length(tickval):-2:1)));
                           else
                               tickval = yticks; yticks(tickval(1:2:length(tickval)));
                           end
    end
    xticklabels(xticks+sign(xticks)+(xticks==0));    ticknames = xticklabels;
    ticknames(xticks<0)  = cellfun(@(d) sprintf('1/%d',-str2double(d)),ticknames(xticks<0),'UniformOutput',false);
    xticklabels(ticknames);
    yticklabels(yticks+sign(yticks)+(yticks==0));    ticknames = yticklabels;
    ticknames(yticks<0)  = cellfun(@(d) sprintf('1/%d',-str2double(d)),ticknames(yticks<0),'UniformOutput',false);
    yticklabels(ticknames);
set(findall(gcf,'-property','FontSize'),'FontSize',FntSiz);
title(hT,'Response Gain (fold)','FontSize',FntSiz);
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-compare_ch-%s_%s',chmode,tarparam)));  end

end

%% pRF parameters with ROI
plcol_wang = [ones(3,1)*plcol(1,:);...
              ones(6,1)*plcol(4,:);];
plshp_wang = [repmat('o',1,3) repmat('^',1,6)];
roilist = {'V1','V2','V3','V3a','V3b','LO1','LO2','TO','IPS'}';

if ispltall
    selch_bb = true(size(okch1bb));  selch_a = true(size(okch1a));
    chmode = 'all';
else
    selch_bb = okchORbb;   selch_a = okchORa;
    chmode = 'thrsh';
end

tarparam = 'xval';
hF(end+1) = figure('Menubar','none','Position',[200 200 900 440],'defaultAxesColorOrder',plcol_wang);
hT=tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
nexttile; ii=1;
    for iroi = roilist'
    roich = ismember(prf_all(modelsets(1)).bb.channels.wangarea,iroi);
    okch  = selch_bb & roich;
    scatter(prf_all(modelsets(1)).bb.(tarparam)(okch),prf_all(modelsets(2)).bb.(tarparam)(okch),plshp_wang(ii),'LineWidth',1.6);
    hold on;    ii=ii+1;
    end
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom'); % diagonal
  hl = plot(xlim,[0,0],'k-',[0,0],ylim,'k-','LineWidth',1.5); uistack(hl,'bottom');  % axis
  hl = plot([getelement(xlim,1) threshold(modelsets(1)).bb],[1 1].*threshold(modelsets(2)).bb,'k-.',...
            [1 1].*threshold(modelsets(1)).bb,[getelement(ylim,1) threshold(modelsets(2)).bb],'k-.');  uistack(hl,'bottom'); % threshold
  axis([max(min(xlim),-100),min(max(xlim),100), max(min(xlim),-100),min(max(xlim),100)],'square');
  title('Broadband'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
nexttile; ii=1;
    for iroi = roilist'
    roich = ismember(prf_all(modelsets(1)).a.channels.wangarea,iroi);
    okch  = selch_a & roich;
    scatter(prf_all(modelsets(1)).a.(tarparam)(okch),prf_all(modelsets(2)).a.(tarparam)(okch),plshp_wang(ii),'LineWidth',1.6);
    hold on;    ii=ii+1;
    end
  hl = plot([min([xlim,ylim]),max([xlim,ylim])],[min([xlim,ylim]),max([xlim,ylim])],'k--','LineWidth',1.5); uistack(hl,'bottom'); % diagonal
  hl = plot(xlim,[0,0],'k-',[0,0],ylim,'k-','LineWidth',1.5); uistack(hl,'bottom');  % axis
  hl = plot([getelement(xlim,1) threshold(modelsets(1)).a],[1 1].*threshold(modelsets(2)).a,'k-.',...
            [1 1].*threshold(modelsets(1)).a,[getelement(ylim,1) threshold(modelsets(2)).a],'k-.');  uistack(hl,'bottom'); % threshold
  axis([max(min(xlim),-100),min(max(xlim),100), max(min(xlim),-100),min(max(xlim),100)],'square');
  title('alpha'); xlabel(modelNames(modelsets(1))); ylabel(modelNames(modelsets(2)));
hs = flipud(findobj(get(gca,'Children'),'Type','Scatter'));
hl=legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
set(findall(gcf,'-property','FontSize'),'FontSize',FntSiz);
title(hT,'Variance Explained (%)','FontSize',FntSiz);
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-compare_ch-%s_%s-rois',chmode,tarparam)));  end

%%
% close all;