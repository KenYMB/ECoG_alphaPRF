% Compare Gaussian Models, Channel Selection
%   Folk from ecog_APRF_07b_checkmodel_fullts
 
% 20210510 Yuasa
% 20240528 Yuasa - add individual plots
 
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
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth), 'modelselection', 'noACorr_bothgain');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');
 
%-- Plotting Setting
FntSiz    = 24;
SFntSiz   = 14;

%% general parameter
clear alphaType broadbandType

average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';
selectchs      = 'wangprobchs';
    allowlag       = false;
va_area = 'wangarea';

usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

%% load analyzePRF @DOG,OG,LFS
noACorrects = [false, true];
targetBANDs = {'FaCLb','FaLb'};
modelNames    = ["Model-based","Frequency-band"];
modelNames2l  = string({sprintf('Model-based\nAlpha pRF'),sprintf('Freq-band\nAlpha pRF')});
modelNamesSrt = ["Model-based","Freq-band"];

modeldata  = struct();
prf_params = struct();
model_all  = struct();
prf_all    = struct();
threshold  = struct();
for ii=length(noACorrects):-1:1
noACorrect   = noACorrects(ii);
if noACorrect
    allowbeta      = false;
    allowwide      = false;
    allowmixbeta   = false;
else
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;
end
 
clear alphaType broadbandType prfID modelID
INITload = 'bb';
ecog_APRFF_INITa_loaddata;
INITload = 'a'; prfID = 'prfbidirgain';         % allow negative gain for prf fitting
ecog_APRFF_INITa_loaddata;
INITload = [];
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold
 
%-- merge loaded data
% modeldata(ii).a   = modeldata_a;
% modeldata(ii).bb  = modeldata_bb;
prf_params(ii).a  = prf_params_a;
prf_params(ii).bb = prf_params_bb;
 
% model_all(ii).a   = model_all_a;
% model_all(ii).bb  = model_all_bb;
prf_all(ii).a     = prf_all_a;
prf_all(ii).bb    = prf_all_bb;

threshold(ii).a   = threshold_a;
threshold(ii).bb  = threshold_bb;

%-- flip gain for alpha (convert from ratio of suppression to power change)
% model_all(ii).a.datats = cellfun(@(C)-C,model_all(ii).a.datats,'UniformOutput',false);
prf_all(ii).a.gain     = -prf_all(ii).a.gain;
prf_all(ii).a.params   = prf_all(ii).a.params.*[1,1,1,-1,1,1,1];

end

%% %%%%%%%%%%%%%%
%% Visualize
%% %%%%%%%%%%%%%%
%%% prepare visualize
res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;
 
alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');
hF = gobjects(0);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FaCLb(BW) vs FaLb
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelsets = [1 2];      % FaCLb(BW) vs FaLb
targetBAND = strjoin(targetBANDs,'VS');

idat = 1;
okch1bb = ~(prf_all(modelsets(idat)).bb.xval<=threshold(modelsets(idat)).bb | prf_all(modelsets(idat)).bb.ecc >= eclimit);  % use tilde for nan
% okch1a  = ~(prf_all(modelsets(idat)).a.xval<=threshold(modelsets(idat)).a | prf_all(modelsets(idat)).a.ecc >= eclimit);  % use tilde for nan
okch1a  = ~(prf_all(modelsets(idat)).a.xval<=threshold(modelsets(idat)).a);  % use tilde for nan
idat = 2;
okch2bb = ~(prf_all(modelsets(idat)).bb.xval<=threshold(modelsets(idat)).bb | prf_all(modelsets(idat)).bb.ecc >= eclimit);  % use tilde for nan
% okch2a  = ~(prf_all(modelsets(idat)).a.xval<=threshold(modelsets(idat)).a | prf_all(modelsets(idat)).a.ecc >= eclimit);  % use tilde for nan
okch2a  = ~(prf_all(modelsets(idat)).a.xval<=threshold(modelsets(idat)).a);  % use tilde for nan

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Model Responses (as N-1 fold increase or decrease: 100% == 0 == no change)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Errorbar plot
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;
%% Errorbar plot of xval
if issaveplot

% nboot = 5000;
nboot = 1;

roilist = {{'V1','V2','V3'}};
roiname = {'V1–V3'};
conname = {'No Gain'};      % {'Neg Gain'},{'Pos Gain'}
nsel = length(roilist);
showelecnum = ' n=%.0f';

%-- plot parameter
shiftframe = 0.0;
if nsel==1,  shiftbar   = 0.25;
else,        shiftbar   = 0.15;
end
plcol_wang = [ones(1,1)*plcol(2,:);...
              ones(1,1)*plcol(7,:)];

nmdl = 2;
n_keep_roi     = repmat({{}},1,nmdl);
mn_keep_roi    = repmat({{}},1,nmdl);
for iroi = 1:length(roilist)
for imdl = 1:nmdl
roich  = ismember(prf_all(modelsets(imdl)).a.channels.wangarea,roilist{iroi});
% roicol = false(1,size(prf_all(modelsets(imdl)).a.channels,2));
% for ii = 1:length(roilist{iroi})
% roicol = roicol | (...
%          contains(prf_all(modelsets(imdl)).a.channels.Properties.VariableNames,'wang') &...
%          contains(prf_all(modelsets(imdl)).a.channels.Properties.VariableNames,roilist{iroi}));
% end
% roich  = sum(prf_all(modelsets(imdl)).a.channels{:,roicol},2)>0.05;

if startsWith(lower(conname),{'neg'})
    gainch = (prf_all(modelsets(imdl)).a.gain <= 0);
elseif startsWith(lower(conname),{'pos'})
    gainch = (prf_all(modelsets(imdl)).a.gain >  0);
else
    gainch = true(size(prf_all(modelsets(imdl)).a.gain));
end

okch   = (okch1bb&okch2bb)   & roich & gainch;
if sum(okch) < 5, okch(:) = false;  end

selch = find(okch)';
if sum(okch) > 4,   if nboot == 1 
                    else,          bootch = randi(sum(okch),nboot,sum(okch));
                                   selch = selch(bootch);
                    end
else,               selch  = nan(nboot,0);
end
n_keep_roi{imdl}{end+1}     = sum(diff(sort(selch,2),1,2)~=0,2)+1;
mn_keep_roi{imdl}{end+1}    = mean(prf_all(modelsets(imdl)).a.xval(selch),2,'omitnan');
end
end

%-- figure out
if nboot == 1, errlfunc = @(d)nanstd(d)./sqrt(sum(~isnan(d)));
               errufunc = @(d)nanstd(d)./sqrt(sum(~isnan(d)));
else,          errlfunc = @(d)nanmean(d)-prctile(d,alpha/2*100);
               errufunc = @(d)prctile(d,(1-alpha/2)*100)-nanmean(d);
end
hF(end+1) = figure;
set(gcf,'Position',get(gcf,'Position').*[.3 .5 0 0]+[0 0 nsel*180+50 420]);
colororder(plcol_wang);
hold on;
heb = gobjects(0);
for imdl = 1:nmdl
heb(end+1) = errorbar((1:nsel)+shiftbar.*((imdl-1)./(nmdl-1)*2-1),...
    cellfun(@nanmean,mn_keep_roi{imdl}),...
    cellfun(@(C) errlfunc(C),mn_keep_roi{imdl}),...
    cellfun(@(C) errufunc(C),mn_keep_roi{imdl}),...
    'o','LineStyle','none',...
    'LineWidth',2,'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
end
if nsel==1
    ylim([floor(min(ylim)/10)*10,ceil(max(ylim)/10)*10]);
    xticks(1+shiftbar.*(((1:nmdl)-1)./(nmdl-1)*2-1)); 
    xticklabels(strrep(modelNamesSrt,'-','-\newline')); xtickangle(0);
else
    ylim([floor(min(ylim)/10)*10,ceil((max(ylim)+0.4*diff(ylim))/10)*10]);
    xticks(1:nsel); xticklabels(roiname); xtickangle(0);
    legend(modelNamesSrt,'FontSize',FntSiz*0.8);
end
xlim([0.5,nsel+0.5]+shiftframe); if prod(ylim)<0, plot(xlim,[0,0],'k--'); end
set(gca,'FontSize',FntSiz);

% shiftframe = 0.05;
% for ii=1:nsel
% for imdl = 1:nmdl
%     mn_shift = 0; 
%     text(ii+0.05+shiftbar.*((imdl-1)./(nmdl-1)*2-1),...
%         nanmean(mn_keep_roi{imdl}{ii})+mn_shift, sprintf(showelecnum,nanmean(n_keep_roi{imdl}{ii})),...
%         'FontSize',FntSiz,'Color',plcol_wang(imdl,:));
% end
% end

figname =sprintf('xR2-%s%s_%s_errorbar',targetBAND,R2mode,selectchs);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Additional Plots in V1–V3 (classify pos & neg gain in Power Change in Alpha)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;

% ispltall = false;
ispltall = true;
mkSiz = 60;
%% Difference of Gaussian Model in training accuracy (scatter) w/ rough ROI
if issaveplot
npanel = 2;

% params = 'gain';
% params = 'xval';
% params = 'ecc';
% params = 'rfsize';
for params=["peakresp"]

plcol_wang = [ones(1,1)*plcol(6,:);...
              ones(1,1)*plcol(4,:)];
plshp_wang = [repmat('o',1,1) repmat('^',1,1)];
if ispltall
    clrate     = [1 9];
    plcol_gray = [plcol_wang .* clrate(1)  +  ones(1,3).*0.7 .* clrate(2)]./sum(clrate);
end
          
% hF(end+1) = figure('Menubar','none','Position',[200 200 860 420],'defaultAxesColorOrder',plcol_wang);
hF(end+1) = figure('Menubar','none','Position',[200 200 450*npanel 440],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(1,npanel,'Padding','compact','TileSpacing','compact');
roilist = {'V1','V2','V3'}';
nexttile;
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    otherwise,              ll = 1e3; ul = -1e3;  roundunit = 5;
end
    %--- get data range & set axis
    roich  = ismember(prf_all(modelsets(2)).a.channels.wangarea,roilist);
    okch = ((okch1a)|(okch2a)) & roich;
    datrange = prctile([prf_all(modelsets(2)).a.(params)(okch);prf_all(modelsets(1)).a.(params)(okch)],[15 85]);
    datrange = fliplr(datrange) + diff(datrange)/.7.*[-1 1];
    ll = min(ll,floor(min(datrange)/roundunit)*roundunit);
    ul = max(ul,ceil(max(datrange)/roundunit)*roundunit);
    switch params
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end
    axis([ll ul ll ul],'square');
    hold on;
    %--- plot scatter
    hs = gobjects(0);
    for ii=1:2
    scatterprop = {mkSiz,plshp_wang(ii),'LineWidth',1.6,'MarkerEdgeColor',plcol_wang(ii,:),'MarkerFaceColor',plcol_wang(ii,:)};
    if ii==1,    gainch = (prf_all(modelsets(2)).a.gain <= 0);
    else,        gainch = (prf_all(modelsets(2)).a.gain > 0);
    end
    okch   = ((okch1a)|(okch2a)) & roich & gainch;
    hs(end+1) = scatter(prf_all(modelsets(2)).a.(params)(okch),prf_all(modelsets(1)).a.(params)(okch),scatterprop{:});
    if ispltall
        okch   = roich & gainch;
        scatter(prf_all(modelsets(2)).a.(params)(okch),prf_all(modelsets(1)).a.(params)(okch),scatterprop{:},'MarkerEdgeColor',plcol_gray(ii,:),'MarkerFaceColor',plcol_gray(ii,:));
    end
    end
    uistack(hs,'top');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    if ismember(params,'xval')
    plot([getelement(xlim,1) threshold(modelsets(2)).a],[1 1].*threshold(modelsets(1)).a,'k-.',...
         [1 1].*threshold(modelsets(2)).a,[getelement(ylim,1) threshold(modelsets(1)).a],'k-.');  % threshold
%     legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
    end
    
    set(gca,'FontSize',FntSiz);
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
    title(params);
    xticks(yticks); xtickangle(0);
    
if npanel > 1
nexttile; hold on;

    %--- set axis
    axis([ll ul ll ul],'square');
    hold on;
    %--- plot scatter
    hs = gobjects(0);
    for ii=1:2
    scatterprop = {mkSiz,plshp_wang(ii),'LineWidth',1.6,'MarkerEdgeColor',plcol_wang(ii,:),'MarkerFaceColor',plcol_wang(ii,:)};
    if ii==1,    gainch = (prf_all(modelsets(2)).a.gain <= 0);
    else,        gainch = (prf_all(modelsets(2)).a.gain > 0);
    end
    okch   = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich & gainch;
    hs(end+1) = scatter(prf_all(modelsets(2)).a.(params)(okch),prf_all(modelsets(1)).a.(params)(okch),scatterprop{:});
    if ispltall
        okch   = (okch1bb&okch2bb) & roich & gainch;
        scatter(prf_all(modelsets(2)).a.(params)(okch),prf_all(modelsets(1)).a.(params)(okch),scatterprop{:},'MarkerEdgeColor',plcol_gray(ii,:),'MarkerFaceColor',plcol_gray(ii,:));
    end
    end
    uistack(hs,'top');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    if ismember(params,'xval')
    plot([getelement(xlim,1) threshold(modelsets(2)).a],[1 1].*threshold(modelsets(1)).a,'k-.',...
         [1 1].*threshold(modelsets(2)).a,[getelement(ylim,1) threshold(modelsets(1)).a],'k-.');  % threshold
%     legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
    end
    
    set(gca,'FontSize',FntSiz);
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
    title([params ' w/ Broadband threshold']);
    xticks(yticks); xtickangle(0);
end
 
figname =sprintf('Scatter-%s_%s%s_%s_ROI-low',params,targetBAND,R2mode,selectchs);
if ispltall,     figname = [figname '_withrejected'];   end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

end
end

%% Cross Relations
if issaveplot
npanel = 2;

params=["gain","peakresp"];

plcol_wang = [ones(1,1)*plcol(6,:);...
              ones(1,1)*plcol(4,:)];
plshp_wang = [repmat('o',1,1) repmat('^',1,1)];
if ispltall
    clrate     = [1 9];
    plcol_gray = [plcol_wang .* clrate(1)  +  ones(1,3).*0.7 .* clrate(2)]./sum(clrate);
end
          
% hF(end+1) = figure('Menubar','none','Position',[200 200 860 420],'defaultAxesColorOrder',plcol_wang);
hF(end+1) = figure('Menubar','none','Position',[200 200 450*npanel 440*2],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(2,npanel,'Padding','compact','TileSpacing','compact');
roilist = {'V1','V2','V3'}';

for mdlidx=1:2
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    otherwise,              ll = 1e3; ul = -1e3;  roundunit = 5;
end
    %--- get data range
    roich  = ismember(prf_all(modelsets(mdlidx)).a.channels.wangarea,roilist);
    okch = ((okch1a)|(okch2a)) & roich;
    datrange = prctile([prf_all(modelsets(mdlidx)).a.(params(2))(okch);prf_all(modelsets(mdlidx)).a.(params(1))(okch)],[15 85]);
    datrange = fliplr(datrange) + diff(datrange)/.7.*[-1 1];
    ll = min(ll,floor(min(datrange)/roundunit)*roundunit);
    ul = max(ul,ceil(max(datrange)/roundunit)*roundunit);
    switch params(2)
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end

nexttile;
    %--- set axis
    axis([ll ul ll ul],'square');
    hold on;
    %--- plot scatter
    hs = gobjects(0);
    for ii=1:2
    scatterprop = {mkSiz,plshp_wang(ii),'LineWidth',1.6,'MarkerEdgeColor',plcol_wang(ii,:),'MarkerFaceColor',plcol_wang(ii,:)};
    if ii==1,    gainch = (prf_all(modelsets(2)).a.gain <= 0);
    else,        gainch = (prf_all(modelsets(2)).a.gain > 0);
    end
    okch   = ((okch1a)|(okch2a)) & roich & gainch;
    hs(end+1) = scatter(prf_all(modelsets(mdlidx)).a.(params(2))(okch),prf_all(modelsets(mdlidx)).a.(params(1))(okch),scatterprop{:});
    if ispltall
        okch   = roich & gainch;
        scatter(prf_all(modelsets(mdlidx)).a.(params(2))(okch),prf_all(modelsets(mdlidx)).a.(params(1))(okch),scatterprop{:},'MarkerEdgeColor',plcol_gray(ii,:),'MarkerFaceColor',plcol_gray(ii,:));
    end
    end
    uistack(hs,'top');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    
    set(gca,'FontSize',FntSiz);
    xlabel(params(2)); ylabel(params(1));
    title(modelNames(modelsets(mdlidx)));
    xticks(yticks); xtickangle(0);
    
if npanel > 1
nexttile;
    %--- set axis
    axis([ll ul ll ul],'square');
    hold on;
    %--- plot scatter
    hs = gobjects(0);
    for ii=1:2
    scatterprop = {mkSiz,plshp_wang(ii),'LineWidth',1.6,'MarkerEdgeColor',plcol_wang(ii,:),'MarkerFaceColor',plcol_wang(ii,:)};
    if ii==1,    gainch = (prf_all(modelsets(2)).a.gain <= 0);
    else,        gainch = (prf_all(modelsets(2)).a.gain > 0);
    end
    okch   = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich & gainch;
    hs(end+1) = scatter(prf_all(modelsets(mdlidx)).a.(params(2))(okch),prf_all(modelsets(mdlidx)).a.(params(1))(okch),scatterprop{:});
    if ispltall
        okch   = (okch1bb&okch2bb) & roich & gainch;
        scatter(prf_all(modelsets(mdlidx)).a.(params(2))(okch),prf_all(modelsets(mdlidx)).a.(params(1))(okch),scatterprop{:},'MarkerEdgeColor',plcol_gray(ii,:),'MarkerFaceColor',plcol_gray(ii,:));
    end
    end
    uistack(hs,'top');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    
    set(gca,'FontSize',FntSiz);
    xlabel(params(2)); ylabel(params(1));
    title([modelNames(modelsets(mdlidx)) ' w/ Broadband threshold']);
    xticks(yticks); xtickangle(0);
end
end
 
figname =sprintf('xScatter-%s-%s_%s%s_%s_ROI-low',params(1),params(2),targetBAND,R2mode,selectchs);
if ispltall,     figname = [figname '_withrejected'];   end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Histogram of Peak response
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if issaveplot
npanel = 2;

if issaveplot,  plotEach = false;
else,           plotEach = true;
end

nBins = 36;
params = 'peakresp';

plcol_wang = [ones(1,1)*plcol(2,:);...
              ones(1,1)*plcol(7,:);];
          
% hF(end+1) = figure('Menubar','none','Position',[200 200 860 420],'defaultAxesColorOrder',plcol_wang);
hF(end+1) = figure('Menubar','none','Position',[200 200 600 350*npanel+25],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(npanel,1,'Padding','compact','TileSpacing','compact');
roilist = {'V1','V2','V3'}';
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    otherwise,              ll = 1e3; ul = -1e3;  roundunit = 5;
end
    %--- get data range
    roich  = ismember(prf_all(modelsets(2)).a.channels.wangarea,roilist);
    okch = ((okch1a)|(okch2a)) & roich;
    datrange = prctile([prf_all(modelsets(2)).a.(params)(okch);prf_all(modelsets(1)).a.(params)(okch)],[15 85]);
    datrange = fliplr(datrange) + diff(datrange)/.7.*[-1 1];
    ll = min(ll,floor(min(datrange)/roundunit)*roundunit);
    ul = max(ul,ceil(max(datrange)/roundunit)*roundunit);
    switch params
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end

if npanel > 1
nexttile;
    %--- set axis
    xlim([ll ul]);
    hold on;
    %--- plot histogram
    hs = gobjects(0);
    if plotEach, mdlidx = 1;
    else,        mdlidx = 1:2;
    end
    for modelidx=mdlidx
      if plotEach
        if ispltall
            okch   = (okch1bb&okch2bb) & roich;
        else
            okch   = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich;
        end
      else
        if ispltall
            okch   = roich;
        else
            okch   = ((okch1a)|(okch2a)) & roich;
        end
      end
    hs(end+1) = histogram(prf_all(modelsets(modelidx)).a.(params)(okch),'BinLimits',xlim,'NumBins',nBins);
    end
    legend(modelNames2l(mdlidx),'AutoUpdate','off');
    set(gca,'FontSize',FntSiz);
    plot([0 0],ylim,'k-','LineWidth',4);  % axis
    if plotEach
        title([params ' w/ Broadband threshold']);
    else
        title(params);
    end
end
    
nexttile;
    %--- set axis
    xlim([ll ul]);
    hold on;
    %--- plot histogram
    hs = gobjects(0);
    if plotEach, mdlidx = 2;     colororder(gca,flipud(plcol_wang)); 
    else,        mdlidx = 1:2;
    end
    for modelidx=mdlidx
        if ispltall
            okch   = (okch1bb&okch2bb) & roich;
        else
            okch   = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich;
        end
    hs(end+1) = histogram(prf_all(modelsets(modelidx)).a.(params)(okch),'BinLimits',xlim,'NumBins',nBins);
    end
    if plotEach
        legend(modelNames2l(mdlidx),'AutoUpdate','off');
    else
        title([params ' w/ Broadband threshold']);
    end
    set(gca,'FontSize',FntSiz);
    plot([0 0],ylim,'k-','LineWidth',4);  % axis
    xlabel(ht,  'Peak Response','FontSize',FntSiz*1.2);
    ylabel(ht,  '# of electrodes','FontSize',FntSiz*1.2);

figname =sprintf('Hist-%s_%s%s_%s_ROI-low',params,targetBAND,R2mode,selectchs);
if ispltall,     figname = [figname '_withrejected'];   end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Histogram (Violin) of R2 / Gain &
%       Ecc vs Size with fitting
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if issaveplot
npanel = 2;
colpanel = 3;
nBins = 28;

doviolin = [0 0];
% doviolin = [1 1];

plcol_wang = [ones(1,1)*plcol(2,:);...
              ones(1,1)*plcol(7,:);];
plshp_wang = [repmat('o',1,1) repmat('o',1,1)];
          
% hF(end+1) = figure('Menubar','none','Position',[200 200 600*colpanel 350*npanel+25],'defaultAxesColorOrder',plcol_wang);
hF(end+1) = figure('Menubar','none','Position',[200 200 400*colpanel+200 400*npanel],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(1,colpanel,'Padding','compact','TileSpacing','loose');
htIn = gobjects(0);
roilist = {'V1','V2','V3'}';


%%% R2 Histogram / Violin
params = 'xval';
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    case {'peakresp'},      ll = 10;  ul = -10;     roundunit = 3;
    otherwise,              ll = 1e3; ul = -1e3;    roundunit = 5;
end
    %--- get data range
    rangethresh = 70;
    roich  = ismember(prf_all(modelsets(2)).a.channels.wangarea,roilist);
    if ispltall,    okch   = (okch1bb&okch2bb) & roich;
    else,           okch   = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich;
    end
    if doviolin(length(htIn)+1)<0
    errfunc = @(d)std(d,'omitnan')./sqrt(sum(~isnan(d)));
    ll = min(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    ul = max(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    else
    datrange = prctile([prf_all(modelsets(2)).a.(params)(okch);prf_all(modelsets(1)).a.(params)(okch)],50+[-0.5 0.5].*rangethresh);
    datrange = fliplr(datrange) + diff(datrange)/rangethresh*100.*[-1 1];
    ll = min(ll,min(datrange));
    ul = max(ul,max(datrange));
    end
    ll = floor(ll/roundunit)*roundunit;
    ul = ceil(ul/roundunit)*roundunit;
    switch params
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end

htIn(end+1) = tiledlayout(ht,npanel,1,'Padding','tight','TileSpacing','compact');
htIn(end).Layout.Tile=length(htIn);
for mdlidx=1:2
nexttile(htIn(end));    hold on;
    %--- set axis
    hl = gobjects(0);
    if doviolin(length(htIn))>0,       ylim([ll ul]);
    else,                              xlim([ll ul]);
    end
    %--- plot histogram
    if doviolin(length(htIn))>0
        colororder(plcol_wang(mdlidx,:)); 
        violinplot(prf_all(modelsets(mdlidx)).a.(params)(okch));
    elseif doviolin(length(htIn))<0
        colororder(gca,plcol_wang(mdlidx,:)); 
        errorbar(mean(prf_all(modelsets(mdlidx)).a.(params)(okch),'omitnan'),1,...
                             errfunc(prf_all(modelsets(mdlidx)).a.(params)(okch)),...
                            'horizontal','o','LineStyle','none','LineWidth',2,...
                            'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
    else
        colororder(gca,plcol_wang(mdlidx,:)); 
        histogram(prf_all(modelsets(mdlidx)).a.(params)(okch),'BinLimits',xlim,'NumBins',nBins);
    end
    set(gca,'FontSize',FntSiz);
%     legend(modelNames2l(mdlidx),'AutoUpdate','off');
if doviolin(length(htIn))>0
    view([90 90]); xticks([]);
    if prod(ylim)<0, hl(end+1)=plot(xlim,[0,0],'k--','LineWidth',2); end  % axis
    ylabel(gca,  'Cross-validated R^2 (%)','FontSize',FntSiz*1.2,'VerticalAlignment','cap');
elseif doviolin(length(htIn))<0
    yticks([]);
    if prod(xlim)<0, hl(end+1)=plot(ylim,[0,0],'k--','LineWidth',2); end  % axis
    xlabel(gca,  'Cross-validated R^2 (%)','FontSize',FntSiz*1.2,'VerticalAlignment','cap');
else
    ylim(ylim+[0 0.2]);
    if prod(xlim)<0, hl(end+1)=plot([0,0], ylim,'k-','LineWidth',4); end  % axis
    xlabel(gca,  'Cross-validated R^2 (%)','FontSize',FntSiz*1.2,'VerticalAlignment','cap');
    ylabel(gca,  '# of electrodes','FontSize',FntSiz*1.2);
end
uistack(hl,'bottom');

end
title(htIn(end),  'Variance Explained','FontSize',FntSiz*1.1,'FontWeight','bold');
    

%%% Gain Histogram / Violin
params = 'peakresp';
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    case {'peakresp'},      ll = 10;  ul = -10;     roundunit = 3;
    otherwise,              ll = 1e3; ul = -1e3;    roundunit = 5;
end
    %--- get data range
    rangethresh = 95;
    roich  = ismember(prf_all(modelsets(2)).a.channels.wangarea,roilist);
    if ispltall,    okch   = (okch1bb&okch2bb) & roich;
    else,           okch   = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich;
    end
    if doviolin(length(htIn)+1)<0
    errfunc = @(d)std(d,'omitnan')./sqrt(sum(~isnan(d)));
    ll = min(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    ul = max(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    else
    datrange = prctile([prf_all(modelsets(2)).a.(params)(okch);prf_all(modelsets(1)).a.(params)(okch)],50+[-0.5 0.5].*rangethresh);
    datrange = fliplr(datrange) + diff(datrange)/rangethresh*100.*[-1 1];
    ll = min(ll,min(datrange));
    ul = max(ul,max(datrange));
    end
    ll = floor(ll/roundunit)*roundunit;
    ul = ceil(ul/roundunit)*roundunit;
    switch params
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end

htIn(end+1) = tiledlayout(ht,npanel,1,'Padding','tight','TileSpacing','compact');
htIn(end).Layout.Tile=length(htIn);
for mdlidx=1:2
nexttile(htIn(end));    hold on;
    %--- set axis
    hl = gobjects(0);
    if doviolin(length(htIn))>0,       ylim([ll ul]);
    else,                              xlim([ll ul]);
    end
    %--- plot histogram
    if doviolin(length(htIn))>0
        colororder(plcol_wang(mdlidx,:)); 
        violinplot(prf_all(modelsets(mdlidx)).a.(params)(okch));
    elseif doviolin(length(htIn))<0
        colororder(gca,plcol_wang(mdlidx,:)); 
        errorbar(mean(prf_all(modelsets(mdlidx)).a.(params)(okch),'omitnan'),1,...
                             errfunc(prf_all(modelsets(mdlidx)).a.(params)(okch)),...
                            'horizontal','o','LineStyle','none','LineWidth',2,...
                            'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
    else
        colororder(gca,plcol_wang(mdlidx,:)); 
        histogram(prf_all(modelsets(mdlidx)).a.(params)(okch),'BinLimits',xlim,'NumBins',nBins);
    end
    set(gca,'FontSize',FntSiz);
    legend(modelNames2l(mdlidx),'AutoUpdate','off');
if doviolin(length(htIn))>0
    view([90 90]); xticks([]);
    if prod(ylim)<0, hl(end+1)=plot(xlim,[0,0],'k--','LineWidth',2); end  % axis
    ylabel(gca,  'Ratio of Response','FontSize',FntSiz*1.2,'VerticalAlignment','cap');
elseif doviolin(length(htIn))<0
    yticks([]);
    if prod(xlim)<0, hl(end+1)=plot(ylim,[0,0],'k--','LineWidth',2); end  % axis
    xlabel(gca,  'Ratio of Response','FontSize',FntSiz*1.2,'VerticalAlignment','cap');
else
    ylim(ylim+[0 0.2]);
    if prod(xlim)<0, hl(end+1)=plot([0,0], ylim,'k-','LineWidth',4); end  % axis
    xlabel(gca,  'Ratio of Response','FontSize',FntSiz*1.2,'VerticalAlignment','cap');
    ylabel(gca,  '# of electrodes','FontSize',FntSiz*1.2);
end
%-- Update TickLabel
if doviolin(length(htIn))>0
    if length(yticks)<5
        yticks(linspace(floor(min(ylim)),ceil(max(ylim)),5));
    end
    yticklabels(yticks+sign(yticks)+(yticks==0));
    ticknames = yticklabels;
    ticknames(yticks<0)  = cellfun(@(d) sprintf('1/%d',-str2double(d)),ticknames(yticks<0),'UniformOutput',false);
    yticklabels(ticknames);
else
    if length(xticks)<5
        xticks(linspace(floor(min(xlim)),ceil(max(xlim)),5));
    end
    xticklabels(xticks+sign(xticks)+(xticks==0));
    ticknames = xticklabels;
    ticknames(xticks<0)  = cellfun(@(d) sprintf('1/%d',-str2double(d)),ticknames(xticks<0),'UniformOutput',false);
    xticklabels(ticknames);
end
uistack(hl,'bottom');

end
title(htIn(end), 'Response Gain','FontSize',FntSiz*1.1,'FontWeight','bold');
    

%%% Ecc vs Size
params=["rfsize","ecc"];

%--- get data range
roich  = ismember(prf_all(modelsets(1)).a.channels.wangarea,roilist);
okch   = (okch1bb&okch2bb) & roich;
if ~ispltall,    okch   = ((okch1a)|(okch2a)) & okch;    end

datrange = [0 prctile([prf_all(modelsets(1)).a.(params(2))(okch);prf_all(modelsets(2)).a.(params(2))(okch)],98)*cfactor];
datrange = [0 ceil(diff(datrange)/.98*1.05)];
datrange = [datrange datrange.*[1 10/8.5]];    % in degree [Ecc, Size]

htIn(end+1) = tiledlayout(ht,npanel,1,'Padding','tight','TileSpacing','compact');
htIn(end).Layout.Tile=length(htIn);
for mdlidx=1:2
nexttile(htIn(end));    hold on;
    %--- plot scatter
    hs = gobjects(0);
    for ii=1:2
    scatterprop = {mkSiz,plshp_wang(ii),'LineWidth',2.0,'MarkerEdgeColor',plcol_wang(mdlidx,:)};
    if ii==1,    gainch = (prf_all(modelsets(mdlidx)).a.gain <= 0);  scatterprop = [scatterprop, {'MarkerFaceColor',plcol_wang(mdlidx,:)}];
    else,        gainch = (prf_all(modelsets(mdlidx)).a.gain > 0);
    end
    okch   = (okch1bb&okch2bb) & roich & gainch;
    if ~ispltall,    okch   = ((okch1a)|(okch2a)) & okch;    end
    dat1 = prf_all(modelsets(mdlidx)).a.(params(2))(okch).*cfactor;
    dat2 = prf_all(modelsets(mdlidx)).a.(params(1))(okch).*cfactor;
    hs(end+1) = scatter(dat1,dat2,scatterprop{:});
    end
    %-- fit line
    ii=1;
    okch   = (okch1bb&okch2bb) & roich;
    if ~ispltall,    okch   = ((okch1a)|(okch2a)) & okch;    end
    nchan = sum(okch);
    if nchan>3
    lineprop = {'Color',plcol_wang(mdlidx,:),'LineWidth',2.0};
    if ii==1, lineprop = [{'--'},lineprop];
    else,     lineprop = [{':'},lineprop];
    end
    warnstat = warning; warning('off');
    dat1 = prf_all(modelsets(mdlidx)).a.(params(2))(okch).*cfactor;
    dat2 = prf_all(modelsets(mdlidx)).a.(params(1))(okch).*cfactor;
    
    C = cov(dat1,dat2);
    sl = C(1,2)./C(1,1); intc = median(dat2) - sl*median(dat1);
    fitcoef = [intc;sl];
    warning(warnstat);
    datX = minmax([xlim, datrange(1:2)]); datX = min(datX):range(datX)/100:max(datX);
    datY = fitcoef(1,:)'+datX.*fitcoef(2,:)'; 
    plot(datX,mean(datY,1,'omitnan'),lineprop{:});
    end
    %-- set axis
    uistack(hs,'top');
%     plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    axis(datrange);
    set(gca,'FontSize',FntSiz);

xlabel(gca,  'Eccentricity (deg)','FontSize',FntSiz*1.2,'VerticalAlignment','cap');
ylabel(gca,  'Size (deg)','FontSize',FntSiz*1.2);
end
legend(hs,["Negative Gain","Positive Gain"],'AutoUpdate','off');
title(htIn(end),  'Eccentricity VS Size','FontSize',FntSiz*1.1,'FontWeight','bold');


figname =sprintf('Summary-%s_%s%s_%s_ROI-low',params,targetBAND,R2mode,selectchs);
if ispltall,     figname = [figname '_withrejected'];   end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');   end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Histogram (Violin) of R2 / Gain &
%       Ecc vs Size with fitting <Simple Ver>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 % ~issaveplot
npanel = 1;
colpanel = [5 5 6];
nBins = 28;

doviolin = [1 1];       % for R2, gain
% doviolin = [-1 -1];

showerror = true;       % for ecc vs size
nboot     = 5000;       % if showerror

colpanel(doviolin<0) = colpanel(doviolin<0) - 1;

plcol_wang = [ones(1,1)*plcol(2,:);...
              ones(1,1)*plcol(4,:);];
plshp_wang = [repmat('o',1,1) repmat('o',1,1)];
          
% hF(end+1) = figure('Menubar','none','Position',[200 200 400*length(colpanel)+200 400*npanel],'defaultAxesColorOrder',plcol_wang);
hF(end+1) = figure('Menubar','none','Position',[200 200 round(400*sum(colpanel)./colpanel(end))+200 400*npanel],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(1,sum(colpanel),'Padding','compact','TileSpacing','loose');
htIn = gobjects(0);
roilist = {'V1','V2','V3'}';
nmdl = 2;


%%% R2 Histogram / Violin
params = 'xval';
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    case {'peakresp'},      ll = 10;  ul = -10;     roundunit = 3;
    otherwise,              ll = 1e3; ul = -1e3;    roundunit = 5;
end
    %--- get data range
    rangethresh = 70;
    roich  = ismember(prf_all(modelsets(2)).a.channels.wangarea,roilist);
    if ispltall,    okch   = (okch1bb&okch2bb) & roich;
    else,           okch   = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich;
    end
    if doviolin(length(htIn)+1)<0
    errfunc = @(d)std(d,'omitnan')./sqrt(sum(~isnan(d)));
    ll = min(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    ul = max(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    else
    datrange = prctile([prf_all(modelsets(2)).a.(params)(okch);prf_all(modelsets(1)).a.(params)(okch)],50+[-0.5 0.5].*rangethresh);
    datrange = fliplr(datrange) + diff(datrange)/rangethresh*100.*[-1 1];
    ll = min(ll,min(datrange));
    ul = max(ul,max(datrange));
    end
    ll = floor(ll/roundunit)*roundunit;
    ul = ceil(ul/roundunit)*roundunit;
    switch params
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end

htIn(end+1) = tiledlayout(ht,npanel,1,'Padding','tight','TileSpacing','compact');
htIn(end).Layout.Tile=sum(colpanel(1:(length(htIn)-1)))+1;
htIn(end).Layout.TileSpan=[1,colpanel(length(htIn))];
nexttile(htIn(end));    hold on;
    %--- set axis
    hl = gobjects(0);
    if doviolin(length(htIn)),       ylim([ll ul]); xlim([1 nmdl]+[-1 1]*0.6)
    else,                            xlim([ll ul]);
    end
% for mdlidx=1:2
    %--- plot histogram
    if doviolin(length(htIn))>0
        colororder(plcol_wang); 
        violinplot([prf_all(modelsets(1)).a.(params)(okch),...
                    prf_all(modelsets(2)).a.(params)(okch)]);
    elseif doviolin(length(htIn))<0
        colororder(gca,plcol_wang); 
        imdl = 1; errorbar(imdl,...
                 [mean(prf_all(modelsets(imdl)).a.(params)(okch),'omitnan')],...
                 [errfunc(prf_all(modelsets(imdl)).a.(params)(okch))],...
                            'vertical','o','LineStyle','none','LineWidth',2,...
                            'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
        imdl = 2; errorbar(imdl,...
                 [mean(prf_all(modelsets(imdl)).a.(params)(okch),'omitnan')],...
                 [errfunc(prf_all(modelsets(imdl)).a.(params)(okch))],...
                            'vertical','o','LineStyle','none','LineWidth',2,...
                            'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
    else
        colororder(gca,plcol_wang); 
        histogram(prf_all(modelsets(mdlidx)).a.(params)(okch),'BinLimits',xlim,'NumBins',nBins);
    end
% end
    set(gca,'FontSize',FntSiz);
%     legend(modelNames2l,'AutoUpdate','off');
if doviolin(length(htIn))
%     xticks([]);
    xticklabels(strtok(modelNamesSrt,'-'));
    if prod(ylim)<0, hl(end+1)=plot(xlim,[0,0],'k--','LineWidth',2); end  % axis
    ylabel(gca,  'Cross-validated R^2 (%)','FontSize',FntSiz*1.2,'VerticalAlignment','bottom');
    ylim(ylim+[0 1].*range(ylim)./20);
else
    ylim(ylim+[0 0.2]);
    if prod(xlim)<0, hl(end+1)=plot([0,0], ylim,'k-','LineWidth',4); end  % axis
    xlabel(gca,  'Cross-validated R^2 (%)','FontSize',FntSiz*1.2,'VerticalAlignment','bottom');
    ylabel(gca,  '# of electrodes','FontSize',FntSiz*1.2);
end
uistack(hl,'bottom');
title(htIn(end),  'Variance Explained','FontSize',FntSiz*1.1,'FontWeight','bold');
    


%%% Gain Histogram / Violin
params = 'peakresp';
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    case {'peakresp'},      ll = 10;  ul = -10;     roundunit = 3;
                      if doviolin(length(htIn))<0,  roundunit = 2;  end
    otherwise,              ll = 1e3; ul = -1e3;    roundunit = 5;
end
    %--- get data range
    rangethresh = 95;
    roich  = ismember(prf_all(modelsets(2)).a.channels.wangarea,roilist);
    if ispltall,    okch   = (okch1bb&okch2bb) & roich;
    else,           okch   = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich;
    end
    if doviolin(length(htIn)+1)<0
    errfunc = @(d)std(d,'omitnan')./sqrt(sum(~isnan(d)));
    ll = min(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    ul = max(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    else
    datrange = prctile([prf_all(modelsets(2)).a.(params)(okch);prf_all(modelsets(1)).a.(params)(okch)],50+[-0.5 0.5].*rangethresh);
    datrange = fliplr(datrange) + diff(datrange)/rangethresh*100.*[-1 1];
    ll = min(ll,min(datrange));
    ul = max(ul,max(datrange));
    end
    ll = floor(ll/roundunit)*roundunit;
    ul = ceil(ul/roundunit)*roundunit;
    switch params
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end

htIn(end+1) = tiledlayout(ht,npanel,1,'Padding','tight','TileSpacing','compact');
htIn(end).Layout.Tile=sum(colpanel(1:(length(htIn)-1)))+1;
htIn(end).Layout.TileSpan=[1,colpanel(length(htIn))];
nexttile(htIn(end));    hold on;
    %--- set axis
    hl = gobjects(0);
    if doviolin(length(htIn)),       ylim([ll ul]); xlim([1 nmdl]+[-1 1]*0.6)
    else,                            xlim([ll ul]);
    end
% for mdlidx=1:2
    %--- plot histogram
    if doviolin(length(htIn))>0
        colororder(plcol_wang); 
        violinplot([prf_all(modelsets(1)).a.(params)(okch),...
                    prf_all(modelsets(2)).a.(params)(okch)]);
    elseif doviolin(length(htIn))<0
        colororder(gca,plcol_wang); 
        imdl = 1; errorbar(imdl,...
                 [mean(prf_all(modelsets(imdl)).a.(params)(okch),'omitnan')],...
                 [errfunc(prf_all(modelsets(imdl)).a.(params)(okch))],...
                            'vertical','o','LineStyle','none','LineWidth',2,...
                            'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
        imdl = 2; errorbar(imdl,...
                 [mean(prf_all(modelsets(imdl)).a.(params)(okch),'omitnan')],...
                 [errfunc(prf_all(modelsets(imdl)).a.(params)(okch))],...
                            'vertical','o','LineStyle','none','LineWidth',2,...
                            'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
    else
        colororder(gca,plcol_wang); 
        histogram(prf_all(modelsets(mdlidx)).a.(params)(okch),'BinLimits',xlim,'NumBins',nBins);
    end
% end
    set(gca,'FontSize',FntSiz);
%     legend(modelNamesSrt,'AutoUpdate','off');
if doviolin(length(htIn))
%     xticks([]);
    xticklabels(strtok(modelNamesSrt,'-'));
    if prod(ylim)<0, hl(end+1)=plot(xlim,[0,0],'k--','LineWidth',2); end  % axis
    ylabel(gca,  'Ratio of Response','FontSize',FntSiz*1.2,'VerticalAlignment','bottom');
    ylim(ylim+[-1 1].*range(ylim)./20);
else
    ylim(ylim+[0 0.2]);
    if prod(xlim)<0, hl(end+1)=plot([0,0], ylim,'k-','LineWidth',4); end  % axis
    xlabel(gca,  'Ratio of Response','FontSize',FntSiz*1.2,'VerticalAlignment','bottom');
    ylabel(gca,  '# of electrodes','FontSize',FntSiz*1.2);
end
%-- Update TickLabel
if doviolin(length(htIn))
    if length(yticks)<5
        yticks(linspace(floor(min(ylim)),ceil(max(ylim)),5));
    end
    yticklabels(yticks+sign(yticks)+(yticks==0));
    ticknames = yticklabels;
    ticknames(yticks<0)  = cellfun(@(d) sprintf('1/%d',-str2double(d)),ticknames(yticks<0),'UniformOutput',false);
    yticklabels(ticknames);
else
    if length(xticks)<5
        xticks(linspace(floor(min(xlim)),ceil(max(xlim)),5));
    end
    xticklabels(xticks+sign(xticks)+(xticks==0));
    ticknames = xticklabels;
    ticknames(xticks<0)  = cellfun(@(d) sprintf('1/%d',-str2double(d)),ticknames(xticks<0),'UniformOutput',false);
    xticklabels(ticknames);
end
uistack(hl,'bottom');
title(htIn(end), 'Response Gain','FontSize',FntSiz*1.1,'FontWeight','bold');
    

%%% Ecc vs Size
params=["rfsize","ecc"];

%--- get data range
roich  = ismember(prf_all(modelsets(1)).a.channels.wangarea,roilist);
okch   = (okch1bb&okch2bb) & roich;
if ~ispltall,    okch   = ((okch1a)|(okch2a)) & okch;    end

datrange = [0 prctile([prf_all(modelsets(1)).a.(params(2))(okch);prf_all(modelsets(2)).a.(params(2))(okch)],98)*cfactor];
datrange = [0 ceil(diff(datrange)/.98*1.05)];
datrange = [datrange datrange.*[1 10/8.5]];    % in degree [Ecc, Size]

htIn(end+1) = tiledlayout(ht,npanel,1,'Padding','tight','TileSpacing','compact');
htIn(end).Layout.Tile=sum(colpanel(1:(length(htIn)-1)))+1;
htIn(end).Layout.TileSpan=[1,colpanel(length(htIn))];
nexttile(htIn(end));    hold on;
hs = gobjects(0);
for mdlidx=1:2
    %--- plot scatter (for each gain)
    for ii=1:2
    scatterprop = {mkSiz,plshp_wang(ii),'LineWidth',2.0,'MarkerEdgeColor',plcol_wang(mdlidx,:)};
    if ii==1,    gainch = (prf_all(modelsets(mdlidx)).a.gain <= 0);  scatterprop = [scatterprop, {'MarkerFaceColor',plcol_wang(mdlidx,:)}];
    else,        gainch = (prf_all(modelsets(mdlidx)).a.gain > 0);
    end
    okch   = (okch1bb&okch2bb) & roich & gainch;
    if ~ispltall,    okch   = ((okch1a)|(okch2a)) & okch;    end
    dat1 = prf_all(modelsets(mdlidx)).a.(params(2))(okch).*cfactor;
    dat2 = prf_all(modelsets(mdlidx)).a.(params(1))(okch).*cfactor;
    hs(end+1) = scatter(dat1,dat2,scatterprop{:});
    end
    %-- fit line (Robust method)
    %%-- collect data (ignore gains)
    ii=1;
    okch   = (okch1bb&okch2bb) & roich;
    if ~ispltall,    okch   = ((okch1a)|(okch2a)) & okch;    end
    nchan = sum(okch);
    if nchan>3
    lineprop = {'Color',plcol_wang(mdlidx,:),'LineWidth',2.0};
    if ii==1, lineprop = [{'--'},lineprop];
    else,     lineprop = [{':'},lineprop];
    end
    warnstat = warning; warning('off');
    dat1 = prf_all(modelsets(mdlidx)).a.(params(2))(okch).*cfactor;
    dat2 = prf_all(modelsets(mdlidx)).a.(params(1))(okch).*cfactor;
    %%-- exclude outliers (exclude data which are out of more than 2 sigma on the major and minor axes)
    C = cov(dat1,dat2);  sl = C(1,2)./C(1,1); 
    rot2d = @(D,deg) D*[cos(deg),-sin(deg);sin(deg),cos(deg)];
    rotD = rot2d([dat1,dat2],atan(sl));
    dat1r = rotD(:,1); dat2r = rotD(:,2);
    okelec = dat1r < mean(dat1r)+std(dat1r)*2 & dat2r < mean(dat2r)+std(dat2r)*2 &...
             dat1r > mean(dat1r)-std(dat1r)*2 & dat2r > mean(dat2r)-std(dat2r)*2;
    dat1(~okelec) = [];    dat2(~okelec) = [];
    %%-- compute slope of the major axis
    C = cov(dat1,dat2);
    sl = C(1,2)./C(1,1); intc = median(dat2) - sl*median(dat1);
    fitcoef = [intc;sl];
        if showerror
            fitcoefbooot = nan(2,nboot);
            for iboot=1:nboot
                okchboot = randsample(find(okch),sum(okch),true);
                dat1 = prf_all(modelsets(mdlidx)).a.(params(2))(okchboot).*cfactor;
                dat2 = prf_all(modelsets(mdlidx)).a.(params(1))(okchboot).*cfactor;
                %%-- exclude outliers (exclude data which are out of more than 2 sigma on the major and minor axes)
                C = cov(dat1,dat2);  sl = C(1,2)./C(1,1); 
                rotD = rot2d([dat1,dat2],atan(sl));
                dat1r = rotD(:,1); dat2r = rotD(:,2);
                okelec = dat1r < mean(dat1r)+std(dat1r)*2 & dat2r < mean(dat2r)+std(dat2r)*2 &...
                         dat1r > mean(dat1r)-std(dat1r)*2 & dat2r > mean(dat2r)-std(dat2r)*2;
                dat1(~okelec) = [];    dat2(~okelec) = [];
                %%-- compute slope of the major axis
                C = cov(dat1,dat2);
                sl = C(1,2)./C(1,1); intc = median(dat2) - sl*median(dat1);
                fitcoefboot(:,iboot) = [intc;sl];
            end
        end
    warning(warnstat);
    datX = minmax([xlim, datrange(1:2)]); datX = min(datX):range(datX)/100:max(datX);
    datY = fitcoef(1,:)'+datX.*fitcoef(2,:)'; 
    hl = plot(datX,mean(datY,1,'omitnan'),lineprop{:});
        if showerror
            datYboot = fitcoefboot(1,:)'+datX.*fitcoefboot(2,:)'; 
            patch([datX, fliplr(datX)], [quantile(datYboot,alpha/2,1), fliplr(quantile(datYboot,1-alpha/2,1))],...
                plcol_wang(mdlidx,:), FaceAlpha=0.3, Edgecolor='none');
        end
    end
    %-- set axis
%     plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    uistack(hl,'top');
    axis(datrange);
    set(gca,'FontSize',FntSiz);

xlabel(gca,  'Eccentricity (deg)','FontSize',FntSiz*1.2,'VerticalAlignment','top');
ylabel(gca,  'Size (deg)','FontSize',FntSiz*1.2);
end
% uistack(flip(hs),'top');
% legend(hs(end-1:end),["Negative Gain","Positive Gain"],'AutoUpdate','off');
legend(hs(1:2:end),modelNames,'AutoUpdate','off');
title(htIn(end),  'pRF size','FontSize',FntSiz*1.1,'FontWeight','bold');


figname =sprintf('Summary1R-%s_%s%s_%s_ROI-low',params,targetBAND,R2mode,selectchs);
if ispltall,     figname = [figname '_withrejected'];   end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');   end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Histogram (Violin) of R2 / Gain &
%       Ecc vs Size with fitting <Individual Plots>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1 % ~issaveplot
npanel = 1;
colpanel = [5 5 6];
nBins = 28;

doviolin = [1 1];
% doviolin = [-1 -1];

colpanel(doviolin<0) = colpanel(doviolin<0) - 1;

plcol_wang = [ones(1,1)*plcol(2,:);...
              ones(1,1)*plcol(4,:);];
plshp_wang = [repmat('o',1,1) repmat('o',1,1)];

for isbj=1:length(prf_all(1).bb.subjects)
subject = prf_all(1).bb.subjects{isbj};
          
% hF(end+1) = figure('Menubar','none','Position',[200 200 400*length(colpanel)+200 400*npanel],'defaultAxesColorOrder',plcol_wang);
hF(end+1) = figure('Menubar','none','Position',[200 200 round(400*sum(colpanel)./colpanel(end))+200 400*npanel],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(1,sum(colpanel),'Padding','compact','TileSpacing','loose');
htIn = gobjects(0);
roilist = {'V1','V2','V3'}';
nmdl = 2;


%%% R2 Histogram / Violin
params = 'xval';
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    case {'peakresp'},      ll = 10;  ul = -10;     roundunit = 3;
    otherwise,              ll = 1e3; ul = -1e3;    roundunit = 5;
end
    %--- get data range
    rangethresh = 70;
    roich  = ismember(prf_all(modelsets(2)).a.channels.wangarea,roilist);
    sbjch  = ismember(prf_all(modelsets(2)).a.channels.subject_name,prf_all(modelsets(2)).a.subjects{isbj});
    if ispltall,    okch   = (okch1bb&okch2bb) & roich & sbjch;
    else,           okch   = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich & sbjch;
    end
    if ~any(okch), delete(hF(end)); continue;  end

    if doviolin(length(htIn)+1)<0
    errfunc = @(d)std(d,'omitnan')./sqrt(sum(~isnan(d)));
    ll = min(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    ul = max(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    else
    datrange = prctile([prf_all(modelsets(2)).a.(params)(okch);prf_all(modelsets(1)).a.(params)(okch)],50+[-0.5 0.5].*rangethresh);
    datrange = fliplr(datrange) + diff(datrange)/rangethresh*100.*[-1 1];
    ll = min(ll,min(datrange));
    ul = max(ul,max(datrange));
    end
    ll = floor(ll/roundunit)*roundunit;
    ul = ceil(ul/roundunit)*roundunit;
    switch params
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end
    if ll > ul, ll=0; ul=1; end

htIn(end+1) = tiledlayout(ht,npanel,1,'Padding','tight','TileSpacing','compact');
htIn(end).Layout.Tile=sum(colpanel(1:(length(htIn)-1)))+1;
htIn(end).Layout.TileSpan=[1,colpanel(length(htIn))];
nexttile(htIn(end));    hold on;
    %--- set axis
    hl = gobjects(0);
    if doviolin(length(htIn)),       ylim([ll ul]); xlim([1 nmdl]+[-1 1]*0.6)
    else,                            xlim([ll ul]);
    end
% for mdlidx=1:2
    %--- plot histogram
    if doviolin(length(htIn))>0
        colororder(plcol_wang); 
        violinplot([prf_all(modelsets(1)).a.(params)(okch),...
                    prf_all(modelsets(2)).a.(params)(okch)]);
    elseif doviolin(length(htIn))<0
        colororder(gca,plcol_wang); 
        imdl = 1; errorbar(imdl,...
                 [mean(prf_all(modelsets(imdl)).a.(params)(okch),'omitnan')],...
                 [errfunc(prf_all(modelsets(imdl)).a.(params)(okch))],...
                            'vertical','o','LineStyle','none','LineWidth',2,...
                            'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
        imdl = 2; errorbar(imdl,...
                 [mean(prf_all(modelsets(imdl)).a.(params)(okch),'omitnan')],...
                 [errfunc(prf_all(modelsets(imdl)).a.(params)(okch))],...
                            'vertical','o','LineStyle','none','LineWidth',2,...
                            'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
    else
        colororder(gca,plcol_wang); 
        histogram(prf_all(modelsets(mdlidx)).a.(params)(okch),'BinLimits',xlim,'NumBins',nBins);
    end
% end
    set(gca,'FontSize',FntSiz);
%     legend(modelNames2l,'AutoUpdate','off');
if doviolin(length(htIn))
%     xticks([]);
    xticklabels(strtok(modelNamesSrt,'-'));
    if prod(ylim)<0, hl(end+1)=plot(xlim,[0,0],'k--','LineWidth',2); end  % axis
    ylabel(gca,  'Cross-validated R^2 (%)','FontSize',FntSiz*1.2,'VerticalAlignment','bottom');
    ylim(ylim+[0 1].*range(ylim)./20);
else
    ylim(ylim+[0 0.2]);
    if prod(xlim)<0, hl(end+1)=plot([0,0], ylim,'k-','LineWidth',4); end  % axis
    xlabel(gca,  'Cross-validated R^2 (%)','FontSize',FntSiz*1.2,'VerticalAlignment','bottom');
    ylabel(gca,  '# of electrodes','FontSize',FntSiz*1.2);
end
uistack(hl,'bottom');
title(htIn(end),  'Variance Explained','FontSize',FntSiz*1.1,'FontWeight','bold');
    


%%% Gain Histogram / Violin
params = 'peakresp';
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    case {'peakresp'},      ll = 10;  ul = -10;     roundunit = 3;
                      if doviolin(length(htIn))<0,  roundunit = 2;  end
    otherwise,              ll = 1e3; ul = -1e3;    roundunit = 5;
end
    %--- get data range
    rangethresh = 95;
    roich  = ismember(prf_all(modelsets(2)).a.channels.wangarea,roilist);
    sbjch  = ismember(prf_all(modelsets(2)).a.channels.subject_name,prf_all(modelsets(2)).a.subjects{isbj});
    if ispltall,    okch   = (okch1bb&okch2bb) & roich & sbjch;
    else,           okch   = ((okch1a&okch1bb)|(okch2a&okch2bb)) & roich & sbjch;
    end
    if doviolin(length(htIn)+1)<0
    errfunc = @(d)std(d,'omitnan')./sqrt(sum(~isnan(d)));
    ll = min(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') - errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    ul = max(mean(prf_all(modelsets(1)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(1)).a.(params)(okch)*1.05),...
             mean(prf_all(modelsets(2)).a.(params)(okch),'omitnan') + errfunc(prf_all(modelsets(2)).a.(params)(okch)*1.05));
    else
    datrange = prctile([prf_all(modelsets(2)).a.(params)(okch);prf_all(modelsets(1)).a.(params)(okch)],50+[-0.5 0.5].*rangethresh);
    datrange = fliplr(datrange) + diff(datrange)/rangethresh*100.*[-1 1];
    ll = min(ll,min(datrange));
    ul = max(ul,max(datrange));
    end
    ll = floor(ll/roundunit)*roundunit;
    ul = ceil(ul/roundunit)*roundunit;
    switch params
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end
    if ll > ul, ll=0; ul=1; end

htIn(end+1) = tiledlayout(ht,npanel,1,'Padding','tight','TileSpacing','compact');
htIn(end).Layout.Tile=sum(colpanel(1:(length(htIn)-1)))+1;
htIn(end).Layout.TileSpan=[1,colpanel(length(htIn))];
nexttile(htIn(end));    hold on;
    %--- set axis
    hl = gobjects(0);
    if doviolin(length(htIn)),       ylim([ll ul]); xlim([1 nmdl]+[-1 1]*0.6)
    else,                            xlim([ll ul]);
    end
% for mdlidx=1:2
    %--- plot histogram
    if doviolin(length(htIn))>0
        colororder(plcol_wang); 
        violinplot([prf_all(modelsets(1)).a.(params)(okch),...
                    prf_all(modelsets(2)).a.(params)(okch)]);
    elseif doviolin(length(htIn))<0
        colororder(gca,plcol_wang); 
        imdl = 1; errorbar(imdl,...
                 [mean(prf_all(modelsets(imdl)).a.(params)(okch),'omitnan')],...
                 [errfunc(prf_all(modelsets(imdl)).a.(params)(okch))],...
                            'vertical','o','LineStyle','none','LineWidth',2,...
                            'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
        imdl = 2; errorbar(imdl,...
                 [mean(prf_all(modelsets(imdl)).a.(params)(okch),'omitnan')],...
                 [errfunc(prf_all(modelsets(imdl)).a.(params)(okch))],...
                            'vertical','o','LineStyle','none','LineWidth',2,...
                            'CapSize',6,'MarkerFaceColor','auto','MarkerSize',8);
    else
        colororder(gca,plcol_wang); 
        histogram(prf_all(modelsets(mdlidx)).a.(params)(okch),'BinLimits',xlim,'NumBins',nBins);
    end
% end
    set(gca,'FontSize',FntSiz);
%     legend(modelNamesSrt,'AutoUpdate','off');
if doviolin(length(htIn))
%     xticks([]);
    xticklabels(strtok(modelNamesSrt,'-'));
    if prod(ylim)<0, hl(end+1)=plot(xlim,[0,0],'k--','LineWidth',2); end  % axis
    ylabel(gca,  'Ratio of Response','FontSize',FntSiz*1.2,'VerticalAlignment','bottom');
    ylim(ylim+[-1 1].*range(ylim)./20);
else
    ylim(ylim+[0 0.2]);
    if prod(xlim)<0, hl(end+1)=plot([0,0], ylim,'k-','LineWidth',4); end  % axis
    xlabel(gca,  'Ratio of Response','FontSize',FntSiz*1.2,'VerticalAlignment','bottom');
    ylabel(gca,  '# of electrodes','FontSize',FntSiz*1.2);
end
%-- Update TickLabel
if doviolin(length(htIn))
    if length(yticks)<5
        yticks(linspace(floor(min(ylim)),ceil(max(ylim)),5));
    end
    yticklabels(yticks+sign(yticks)+(yticks==0));
    ticknames = yticklabels;
    ticknames(yticks<0)  = cellfun(@(d) sprintf('1/%d',-str2double(d)),ticknames(yticks<0),'UniformOutput',false);
    yticklabels(ticknames);
else
    if length(xticks)<5
        xticks(linspace(floor(min(xlim)),ceil(max(xlim)),5));
    end
    xticklabels(xticks+sign(xticks)+(xticks==0));
    ticknames = xticklabels;
    ticknames(xticks<0)  = cellfun(@(d) sprintf('1/%d',-str2double(d)),ticknames(xticks<0),'UniformOutput',false);
    xticklabels(ticknames);
end
uistack(hl,'bottom');
title(htIn(end), 'Response Gain','FontSize',FntSiz*1.1,'FontWeight','bold');
    

%%% Ecc vs Size
params=["rfsize","ecc"];

%--- get data range
roich  = ismember(prf_all(modelsets(1)).a.channels.wangarea,roilist);
sbjch  = ismember(prf_all(modelsets(1)).a.channels.subject_name,prf_all(modelsets(1)).a.subjects{isbj});
okch   = (okch1bb&okch2bb) & roich & sbjch;
if ~ispltall,    okch   = ((okch1a)|(okch2a)) & okch;    end

datrange = [0 prctile([prf_all(modelsets(1)).a.(params(2))(okch);prf_all(modelsets(2)).a.(params(2))(okch)],98)*cfactor];
datrange = [0 ceil(diff(datrange)/.98*1.05)];
datrange = [datrange datrange.*[1 10/8.5]];    % in degree [Ecc, Size]

htIn(end+1) = tiledlayout(ht,npanel,1,'Padding','tight','TileSpacing','compact');
htIn(end).Layout.Tile=sum(colpanel(1:(length(htIn)-1)))+1;
htIn(end).Layout.TileSpan=[1,colpanel(length(htIn))];
nexttile(htIn(end));    hold on;
hs = gobjects(0);
for mdlidx=1:2
    %--- plot scatter (for each gain)
    for ii=1:2
    scatterprop = {mkSiz,plshp_wang(ii),'LineWidth',2.0,'MarkerEdgeColor',plcol_wang(mdlidx,:)};
    if ii==1,    gainch = (prf_all(modelsets(mdlidx)).a.gain <= 0);  scatterprop = [scatterprop, {'MarkerFaceColor',plcol_wang(mdlidx,:)}];
    else,        gainch = (prf_all(modelsets(mdlidx)).a.gain > 0);
    end
    okch   = (okch1bb&okch2bb) & roich & gainch;
    if ~ispltall,    okch   = ((okch1a)|(okch2a)) & okch;    end
    dat1 = prf_all(modelsets(mdlidx)).a.(params(2))(okch).*cfactor;
    dat2 = prf_all(modelsets(mdlidx)).a.(params(1))(okch).*cfactor;
    hs(end+1) = scatter(dat1,dat2,scatterprop{:});
    end
    %-- fit line (Robust method)
    %%-- collect data (ignore gains)
    ii=1;
    okch   = (okch1bb&okch2bb) & roich;
    if ~ispltall,    okch   = ((okch1a)|(okch2a)) & okch;    end
    nchan = sum(okch);
    if nchan>3
    lineprop = {'Color',plcol_wang(mdlidx,:),'LineWidth',2.0};
    if ii==1, lineprop = [{'--'},lineprop];
    else,     lineprop = [{':'},lineprop];
    end
    warnstat = warning; warning('off');
    dat1 = prf_all(modelsets(mdlidx)).a.(params(2))(okch).*cfactor;
    dat2 = prf_all(modelsets(mdlidx)).a.(params(1))(okch).*cfactor;
    %%-- exclude outliers (exclude data which are out of more than 2 sigma on the major and minor axes)
    C = cov(dat1,dat2);  sl = C(1,2)./C(1,1); 
    rot2d = @(D,deg) D*[cos(deg),-sin(deg);sin(deg),cos(deg)];
    rotD = rot2d([dat1,dat2],atan(sl));
    dat1r = rotD(:,1); dat2r = rotD(:,2);
    okelec = dat1r < mean(dat1r)+std(dat1r)*2 & dat2r < mean(dat2r)+std(dat2r)*2 &...
             dat1r > mean(dat1r)-std(dat1r)*2 & dat2r > mean(dat2r)-std(dat2r)*2;
    dat1(~okelec) = [];    dat2(~okelec) = [];
    %%-- compute slope of the major axis
    C = cov(dat1,dat2);
    sl = C(1,2)./C(1,1); intc = median(dat2) - sl*median(dat1);
    fitcoef = [intc;sl];
    warning(warnstat);
    datX = minmax([xlim, datrange(1:2)]); datX = min(datX):range(datX)/100:max(datX);
    datY = fitcoef(1,:)'+datX.*fitcoef(2,:)'; 
    plot(datX,mean(datY,1,'omitnan'),lineprop{:});
    end
    %-- set axis
%     plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    axis(datrange);
    set(gca,'FontSize',FntSiz);

xlabel(gca,  'Eccentricity (deg)','FontSize',FntSiz*1.2,'VerticalAlignment','top');
ylabel(gca,  'Size (deg)','FontSize',FntSiz*1.2);
end
% uistack(flip(hs),'top');
% legend(hs(end-1:end),["Negative Gain","Positive Gain"],'AutoUpdate','off');
legend(hs(1:2:end),modelNames,'AutoUpdate','off');
title(htIn(end),  'pRF size','FontSize',FntSiz*1.1,'FontWeight','bold');


figname =sprintf('Summary1R-%s-%s_%s%s_%s_ROI-low',params,subject,targetBAND,R2mode,selectchs);
if ispltall,     figname = [figname '_withrejected'];   end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');   end
end
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%     Additional Plots in V1–V3 w/ histogram (classify pos & neg gain in Power Change in Alpha)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;

% ispltall  = false;
ispltall  = true;
dobbthrsh = true;
mkSiz = 60;
fithist = false;
%% Difference of Gaussian Model in training accuracy (scatter) w/ rough ROI
if issaveplot
npanel = 1;

nBins = 36;
params = 'peakresp';

plcol_wang = [ones(1,1)*plcol(6,:);...
              ones(1,1)*plcol(4,:);];
plshp_wang = [repmat('o',1,1) repmat('^',1,1)];
          
% hF(end+1) = figure('Menubar','none','Position',[200 200 860 420],'defaultAxesColorOrder',plcol_wang);
hF(end+1) = figure('Menubar','none','Position',[200 200 550*npanel 540],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(1,npanel,'Padding','compact','TileSpacing','compact');
roilist = {'V1','V2','V3'}';
nexttile;
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    otherwise,              ll = 1e3; ul = -1e3;  roundunit = 5;
end
    %--- get data range & set axis
    roich  = ismember(prf_all(modelsets(2)).a.channels.wangarea,roilist);
    if ispltall, okch = roich;
    else,        okch = ((okch1a)|(okch2a)) & roich;
    end
    if dobbthrsh,  okch = okch & okch1bb & okch2bb;  end
    datrange = prctile([prf_all(modelsets(2)).a.(params)(okch);prf_all(modelsets(1)).a.(params)(okch)],[15 85]);
    datrange = fliplr(datrange) + diff(datrange)/.7.*[-1 1];
    ll = min(ll,floor(min(datrange)/roundunit)*roundunit);
    ul = max(ul,ceil(max(datrange)/roundunit)*roundunit);
    switch params
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end
    %--- make table
    dattbl = table(prf_all(modelsets(1)).a.(params),prf_all(modelsets(2)).a.(params),prf_all(modelsets(2)).a.gain <= 0);
    dattbl(~okch,:) = [];
    dattbl.Properties.VariableNames = [modelNames "Neg Gain"];
    scatterprop = {'MarkerSize',1,'Marker','.','LineWidth',1.6,'Color',ones(1,3).*0.3,'Location','NorthEast','Direction','out'};
    if fithist, scatterprop = [scatterprop {'Kernel','on'}]; end
    hh = scatterhist(dattbl.(modelNames(2)),dattbl.(modelNames(1)),'Group',ones(height(dattbl),1),scatterprop{:});
    axis([ll ul ll ul],'square');
    hold on;
    %--- plot scatter
    hs = gobjects(0);
    for ii=1:2
    scatterprop = {mkSiz,plshp_wang(ii),'LineWidth',1.6,'MarkerEdgeColor',plcol_wang(ii,:),'MarkerFaceColor',plcol_wang(ii,:)};
    if ii==1,    gainch = (prf_all(modelsets(2)).a.gain <= 0);
    else,        gainch = (prf_all(modelsets(2)).a.gain > 0);
    end
    okch   = ((okch1a)|(okch2a)) & roich & gainch;
    if dobbthrsh,  okch = okch & okch1bb & okch2bb;  end
    hs(end+1) = scatter(prf_all(modelsets(2)).a.(params)(okch),prf_all(modelsets(1)).a.(params)(okch),scatterprop{:});
    if ispltall
        okch   = roich & gainch;
        if dobbthrsh,  okch = okch & okch1bb & okch2bb;  end
        scatter(prf_all(modelsets(2)).a.(params)(okch),prf_all(modelsets(1)).a.(params)(okch),scatterprop{:},'MarkerEdgeColor',plcol_gray(ii,:),'MarkerFaceColor',plcol_gray(ii,:));
    end
    end
    uistack(hs,'top');
    plot(xlim,ylim,'k--');  % diagonal
    plot(xlim,[0 0],'k:',[0 0],ylim,'k:');  % axis
    if ismember(params,'xval')
    plot([getelement(xlim,1) threshold(modelsets(2)).a],[1 1].*threshold(modelsets(1)).a,'k-.',...
         [1 1].*threshold(modelsets(2)).a,[getelement(ylim,1) threshold(modelsets(1)).a],'k-.');  % threshold
%     legend(hs([1,4]),["V1-V3","Dorsolateral"],'Location','southeast');
    end
    legend(hs,{'Neg Resp','Pos Resp'});
    
    set(gca,'FontSize',FntSiz);
    xlabel(modelNames(modelsets(2))); ylabel(modelNames(modelsets(1)));
%     title(params);
    xticks(yticks); xtickangle(0);

    hh(1).XAxisLocation = 'bottom';
    hh(1).YAxisLocation = 'left';
    if ~fithist
        hh(2).Children.BinLimits = xlim;
        hh(2).Children.NumBins   = nBins;
        ylim(hh(2),'auto');
        hh(3).Children.BinLimits = xlim;
        hh(3).Children.NumBins   = nBins;
        ylim(hh(3),'auto');

%         drawnow;            
%         hh(2).OuterPosition = hh(2).OuterPosition.*[1 1 1 0.8];
%         hh(3).OuterPosition = hh(3).OuterPosition.*[1 1 0.8 1];
    end
    hh(1).Position(1:2) = [1 1].*0.15;
 
figname =sprintf('ScatterHist-%s_%s%s_%s_ROI-low',params,targetBAND,R2mode,selectchs);
if ispltall,     figname = [figname '_withrejected'];   end
if dobbthrsh,    figname = [figname '_bbthresholded'];   end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end

