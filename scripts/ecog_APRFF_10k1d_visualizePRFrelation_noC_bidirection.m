% visualize pRF correlations computed in 06c (equivalant 03c2 + 03d2)
%   Update xR2 with full time-series (no decimation)
%   Folk from ecog_APRF_06d_visualizePRFparamsFPM_bootall_fullts.m
%         and ecog_APRF_06e_visualizePRFsFPM_bootall_fullts.m

% 20210107 Yuasa - Update from 06d,06e
% 20210809 Yuasa - modify for paper
% 20210916 Yuasa - for APRFF

%% prefix
% close all; clearvars;
% startupToolboxToolbox;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'pRFrelations-representative');
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
else
plotsavePth    = 'pRFrelations-representative';
prfPth         = 'pRFmodel';
prfstatPth     = 'pRFanalysis';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'noACorr_bothgain');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%-- Plotting Setting
FntSiz  = 18;
LFntSiz = 22;

%% Load dataset
clear alphaType broadbandType

%%% load analyzePRF & recompute full-ts R2
average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';
selectchs      = 'wangprobchs';
    allowlag       = false;
    allowbeta      = false;
    allowwide      = false;
    allowmixbeta   = true;
    noACorrect     = false;

va_area = 'wangarea';
usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

prfID = 'prfbidirgain';
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
prf_all_bbU = prf_all_bb;
prf_all_aU  = prf_all_a;

clear alphaType broadbandType modeldataID prfID 
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

%% Load correlations
%-- Load main data
% nboot         = 5000;
boottype      = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling
robustfit     = true;
noalphathresh = true;
noeccthresh   = true;          % if true, avoid to apply eccentricity restriction

if noalphathresh,       threshold_a = nan;  end   % just apply broadband threshold

[prfs,prfs_corr,prfs_corrx,rois,nroi] = ...
    ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit,1,1);

%%% prepare visualize
% nsubjects = length(prf_params_bb);
% nroi  = length(rois);
% selroi = (nroi-1):nroi;
% nvisarea  = length(prf_va_bb);
% areaList  = unique(prf_all_bb.channels.(va_area));

res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');
hF = gobjects(0);

%% Add parameters on prf-bootstrapping

for iroi=1:nroi
    
    %-- get bootstrap channel index
    if iscell(prfs.chanidx)
      ichx = [prfs.chanidx{iroi,:}];
    else
      chidx = [];
      chidx(:,1) = nearlyeq(prf_all_bb.ecc,prfs.ecc_bb{iroi}');
      chidx(:,2) = nearlyeq(prf_all_bb.ang,prfs.ang_bb{iroi}');
      chidx(:,3) = nearlyeq(prf_all_bb.rfsize,prfs.rfsize_bb{iroi}');
      chidx(:,4) = nearlyeq(prf_all_a.ecc,prfs.ecc_a{iroi}');
      chidx(:,5) = nearlyeq(prf_all_a.ang,prfs.ang_a{iroi}');
      chidx(:,6) = nearlyeq(prf_all_a.rfsize,prfs.rfsize_a{iroi}');
      ichx = mode(chidx,2);
    end
    
    %-- collect data
    prfs.R2_aU{iroi,1}       = prf_all_aU.xval(ichx)';
    prfs.ecc_aU{iroi,1}      = prf_all_aU.ecc(ichx)';
    prfs.ang_aU{iroi,1}      = prf_all_aU.ang(ichx)';
    prfs.rfsize_aU{iroi,1}   = prf_all_aU.rfsize(ichx)';
    prfs.R2_bbU{iroi,1}     = prf_all_bbU.xval(ichx)';
    prfs.ecc_bbU{iroi,1}    = prf_all_bbU.ecc(ichx)';
    prfs.ang_bbU{iroi,1}    = prf_all_bbU.ang(ichx)';
    prfs.rfsize_bbU{iroi,1} = prf_all_bbU.rfsize(ichx)';

    %-- collect gain
    prfs.gain_a{iroi,1}     = -prf_all_a.gain(ichx)';
    prfs.gain_bb{iroi,1}    = prf_all_bb.gain(ichx)';
    prfs.gain_aU{iroi,1}    = -prf_all_aU.gain(ichx)';
    prfs.gain_bbU{iroi,1}   = prf_all_bbU.gain(ichx)';

end

%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations
% close all;

% fitting = '-ellipse';
% fitting = '-robustellipse';     % exclude outliers
fitting = '-none';     % no fitting line

%% plot 2D histogram (replication around 90 deg for angle)
%% %%%%%%%%%%%%%%%%%%%%
%%% arrange in one figure
if issaveplot

relt    = {'bb_bbU','a_aU','bb_a','bbU_aU'};
flds    = {'R2','gain','ang','ecc','rfsize'};
selroi  = {'V1-V3','Dorsolateral'};
opts = [];
opts.mothod      = 'hist';
opts.ang_limpad  = 90;       % degree for angle
opts.pix2deg     = cfactor;
opts.FontSize    = FntSiz;
opts.showfitted  = strtok(fitting,'-');
opts.fit_alpha   = alpha;
opts.condNames   = {'Broadband','Alpha','Unrestricted Broadband','Unrestricted Alpha'};
h = ecog_prf_plotPRFrelations(prfs,prfs_corr,relt,flds,selroi,rois,opts);
hF = [hF h];
    
for ipl = 1:length(relt)
figname = sprintf('prf-gaincomparison_2Dhist%s_%s',fitting,relt{ipl});
if issaveplot,   savefigauto(h(ipl),fullfile(plotsavedir, figname),'-vector');   end
end

end

%% %%%%%%%%%%%%%%%%%%%%
%% Compare models
%% %%%%%%%%%%%%%%%%%%%%

prf_all_bbO = prf_all_bb;
prf_all_aO  = prf_all_a;

clear alphaType broadbandType modeldataID prfID 
prfID = 'prfbidirgain';
    allowmixbeta   = false;
    noACorrect     = true;
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;

prf_all_aNCU  = prf_all_a;
prf_all_bb   = prf_all_bbO;
prf_all_a    = prf_all_aO;


%% get parameters

for iroi=1:nroi
    
    %-- get bootstrap channel index
    if iscell(prfs.chanidx)
      ichx = [prfs.chanidx{iroi,:}];
    else
      chidx = [];
      chidx(:,1) = nearlyeq(prf_all_bb.ecc,prfs.ecc_bb{iroi}');
      chidx(:,2) = nearlyeq(prf_all_bb.ang,prfs.ang_bb{iroi}');
      chidx(:,3) = nearlyeq(prf_all_bb.rfsize,prfs.rfsize_bb{iroi}');
      chidx(:,4) = nearlyeq(prf_all_a.ecc,prfs.ecc_a{iroi}');
      chidx(:,5) = nearlyeq(prf_all_a.ang,prfs.ang_a{iroi}');
      chidx(:,6) = nearlyeq(prf_all_a.rfsize,prfs.rfsize_a{iroi}');
      ichx = mode(chidx,2);
    end
    
    %-- collect data
    prfs.R2_aNCU{iroi,1}       = prf_all_aNCU.xval(ichx)';
    prfs.ecc_aNCU{iroi,1}      = prf_all_aNCU.ecc(ichx)';
    prfs.ang_aNCU{iroi,1}      = prf_all_aNCU.ang(ichx)';
    prfs.rfsize_aNCU{iroi,1}   = prf_all_aNCU.rfsize(ichx)';
    prfs.gain_aNCU{iroi,1}    = -prf_all_aNCU.gain(ichx)';


end

prfsC = [];
roisC = {'BroadbandVSModelAlpha','BroadbandVSAlphaPower','ModelAlphaVSAlphaPower'};
reltC  = {'bb_aU','bb_aNCU','aU_aNCU'};
selroiC  = {'V1-V3'};
for iroi = 1:length(roisC)
    tarroi = find(ismember(rois,selroiC));
    [cond1,cond2] = strtok(reltC{iroi},'_'); cond2(1) = [];

    prfsC.R2_base{iroi,1}       = prfs.(['R2_' cond1]){tarroi,1};
    prfsC.ecc_base{iroi,1}      = prfs.(['ecc_' cond1]){tarroi,1};
    prfsC.ang_base{iroi,1}      = prfs.(['ang_' cond1]){tarroi,1};
    prfsC.rfsize_base{iroi,1}   = prfs.(['rfsize_' cond1]){tarroi,1};
    prfsC.gain_base{iroi,1}     = prfs.(['gain_' cond1]){tarroi,1};

    prfsC.R2_tar{iroi,1}       = prfs.(['R2_' cond2]){tarroi,1};
    prfsC.ecc_tar{iroi,1}      = prfs.(['ecc_' cond2]){tarroi,1};
    prfsC.ang_tar{iroi,1}      = prfs.(['ang_' cond2]){tarroi,1};
    prfsC.rfsize_tar{iroi,1}   = prfs.(['rfsize_' cond2]){tarroi,1};
    prfsC.gain_tar{iroi,1}     = prfs.(['gain_' cond2]){tarroi,1};
end
    

%%
if issaveplot
relt    = 'base_tar';
flds    = {'R2','gain','ang','ecc','rfsize'};
opts = [];
opts.mothod      = 'hist';
opts.ang_limpad  = 90;       % degree for angle
opts.pix2deg     = cfactor;
opts.FontSize    = FntSiz;
opts.showfitted  = strtok(fitting,'-');
opts.fit_alpha   = alpha;
h = ecog_prf_plotPRFrelations(prfsC,[],relt,flds,roisC,roisC,opts);
hF = [hF h];
    
figname = sprintf('prf-gaincomparison_2Dhist%s_BiDir_%s',fitting,selroiC{:});
if issaveplot,   savefigauto(h,fullfile(plotsavedir, figname),'-vector');   end
end

%%
roisC2 = {sprintf('Model-based\nAlpha Suppression');sprintf('Power Change\nin Alpha')};
fitting = '-robustellipse';     % exclude outliers
relt    = 'base_tar';
flds    = {'ang','ecc','rfsize'};
opts = [];
opts.mothod      = 'hist';
opts.ang_limpad  = 90;       % degree for angle
opts.pix2deg     = cfactor;
opts.FontSize    = FntSiz;
opts.showfitted  = strtok(fitting,'-');
opts.fit_alpha   = alpha;
opts.condNames   = {'Broadband','Alpha'};
h = ecog_prf_plotPRFrelations(prfsC,[],relt,flds,roisC2,roisC2,opts);
hF = [hF h];
    
figname = sprintf('prf-modelcomparison_2Dhist%s_%s',fitting,'all');
if issaveplot,   savefigauto(h,fullfile(plotsavedir, figname),'-vector');   end



%%




%%
%%
%%
%%
if ~issaveplot
%% Load results
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
prf_all_a_noAC  = prf_all_a;
%% Load prfs
%-- Set parameter for main data
clear alphaType broadbandType
    allowmixbeta   = true;
    noACorrect     = false;
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

%-- Load main data
% nboot         = 5000;
boottype      = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling
robustfit     = true;
noalphathresh = true;
noeccthresh   = true;          % if true, avoid to apply eccentricity restriction

if noalphathresh,       threshold_a = nan;  end   % just apply broadband threshold
[prfsO,prfsO_corr,prfsO_corrx,rois,nroi] = ...
    ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit,1,1);

%% Update prfs to plot same electrodes

prfsO.gain_a  = cell(size(prfsO.chanidx,1),1);
prfsO.gain_bb = cell(size(prfsO.chanidx,1),1);

prfsU = prfsO;
prfsU_corr = prfsO_corr;
for iroi=1:nroi
    
    %-- get bootstrap channel index
    if iscell(prfsO.chanidx)
      ichx = [prfsO.chanidx{iroi,:}];
    else
      chidx = [];
      chidx(:,1) = nearlyeq(prf_all_bb.ecc,prfsO.ecc_bb{iroi}');
      chidx(:,2) = nearlyeq(prf_all_bb.ang,prfsO.ang_bb{iroi}');
      chidx(:,3) = nearlyeq(prf_all_bb.rfsize,prfsO.rfsize_bb{iroi}');
      chidx(:,4) = nearlyeq(prf_all_a.ecc,prfsO.ecc_a{iroi}');
      chidx(:,5) = nearlyeq(prf_all_a.ang,prfsO.ang_a{iroi}');
      chidx(:,6) = nearlyeq(prf_all_a.rfsize,prfsO.rfsize_a{iroi}');
      ichx = mode(chidx,2);
    end
    
    %-- collect data
    prfsU.R2_a{iroi,1}       = prf_all_a_noAC.xval(ichx)';
    prfsU.ecc_a{iroi,1}      = prf_all_a_noAC.ecc(ichx)';
    prfsU.ang_a{iroi,1}      = prf_all_a_noAC.ang(ichx)';
    prfsU.rfsize_a{iroi,1}   = prf_all_a_noAC.rfsize(ichx)';

    %-- collect gain
    prfsO.gain_a{iroi,1}     = prf_all_a.gain(ichx)';
    prfsO.gain_bb{iroi,1}    = prf_all_bb.gain(ichx)';
    prfsU.gain_a{iroi,1}     = prf_all_a_noAC.gain(ichx)';
    prfsU.gain_bb{iroi,1}    = prfsO.gain_bb{iroi,1};

end

%% Merge prfs
tarroi = 'V1-V3';
iroi   = find(ismember(rois,tarroi),1);

prfsC = [];
paramnames = fieldnames(prfsO);
for iprm = reshape(paramnames,1,[])
    if size(prfsO.(iprm{:}),1) == nroi
        prfsC.(iprm{:}) = vertcat(prfsO.(iprm{:})(iroi,:),prfsU.(iprm{:})(iroi,:));
    else
        prfsC.(iprm{:}) = prfsO.(iprm{:});
    end
end

prfsC_corr = [];
% paramnames = fieldnames(prfsO_corr);
% for iprm = reshape(paramnames,1,[])
%     if size(prfsO_corr.(iprm{:}),1) == nroi
%         prfsC_corr.(iprm{:}) = vertcat(prfsO_corr.(iprm{:})(iroi,:),prfsU_corr.(iprm{:})(iroi,:));
%     else
%         prfsC_corr.(iprm{:}) = prfsO_corr.(iprm{:});
%     end
% end

%-- Exclude electrodes with too large eccentricity
if ~noeccthresh
% badch = (prfsC.ecc_a{1}>eclimit)|(prfsC.ecc_a{2}>eclimit);
for iroi = 1:2
badch = (prfsC.ecc_a{iroi}>eclimit);
for iprm = reshape(paramnames,1,[])
    if iscell(prfsC.(iprm{:})) && size(prfsC.(iprm{:}){iroi},2) == length(badch)
        prfsC.(iprm{:}){iroi}(:,badch) = [];
    end
end
end
end

%-- Set ROI names
roisC = {sprintf('Model-based\nAlpha Suppression');sprintf('Power Change\nin Alpha')};

%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations
% close all;

% fitting = '-ellipse';
fitting = '-robustellipse';     % exclude outliers

%% plot 2D histogram (replication around 90 deg for angle)
%% %%%%%%%%%%%%%%%%%%%%
%%% arrange in one figure

relt    = 'bb_a';
flds    = {'ang','ecc','rfsize'};
opts = [];
opts.mothod      = 'hist';
opts.ang_limpad  = 90;       % degree for angle
opts.pix2deg     = cfactor;
opts.FontSize    = FntSiz;
opts.showfitted  = strtok(fitting,'-');
opts.fit_alpha   = alpha;
h = ecog_prf_plotPRFrelations(prfsC,prfsC_corr,relt,flds,roisC,roisC,opts);
hF = [hF h];
    
figname = sprintf('prf-modelcomparison_2Dhist%s_%s',fitting,'all');
if noeccthresh, figname = sprintf('%s-noecc',figname);  end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');   end


%% plot gain in 2D histogram
%%% arrange in one figure

relt    = 'bb_a';
flds    = {'gain'};
opts = [];
opts.mothod      = 'hist';
opts.ang_limpad  = 90;       % degree for angle
opts.pix2deg     = cfactor;
opts.FontSize    = FntSiz;
opts.showfitted  = strtok(fitting,'-');
opts.fit_alpha   = alpha;
opts.showdiag = false;
h = ecog_prf_plotPRFrelations(prfsC,prfsC_corr,relt,flds,roisC,roisC,opts);
hF = [hF h];
    
figname = sprintf('prf-modelcomparison_gain2Dhist%s_%s',fitting,'all');
if noeccthresh, figname = sprintf('%s-noecc',figname);  end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');   end

end


%% plot gain in 1D histogram
if issaveplot,      npanel = 2;
else,               npanel = 1;
end

nBins = 36;
params = 'gain';

plcol_wang = [ones(1,1)*plcol(2,:);...
              ones(1,1)*plcol(4,:);];
          
hF(end+1) = figure('Menubar','none','Position',[200 200 600 350*npanel],'defaultAxesColorOrder',plcol_wang);
ht = tiledlayout(npanel,1,'Padding','compact','TileSpacing','compact');
selroi  = {'V1-V3','Dorsolateral'}; selroi(npanel+1:end) = [];
for ii=1:npanel
nexttile;
switch params
    case {'xval'},          ll = 100; ul = -100;    roundunit = 20;
    case {'ecc','rfsize'},  ll = 0;   ul = 50;      roundunit = 5;
    otherwise,              ll = 1e3; ul = -1e3;  roundunit = 5;
end
    %--- get data range & set axis
    iroi = find(ismember(rois,selroi(ii)));
    datrange = prctile([prfs.([params '_aU']){iroi},prfs.([params '_aNCU']){iroi}],[15 85]);
    datrange = fliplr(datrange) + diff(datrange)/.7.*[-1 1];
    ll = min(ll,floor(min(datrange)/roundunit)*roundunit);
    ul = max(ul,ceil(max(datrange)/roundunit)*roundunit);
    switch params
        case {'xval'},          ll = max(-100,ll);    ul = min(100,ul);
        case {'ecc','rfsize'},  ll = max(0,ll);       ul = min(60,ul);
    end
    xlim([ll ul]);
    hold on;
    %--- plot histogram
    hs = gobjects(0);
    hs(end+1) = histogram(prfs.([params '_aU']){iroi},'BinLimits',xlim,'NumBins',nBins);
    hs(end+1) = histogram(prfs.([params '_aNCU']){iroi},'BinLimits',xlim,'NumBins',nBins);
    if ii==1, legend(roisC,'AutoUpdate','off'); end
    set(gca,'FontSize',FntSiz);
    plot([0 0],ylim,'k-','LineWidth',4);  % axis
    title(sprintf('%s - %s',params,selroi{ii}));
end
set(gcf,'Name','pRF gain histogram');

figname = sprintf('prf-modelcomparison_gain1Dhist%s_%s',fitting,'all');
if noeccthresh, figname = sprintf('%s-noecc',figname);  end
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');   end
    