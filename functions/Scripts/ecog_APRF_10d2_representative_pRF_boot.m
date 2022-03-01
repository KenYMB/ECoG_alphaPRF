%% ECoG Alpha pRF (simple broadband computation)
% plot bootstrap results
% for figure 4

% 20201215
% 20210224 - update for finalize figure 4


%% Define paths and dataset
% close all; clear all;
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
SetDefault('issaveplot',true);
if issaveplot
    plotsavedir    = fullfile(figPth, 'pRF-representative','boot');
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%% plot parameter
repelec = {'name',{'GB103'},'subject_name',{'som726'}};
tarBAND = 'both';

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

modeldataID  = 'freq_spectra-timeseries-boot';
prfID        = 'prf-boot';
usefulltsxR2 = false;
ecog_APRF_INITa_loaddata;
ecog_APRF_INITb_mergedata;

model_all_bb_boot = model_all_bb;
model_all_a_boot  = model_all_a;
prf_all_bb_boot   = prf_all_bb;
prf_all_a_boot    = prf_all_a;

%%% load original pRF data
modeldataID = [];
prfID       = [];
usefulltsxR2 = false;
ecog_APRF_INITa_loaddata;
ecog_APRF_INITb_mergedata;

model_all_bb_orig = model_all_bb;
model_all_a_orig  = model_all_a;
prf_all_bb_orig   = prf_all_bb;
prf_all_a_orig    = prf_all_a;

%% set threshold
ecog_APRF_INITc_postfix;
ecog_APRF_INITd_threshold;
thresh_ecc = eclimit;
% thresh_ecc = nan;

elec_ok = ~(prf_all_bb.xval<=threshold_bb | prf_all_a.xval<=threshold_a ...
         | prf_all_bb.ecc >= thresh_ecc | prf_all_a.ecc >= thresh_ecc); 

%% visualize parameter
res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

plcol = get(groot,'defaultAxesColorOrder');
     
% close all;

%%
%%
%% pick up bootstrap trials

Nboot = 1000;
pickidx = (1:100)+100;
prf_all_bb_bootP = prf_all_bb_boot;
prf_all_a_bootP  = prf_all_a_boot;

prffields = fieldnames(prf_all_bb_boot);
for iflds = prffields'
    if size(prf_all_bb_bootP.(iflds{:}),3)==Nboot
        prf_all_bb_bootP.(iflds{:}) = prf_all_bb_bootP.(iflds{:})(:,:,pickidx);
        prf_all_a_bootP.(iflds{:}) = prf_all_a_bootP.(iflds{:})(:,:,pickidx);
    elseif size(prf_all_bb_bootP.(iflds{:}),4)==Nboot
        prf_all_bb_bootP.(iflds{:}) = prf_all_bb_bootP.(iflds{:})(:,:,:,pickidx);
        prf_all_a_bootP.(iflds{:}) = prf_all_a_bootP.(iflds{:})(:,:,:,pickidx);
    end
end

%% plot averaged pRF location overlay
% close all

cntcol      = plcol([2,1],:);

el = findtable(prf_all_bb_boot.channels,repelec{:});
for ii=1:numel(repelec{2})
figureTitle = sprintf('pRFplot-%s-%s_loc',prf_all_bb_boot.channels.subject_name(el(ii)),prf_all_bb_boot.channels.name{el(ii)});

opt = [];
opt.plot.fontSize   = 20;
opt.plot.resolution = resmx;
opt.plot.pix2deg    = cfactor;
opt.plot.XLim       = [-10 10];
opt.plot.YLim       = [-10 10];
opt.plot.addWangToTitle   = 'no';
opt.plot.addBensonToTitle = 'no';
opt.plot.addSbjToTitle    = 'yes';
opt.plot.Type       = 'contour';
opt.plot.options    = {'LineWidth',2.6};

%-- plot bootstrap each
opt.plot.sigma  = [1];
opt.plot.isavg  = false;
opt.plot.colors = cntcol;
% opt.plot.colors = [cntcol [0.3;0.3]];
opt.plot.addCenter = 'no';
opt.plot.addCenter = 'yes';
opt.plot.options    = {'LineWidth',0.5};
hF = ecog_plotGridPRF(el(ii), opt, prf_all_a_bootP, prf_all_bb_bootP);
  set(gcf,'MenuBar','none','Position',get(gcf,'Position').*[1,1,2.2,1.8]);
  title('');
  adjustplot(gcf);

    figureName = sprintf('%s_%s-%s',figureTitle,'booteach-1sd-center',tarBAND);
    if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName));  end
    
end
    