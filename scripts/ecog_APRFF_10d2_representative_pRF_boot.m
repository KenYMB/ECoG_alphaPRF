%% ECoG Alpha pRF (simple broadband computation)
% plot bootstrap results
% for figure 4

% 20201215
% 20210224 - update for finalize figure 4
% %% without ERP %%

%%
% close all; clear all;
% startupToolboxToolbox;
run_checkPath;

%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'pRF-representative');
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
else
plotsavePth    = 'pRF-representative';
prfPth         = 'pRFmodel';
prfstatPth     = 'pRFanalysis';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'boot');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%-- Plotting Setting
FntSiz  = 20;

%% Representative channels
if issaveplot
repelec = {'name',{'Oc18','GB103','Oc17','GB102'},...
            'subject_name',{'p02','p10','p02','p10'}};
else
% repelec = {'name',{'Oc17'},'subject_name',{'p02'}};
repelec = {'name',{'GB103'},'subject_name',{'p10'}};
end

selsbj      = ismember(subjectList,repelec{4});
subjectList = subjectList(selsbj);

selset      = ismember(repelec{4},subjectList);
repelec{2}  = repelec{2}(selset);
repelec{4}  = repelec{4}(selset);

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

try
% modeldataID  = 'freq_spectra-timeseries-bootrep';
prfID        = 'prf-bootrep';
usefulltsxR2 = false;
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
catch exrep
try
% modeldataID  = 'freq_spectra-timeseries-boot';
prfID        = 'prf-boot';
usefulltsxR2 = false;
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
catch
dbstack('-completenames');
error(exrep.identifier, ...
            '%s\nNeed to run %s in advance.',exrep.message,'ecog_APRFF_10d0_analyzePRFboot_REP');
end
end

% model_all_bb_boot = model_all_bb;
% model_all_a_boot  = model_all_a;
prf_all_bb_boot   = prf_all_bb;
prf_all_a_boot    = prf_all_a;

%%% load original pRF data
modeldataID = [];
prfID       = [];
usefulltsxR2 = false;
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;

model_all_bb_orig = model_all_bb;
model_all_a_orig  = model_all_a;
prf_all_bb_orig   = prf_all_bb;
prf_all_a_orig    = prf_all_a;

%% set threshold
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;
thresh_ecc = eclimit;
% thresh_ecc = nan;

elec_ok = ~(prf_all_bb.xval<=threshold_bb | prf_all_a.xval<=threshold_a ...
         | prf_all_bb.ecc >= thresh_ecc | prf_all_a.ecc >= thresh_ecc); 
     
%% visualize parameter
res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

plcol = get(groot,'defaultAxesColorOrder');
     
hF = gobjects(0);

%%
%%
%% pick up bootstrap trials
showboot = 100;
Nboot    = size(prf_all_bb_boot.params,4);

prf_all_bb_bootP = prf_all_bb_boot;
prf_all_a_bootP  = prf_all_a_boot;
if Nboot > showboot
pickidx = (1:showboot)+min(fix((Nboot-showboot)./2),100);

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
end

%% plot averaged pRF location overlay

tarBAND = 'both';

cntcol      = plcol([2,1],:);

el = findtable(prf_all_bb_boot.channels,repelec{:});
for ii=1:numel(repelec{2})
figureTitle = sprintf('pRFplot-%s-%s_loc',prf_all_bb_boot.channels.subject_name(el(ii)),prf_all_bb_boot.channels.name{el(ii)});

opt = [];
opt.plot.fontSize   = FntSiz;
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
p = ecog_plotGridPRF(el(ii), opt, prf_all_a_bootP, prf_all_bb_bootP);
  set(gcf,'MenuBar','none','Position',get(gcf,'Position').*[1,1,2.2,1.8]);
  title('');
  adjustplot(gcf);
  hF = [hF p{:}];

    figureName = sprintf('%s_%s-%s',figureTitle,'booteach-1sd-center',tarBAND);
    if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName),'-vector');    end
 
if issaveplot
opt.plot.options    = {'LineWidth',2.6};
%-- plot bootstrap 1sd w/ center
opt.plot.sigma  = 1;
opt.plot.isavg  = true;
opt.plot.colors = cntcol;
opt.plot.addCenter = 'yes';
p = ecog_plotGridPRF(el(ii), opt, prf_all_a_boot, prf_all_bb_boot);
  set(gcf,'MenuBar','none','Position',get(gcf,'Position').*[1,1,2.2,1.8]);
  title('');
  adjustplot(gcf);
  hF = [hF p{:}];

    figureName = sprintf('%s_%s-%s',figureTitle,'bootavg-1sd-center',tarBAND);
    if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName),'-vector');    end
    
%-- plot original 1sd w/ center
opt.plot.sigma  = 1;
opt.plot.isavg  = true;
opt.plot.colors = cntcol;
opt.plot.addCenter = 'yes';
p = ecog_plotGridPRF(el(ii), opt, prf_all_a_orig, prf_all_bb_orig);
  set(gcf,'MenuBar','none','Position',get(gcf,'Position').*[1,1,2.2,1.8]);
  title('');
  adjustplot(gcf);
  hF = [hF p{:}];

    figureName = sprintf('%s_%s-%s',figureTitle,'orig-1sd-center',tarBAND);
    if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName),'-vector');    end
end

end
    