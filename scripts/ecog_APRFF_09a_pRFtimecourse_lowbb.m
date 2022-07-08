%% ECoG Alpha pRF (simple broadband computation)
% plot bootstrap results
% for figure 4

% 20201215
% 20210224 - update for finalize figure 4
% 20210508 - use another y axis
% 20211118 - for low broadband
% %% without ERP %%

%% prefix
close all; clear all;
% startupToolboxToolbox;
%% Define paths and dataset
run_checkPath;
%-- Input & Output path
plotsavePth    = 'pRF-lowbb';
prfPth         = 'pRFmodel';
%-- Set save figure dirctory
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth));
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

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

%%% load original pRF data
modeldataID = [];
prfID       = [];
usefulltsxR2 = false;
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;

[modeldata_bbl, prf_params_bbl] = ecog_prf_loadprfs(subjectList,'bbL',prfPth,modeldataID,prfID,average,smoothingMode,smoothingN,prfmodel,gaussianmode,selectchs,selectch_exFEF,selectch_thresh,usefulltsR2,usefulltsxR2,skipsummarizeROIs);
[model_all_bbl, prf_all_bbl]    = ecog_prf_mergeprfs(modeldata_bbl,prf_params_bbl,va_area,usexvalparams,rearropt);

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
     
close all;

%% plot averaged pRF time-series

repelec = {'name',{'Oc18','GB103','Oc17','GB102'},'subject_name',{'p02','p10','p02','p10'}};
% tarBAND = 'broadband';
% tarBAND = 'alpha';

% for tarBANDs = {'broadband', 'alpha-noflip','low-broadband'}
for tarBANDs = {'low-broadband'}
tarBAND = tarBANDs{:};

switch tarBAND
    case {'broadband'}
        mdlall      = model_all_bb;
        prfall      = prf_all_bb;
        cntcol      = plcol(1,:);
        tname  = 'Broadband';
    case {'alpha','alpha-noflip'}
        mdlall      = model_all_a;
        prfall      = prf_all_a;
        cntcol      = plcol(2,:);
        tname  = 'Alpha';
    case {'low-broadband'}
        mdlall      = model_all_bbl;
        prfall      = prf_all_bbl;
        cntcol      = plcol(3,:);
        tname  = 'Low Broadband';
end

el = findtable(prfall.channels,repelec{:});
for ii=1:numel(repelec{2})
figureTitle = sprintf('pRFplot-%s-%s_ts',prf_all_bb.channels.subject_name(el(ii)),prf_all_bb.channels.name{el(ii)});

opt = [];
opt.plot.fontSize   = 18;
opt.plot.addWangToTitle   = 'no';
opt.plot.addBensonToTitle = 'no';
opt.plot.addSbjToTitle    = 'yes';
opt.plot.addR2ToTitle     = 'yes';
opt.plot.model_options    = {'Color',cntcol};

if contains(tarBAND,'-noflip')
    opt.plot.reverseYaxis = 'no';
end

% boundary = [0 28 40 56 84 96 112 140 152 168 196 208 224]+.5;     % for all event change
boundary = [0 40 56 96 112 152 168 208 224]+.5;         % for BLANK

%-- plot w/o polynomial projection
bndfactor = 1/224*75;
opt.plot.boundaries = boundary.*bndfactor;
opt.plot.options = {'XTick', boundary.*bndfactor,'XTickLabel',{}};
opt.plot.XLim    = [floor(boundary(1)) ceil(boundary(end))].*bndfactor;
    
opt.skipprojection = 'yes';
ecog_plotGridPRFts(mdlall.datats,mdlall.stimulus,prfall,el(ii), opt);
  set(gcf,'MenuBar','none','Position',get(gcf,'Position').*[1,1,2.1,1.6]);
  R2val = get(gca,'Title'); R2val = R2val.String(end-6:end-2);
  title(sprintf('\\fontsize{22}\\color[rgb]{%g %g %g}%s (%s)',cntcol,tname,R2val));
  adjustplot(gcf);
  
switch tarBAND
    case {'broadband','low-broadband'}
        set(gca,'YTickLabel',num2str(get(gca,'YTick')'./100+1));
    case {'alpha','alpha-noflip'}
%         set(gca,'YTick',log10(round(10.^get(gca,'YTick'),1,'significant')));
        if min(yticks)>-1
            set(gca,'YTick',log10([1/8 1/4 1/2 1 2 4 8]));
            set(gca,'YTickLabel',{'1/8','1/4','1/2','1','2','4','8'});
%             set(gcf,'Position',get(gcf,'Position').*[1,1,1.025,1]);
        else
            set(gca,'YTick',log10([.1 .3 1 3 9]));
            set(gca,'YTickLabel',num2str(10.^get(gca,'YTick')'));
        end
end
    
    %-- add gray-box for blank
    pb=plotbox(gca,mat2cell(boundary(2:end).*bndfactor,1,2*ones(1,4)),ylim,'gray','FaceAlpha',0.5);
    uistack(pb,'bottom')
    
      figureName = sprintf('%s_%s-%s',figureTitle,'orig-noproj-showblank-ratio',tarBAND);
      savefigauto(gcf, fullfile(plotsavedir, figureName));
      
        %-- show boundaries
        boundaries = ([28 84 140 196]+0.5).*bndfactor;
        hl = straightline(boundaries,'v','k--');
        uistack(hl,'bottom')

          figureName = sprintf('%s_%s-%s',figureTitle,'orig-noproj-showbound-ratio',tarBAND);
          savefigauto(gcf, fullfile(plotsavedir, figureName));
    
end
end

