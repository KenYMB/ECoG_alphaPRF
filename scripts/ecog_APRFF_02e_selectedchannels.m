% Compare decimation & full-ts, Channel Selection

% 20200526 Yuasa: Separate from ecog_APRFF_10f_channelSelection

%%
% close all; clear all;
% startupToolboxToolbox;
checkPath;

%-- Input & Output path
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
else
prfPth         = 'pRFmodel';
prfstatPth     = 'pRFanalysis';
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%% load analyzePRF & recompute full-ts R2 & set data for bootstrap with half-trials
clear alphaType broadbandType
    
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

usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

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

%% Channel selection with FPM - reverse

threshold = threshold_a;
nboot     = 5000;
%-- compute xR2_a in ROIs (wang_prob)
filename = sprintf('all_chansel-wangprob%s-%s%s_thresh%02d-rev',R2mode,selectchs,postfix,threshold);
filepath = fullfile(SetDefaultAnalysisPath('DAT',prfstatPth),filename);

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