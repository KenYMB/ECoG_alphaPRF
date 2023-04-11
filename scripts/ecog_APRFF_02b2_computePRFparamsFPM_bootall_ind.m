% Just combine 03c & 03cP with new threshold estimated in 06b
%   without visualization -> 06d
%   Update xR2 with full time-series (no decimation)

% 20201008 Yuasa - Update from 05c
% 20210824 Yuasa - for other datatypes
% 20211118 Yuasa - fix unintended script for prfs.chanidx

%% Define paths and dataset
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
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

%% load analyzePRF & recompute full-ts R2
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
    SetDefault('noecclimit',false);

usefulltsR2   = false;
usefulltsxR2  = false;

ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

%% for bootstrap
va_area = 'wangarea';
if noecclimit,  eclimit = nan;  end       % defult=50 pixel (=radius)
nboot     = 5000;
boottype  = 'bootind';         % 'boot': w/ resampling; 'iter': w/o resampling
robustfit = true;

filename = sprintf('all_prfparams-%sall%s-%s%s_thresh%02d',boottype,R2mode,selectchs,postfix,threshold_bb);
if exist('threshold_a','var')&&~isnan(threshold_a)
    filename = sprintf('%s-%02d',filename,threshold_a); end
if exist('eclimit','var')&&~isnan(eclimit)
    filename = sprintf('%s-ecc%02d',filename,eclimit);     end
if ~robustfit
    filename = sprintf('%s-exactfit',filename,eclimit);     end
filepath = fullfile(SetDefaultAnalysisPath('DAT',prfstatPth),filename);

if exist([filepath '.mat'],'file')
 fprintf('loading %s...',filename);
 load(filepath);
 fprintf('\n');
 nroi = length(rois);
else
 fprintf('computing %s\n',filename);
    
 [~,~,rois]   = ecog_rearrangePRF(prf_params_a,va_area);
 roiidx   = ~ismember(rois,{'none','PHC','SPL1','FEF'});
 rois  = rois(roiidx);
 rois  = [rois;{'V1-V3';'dorsolateral'}];      % add low/high visual area
 nroi  = length(rois);
 nchan = sum(cellfun(@(C) height(C.channels),prf_params_bb));
 
 prfs = [];
 prfs.R2_bb       = cell(nroi,1);
 prfs.R2_a        = cell(nroi,1);
 prfs.ecc_bb      = cell(nroi,1);
 prfs.ecc_a       = cell(nroi,1);
 prfs.ang_bb      = cell(nroi,1);
 prfs.ang_a       = cell(nroi,1);
 prfs.rfsize_bb   = cell(nroi,1);
 prfs.rfsize_a    = cell(nroi,1);
%  prfs.cntr_dist   = cell(nroi,1);
%  prfs.num         = nan(nroi,nboot);
%  prfs.chanidx     = cell(nroi,nboot);
 prfs.num_bb      = nan(nroi,nboot);
 prfs.num_a       = nan(nroi,nboot);
 prfs.chanidx_bb  = cell(nroi,nboot);
 prfs.chanidx_a   = cell(nroi,nboot);
 
 tic;
 for iboot=1:nboot       % 4m for 2,000 iteration
     %% classify pRF results into visual areas (Wang atlas using probabilities)
     [~,prf_all_a]       = ecog_rearrangePRF(prf_params_a,va_area,'prob');
     channels  = prf_all_a.channels;
     [~,prf_all_bb]     = ecog_rearrangePRF(prf_params_bb,va_area,channels);
     
     %% save data
     %-- resampling
     switch boottype
         case {'boot','bootind'}
             boot_keep = randi(nchan,1,nchan);
         case {'iter','iterind'}
             boot_keep = 1:nchan;
         otherwise
             error('''%s'' is unknown parameter',boottype);
     end
     %-- channel selection
     elec_ok_bb = ~(prf_all_bb.xval(boot_keep)<=threshold_bb | prf_all_bb.ecc(boot_keep) >= eclimit);  % use tilde for nan
     elec_ok_a  = ~(prf_all_a.xval(boot_keep)<=threshold_a   | prf_all_a.ecc(boot_keep) >= eclimit);   % use tilde for nan
     %-- computation in roi
     for iroi = 1:nroi
         %-- roi selection
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
                 elec_roi = elec_roi | (channels(boot_keep,:).(va_area)==catrois{jroi});
             end
         else
             elec_roi = channels(boot_keep,:).(va_area)==rois{iroi};
         end
         elec_keep_bb = elec_ok_bb & elec_roi;
         elec_keep_a  = elec_ok_a  & elec_roi;
             
         if any(elec_keep_bb)
             %-- compute parameters
             prfs.R2_bb{iroi}       = cat(2,prfs.R2_bb{iroi},prf_all_bb.xval(boot_keep(elec_keep_bb))');
             prfs.ecc_bb{iroi}      = cat(2,prfs.ecc_bb{iroi},prf_all_bb.ecc(boot_keep(elec_keep_bb))');
             prfs.ang_bb{iroi}      = cat(2,prfs.ang_bb{iroi},prf_all_bb.ang(boot_keep(elec_keep_bb))');
             prfs.rfsize_bb{iroi}   = cat(2,prfs.rfsize_bb{iroi},prf_all_bb.rfsize(boot_keep(elec_keep_bb))');
         end
         if any(elec_keep_a)
             %-- compute parameters
             prfs.R2_a{iroi}        = cat(2,prfs.R2_a{iroi},prf_all_a.xval(boot_keep(elec_keep_a))');
             prfs.ecc_a{iroi}       = cat(2,prfs.ecc_a{iroi},prf_all_a.ecc(boot_keep(elec_keep_a))');
             prfs.ang_a{iroi}       = cat(2,prfs.ang_a{iroi},prf_all_a.ang(boot_keep(elec_keep_a))');
             prfs.rfsize_a{iroi}    = cat(2,prfs.rfsize_a{iroi},prf_all_a.rfsize(boot_keep(elec_keep_a))');
         end
         %-- electrode numbers (save even for invalid condition)
         prfs.num_bb(iroi,iboot)         = sum(elec_keep_bb);
         prfs.num_a(iroi,iboot)          = sum(elec_keep_a);
         prfs.chanidx_bb{iroi,iboot}     = boot_keep(elec_keep_bb);
         prfs.chanidx_a{iroi,iboot}      = boot_keep(elec_keep_a);
     end
 end
 toc;
 
 saveauto(filepath,'rois','prfs','nboot','threshold_bb','threshold_a','eclimit');
end
selroi = 1:nroi;

if robustfit,  fitopt = {'RobustOpts','on'};
else,          fitopt = {'RobustOpts','off'};
end
    
%% for cross-correlations
if ~ismember('prfs_corrx',who('-file',filepath))
    flds  = {'R2_ecc','R2_rfsize','ecc_rfsize'};
    bnds  = {'_bb','_a'};
 
    prfs_corrx =[];
    
    for sbbnd = bnds
         for sbfld = flds
           prfs_corrx.([sbfld{:} sbbnd{:} '_fit'])    = cell(nroi,1);
           prfs_corrx.([sbfld{:} sbbnd{:} '_corr'])   = cell(nroi,1);
           prfs_corrx.([sbfld{:} sbbnd{:} '_R2'])     = cell(nroi,1);
         end

        for iroi = selroi
            for iboot = 1:nboot
                dnum = prfs.(['num' sbbnd{:}])(iroi,iboot);
                doffset = sum(prfs.(['num' sbbnd{:}])(iroi,1:(iboot-1)));
                if dnum > 3
                    for sbfld = flds
                        sbflds = strsplit(sbfld{:},'_');
                        dat1 = prfs.([sbflds{1} sbbnd{:}]){iroi}((1:dnum)+doffset);
                        dat2 = prfs.([sbflds{2} sbbnd{:}]){iroi}((1:dnum)+doffset);

                        fitparms = fitlm(dat1,dat2,fitopt{:});
                        prfs_corrx.([sbfld{:} sbbnd{:} '_fit']){iroi}  = [prfs_corrx.([sbfld{:} sbbnd{:} '_fit']){iroi} fitparms.Coefficients.Estimate];
                        prfs_corrx.([sbfld{:} sbbnd{:} '_corr']){iroi} = [prfs_corrx.([sbfld{:} sbbnd{:} '_corr']){iroi} sqrt(fitparms.Rsquared.Ordinary)];
                        prfs_corrx.([sbfld{:} sbbnd{:} '_R2']){iroi}   = [prfs_corrx.([sbfld{:} sbbnd{:} '_R2']){iroi} fitparms.Rsquared.Adjusted];
                    end
                end
            end
        end
    end
    
    saveauto(filepath,'prfs_corrx','-APPEND');
end
