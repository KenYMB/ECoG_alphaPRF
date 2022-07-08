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
boottype  = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling
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
 prfs.cntr_dist   = cell(nroi,1);
 prfs.num         = nan(nroi,nboot);
 prfs.chanidx     = cell(nroi,nboot);
 
 tic;
 for iboot=1:nboot       % 4m for 2,000 iteration
     %% classify pRF results into visual areas (Wang atlas using probabilities)
     [~,prf_all_a]       = ecog_rearrangePRF(prf_params_a,va_area,'prob');
     channels  = prf_all_a.channels;
     [~,prf_all_bb]     = ecog_rearrangePRF(prf_params_bb,va_area,channels);
     
     %% save data
     %-- resampling
     switch boottype
         case 'boot'
             boot_keep = randi(nchan,1,nchan);
         case 'iter'
             boot_keep = 1:nchan;
         otherwise
             error('''%s'' is unknown parameter',boottype);
     end
     %-- channel selection
     elec_ok = ~(prf_all_bb.xval(boot_keep)<=threshold_bb | prf_all_a.xval(boot_keep)<=threshold_a ...
         | prf_all_bb.ecc(boot_keep) >= eclimit | prf_all_a.ecc(boot_keep) >= eclimit);  % use tilde for nan
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
         elec_keep = elec_ok & elec_roi;
             
         if any(elec_keep)
             %-- compute parameters
             prfs.R2_bb{iroi}       = cat(2,prfs.R2_bb{iroi},prf_all_bb.xval(boot_keep(elec_keep))');
             prfs.R2_a{iroi}        = cat(2,prfs.R2_a{iroi},prf_all_a.xval(boot_keep(elec_keep))');
             prfs.ecc_bb{iroi}      = cat(2,prfs.ecc_bb{iroi},prf_all_bb.ecc(boot_keep(elec_keep))');
             prfs.ecc_a{iroi}       = cat(2,prfs.ecc_a{iroi},prf_all_a.ecc(boot_keep(elec_keep))');
             prfs.ang_bb{iroi}      = cat(2,prfs.ang_bb{iroi},prf_all_bb.ang(boot_keep(elec_keep))');
             prfs.ang_a{iroi}       = cat(2,prfs.ang_a{iroi},prf_all_a.ang(boot_keep(elec_keep))');
             prfs.rfsize_bb{iroi}   = cat(2,prfs.rfsize_bb{iroi},prf_all_bb.rfsize(boot_keep(elec_keep))');
             prfs.rfsize_a{iroi}    = cat(2,prfs.rfsize_a{iroi},prf_all_a.rfsize(boot_keep(elec_keep))');
         end
         %-- electrode numbers (save even for invalid condition)
         prfs.num(iroi,iboot)         = sum(elec_keep);
         prfs.chanidx{iroi,iboot}     = boot_keep(elec_keep);
     end
 end
 for iroi = 1:nroi
     %-- distance
     cntr_bb = prfs.ecc_bb{iroi}.*exp(1i*prfs.ang_bb{iroi}*pi/180);
     cntr_a  = prfs.ecc_a{iroi}.*exp(1i*prfs.ang_a{iroi}*pi/180);
     prfs.cntr_dist{iroi} = abs(cntr_bb-cntr_a);
 end
 toc;
 
 saveauto(filepath,'rois','prfs','nboot','threshold_bb','threshold_a','eclimit');
end
selroi = 1:nroi;

if robustfit,  fitopt = {'RobustOpts','on'};
else,          fitopt = {'RobustOpts','off'};
end
    
%% for correlations
if ~ismember('prfs_corr',who('-file',filepath))
    flds = {'R2','ecc','rfsize','ang'};
 
    prfs_corr =[];
     for sbfld = flds
       prfs_corr.([sbfld{:} '_fit'])      = cell(nroi,1);
       prfs_corr.([sbfld{:} '_corr'])      = cell(nroi,1);
       prfs_corr.([sbfld{:} '_R2'])      = cell(nroi,1);     
     end

    for iroi = selroi
        for iboot = 1:nboot
            dnum = prfs.num(iroi,iboot);
            doffset = sum(prfs.num(iroi,1:(iboot-1)));
            if dnum > 3
                for sbfld = flds
                    dat1 = prfs.([sbfld{:} '_bb']){iroi}((1:dnum)+doffset);
                    dat2 = prfs.([sbfld{:} '_a']){iroi}((1:dnum)+doffset);

                    fitparms = fitlm(dat1,dat2,fitopt{:});
                    if ismember(sbfld, 'ang'),  Rval = circ_corrcc(deg2rad(dat1'),deg2rad(dat2'));
                    else,                       Rval = sqrt(fitparms.Rsquared.Ordinary);
                    end
                    prfs_corr.([sbfld{:} '_fit']){iroi}  = [prfs_corr.([sbfld{:} '_fit']){iroi} fitparms.Coefficients.Estimate];
                    prfs_corr.([sbfld{:} '_corr']){iroi} = [prfs_corr.([sbfld{:} '_corr']){iroi} sqrt(fitparms.Rsquared.Ordinary)];
                    prfs_corr.([sbfld{:} '_R2']){iroi}   = [prfs_corr.([sbfld{:} '_R2']){iroi} fitparms.Rsquared.Adjusted];
                end
            end
        end
    end
    
    saveauto(filepath,'prfs_corr','-APPEND');
end

%-- for distance
if ismember('prfs_corr',who('-file',filepath)) && ~isfield(prfs_corr,'cntr_dist')
    prfs_corr.cntr_dist      = cell(nroi,1);
    for iroi = selroi
        for iboot = 1:nboot
            dnum = prfs.num(iroi,iboot);
            doffset = sum(prfs.num(iroi,1:(iboot-1)));
            
            dat = prfs.cntr_dist{iroi}((1:dnum)+doffset);
            prfs_corr.cntr_dist{iroi} = [prfs_corr.cntr_dist{iroi} nanmean(dat)];
        end
    end
    saveauto(filepath,'prfs_corr','-APPEND');
end

%% for cross-correlations
if ~ismember('prfs_corrx',who('-file',filepath))
    flds  = {'R2_ecc','R2_rfsize','ecc_rfsize','R2_dist','ecc_dist','rfsize_dist'};
 
    prfs_corrx =[];
     for sbfld = flds
       prfs_corrx.([sbfld{:} '_bb' '_fit'])    = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_bb' '_corr'])   = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_bb' '_R2'])     = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_a' '_fit'])     = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_a' '_corr'])    = cell(nroi,1);
       prfs_corrx.([sbfld{:} '_a' '_R2'])      = cell(nroi,1);
     end

    for iroi = selroi
        for iboot = 1:nboot
            dnum = prfs.num(iroi,iboot);
            doffset = sum(prfs.num(iroi,1:(iboot-1)));
            if dnum > 3
                for sbfld = flds
                    sbflds = strsplit(sbfld{:},'_');
                    dat1 = prfs.([sbflds{1} '_bb']){iroi}((1:dnum)+doffset);
                    if ~strcmp(sbflds{2},'dist')
                      dat2 = prfs.([sbflds{2} '_bb']){iroi}((1:dnum)+doffset);
                    else
                      dat2 = prfs.('cntr_dist'){iroi}((1:dnum)+doffset);
                    end

                    fitparms = fitlm(dat1,dat2,fitopt{:});
                    prfs_corrx.([sbfld{:} '_bb' '_fit']){iroi}  = [prfs_corrx.([sbfld{:} '_bb' '_fit']){iroi} fitparms.Coefficients.Estimate];
                    prfs_corrx.([sbfld{:} '_bb' '_corr']){iroi} = [prfs_corrx.([sbfld{:} '_bb' '_corr']){iroi} sqrt(fitparms.Rsquared.Ordinary)];
                    prfs_corrx.([sbfld{:} '_bb' '_R2']){iroi}   = [prfs_corrx.([sbfld{:} '_bb' '_R2']){iroi} fitparms.Rsquared.Adjusted];
                    
                    dat1 = prfs.([sbflds{1} '_a']){iroi}((1:dnum)+doffset);
                    if ~strcmp(sbflds{2},'dist')
                      dat2 = prfs.([sbflds{2} '_a']){iroi}((1:dnum)+doffset);
                    else
                      dat2 = prfs.('cntr_dist'){iroi}((1:dnum)+doffset);
                    end

                    fitparms = fitlm(dat1,dat2,fitopt{:});
                    prfs_corrx.([sbfld{:} '_a' '_fit']){iroi}  = [prfs_corrx.([sbfld{:} '_a' '_fit']){iroi} fitparms.Coefficients.Estimate];
                    prfs_corrx.([sbfld{:} '_a' '_corr']){iroi} = [prfs_corrx.([sbfld{:} '_a' '_corr']){iroi} sqrt(fitparms.Rsquared.Ordinary)];
                    prfs_corrx.([sbfld{:} '_a' '_R2']){iroi}   = [prfs_corrx.([sbfld{:} '_a' '_R2']){iroi} fitparms.Rsquared.Adjusted];
                end
            end
        end
    end
    
    saveauto(filepath,'prfs_corrx','-APPEND');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% permutation (from 03cP)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Permutation
flds = {'R2','ecc','rfsize','ang'};
for iroi = selroi
    for iboot = 1:nboot
        dnum = prfs.num(iroi,iboot);
        doffset = sum(prfs.num(iroi,1:(iboot-1)));
        dsel = (1:dnum)+doffset;
        dsel1 = dsel(randperm(dnum));
        dsel2 = dsel(randperm(dnum));
        for sbfld = flds
            prfs.([sbfld{:} '_bb']){iroi}(dsel) = prfs.([sbfld{:} '_bb']){iroi}(dsel1);
            prfs.([sbfld{:} '_a']){iroi}(dsel)  = prfs.([sbfld{:} '_a']){iroi}(dsel1);
        end
    end
end


%% for correlations
if ~ismember('prfs_perm',who('-file',filepath))
    flds = {'R2','ecc','rfsize','ang'};
 
    prfs_perm =[];
     for sbfld = flds
       prfs_perm.([sbfld{:} '_fit'])      = cell(nroi,1);
       prfs_perm.([sbfld{:} '_corr'])      = cell(nroi,1);
       prfs_perm.([sbfld{:} '_R2'])      = cell(nroi,1);     
     end

    for iroi = selroi
        for iboot = 1:nboot
            dnum = prfs.num(iroi,iboot);
            doffset = sum(prfs.num(iroi,1:(iboot-1)));
            if dnum > 3
                for sbfld = flds
                    dat1 = prfs.([sbfld{:} '_bb']){iroi}((1:dnum)+doffset);
                    dat2 = prfs.([sbfld{:} '_a']){iroi}((1:dnum)+doffset);

                    fitparms = fitlm(dat1,dat2,fitopt{:});
                    if ismember(sbfld, 'ang'),  Rval = circ_corrcc(deg2rad(dat1'),deg2rad(dat2'));
                    else,                       Rval = sqrt(fitparms.Rsquared.Ordinary);
                    end
                    prfs_perm.([sbfld{:} '_fit']){iroi}  = [prfs_perm.([sbfld{:} '_fit']){iroi} fitparms.Coefficients.Estimate];
                    prfs_perm.([sbfld{:} '_corr']){iroi} = [prfs_perm.([sbfld{:} '_corr']){iroi} sqrt(fitparms.Rsquared.Ordinary)];
                    prfs_perm.([sbfld{:} '_R2']){iroi}   = [prfs_perm.([sbfld{:} '_R2']){iroi} fitparms.Rsquared.Adjusted];
                end
            end
        end
    end
    
    saveauto(filepath,'prfs_perm','-APPEND');
end

%-- for distance
if ~isfield(prfs_perm,'cntr_dist')
    prfs_perm.cntr_dist      = cell(nroi,1);
    for iroi = selroi
        for iboot = 1:nboot
            dnum = prfs.num(iroi,iboot);
            doffset = sum(prfs.num(iroi,1:(iboot-1)));
            %-- distance
            dat = prfs.cntr_dist{iroi}((1:dnum)+doffset);
            prfs_perm.cntr_dist{iroi} = [prfs_perm.cntr_dist{iroi} nanmean(dat)];
        end
    end
    saveauto(filepath,'prfs_perm','-APPEND');
end

%% for cross-correlations
if ~ismember('prfs_permx',who('-file',filepath))
    flds  = {'R2_ecc','R2_rfsize','ecc_rfsize'};
 
    prfs_permx =[];
     for sbfld = flds
       prfs_permx.([sbfld{:} '_bb' '_fit'])    = cell(nroi,1);
       prfs_permx.([sbfld{:} '_bb' '_corr'])   = cell(nroi,1);
       prfs_permx.([sbfld{:} '_bb' '_R2'])     = cell(nroi,1);
       prfs_permx.([sbfld{:} '_a' '_fit'])     = cell(nroi,1);
       prfs_permx.([sbfld{:} '_a' '_corr'])    = cell(nroi,1);
       prfs_permx.([sbfld{:} '_a' '_R2'])      = cell(nroi,1);
     end

    for iroi = selroi
        for iboot = 1:nboot
            dnum = prfs.num(iroi,iboot);
            doffset = sum(prfs.num(iroi,1:(iboot-1)));
            if dnum > 3
                for sbfld = flds
                    sbflds = strsplit(sbfld{:},'_');
                    dat1 = prfs.([sbflds{1} '_bb']){iroi}((1:dnum)+doffset);
                    if ~strcmp(sbflds{2},'dist')
                      dat2 = prfs.([sbflds{2} '_bb']){iroi}((1:dnum)+doffset);
                    else
                      dat2 = prfs.('cntr_dist'){iroi}((1:dnum)+doffset);
                    end

                    fitparms = fitlm(dat1,dat2,fitopt{:});
                    prfs_permx.([sbfld{:} '_bb' '_fit']){iroi}  = [prfs_permx.([sbfld{:} '_bb' '_fit']){iroi} fitparms.Coefficients.Estimate];
                    prfs_permx.([sbfld{:} '_bb' '_corr']){iroi} = [prfs_permx.([sbfld{:} '_bb' '_corr']){iroi} sqrt(fitparms.Rsquared.Ordinary)];
                    prfs_permx.([sbfld{:} '_bb' '_R2']){iroi}   = [prfs_permx.([sbfld{:} '_bb' '_R2']){iroi} fitparms.Rsquared.Adjusted];
                    
                    dat1 = prfs.([sbflds{1} '_a']){iroi}((1:dnum)+doffset);
                    if ~strcmp(sbflds{2},'dist')
                      dat2 = prfs.([sbflds{2} '_a']){iroi}((1:dnum)+doffset);
                    else
                      dat2 = prfs.('cntr_dist'){iroi}((1:dnum)+doffset);
                    end

                    fitparms = fitlm(dat1,dat2,fitopt{:});
                    prfs_permx.([sbfld{:} '_a' '_fit']){iroi}  = [prfs_permx.([sbfld{:} '_a' '_fit']){iroi} fitparms.Coefficients.Estimate];
                    prfs_permx.([sbfld{:} '_a' '_corr']){iroi} = [prfs_permx.([sbfld{:} '_a' '_corr']){iroi} sqrt(fitparms.Rsquared.Ordinary)];
                    prfs_permx.([sbfld{:} '_a' '_R2']){iroi}   = [prfs_permx.([sbfld{:} '_a' '_R2']){iroi} fitparms.Rsquared.Adjusted];
                end
            end
        end
    end
    
    saveauto(filepath,'prfs_permx','-APPEND');
end
