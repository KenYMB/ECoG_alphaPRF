function [prfs,prfs_corr,prfs_corrx,rois,nroi,nboot] = ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit,isupdaterois,excl_inaccROIs)
% [prfs,prfs_corr,prfs_corrx,rois,nroi,nboot] = ECOG_PRF_LOADPRFRELATIONS(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,[threshold_a,eclimit,robustfit,isupdaterois,excl_inaccROIs])
%    loads prf relations data
% 
% Outputs:
%   prfs
%   prfs_corr
%   prfs_corrx
%   rois
% 
% Inputs:
%   [fine name related]
%   prfstatPth      = Directory where mat files are saved 
%   R2mode          = 
%   selectchs       = 
%   postfix         = 
% 
%   [prf relations parameter]
%   boottype        = 'boot','iter'
%   threshold_bb   	= threshold for broadband [%]
%   threshold_a     = threshold for alpha [%]
%   eclimit         = threshold in eccentricity [pixel]
%   robustfit       = if treu, applies robust fit (default: true)
% 
%   [load options]
%  isupdaterois     = if true, updates ROI labels (default: true)
%  excl_inaccROIs   = if true, rejects inaccurate ROIs (default: true)

% 20220622 Yuasa

% Dependency: SetDefault (ky_utils)

%% Load pRF relations data
%-- Set parameter
SetDefault('threshold_bb',nan);
SetDefault('threshold_a',nan);
SetDefault('eclimit',nan);
SetDefault('robustfit',true);
SetDefault('isupdaterois',true);
SetDefault('excl_inaccROIs',true);

%-- Get 1st value if selectchs is cell
if iscell(selectchs),   selectchs = selectchs{1};   end

%-- Set file names
filename = sprintf('all_prfparams-%sall%s-%s%s_thresh%02d',boottype,R2mode,selectchs,postfix,threshold_bb);
if ~isnan(threshold_a)
    filename = sprintf('%s-%02d',filename,threshold_a); end
if ~isnan(eclimit)
    filename = sprintf('%s-ecc%02d',filename,eclimit);     end
if ~robustfit
    filename = sprintf('%s-exactfit',filename);     end
filepath = fullfile(SetDefaultAnalysisPath('DAT',prfstatPth),filename);

%-- Load
if exist([filepath '.mat'],'file')
 fprintf('loading %s...',filename);
 load(filepath,'prfs','prfs_corr','prfs_corrx','rois');
 fprintf('\n');
 nroi = length(rois);
 nboot = size(prfs.num,2);
 SetDefault('threshold_a',nan); SetDefault('eclimit',nan);
else
 error('Failed to find %s',filename);
end

%-- update ROI labels
if isupdaterois
rois(ismember(rois,'low'))={'V1-V3'};
rois(ismember(rois,'high'))={'Dorsolateral'};
rois(ismember(rois,'dorsolateral'))={'Dorsolateral'};
end

%-- reject inaccurate ROIs
if excl_inaccROIs
rejthr = nboot .* 2;
for ifld = fieldnames(prfs)'
    for iroi = 1:nroi
        if ~ismember(ifld{:},{'num','chanidx'}) && length(prfs.(ifld{:}){iroi})<rejthr
        prfs.(ifld{:}){iroi} = nan;
        end
    end
end
end

end