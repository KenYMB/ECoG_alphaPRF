function [prf_va, prf_all, outarea] = ecog_rearrangePRF(prf_subj,va_atlas,option,option2)

% Description: 
% ECOG_REARRANGEPRF combine cell-array of pRF information structure for 
% each subject (prf_subj) into one pRF information structure (prf_all), and
% then categorize into each visual area (prf_va) based on refed atlas label.
%
% [prf_all] = ecog_rearrangePRF(prf_subj)
%   Without va_atlas, it outputs the concatenated data across subjects
% 
% [prf_va] = ecog_rearrangePRF(prf_subj, refer_atlas)
% [prf_va, prf_all] = ecog_rearrangePRF(prf_subj, refer_atlas)
% [prf_va, prf_all, area_list] = ecog_rearrangePRF(prf_subj, refer_atlas)
%   With va_atlas, it outputs the concatenated data and ROI-based data in which
%   each electrode is assigned to a visual area specified by the va_atlas.
% 
% [prf_va, prf_all] = ecog_rearrangePRF(prf_subj, refer_atlas_fpm, [threshold=0])
%   If refer_atlas is a full-probability map name, each electrode is assigned to
%   ALL visual areas where the probability is over the threshold.
%   The threshold can be specified in [0,1] and the default value is 0.
% 
% [prf_va, prf_all] = ecog_rearrangePRF(prf_subj, refer_atlas, 'prob', [threshold=0])
%   If refer_atlas is not a probabilistic map and data includes full-probability map,
%   each electrode assign to a visual area or 'none' based on the probability.
%   Since it works probability, the outputs are different run-by-run.
%   Electrodes are aggined to a visual area excluding 'none', where
%   their probability of 'none' is under the threshold.
%   If threshold is 1, all electrodes are assigned to one of visual areas.
%   Information of channels.(refer_atlas) in prf_va and prf_all is updated.
% 
% [prf_va, prf_all] = ecog_rearrangePRF(prf_subj, refer_atlas, 'norm', [threshold=0])
%   If refer_atlas is not a probabilistic map and data includes full-probability map,
%   each electrode assign to the visual area where the probability is the largest. 
%   The electrode where the total probability across all visual areas is
%   under the threshold is assinged to none.
%   Electrodes are aggined to the 'none', where their probability of 'none'
%   is over (1 - threshold).
%   If threshold is 0, all electrodes are assigned to one of visual areas.
%   If threshold is 1, the visual areas should not be updated.
%   Information of channels.(refer_atlas) in prf_va and prf_all is updated.
% 
% [prf_va, prf_all] = ecog_rearrangePRF(prf_subj, refer_atlas, channles)
%   You can specify channel talbe, then channel information refers the table 
%   instead of prf_subj{i}.channels.
%   It is usefull if you want to re-use probabilistic assigned visual label. 
% 
%
% Input
% - prf_subj     	= Nx1 cell-array of pRF information structure for each subject
% - refer_atlas    	= 'benson' (default), 'wang' or 'wangprob'
% - threshold       = probability threshold (only valid for probability map)
%
% Output
% - prf_va         	= Mx1 cell-array of pRF information structure for each visual area
% - prf_all        	= pRF information structure
% - area_list     	= Mx1 cell-array of visual area names
% 
% Concatenate prf results from N subjects, and segrigate into M visual areas

% Dependency: SetDefault, istablefield, int2ordinal

% 20200313 - Yuasa
% 20200414 - Yuasa: update for fpm
% 20200417 - Yuasa: 'prob' option
% 20200914 - Yuasa: support 'hcpprob'
% 20200922 - Yuasa: enable to merge model time-series
% 20210714 - Yuasa: 'norm' option
% 20210916 - Yuasa: enable to run without refer_atlas
% 20220106 - Yuasa: enable to merge spectral data

%% parameter setting
%-- convert structure-array to cell-array
if isstruct(prf_subj)
    prf_subj = mat2cell(prf_subj,ones(size(prf_subj)));
elseif istable(prf_subj)    % consider prf_subj as channels itself
    prf_subj = {struct('channels',prf_subj)};
elseif ~iscell(prf_subj)
    prf_subj = {prf_subj};
end

%-- ignore empty subjects
emptysbj = cellfun(@(C) isempty(C.channels),prf_subj);
prf_subj(emptysbj) = [];
    
%-- set ROI names (get category names if channel table has categorical array)
bensonarea = ["V1","V2","V3","hV4","VO1","VO2","V3a","V3b","LO1","LO2","TO1","TO2","none"];
wangarea   = ["V1v","V1d","V2v","V2d","V3v","V3d","hV4","VO1","VO2","PHC1","PHC2","V3a","V3b","LO1","LO2","TO1","TO2","IPS0","IPS1","IPS2","IPS3","IPS4","IPS5","SPL1","FEF","none"];
hcparea    = ["V1","V2","V3","none"];
if istablefield(prf_subj{1}.channels,'bensonarea')&&iscategorical(prf_subj{1}.channels.bensonarea)
    bensonarea = categories(prf_subj{1}.channels.bensonarea);
end
if istablefield(prf_subj{1}.channels,'wangarea')&&iscategorical(prf_subj{1}.channels.wangarea)
    wangarea   = categories(prf_subj{1}.channels.wangarea);
end
if istablefield(prf_subj{1}.channels,'hcparea')&&iscategorical(prf_subj{1}.channels.hcparea)
    hcparea   = categories(prf_subj{1}.channels.hcparea);
end

prfflds = {'ang', 'ecc', 'expt', 'rfsize', 'R2', 'gain', 'xR2', 'xval', 'resnorms', 'numiters', 'meanvol', 'params','testperformance','aggregatedtestperformance'}; % exclude 'noisereg'
prfflds = [prfflds, {'datats'}]; % add fields for model time-series
prfflds = [prfflds, {'spectra','spectra_off'}]; % add fields for spectral data
if isfield(prf_subj{1},'stimulus') && iscell(prf_subj{1}.stimulus) && ~isempty(prf_subj{1}.stimulus) && iscell(prf_subj{1}.stimulus{1})
prfflds = [prfflds, {'stimulus'}]; % add fields for model time-series & has double cell-array stimulus
end
indflds = {'subject','results_xval'};

%-- target visual area
SetDefault('va_atlas','');
skipva = nargout <= 1 && isempty(va_atlas);
fullprob = false;
if startsWith(va_atlas,'benson')
    va_atlas = 'bensonarea';
    va_prob  = 'hcpprob';
    outarea  = bensonarea;
elseif startsWith(va_atlas,'wangprob')
    va_atlas = 'wangprob';
    outarea  = wangarea;
    fullprob = true;
elseif startsWith(va_atlas,'wang')
    va_atlas = 'wangarea';
    va_prob  = 'wangprob';
    outarea  = wangarea;
elseif startsWith(va_atlas,'hcpprob')
    va_atlas = 'hcpprob';
    outarea  = hcparea;
    fullprob = true;
elseif startsWith(va_atlas,'hcp')
    va_atlas = 'hcparea';
    va_prob  = 'hcpprob';
    outarea  = hcparea;
elseif startsWith(va_atlas,'subj')
    va_atlas = 'subject';
elseif ~skipva
    error('unknown atlas');
end

threshold = 0;
donorm    = false;
useprob   = false;
haschan   = false;
if nargin >= 3
    if fullprob,            if ~isempty(option), threshold = option; end
    elseif istable(option), haschan   = ismember(va_atlas, option.Properties.VariableNames);
    elseif ischar(option),  useprob   = ismember(option,{'prob','boot'});
                            donorm    = ismember(option,{'norm'});
           if nargin >= 4 && ~isempty(option2), threshold = option2; end
    end
end
threshold = min(max(threshold,0),1);   % threshold = [0,1]

%% concatenate all subjects
prf_all    = rmfield(prf_subj{1},indflds(isfield(prf_subj{1},indflds)));
nsubjects  = length(prf_subj);
fldList    = prfflds(isfield(prf_all,prfflds));

%-- concatenate fields
for ii = 2:nsubjects
    for jj = reshape(fldList,1,[])
        if strcmp(jj{:},'params'),  idim = 3;
        else,                       idim = 1;
        end
        if ~isempty(prf_all.(jj{:}))
            if strcmp(jj{:},'datats')
              for pp = 1:numel(prf_all.(jj{:}))
                assert(size(prf_all.(jj{:}){pp},2)==size(prf_subj{ii}.(jj{:}){pp},2),...
                    'the %s data has different time series to the %s data',...
                    int2ordinal(ii),int2ordinal(1));
                prf_all.(jj{:}){pp} = ...
                    cat(idim,prf_all.(jj{:}){pp},prf_subj{ii}.(jj{:}){pp});
              end
            else
                prf_all.(jj{:}) = ...
                    cat(idim,prf_all.(jj{:}),prf_subj{ii}.(jj{:}));
            end
        end
    end
end

%-- channels
if haschan
    channelList = option;
else
    channelList = [];
    for ii=1:nsubjects
        channels = prf_subj{ii}.channels;
        if iscell(channels.low_cutoff)
            channels.low_cutoff = nan(height(channels),1);
        end
        if isempty(channelList)
            channelList = channels;
        else
            [~,ia,ib] = intersect(fieldnames(summary(channelList)),fieldnames(summary(channels)),'stable');
            channelList = cat(1,channelList(:,ia),channels(:,ib));
        end
    end
    if istablefield(channelList,'bensonarea')
        channelList.bensonarea   = categorical(channelList.bensonarea,bensonarea);
    end
    if istablefield(channelList,'wangarea')
        channelList.wangarea     = categorical(channelList.wangarea,wangarea);
    end
    channelList.subject_name = categorical(channelList.subject_name);
end
prf_all.channels = channelList;
prf_all.subjects = categories(channelList.subject_name);

%% Rearrange based on va_area
if skipva
%% Output for skipva
prf_va = prf_all;
elseif strcmp(va_atlas,'subject')
%% classify pRF results into visual areas
prf_va    = cell(nsubjects,1);
cpflds    = {'subjects'};
fldList   = setdiff(fieldnames(prf_all),cpflds);

for ii=1:nvisarea
    if fullprob     % Full Probability Map
        va_ch_idx = prf_all.channels.([va_atlas '_' outarea{ii}]) > threshold;
    else            % Maximum Probability Map
        va_ch_idx = prf_all.channels.(va_atlas)==outarea{ii};
    end
    if any(va_ch_idx)
        for jj = reshape(fldList,1,[])
            if ismember(jj,[prfflds,{'channels'}])
                if strcmp(jj{:},'params')
                    prf_va{ii}.(jj{:}) = ...
                        prf_all.(jj{:})(:,:,va_ch_idx);
                elseif strcmp(jj{:},'datats')
                  for pp = 1:numel(prf_all.(jj{:}))
                    prf_va{ii}.(jj{:}){pp} = ...
                        prf_all.(jj{:}){pp}(va_ch_idx,:);
                  end
                else
                    prf_va{ii}.(jj{:}) = ...
                        prf_all.(jj{:})(va_ch_idx,:);
                end
            else
                prf_va{ii}.(jj{:}) = prf_all.(jj{:});
            end
        end
        for jj = reshape(cpflds,1,[])
            prf_va{ii}.(jj{:}) = prf_all.(jj{:});
        end
        prf_va{ii}.areaname = outarea{ii};
        if fullprob     % Full Probability Map
            prf_va{ii}.channels.(va_atlas) = prf_va{ii}.channels.([va_atlas '_' outarea{ii}]);
        end
    end
end
    
else
%% Resample visual area label if need
nchans = height(prf_all.channels);
nvisarea  = length(outarea)-1;
if useprob || donorm
    %-- probability map
    probmat = zeros(nchans,nvisarea);
    for ii=1:nvisarea
      if ismember([va_prob '_' outarea{ii}],prf_all.channels.Properties.VariableNames)
        probmat(:,ii) = prf_all.channels.([va_prob '_' outarea{ii}]);
      else
        probmat(:,ii) = ismember(prf_all.channels.(va_atlas),outarea{ii});
      end
    end
    probmat = [probmat, 1-sum(probmat,2)];
    probmat(probmat<0) = 0;
    
    %-- update visual area label 
    if useprob
        probmat(probmat(:,end)<threshold,end) = 0;
        for ii=1:nchans
            va_boot = randsample(outarea,1,true,probmat(ii,:));
            prf_all.channels.(va_atlas)(ii) = va_boot;
        end
    elseif donorm
        probmat(probmat(:,end)<(1-threshold),end) = 0;
        for ii=1:nchans
            [~,maxidx] = max(probmat(ii,:));
            prf_all.channels.(va_atlas)(ii) = outarea(maxidx);
        end
    end
end

%% classify pRF results into visual areas
prf_va    = cell(nvisarea,1);
cpflds    = {'subjects'};
fldList   = setdiff(fieldnames(prf_all),cpflds);

for ii=1:nvisarea
    if fullprob     % Full Probability Map
        va_ch_idx = prf_all.channels.([va_atlas '_' outarea{ii}]) > threshold;
    else            % Maximum Probability Map
        va_ch_idx = prf_all.channels.(va_atlas)==outarea{ii};
    end
    if any(va_ch_idx)
        for jj = reshape(fldList,1,[])
            if ismember(jj,[prfflds,{'channels'}])
                if strcmp(jj{:},'params')
                    prf_va{ii}.(jj{:}) = ...
                        prf_all.(jj{:})(:,:,va_ch_idx);
                elseif strcmp(jj{:},'datats')
                  for pp = 1:numel(prf_all.(jj{:}))
                    prf_va{ii}.(jj{:}){pp} = ...
                        prf_all.(jj{:}){pp}(va_ch_idx,:);
                  end
                else
                    prf_va{ii}.(jj{:}) = ...
                        prf_all.(jj{:})(va_ch_idx,:);
                end
            else
                prf_va{ii}.(jj{:}) = prf_all.(jj{:});
            end
        end
        for jj = reshape(cpflds,1,[])
            prf_va{ii}.(jj{:}) = prf_all.(jj{:});
        end
        prf_va{ii}.areaname = outarea{ii};
        if fullprob     % Full Probability Map
            prf_va{ii}.channels.(va_atlas) = prf_va{ii}.channels.([va_atlas '_' outarea{ii}]);
        end
    end
end
% %-- Reject empty areas
% if ~useprob&&~haschan
%     emptyarea = cellfun(@isempty,prf_va);
%     prf_va(emptyarea) = [];
%     outarea(emptyarea) = [];
% end
%-- Reject "none"
outarea(end) = [];

end
