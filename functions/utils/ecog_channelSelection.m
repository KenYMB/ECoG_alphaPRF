function [data] = ecog_channelSelection(data,channelType,exclude_FEF,threshold)

% Description: 
%
% [modeldata] = ecog_channelSelection(modeldata,channelType,[exclude_FEF=TRUE,threshold=0.05])
% [results] = ecog_channelSelection(results,channelType,[exclude_FEF=TRUE,threshold=0.05])
% 	extracts data of some specific channels.
% 
% Input
% - modeldata   = Nx1 cell-array of time-series structure
% - results     = Nx1 cell-array of pRF information structure
% - channelType = Specify intended channels
%       'allchs'        : keep all channels (output is the same as input)
%       'visualchs'     : channels assinged to any kinds of atlas OR
%                         with probability over threshold in any atlas
%       'atlaschs'      : channels assinged to any kinds of atlas
%       'bensonchs'     : channels assinged to benson atlas
%       'wangchs'       : channels assinged to wang atlas
%       'wangprobchs'   : channels with total probability over threshold in wang atlas
%       'HDgridchs'     : channels which group are marked as 'HDgrid'
%       'GB*'           : channels with names starts with 'GB'
% - exclude_FEF = true (default) or false
%                 if true, output does not include channels where
%                 area='FEF' or probability at FEF is the largest
% - threshold   = [0,1] (default=0.05), output does not include channels
%                  where the total probability to be assigned a visual area
%                  is lower than threshold 
% 

% Dependency: <analyzePRF>, SetDefault

% 20210915 Yuasa
% 20211124 Yuasa - fix threshold computation
% 20220121 Yuasa - enable combination of channelTypes

%% Parameter setting
narginchk(2,4);
if nargin==3 && ~islogical(exclude_FEF)
    threshold   = exclude_FEF;
    exclude_FEF = true;
end
SetDefault('exclude_FEF',true);
SetDefault('threshold',0.05);
if ~iscell(channelType), channelType = {channelType}; end

%% loop for cell data
cellflg = iscell(data);
if ~cellflg
    data = {data};
end
    
for ii = 1:numel(data)
idat = data{ii};
%% Get channels index
assert(isfield(idat,'channels'),'all data must have channels');
channels = idat.channels;
nchan    = height(channels);

%-- get table field index
bensonfld   = ismember(channels.Properties.VariableNames,'bensonarea');
wangfld     = ismember(channels.Properties.VariableNames,'wangarea');
hcpfld      = ismember(channels.Properties.VariableNames,'hcparea');
wangprobfld = contains(channels.Properties.VariableNames,'wangprob_');
if exclude_FEF
    wangprobfld = wangprobfld & ~contains(channels.Properties.VariableNames,'_FEF');
end
hcpprobfld  = contains(channels.Properties.VariableNames,'hcpprob_');

%-- get channel list
bensonchs   = ~all(ismember(channels{:,bensonfld},'none'),2);
if exclude_FEF
    wangchs     = ~all(ismember(channels{:,wangfld},{'none','FEF'}),2);
else
    wangchs     = ~all(ismember(channels{:,wangfld},'none'),2);
end
hcpchs      = ~all(ismember(channels{:,hcpfld},'none'),2);
wangprobchs = sum(channels{:,wangprobfld},2)>threshold;
hcpprobchs  = sum(channels{:,hcpprobfld},2)>threshold;
selectchs = false(nchan,1);
for itype = 1:length(channelType)
    switch channelType{itype}
        case 'allchs'
            selectchs = selectchs | true(nchan,1);
        case 'visualchs'
            selectchs = selectchs | bensonchs|wangchs|hcpchs|wangprobchs|hcpprobchs;
        case 'atlaschs'
            selectchs = selectchs | bensonchs|wangchs|hcpchs;
        case 'bensonchs'
            selectchs = selectchs | bensonchs;
        case 'wangchs'
            selectchs = selectchs | wangchs;
        case 'hcpchs'
            selectchs = selectchs | hcpchs;
        case 'wangprobchs'
            selectchs = selectchs | wangprobchs;
        case 'hcpprobchs'
            selectchs = selectchs | hcpprobchs;
        otherwise
          if endsWith(channelType{itype},'chs')
            %-- search channels.group
            pattern = strrep(channelType{itype},'chs','');
            selectchs = selectchs | ismember(channels.group,pattern);
          else
            %-- specify channel names with wildcard or error
            assert(contains(channelType{itype},'*'),'Unknown channelType')
            pattern = ['^' regexptranslate('wildcard', channelType{itype}) '$'];
            selectchs = selectchs | ~cellfun(@isempty,regexp(channels.name,pattern,'once'));
          end
    end
end

%% Channel selection for each field
targetflds = {'channels','ang','ecc','expt','rfsize','R2','gain','xR2','xval',...
    'resnorms','numiters','meanvol','params','testperformance',...
    'aggregatedtestperformance','datats','stimulus','spectra','spectra_off'};
targetdats = {'results_xval'};
%-- check subsubfileds
subdats = reshape(fieldnames(idat),1,[]);
subdats = subdats(ismember(subdats,targetdats));
nidat   = 1+length(subdats);
%-- loop for data
for jj = 1:nidat
    if jj==1,   iidat = idat;
    else,       iidat = idat.(subdats{jj-1});
    end
    subflds = reshape(fieldnames(iidat),1,[]);
    for fld = subflds
        if ismember(fld,targetflds) && ~isempty(iidat.(fld{:}))
            if strcmp(fld{:},'datats')
                %-- cell, dim=1
                iidat.(fld{:}) = cellfun(@(C) C(selectchs,:),iidat.datats,'UniformOutput',false);
            elseif strcmp(fld{:},'channels')
                %--- table, dim=1
                iidat.(fld{:})(~selectchs,:) = [];
            elseif strcmp(fld{:},'params') && size(iidat.(fld{:}),3)==nchan
                %-- dim=3
                iidat.(fld{:})(:,:,~selectchs,:,:) = [];
            elseif size(iidat.(fld{:}),1)==nchan
                %--- dim=1
                iidat.(fld{:})(~selectchs,:,:) = [];
            end
        end
    end
    if jj==1,   idat = iidat;
    else,       idat.(subdats{jj-1}) = iidat;
    end
end
%% Output
data{ii} = idat;
end

if ~cellflg
    data = data{1};
end
