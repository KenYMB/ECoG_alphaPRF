function [channels] = ecog_summarizeROIs(channels)

% Description: 
%
% [channels] = ecog_summarizeROIs(channels)
% [strctWchannels] = ecog_summarizeROIs(strctWchannels)
% 
% Combine visual areas to reduce ROIs
%
% Input
% - channels     	= channel table or cell-array of channel tables
% - strctWchannels  = structure including channel table or cell-array of
%                     structures such as modeldata and prf
% 
% Concatenate prf results from N subjects, and segrigate into M visual areas

% Dependency: SetDefault, istablefield

% 20200609 - Yuasa
% 20200616 - Yuasa: update for wang FPM
% 20210421 - Yuasa: update for hcp

%% Parameter setting
bensonarea = ["V1","V2","V3","hV4","VO1","VO2","V3a","V3b","LO1","LO2","TO1","TO2","none"];
wangarea = ["V1v","V1d","V2v","V2d","V3v","V3d","hV4","VO1","VO2","PHC1","PHC2","V3a","V3b","LO1","LO2","TO1","TO2","IPS0","IPS1","IPS2","IPS3","IPS4","IPS5","SPL1","FEF","none"];
hcparea  = ["V1","V2","V3","none"];

bensonnew  = ["V1","V2","V3","hV4","VO","V3a","V3b","LO1","LO2","TO","none"];
bensoncomb = [1,2,3,4,[5,5],6,7,8,9,[10,10],11];
wangnew    = ["V1","V2","V3","hV4","VO","PHC","V3a","V3b","LO1","LO2","TO","IPS","SPL1","FEF","none"];
wangcomb   = [[1,1],[2,2],[3,3],4,[5,5],[6,6],7,8,9,10,[11,11],[12,12,12,12,12,12],13,14,15];

%-- check input
cellinp = iscell(channels);
if ~cellinp, channels = {channels}; end
strctinp = isstruct(channels{1});
ninp     = length(channels);
validinp = true;
for ii = 1:ninp
    validinp = validinp && ...
        ((~strctinp && istable(channels{ii})) ||...
         (strctinp && isfield(channels{ii},'channels') && istable(channels{ii}.channels)));
end
assert(validinp,'Unknown data type');

%% Update ROI labels
for ii = 1:ninp
    if strctinp, channel = channels{ii}.channels;
    else,        channel = channels{ii};
    end
    %-- update bensonarea
    if istablefield(channel,'bensonarea')
        channel.bensonarea = categorical(channel.bensonarea);
        for jj = 1:length(bensonarea)
            channel.bensonarea(ismember(channel.bensonarea,bensonarea(jj))) = ...
                bensonnew(bensoncomb(jj));
        end
        channel.bensonarea = categorical(channel.bensonarea,bensonnew);
    end
    %-- update wangarea
    if istablefield(channel,'wangarea')
        channel.wangarea = categorical(channel.wangarea);
        for jj = 1:length(wangarea)
            channel.wangarea(ismember(channel.wangarea,wangarea(jj))) = ...
                wangnew(wangcomb(jj));
        end
        channel.wangarea = categorical(channel.wangarea,wangnew);
    end
    %-- update wangprob
    probidx = contains(channel.Properties.VariableNames,'wangprob');
    if any(probidx)
        probtable = table('Size',[height(channel),0]);
        for jj = 1:(length(wangnew)-1)
            roicomb = find(wangcomb == jj);
            roiidx = ismember(channel.Properties.VariableNames,['wangprob_' wangnew{jj}]);
            for kk = roicomb
                roiidx = roiidx | ismember(channel.Properties.VariableNames,['wangprob_' wangarea{kk}]);
            end
            %%-- sum all probabilities included in the new roi
            probtable.(['wangprob_' wangnew{jj}]) = min(1,sum(channel{:,roiidx},2));
        end
        %%-- replace into new probability table
        channel(:,probidx) = [];
        channel = cat(2,channel(:,1:(find(probidx,1,'first')-1)),probtable,channel(:,find(probidx,1,'first'):end));
    end
    %-- categolical hcp
    if istablefield(channel,'hcparea')
        channel.hcparea = categorical(channel.hcparea,hcparea);
    end
    
    if strctinp, channels{ii}.channels = channel;
    else,        channels{ii} = channel;
    end
end

%-- output
if ~cellinp, channels = channels{1}; end
