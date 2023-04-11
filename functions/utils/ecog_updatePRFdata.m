function    [data,data2,xR2fldname] = ecog_updatePRFdata(data,varargin)

% [modeldata] = ECOG_UPDATEPRFDATA(modeldata [,channelDir])
% [prf] = ECOG_UPDATEPRFDATA(prf [,channelDir])
%   loads channel information from channelDir, and updates channel
%   structure in modeldata or prf.
%   If channelDir is not specified, it trys to load channel data from
%   (analysisRootPath)/Data/Channels.
%   Each channel data should be named as (subject)-channels.mat.
% 
% [modeldata,prf] = ECOG_UPDATEPRFDATA(modeldata,prf [,channelDir][,updatechan=True,updateR2=False,updatexR2=False])
% [modeldata,prf,xR2fldname] = ECOG_UPDATEPRFDATA(modeldata,prf [,channelDir][,updatechan=True,updateR2=False,updatexR2=False])
%   If updatechan is true, it loads channel information from channelDir, and updates channel
%   structure in modeldata and prf.
%   If updateR2 is true, it updates R2 values in prf based on the input
%   modeldata and prf parameters.
%   If updatexR2 is true, it updates cross-validated R2 values in prf based
%   on the input modeldata and prf parameters.
%   
%   modeldata and prf are cell-arrays of time-series structure and pRF information structure
%   xR2fldname is an estimated fieldname of cross-validated R2 values in prf

% Dependency: ecog_computePRFtimeseries

% 20210421 - Yuasa

%% check inputs
%-- set defaults
chanDir = fullfile(SetDefaultAnalysisPath('DAT'),'Channels');
updatechan = true;
updateR2   = false;
updatexR2  = false;

iargin = 1;
%-- check # of data
if nargin>iargin && ( ( iscell(varargin{iargin}) && isstruct(varargin{iargin}{1}) ) || isstruct(varargin{iargin}) )
    ndat   = 2;
    data2  = varargin{1};
    iargin = iargin+1;
    nargoutchk(0,3);
else
    ndat   = 1;
    nargoutchk(0,1);
end
%-- check channel directory
if nargin>iargin
    if ischar(varargin{iargin})||isstring(varargin{iargin})
        chanDir = varargin{iargin};
        iargin = iargin+1;
    elseif iscell(varargin{iargin})||isstruct(varargin{iargin})
        channels = varargin{iargin};
        iargin = iargin+1;
        chanDir = '';
    end
end
%-- check which fileds will be updated
if ndat>1 && nargin>iargin && ~isempty(varargin{iargin})
    updatechan = varargin{iargin};
end
if ndat>1 && nargin>(iargin+1) && ~isempty(varargin{iargin+1})
    updateR2   = varargin{iargin+1};
end
if ndat>1 && nargin>(iargin+2) && ~isempty(varargin{iargin+2})
    updatexR2  = varargin{iargin+2};
end

%-- error check
cellinput = iscell(data);
if ~cellinput,  data = {data};  end
if ndat>1           % modeldata and prf are specified
    if ~cellinput,  data2 = {data2};  end
    assert(length(data)==length(data2), 'modeldata and prf should be a pair');
    for isbj = 1:length(data)
        assert(strcmp(data{isbj}.subject,data2{isbj}.subject), 'modeldata and prf should be a pair');
    end
end
if isempty(chanDir)  % channels is specified as a structure
    if ~cellinput,  channels = {channels};  end
    assert(length(data)==length(channels), 'channels should be a pair with data');
end

%% check cross-validated data
if ndat>1
    if isfield(data2{1},'results_xval')
        xR2fldname = intersect({'aggregatedtestperformance','xval','xR2'},fieldnames(data2{1}.results_xval),'stable');
        if ~isempty(xR2fldname),    xR2fldname = xR2fldname{1};
        elseif updatexR2,           xR2fldname = 'xval';
        end
    elseif updatexR2
        error('pRF results do not include cross-validated results');
    end
end

%% load new channels
if updatechan && ~isempty(chanDir)
    channels = cell(length(data),1);
    for isbj = 1:length(data)
        channels{isbj} = load(fullfile(chanDir,sprintf('%s-channels',data{isbj}.subject)),'channels');
        channels{isbj} = channels{isbj}.channels;
    end
end

%% update data
%-- for data2
if ndat>1 && updatechan
   data2 =  ecog_updatePRFdata(data2,channels);
end

%-- for data
for isbj = 1:length(data)
    %-- update channel labels
    if updatechan
        selchan = ismember(channels{isbj}.name,data{isbj}.channels.name);
        roicontaints = {'benson','wang','hcp','matchednode'};
        
        %-- reset channel fields
        fpmname = fieldnames(data{isbj}.channels);
        data{isbj}.channels(:,startsWith(fpmname,roicontaints)) = [];
        
        %-- copy channel fields
        fpmname = fieldnames(channels{isbj});
        fpmname = fpmname(startsWith(fpmname,roicontaints));
        data{isbj}.channels(:,fpmname)  = channels{isbj}(selchan,fpmname);
    end
    
    %-- update R2 with full time-series
    if updateR2
        [~, ~, cod] = ecog_computePRFtimeseries(data{isbj}.stimulus,data{isbj}.datats,data2{isbj});
        data2{isbj}.R2 = cod;
    end
    %-- update xR2 with full time-series
    if updatexR2
        [~, ~, cod] = ecog_computePRFtimeseries(data{isbj}.stimulus,data{isbj}.datats,data2{isbj}.results_xval);
        data2{isbj}.results_xval.(xR2fldname) = cod;
        data2{isbj}.(xR2fldname) = cod;
    end
end
