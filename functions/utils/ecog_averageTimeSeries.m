function [freq] = ecog_averageTimeSeries(freq,average,opt)

% Description: 
%
% [freq_sbj] = ecog_averageTimeSeries(freq_sbj,average_type)
% [freq_sbj] = ecog_averageTimeSeries(freq_sbj,average_type,'realign')
% 
% ECOG_AVERAGETIMESERIES take average spectra data across events based on
% average_type.
% 
% insert NaNs at the end of each time series to
% have the same time length across input data. Useful to merge data time
% series across subjects.
% If realign_events is true, the 'stimulus' field is replicated to have
% the same length as the numbers of channels.
% 
% 
% Input & Output
% - freq_sbj     = Nx1 cell-array of spectral structure for each subject
% - average_type = 'none', 'seesions', 'runs', 'stimuli' or 'trials'

% Dependency: <analyzePRF>, SetDefault, replmat, int2ordinal

% 20211230 Yuasa

%% parameter setting
narginchk(2,3);
if nargin < 3
    isrealign = false;
elseif strcmpi(opt,'realign')
    isrealign = true;
else
    error('Invalid option');
end
assert(ismember(average,{'none', 'seesions', 'runs', 'stimuli', 'trials'}),...
    'Unknown average_type');
    
%-- convert structure-array to cell-array
if isstruct(freq)
    freq = mat2cell(freq,ones(size(freq)));
elseif ~iscell(freq)
    freq = {freq};
end

%-- check data validity & get maximum data size
ndat = numel(freq);
for ii = 1:ndat
    assert(all(ismember({'spectra','events'},fieldnames(freq{ii}))), ...
        'Invalid input (%s data)',int2ordinal(ii));
    assert(size(freq{ii}.spectra,2)==size(freq{ii}.events,1), ...
        'Invalid input (%s data)',int2ordinal(ii));
end

eventflds = {'task_name','run_name','stim_file_index','trial_type','trial_name'};
event_tmp = table([],[],[],[],[],'VariableNames',eventflds);

%% take average
for ii = 1:ndat
    %% get parameters
    f           = freq{ii}.f;
    channels    = freq{ii}.channels;
    events      = freq{ii}.events;
    spectra     = freq{ii}.spectra;
    
    %% average
    %-- prepare for average
    switch average
        case {'none'},      avg_group = [1:height(events)]';
        case {'sessions'},  avg_group = findgroups(events.task_name,events.run_name,events.stim_file_index);
        case {'runs'},      avg_group = findgroups(events.task_name,events.stim_file_index);
        case {'stimuli'},   avg_group = findgroups(events.task_name,events.trial_type);
        case {'trials'},    avg_group = ones(height(events),1);
            bslIndex  = contains(events.trial_name, 'BLANK');
            avg_group(bslIndex) = 2;
        otherwise,          error('''%s'' is unknown average type',average);
    end
    n_avg     = groupcounts(avg_group);
    if ~all(diff(n_avg)==0),    warning('[%s] %s do not consist with the same time series',mfilename,average);  end
    
    %-- take average across repeats (chan x events x f)
    data_spctr      = zeros(height(channels),length(n_avg),length(f));
    events_idx      = zeros(length(n_avg),1);
    for iavg = 1:length(n_avg)
        data_spctr(:,iavg,:) = geomean(spectra(:,avg_group==iavg,:),2,'omitnan');
        events_idx(iavg)     = find(avg_group==iavg,1);
    end
    events      = events(events_idx,:);
    
    %% put back parameters
    switch average
        case {'none','sessions'}
            events.run_name = cellfun(@(t) sprintf('%02d',str2double(t)),events.run_name,'UniformOutput',false);
        otherwise
            events.run_name = cellfun(@(t) sprintf('%02d',t),num2cell(findgroups(events.run_name)),'UniformOutput',false);
    end
    freq{ii}.events      = events;
    freq{ii}.spectra     = data_spctr;
    
    event_tmp = vertcat(event_tmp, events(:,eventflds));
end
event_tmp     = unique(event_tmp);
ntim          = size(event_tmp,1);

%% realign events (insert NaNs)
if isrealign
  for ii = 1:ndat
    eventflg    = ismember(event_tmp,freq{ii}.events(:,eventflds));
    sizdat      = size(freq{ii}.spectra);
    dummydat    = cast(nan(copymatrix(sizdat,2,ntim)),class(freq{ii}.spectra));
    
    freq{ii}.spectra    = copymatrix(dummydat,eventflg,2,freq{ii}.spectra);
    freq{ii}.events     = event_tmp;
    freq{ii}.events.subject(:)   = {freq{ii}.subject};
  end
end
