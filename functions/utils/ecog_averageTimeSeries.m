function [freq] = ecog_averageTimeSeries(freq,average,opt)

% Description: 
%
% [freq_sbj] = ecog_averageTimeSeries(freq_sbj,average_type)
% [freq_sbj] = ecog_averageTimeSeries(freq_sbj,average_type,'realign')
% 
% ECOG_AVERAGETIMESERIES take average spectra data across events based on
% average_type.
% 
% If 'realign' is specified, lacked event is filled by NaN so that all
% subjects have the same events field.
% 
% 
% Input & Output
% - freq_sbj     = Nx1 cell-array of spectral structure for each subject
% - average_type = 'none', 'seesions', 'runs', 'stimuli' or 'trials'

% Dependency: <analyzePRF>, SetDefault, replmat, int2ordinal, ecog_averageEvents

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
    %% average
    events      = freq{ii}.events;
    spectra     = freq{ii}.spectra;     % channels x events x f
    %-- skip permute(spectra,[3,2,1]) & permute(data_spctr,[3,2,1])
    [data_spctr, events] = ecog_averageEvents(spectra,events,average,@geomean);
    
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
