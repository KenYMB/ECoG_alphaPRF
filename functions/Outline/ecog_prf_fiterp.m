function [regressor] = ecog_prf_fiterp(data, opts)

% Description: 
%
% [regressor] = ecog_prf_fiterp(data, [opts])
%
% Input
% - data            = Nx1 cell-array of data structure with the following fields:
%   - subject
%   - epochs        = t x events x channels
%   - t
%   - events
%   - channels
%   - fsample
% - opts
%   - baseline_time = (default = [-0.2 0]s)
%   - average       = 'none','sessions','runs'(default),'stimuli' or 'trials'
%   - allowlag      = true or false(default), if allow lag in regression
%   - maxlag        = maximum lags (default = 0.01s, valid if allowlag is true)
%   - issave
%   - outputDir
%   - doplots
%   - plotsavedir
%   - plot          % see ecog_plotTimecourses for this option
%   ---------------------
%   - fileid
%
% Output
% - regressor       = Nx1 cell-array of regressor structure with the following fields:
%   - subject
%   - coef          = channels x events x boots
%   - predictor     = t x 1 x channels
%                       not output all predictor with iteration to reduce
%                       memory usage
%   - t
%   - events
%   - channels
%   - fsample

% Hidden options
% - opts
%   - compute       = [];
%   - skipexist     = [];
%   - fileid        = 'data_erp-coef';
%   - bootstrap     = 0;

% Hidden usage
% - opts.compute    = false;        % load or bypass data
% 
% - opts.bootstrap  = number (p) > 0: apply bootstrapping with p iterations
%                     0(default)    : not apply bootstrapping
% 
% [regressor] = ecog_prf_fiterp(data, [opts])
% [regressor] = ecog_prf_fiterp(subjectList, [opts])

% Dependency: ECoG_utils, SetDefault, saveauto

% 20210413 - Yuasa
% 20210425 - Yuasa: add allowshift option
% 20220222 - Yuasa: load HDgrid information from subjectlist.tsv

%% Set options
%--Define inputs 
% <opts>
SetDefault('opts.average','runs');
SetDefault('opts.baseline_time',[-0.2 0]);
SetDefault('opts.allowlag',false);
SetDefault('opts.maxlag',0.01);
SetDefault('opts.issave',false);
SetDefault('opts.outputDir',fullfile(analysisRootPath, 'Data', 'ERPs'));
SetDefault('opts.doplots',false);
SetDefault('opts.plotsavedir',fullfile(analysisRootPath, 'Figures', 'ERPs'));
SetDefault('opts.plot.XLim',[-0.2 0.8]);
SetDefault('opts.plot.fontSize',16);
% <hidden opts>
SetDefault('opts.bootstrap',0);
SetDefault('opts.compute',[]);
SetDefault('opts.skipexist',[]);
SetDefault('opts.fileid','data_erp-coef');
SetDefault('opts.plotpredata',false);       % just plot <regress-pre>
SetDefault('opts.plot.RotGrid',false);
% <modify opts>
if opts.maxlag==0,  opts.allowlag = false;  end

%-- check inputs and outputs
assert(~isempty(data), 'Please provide the data struct');
if isempty(opts.compute)
    SetDefault('opts.skipexist',true);
    opts.compute = true;
elseif ~opts.compute
    opts.skipexist = false;
else
    SetDefault('opts.skipexist',false);
end
if opts.issave && ~exist(opts.outputDir, 'dir'),     mkdir(opts.outputDir); end
if opts.doplots && ~exist(opts.plotsavedir, 'dir'), mkdir(opts.plotsavedir); end

%-- Correct Subject Information
subjectList_fname = 'subjectlist.tsv';
SbjInfo    = loadSbjInfo(subjectList_fname,'all');
hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');

%% Loop across subjects
cellinput = iscell(data);
if ~cellinput,  data = {data};  end

regressor = cell(size(data));

for ii = 1 : length(data)
    %-- Set data
    if ~opts.compute && ~isstruct(data{ii}) && ischar(data{ii})
        subject     = data{ii};
        isloadfile  = true;
        compute     = false;
    else
        epochs      = data{ii}.epochs;
        idat        = rmfield(data{ii},'epochs');
        subject     = idat.subject;
        isSave      = opts.issave;
        isloadfile  = opts.skipexist;
        compute     = opts.compute;
    end
    average     = opts.average;
    allowlag    = opts.allowlag;
    if allowlag,    maxlag = opts.maxlag;
    else,           maxlag = 0;
    end
    %-- Try to load files
    if isloadfile
        postfix = cnstpostfix(opts.fileid,allowlag,average);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject, opts.fileid,postfix));
        %-- load files from directory
        if exist(filename,'file') || ~compute
            fprintf('[%s] Loading PRF coefficients for subject %s from %s ',mfilename, subject, opts.outputDir);
            idat        = load(filename);
            isSave      = false;
            fprintf('\n');
            compute = false;
        end
    end
    
    %-- HD grid flag
    if hasSbjInfo && istablefield(SbjInfo,'hasHDgrid')
        hasHDgrid = any(strcmpi(SbjInfo.hasHDgrid(ismember(SbjInfo.participant_id,subject)),'yes'));
    else
        hasHDgrid = false;
    end

    %-- Prepare arguments
    events      = idat.events;
    fsample     = idat.fsample;
    isboot      = logical(opts.bootstrap);
    trials      = rmfield(idat,intersect(fieldnames(idat),{'epochs','predictor','coef','t'}));
    trials.time = idat.t;

    %-- Main
    if compute
        fprintf('[%s] Applying regression for subject %s \n',mfilename, subject);
        
        %-- prepare for average
        switch average
            case {'none'},      avg_group = 1:height(events);
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
        
        %-- change dims & take average across repeats (t x events x channels -> channels x t x averaged events)
        evoked      = zeros(height(trials.channels),length(trials.time),length(n_avg));
        events_idx      = zeros(length(n_avg),1);
        for iavg = 1:length(n_avg)
            evoked(:,:,iavg) = mean(permute(epochs(:,avg_group==iavg,:),[3,1,2]),3,'omitnan');
            events_idx(iavg)     = find(avg_group==iavg,1);
        end
        events      = events(events_idx,:);
        if isboot,  n_iter = opts.bootstrap;
        else,       n_iter = 1;
        end
        
        trials.evoked   = evoked; % channels x t x events
        trials.events   = events;
        
        %-- set baseline
        bslidx  = trials.time >= opts.baseline_time(1) & trials.time < opts.baseline_time(2);
        %-- set 0 for prf BLANK, 1 for prf stimuli, and 10+trial_type for other stimuli
        stmlist = ismember(trials.events.task_name,'prf')&~ismember(trials.events.trial_name,'BLANK') + ...
                  ~ismember(trials.events.task_name,'prf').*(trials.events.trial_type+10); 
        
        %-- set new params
        n_chan      = size(evoked,1);
        n_stims     = size(evoked,3);
        coef  = NaN(n_chan,n_stims,n_iter);
        for iter = 1:n_iter
          if isboot
            %-- recompute data_spctr
            for iavg = 1:length(n_avg)
                stimIdx     = find(avg_group==iavg);
                trials.evoked(:,:,iavg) = mean(permute(epochs(:,randsample(stimIdx,length(stimIdx),true),:),[3,1,2]),3,'omitnan');
            end
          end
        %-- regress with EPRs for pRF stimuli
        maxlagpts = round(maxlag * fsample);
        [tmp_coef,predictor] = ecog_regressERP(trials.evoked,bslidx,stmlist,1,maxlagpts,false); % target=1 (PRF)
        coef(:,:,iter) = tmp_coef;
        end
    else
        %-- Load parameters
        coef          = idat.coef;
        predictor     = permute(idat.predictor,[3,1,2]);
        average       = idat.average;
        n_avg         = idat.n_avg;
            SetDefault('idat.allowlag',opts.allowlag);
            SetDefault('idat.maxlag',nan);
        allowlag     = idat.allowlag;
        maxlag       = idat.maxlag;
    end
    
    %-- Plot figures
    if opts.doplots && ~isboot
        trials.evoked   = predictor;     % channels x t x events
        %%-- set trial name 'PRF' & set specs
        trials.events  = events(find(ismember(events.task_name,'prf')&~ismember(events.trial_name,'BLANK'),1),:);
        trials.events.trial_name = {'PRF'};
        eventList   = table(trials.events.trial_name,1,'VariableNames',{'trial_name','stmlist'});
        eventList   = sortrows(eventList,'stmlist');
        eventList   = unique(eventList.trial_name,'stable');
        
        specs = [];
        specs.dataTypes = {'evoked'};
        specs.plot      = opts.plot;
        
        if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
        else,         	whichElectrodes = trials.channels.name;
        end
        
        %%-- plot pre regression
        figureName = sprintf('evoked-prf_%s_predictor', subject);
        
        ecog_plotTimecourses(trials, whichElectrodes, eventList, specs);
        hgexport(gcf, fullfile(opts.plotsavedir, figureName), hgexport('factorystyle'), 'Format', 'png'); close;
          if hasHDgrid
          [~,p]=ecog_plotGridTimecourses(trials, 'GB', eventList, specs);
          for pp=1:length(p)
            if length(p)==1,  hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB',figureName)), hgexport('factorystyle'), 'Format', 'png');
            else,             hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB-%02d',figureName,pp)), hgexport('factorystyle'), 'Format', 'png');
            end
          end
          close(p);
          end
    end
    
    %-- Collect into an output struct
    idat.coef      = coef;                          % channels x events x boots
    idat.predictor = permute(predictor,[2,3,1]);    % t x 1 x channels
    idat.average   = average;
    idat.n_avg     = n_avg;
    idat.allowlag  = allowlag;
    idat.maxlag    = maxlag;
    idat.events    = events;
    regressor{ii} = idat;
    %-- Save out the data
    postfix = cnstpostfix(opts.fileid,allowlag,average);
    if isSave
        fprintf('[%s] Saving data for subject %s to %s \n',mfilename, subject, opts.outputDir);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject, opts.fileid,postfix));
        saveauto(filename,'-struct','idat');
    end

end
if ~cellinput,  regressor = regressor{1};  end
fprintf('[%s] Done! \n',mfilename);
end

%% Sub function
function postfix = cnstpostfix(fileid,allowlag,average)
postfix = '';
%%-- lag
addlag = allowlag && ~contains(fileid,'-lag');
if addlag,   postfix = sprintf('%s-lag',postfix);	end
%%-- average 
postfix      = sprintf('%s_avg-%s',postfix,average);
end
