function [data] = ecog_prf_regressData(data, opts)

% Description: 
%
% [regressed] = ecog_prf_regressData(data, [opts])
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
%   - baseline_time
%   - epoch_time    = reference period to apply regression
%   - allowlag      = true or false(default), if allow lag in regression
%    - maxlag        = maximum lags (default = 0.01s, valid if allowlag is true)
%    - lag_maxabs    = true or false(default), if estimate time lag with maxabs or max
%   - issave
%   - outputDir
%   - doplots
%   - plotsavedir
%   - plot          % see ecog_plotGridTimecourses for this option
%   ---------------------
%   - fileid
%
% Output
% - regressed       = Nx1 cell-array of data structure with the following fields:
%   - subject
%   - epochs        = t x events x channels
%   - t
%   - events
%   - channels
%   - fsample

% Hidden options
% - opts
%   - compute       = [];
%   - skipexist     = [];
%   - fileid        = 'data_regress';

% Hidden usage
% - opts.compute = false;        % load or bypass data
% 
% [regressed] = ecog_prf_regressData(regressed, [opts])
% [regressed] = ecog_prf_regressData(subjectList, [opts])

% Dependency: ECoG_utils, ecog_regressout, SetDefault, saveauto

% 20200219 - Yuasa
% 20220222 - Yuasa: load HDgrid information from subjectlist.tsv

%% Set options
%--Define inputs 
% <opts>
SetDefaultAnalysisPath('DATA','Preprocessed','opts.outputDir');
SetDefaultAnalysisPath('FIGURE','regression','opts.plotsavedir');
SetDefault('opts.baseline_time',[-0.2 0]);
SetDefault('opts.epoch_time',[-0.2 0.8]);
SetDefault('opts.allowlag',false);
SetDefault('opts.maxlag',0.01);
SetDefault('opts.lag_maxabs',false);
SetDefault('opts.issave',false);
SetDefault('opts.doplots',false);
SetDefault('opts.plot.XLim',opts.epoch_time);
SetDefault('opts.plot.fontSize',16);
% <hidden opts>
SetDefault('opts.compute',[]);
SetDefault('opts.skipexist',[]);
SetDefault('opts.fileid','data_regress');
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
if opts.doplots && ~exist(opts.plotsavedir, 'dir'),  mkdir(opts.plotsavedir); end

%-- Collect Subject Information
subjectList_fname = 'subjectlist.tsv';
SbjInfo    = loadSbjInfo(subjectList_fname,'all');
hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');

%% Loop across subjects
cellinput = iscell(data);
if ~cellinput,  data = {data};  end

for ii = 1 : length(data)
    %-- Set data
    if ~opts.compute && ~isstruct(data{ii}) && ischar(data{ii})
        subject     = data{ii};
        isloadfile  = true;
        compute     = false;
    else
        idat        = data{ii};
        fsample     = idat.fsample;
        subject     = idat.subject;
        isSave      = opts.issave;
        isloadfile  = opts.skipexist;
        compute     = opts.compute;
    end
    allowlag    = opts.allowlag;
    if allowlag,    maxlag = opts.maxlag;
    else,           maxlag = 0;
    end
    %-- Try to load files
    if isloadfile
        postfix     = cnstpostfix(opts.fileid,allowlag);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject, opts.fileid,postfix));
        %-- load files from directory
        if exist(filename,'file') || ~compute
            fprintf('[%s] Loading regressed data for subject %s from %s ',mfilename, subject, opts.outputDir);
            idat        = load(filename);
            isSave      = false;
            fprintf('\n');
            compute = false;
        end
    end
    if compute,         whichplt    = 2;    % pre & post
    elseif isloadfile,  whichplt    = 1;    % post
    else,               whichplt    = double(~opts.plotpredata); % 0:pre, 1:post
    end
    
    %-- HD grid flag
    if hasSbjInfo && istablefield(SbjInfo,'hasHDgrid')
        hasHDgrid = any(strcmpi(SbjInfo.hasHDgrid(ismember(SbjInfo.participant_id,subject)),'yes'));
    else
        hasHDgrid = false;
    end

    trials          = rmfield(idat,{'epochs','t'});
    trials.evoked   = permute(idat.epochs,[3,1,2]); % channels x t x events
    trials.time     = idat.t;

    bslidx  = trials.time >= opts.baseline_time(1) & trials.time < opts.baseline_time(2);
    %-- set 0 for prf BLANK, 1 for prf stimuli, and 10+trial_type for other stimuli
    stmlist = double(ismember(trials.events.task_name,'prf')&~ismember(trials.events.trial_name,'BLANK')) + ...
              ~ismember(trials.events.task_name,'prf').*(trials.events.trial_type+10); 

    %-- Main
    if compute
        %-- regress parameters
        maxlagpts   = round(maxlag * fsample);
        if isempty(opts.epoch_time)
            t_ref   = [];
        else
            t_ref   = trials.time >= opts.epoch_time(1) & trials.time < opts.epoch_time(2);
        end
        negregress  = opts.lag_maxabs;
        
        fprintf('[%s] Applying regression for subject %s \n',mfilename, subject);
        evoked = ecog_regressout(trials.evoked,bslidx,stmlist,maxlagpts,t_ref,negregress);
    else
        %-- Load parameters
        evoked = trials.evoked;
            SetDefault('idat.allowlag',opts.allowlag);
            SetDefault('idat.maxlag',nan);
        allowlag     = idat.allowlag;
        maxlag       = idat.maxlag;
    end
    
    %-- Collect into an output struct
    idat.epochs = permute(evoked,[2,3,1]);
    idat.allowlag = allowlag;
    idat.maxlag   = maxlag;
    data{ii} = idat;
    %-- Save out the data
    postfix = cnstpostfix(opts.fileid,allowlag);
    if isSave
        fprintf('[%s] Saving data for subject %s to %s \n',mfilename, subject, opts.outputDir);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject, opts.fileid, postfix));
        saveauto(filename,'-struct','idat');
    end
    
    %-- Plot figures
    if opts.doplots
        %%-- set trial name 'PRF' & set specs
        trials.events.trial_name(stmlist==1) = {'PRF'};
        eventList   = table(trials.events.trial_name,stmlist,'VariableNames',{'trial_name','stmlist'});
        eventList   = sortrows(eventList,'stmlist');
        eventList   = unique(eventList.trial_name,'stable');
        
        specs = [];
        specs.dataTypes = {'evoked'};
        specs.plot      = opts.plot;
        
        if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
        else,         	whichElectrodes = trials.channels.name;
        end
        
        %%-- plot pre regression
        if whichplt~=1
        figureName = sprintf('evoked-prf%s_%s_regress-pre', postfix,subject);
        
        ecog_plotGridTimecourses(trials, whichElectrodes, eventList, specs);
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
        
        %%-- plot post regression
        if whichplt~=0
        figureName = sprintf('evoked-prf%s_%s_regress-post', postfix,subject);
        
        trials.evoked   = evoked;
        ecog_plotGridTimecourses(trials, whichElectrodes, eventList, specs);
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
    end

end
if ~cellinput,  data = data{1};  end
fprintf('[%s] Done! \n',mfilename);
end

%% Sub function
function postfix = cnstpostfix(fileid,allowlag)
postfix = '';
%%-- lag
addlag = allowlag && ~contains(fileid,'regresslag');
if addlag,   postfix = sprintf('%s_regresslag',postfix);	end
end