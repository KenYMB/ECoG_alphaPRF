function [freq] = ecog_prf_crossspectra(data, opts)

% Description: 
%
% [xfreq] = ecog_prf_crossspectra(data, [opts])
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
%   - target_time
%   - target_freq   = [fmin fmax], frequency bounds for output
%        or
%     targetBAND    = 'alpha' or 'broadband'
%   - excl_freq     = Mx2 array to exclude frequency bands
%   - stimNames     = cell-strings array    % trial names for analysis (important to save memory)
%                     'SHUFFLE" suffle the event list for each electrode and output the average
%   - channels
%   - issave
%   - outputDir
%   ---------------------
%   - fileid
%
% Output
% - xfreq           = Nx1 cell-array of data structure with the following fields:
%   - subject
%   - crossspectra  = channels x channels x events x f
%   - f
%   - events
%   - channels

% Hidden usage
% - opts.compute = false;
% 
% [xfreq] = ecog_prf_crosssspectra(xfreq, [opts])
% [xfreq] = ecog_prf_crosssspectra(subjectList, [opts])

% Dependency: <ECoG_utils>, ecog_regressout, SetDefault, cellstrfind, saveauto

% 20201111 - Yuasa
% 20210512 - Yuasa: update to specify target band
% 20220222 - Yuasa: load recording site information from subjectlist.tsv
% 20230320 - Yuasa: add hidden options and shuffle mode

%% Set options
%--Define inputs 
% <opts>
SetDefaultAnalysisPath('DATA','xSpectrum','opts.outputDir');
SetDefault('opts.target_time',[0 0.5]);
SetDefault('opts.channels',[]);
SetDefault('opts.issave',false);
SetDefault('opts.target_freq',[]);
SetDefault('opts.excl_freq',[]);
SetDefault('opts.targetBAND','');
SetDefault('opts.stimNames',{},'cell');
% <hidden opts>
SetDefault('opts.gammafit',false);      % just use to set default value in target_freq
SetDefault('opts.f_wintyp',@hann);      % window for cpsd; default is hann window
SetDefault('opts.f_winlng',0.2);        % window length in time; default is 0.2s (=Fs/5)
SetDefault('opts.f_winovl',0.5);        % overlap ratio of window; default is 0.5 (=50%) of window length
SetDefault('opts.inclauto',false);
SetDefault('opts.calcmode','default');  % 'default'=cpsd or 'coh'=mscohere
SetDefault('opts.isavg',[]);            % if true, take average across events
SetDefault('opts.taskNames','prf','cell');
SetDefault('opts.compute',[]);
SetDefault('opts.skipexist',[]);
if ismember(lower(opts.calcmode),{'coh','mscoh','mscohere'})
SetDefault('opts.fileid','freq_coherence');
else
SetDefault('opts.fileid','freq_crossspectra');
end
SetDefault('opts.allowlag',false);                      % just for loading files
SetDefault('opts.maxlag',nan);                          % just for saving and loading files

% <target_freq>
if isempty(opts.target_freq) && isempty(opts.targetBAND)
    opts.target_freq = [0 inf];
    opts.targetBAND  = 'all';
elseif ~isempty(opts.target_freq)
    if isempty(opts.targetBAND),  opts.targetBAND  = 'other';  end
    if numel(opts.target_freq)==1
        opts.target_freq = [opts.target_freq opts.target_freq];
    elseif numel(opts.target_freq)~=2
        error('opts.target_freq is invalid');
    end
end

% <StimNames>
isshuffle = ~isempty(opts.stimNames) && ~isempty(opts.stimNames{1}) && strcmpi(opts.stimNames{1},'shuffle');
if isempty(opts.stimNames) || isempty(opts.stimNames{1}) || strcmpi(opts.stimNames{1},'all')
    opts.stimNames = {'*'};
elseif isshuffle && numel(opts.stimNames)==1
    opts.stimNames(2) = {'*'};
end
if isempty(opts.isavg),  opts.isavg = isshuffle;  end   % take average only if shuffled

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

%-- Collect Subject Information
subjectList_fname = 'subjectlist.tsv';
SbjInfo    = loadSbjInfo(subjectList_fname,'all');
hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');

%% Loop across subjects
cellinput = iscell(data);
if ~cellinput,  data = {data};  end

freq = cell(size(data));

for ii = 1 : numel(data)
    ifreq = [];
    %-- Set data
    if ~opts.compute && ~isstruct(data{ii}) && ischar(data{ii})
        subject     = data{ii};
        isloadfile  = true;
        compute     = false;
        %-- Takeover parameters
        allowlag     = opts.allowlag;
        maxlag       = opts.maxlag;
    else
        idat        = data{ii};
        subject     = idat.subject;
        isSave      = opts.issave;
        isloadfile  = opts.skipexist;
        compute     = opts.compute;
        %-- Takeover parameters
        SetDefault('idat.allowlag',opts.allowlag);
        SetDefault('idat.maxlag',opts.maxlag);
        allowlag     = idat.allowlag;
        maxlag       = idat.maxlag;
    end
    %-- Recording site
    if hasSbjInfo && istablefield(SbjInfo,'site')
        RECsite = SbjInfo.site(ismember(SbjInfo.participant_id,subject));
        RECsite = RECsite{1};
    else
        RECsite = 'unknown';
    end
    %-- Other parameters
    targetBAND  = opts.targetBAND;
    taskNames   = opts.taskNames;
    stimNames   = opts.stimNames;
    target_time = opts.target_time;
    [target_freq,excl_freq] = setfreq(RECsite,targetBAND,opts.target_freq,opts.excl_freq,opts.gammafit);
    %-- Try to load files
    if isloadfile
        postfix = cnstpostfix(opts.fileid,targetBAND,target_freq,stimNames,allowlag);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        %-- load files from directory
        if exist(filename,'file') || ~compute
            fprintf('[%s] Loading cross-spectra for subject %s from %s ',mfilename, subject, opts.outputDir);
            idat        = load(filename);
            isSave      = false;
            fprintf('\n');
            compute = false;
        end
    end
    %-- Main
    if compute
        fprintf('[%s] Computing cross-spectra for subject %s \n',mfilename, subject);

        %-- Prepare arguments
        epochs      = permute(idat.epochs,[3,2,1]); % channels x events x t
        t           = idat.t;
        srate       = idat.fsample;
        events      = idat.events;
        channels    = idat.channels;
        
        %-- Select channels
        if isempty(opts.channels)
            chanidx = true(height(channels),1);
        elseif iscellstr(opts.channels) || ischar(opts.channels)
            chanidx = cellstrfind(channels.name,opts.channels);
        else
            chanidx  = opts.channels;
        end
        epochs   = epochs(chanidx,:,:);
        channels = channels(chanidx,:);

        %-- Select the prf events (baseline is BLANK in prf sessions)
        trialIndex  = cellstrfind(events.task_name, taskNames,1) & cellstrfind(events.trial_name, stimNames,1);
        assert(sum(trialIndex)>0, 'No trial is selected');
        epochs      = epochs(:,trialIndex,:);
        events      = events(trialIndex,:);
        if isshuffle
            epochs  = permute(epochs,[3,2,1]); % t x events x channels (permute again to keep order in t)
            [~, trialShuffle]=sort(rand(size(epochs,2:ndims(epochs))),1);
            epochs  = reshape(epochs(:,trialShuffle+(0:size(epochs,2):(prod(size(epochs,2:ndims(epochs)))-1))),size(epochs));
            epochs  = permute(epochs,[3,2,1]); % channels x events x t
            events.trial_type(:) = 0;
            events.stim_file_index(:) = 0;
            events.trial_name(:) = {'SHUFFLE'};
        end

        %-- Calculate spectra
        fft_t   = t>target_time(1) & t<=target_time(2); % time segment for spectrum
        fft_w   = window(opts.f_wintyp,round(srate.*opts.f_winlng)); % window width   % window(@hann,round(srate/5))
        fft_ov  = round(srate.*opts.f_winlng.*opts.f_winovl); % overlap               % round(srate/10)
        reg_erp = 0; % 1 to regress erp out, 0 not to, if yes, make sure to baseline correct first
        if opts.inclauto,   outtype = 'sparse+';
        else,               outtype = 'sparse';
        end
        %%-- stimuli epochs
        [f,crossspectra] = ...
            ecog_crossspectra(epochs,[],fft_w,fft_t,fft_ov,srate,reg_erp,target_freq,outtype,opts.calcmode);
        isfull = ndims(crossspectra) > 3;
        %%-- exclude frequencies
        f_rej  = false(size(f));
        for kk = 1:size(excl_freq)
            f_rej = f_rej | (f>=excl_freq(kk,1) & f<=excl_freq(kk,2));
        end
        f(f_rej) = [];
        if isfull
            crossspectra(:,:,:,f_rej) = [];
                chancmb = 1:height(channels);
        else
            crossspectra(:,:,f_rej) = [];
                sparsemat  = tril(true(height(channels)),-~opts.inclauto);
                chancmb = [];
                [chancmb(:,2), chancmb(:,1)] = find(sparsemat);
        end

        %-- Fill nan for flat signals
        if isfull
            [nanel1,naneltrl] = find(var(crossspectra,0,4)==0);
            nanel2 = mod(naneltrl-1,size(crossspectra,2))+1;
            nantrl = ceil(naneltrl/size(crossspectra,2));
            for jj = 1:length(nanel1)
                crossspectra(nanel1(jj),nanel2(jj),nantrl(jj),:) = nan;
            end
        else
            [nanel1,nantrl] = find(var(crossspectra,0,3)==0);
            for jj = 1:length(nanel1)
                crossspectra(nanel1(jj),nantrl(jj),:) = nan;
            end
        end

        %-- Take average across events
        if opts.isavg
            if isfull, avgdim = 3;
            else,      avgdim = 2;
            end
            events(2:end,:) = [];
            events.trial_type = 0;
            events.stim_file_index = 0;
            events.trial_name = {'AVERAGE'};
            if all(isreal(crossspectra),'all')
                crossspectra = cat(avgdim,mean(crossspectra,avgdim,'omitnan'),...
                                          prctile(crossspectra,[16 84],avgdim),...
                                          prctile(crossspectra,[2.5 97.5],avgdim),...
                                          prctile(crossspectra,[0.15 99.85],avgdim));
                events = repmat(events,4,1);
                events.trial_name(2:end) = {'68% PERCENTILE','95% PERCENTILE','99.7% PERCENTILE'};
            else
                crossspectra = cat(avgdim,mean(crossspectra,avgdim,'omitnan'));
            end
        end

    else
        %-- Load parameters
        crossspectra    = idat.crossspectra;
        f               = idat.f;
        channels        = idat.channels;
        events          = idat.events;
        SetDefault('idat.targetBAND',opts.targetBAND);
        SetDefault('idat.target_freq',opts.target_freq);
        SetDefault('idat.excl_freq',opts.excl_freq);
        if isempty(idat.targetBAND)
          targetBAND      = 'other';
          target_freq     = [0 inf];
          excl_freq       = [];
        else
          targetBAND      = idat.targetBAND;
          [target_freq,excl_freq] = setfreq(subject,targetBAND,idat.target_freq,idat.excl_freq,opts.gammafit);
        end
        if isfield(idat,'chancmb')
          chancmb         = idat.chancmb;
        else
          if ndims(crossspectra) > 3
              chancmb = 1:height(channels);
          else
              [chanx,chany]=meshgrid(1:height(channels),1:height(channels));
              chanl = true(size(chanx));
              chancmb = [chanx(tril(chanl,-~opts.inclauto)), chany(tril(chanl,-~opts.inclauto))];
          end
        end
            SetDefault('idat.target_time',opts.target_time);
            SetDefault('idat.taskNames',opts.taskNames);
            SetDefault('idat.stimNames',opts.stimNames);
        target_time  = idat.target_time;
        taskNames    = idat.taskNames;
        stimNames    = idat.stimNames;
    end
    
    %-- Collect into an output struct
    ifreq.subject        = subject;
    ifreq.crossspectra   = crossspectra;
    ifreq.f              = f;
    ifreq.events         = events;
    ifreq.channels       = channels;
    ifreq.target_time    = target_time;
    ifreq.targetBAND     = targetBAND;
    ifreq.target_freq    = target_freq;
    ifreq.excl_freq      = excl_freq;
    ifreq.chancmb        = chancmb;
    ifreq.taskNames      = taskNames;
    ifreq.stimNames      = stimNames;
    ifreq.allowlag       = allowlag;
    ifreq.maxlag         = maxlag;
    freq{ii} = ifreq;
    
    %-- Save out the data
    postfix = cnstpostfix(opts.fileid,targetBAND,target_freq,stimNames,allowlag);
    if isSave
        fprintf('[%s] Saving cross-spectra for subject %s to %s \n',mfilename, subject, opts.outputDir);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        saveauto(filename,'-struct','ifreq');
    end
    
end
if ~cellinput,  freq = freq{1};  end
fprintf('[%s] Done! \n',mfilename);
end

%% Sub function
function postfix = cnstpostfix(fileid,targetBAND,target_freq,stimNames,allowlag)
postfix = '';
%%-- lag
addlag = allowlag && ~contains(fileid,'regresslag');
if addlag,   postfix = sprintf('%s_regresslag',postfix);	end
%%-- freq
if strcmpi(targetBAND,'other')
    postfix = sprintf('%s-%g-%gHz',postfix,min(target_freq),max(target_freq));
else
    postfix = sprintf('%s-%s',postfix,targetBAND);
end
%-- stim
isshuffle = ~isempty(stimNames) && ~isempty(stimNames{1}) && strcmpi(stimNames{1},'shuffle');
if isshuffle
    postfix     = [postfix '_' stimNames{1}];
    stimNames(1) = [];
end
if ~isempty(stimNames) && ~ismember('*',stimNames)
    if isshuffle,   postfix     = [postfix '-'];
    else,           postfix     = [postfix '_'];
    end
    if length(stimNames)==1
        postfix     = [postfix strrep(stimNames{1},'*','')];
    else
        postfix     = [postfix strrep(stimNames{1},'*','') '-' strrep(stimNames{end},'*','')];
    end
end
end

function [target_freq,excl_freq] = setfreq(RECsite,targetBAND,target_freq,excl_freq,gammafit)
if isempty(target_freq)
    switch lower(targetBAND)
        case {'alpha'}
            target_freq = [6 17];
        case {'broadband'}
            if gammafit
                target_freq = [30 194];
            else
                target_freq = [70 180];
            end
        case {'lowbroadband'}
            target_freq = [3 26];   % lowbb
%             target_freq = [3 34];   % lowbb-wide
            excl_freq   = [excl_freq; [5.1 22]];   % exclude alpha-bump
%             excl_freq   = [excl_freq; [15 30]];   % exclude beta-bump
        case {'all'}
%             if gammafit
%                 target_freq = [1 194];
%             else
%                 target_freq = [1 180];
%             end
            target_freq = [1 230];
        otherwise
            error('%s is unknown',targetBAND);
    end
    switch upper(RECsite)
        case {'NYU'}    % NYU subjects: 60Hz line noise
%             excl_freq = [excl_freq; [53 69; 116 125; 176 185]];
            excl_freq = [excl_freq; [50 69; 110 129; 170 189]];
        case {'UMCU'}   % Utrecht subjects: 50Hz line noise
%             excl_freq = [excl_freq; [43 59; 96 105; 146 155]];
            excl_freq = [excl_freq; [40 59; 90 109; 140 159]];
    end
end
end
