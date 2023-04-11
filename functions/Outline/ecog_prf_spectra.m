function [freq] = ecog_prf_spectra(data, opts)

% Description: 
%
% [freq] = ecog_prf_spectra(data, [opts])
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
%   - issave
%   - outputDir
%   - doplots
%    - doplots_alpha
%    - doplots_broadband
%    - doplots_whole
%    - doplots_wholelog
%   - plotsavedir
%   - plot          % see ecog_plotGridSpectra for this option
%   ---------------------
%   - fileid
%
% Output
% - freq            = Nx1 cell-array of data structure with the following fields:
%   - subject
%   - spectra       = channels x events x f
%   - spectra_off   = channels x 1 x f
%   - f
%   - events
%   - channels

% Hidden options
% - opts
%   - compute        = [];
%   - skipexist      = [];
%   - fileid         = 'freq_spectra';
%   - taskNames      = {'prf'};
%   - freq_alpha     = [3 25];
%   - freq_broadband = [20 220];

% Hidden usage
% - opts.compute = false;        % load or bypass data
%   - opts.allowlag = false;
% 
% [freq] = ecog_prf_spectra(freq, [opts])
% [freq] = ecog_prf_spectra(subjectList, [opts])

% Dependency: <ECoG_utils>, ecog_regressout, SetDefault, cellstrfind, saveauto

% 20200219 - Yuasa
% 20220222 - Yuasa: load HDgrid information from subjectlist.tsv
% 20220810 - Yuasa: update for dataset without BLANK
% 20230320 - Yuasa: add hidden options

%% Set options
%--Define inputs 
% <opts>
SetDefaultAnalysisPath('DATA','Spectrum','opts.outputDir');
SetDefaultAnalysisPath('FIGURE','spectrum','opts.plotsavedir');
SetDefault('opts.target_time',[0 0.5]);
SetDefault('opts.issave',false);
SetDefault('opts.doplots',false);
SetDefault('opts.doplots_alpha',opts.doplots);
SetDefault('opts.doplots_broadband',opts.doplots);
SetDefault('opts.doplots_whole',opts.doplots);
SetDefault('opts.doplots_wholelog',opts.doplots);
SetDefault('opts.freq_alpha',[3 25]);                   % frequency range for alpha plotting
SetDefault('opts.freq_broadband',[20 220]);             % frequency range for braodband plotting
SetDefault('opts.plot.XLim',[]);                        % frequency range for whole plotting
SetDefault('opts.plot.XScale','linear');
SetDefault('opts.plot.fontSize',16);
% <hidden opts>
SetDefault('opts.f_wintyp',@hann);      % window for cpsd; default is hann window
SetDefault('opts.f_winlng',0.2);        % window length in time; default is 0.2s (=Fs/5)
SetDefault('opts.f_winovl',0.5);        % overlap ratio of window; default is 0.5 (=50%) of window length
SetDefault('opts.taskNames','*','cell');
SetDefault('opts.compute',[]);
SetDefault('opts.skipexist',[]);
SetDefault('opts.fileid','freq_spectra');
SetDefault('opts.allowlag',false);                      % just for loading files
SetDefault('opts.maxlag',nan);                          % just for saving and loading files

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
if (opts.doplots||opts.doplots_alpha) && ~exist(opts.plotsavedir, 'dir'), mkdir(opts.plotsavedir); end

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
    %-- Try to load files
    if isloadfile
        postfix = cnstpostfix(opts.fileid,allowlag);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        %-- load files from directory
        if exist(filename,'file') || ~compute
            fprintf('[%s] Loading spectra for subject %s from %s ',mfilename, subject, opts.outputDir);
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
    
    %-- Main
    if compute
        fprintf('[%s] Computing spectra for subject %s \n',mfilename, subject);

        %-- Prepare arguments
        epochs      = permute(idat.epochs,[3,2,1]); % channels x events x t
        t           = idat.t;
        srate       = idat.fsample;
        events      = idat.events;
        channels    = idat.channels;

        %-- Select the prf events (baseline is BLANK in prf sessions)
        trialIndex  = cellstrfind(events.task_name, opts.taskNames,1);
        events      = events(trialIndex,:);
        bslIndex    = contains(events.trial_name, 'BLANK');

        %-- Calculate spectra
        fft_t   = t>opts.target_time(1) & t<=opts.target_time(2); % time segment for spectrum
        fft_w   = window(opts.f_wintyp,round(srate.*opts.f_winlng)); % window width   % window(@hann,round(srate/5))
        fft_ov  = round(srate.*opts.f_winlng.*opts.f_winovl); % overlap               % round(srate/10)
        reg_erp = 0; % 1 to regress erp out, 0 not to, if yes, make sure to baseline correct first
        %%-- stimuli epochs
        [f,spectra] = ...
            ecog_spectra(epochs(:,trialIndex,:),[],fft_w,fft_t,fft_ov,srate,reg_erp);

        %-- Fill nan for flat signals
        [nanel,nantrl] = find(var(spectra,0,3)==0);
        for jj = 1:length(nanel)
            spectra(nanel(jj),nantrl(jj),:) = nan;
        end

        %-- Baseline
        if ~any(bslIndex)    % compute baseline based on pre-stimulus period
            %-- Calculate spectra
            fft_t   = t>diff(opts.target_time([2,1])) & t<=0; % time segment for spectrum
            [~,spectra_off] = ...
                ecog_spectra(epochs(:,trialIndex,:),[],fft_w,fft_t,fft_ov,srate,reg_erp);
            %-- Fill nan for flat signals
            [nanel,nantrl] = find(var(spectra_off,0,3)==0);
            for jj = 1:length(nanel)
                spectra_off(nanel(jj),nantrl(jj),:) = nan;
            end
        else                % pick up BLANK trials & apply geomean across trials
            spectra_off     = geomean(spectra(:,bslIndex,:),2,'omitnan');
        end
    else
        %-- Load parameters
        spectra     = idat.spectra;
        spectra_off = idat.spectra_off;
        f           = idat.f;
        channels    = idat.channels;
        events      = idat.events;
        bslIndex    = contains(events.trial_name, 'BLANK');
    end
    
    %-- Collect into an output struct
    ifreq.subject     = subject;
    ifreq.spectra     = spectra;
    ifreq.spectra_off = spectra_off;
    ifreq.f           = f;
    ifreq.events      = events;
    ifreq.channels    = channels;
    ifreq.allowlag    = allowlag;
    ifreq.maxlag      = maxlag;
    freq{ii} = ifreq;
    
    %-- Save out the data
    postfix = cnstpostfix(opts.fileid,allowlag);
    if isSave
        fprintf('[%s] Saving spectra for subject %s to %s \n',mfilename, subject, opts.outputDir);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        saveauto(filename,'-struct','ifreq');
    end
    
    %%-- set 0 for prf BLANK, 1 for prf stimuli, and 10+trial_type for other stimuli
    stmlist = double(ismember(events.task_name,'prf')&~ismember(events.trial_name,'BLANK')) + ...
              ~ismember(events.task_name,'prf').*(events.trial_type+10); 
    
    %-- Plot figures          
    if opts.doplots_alpha
        foi  = f>=opts.freq_alpha(1) & f<=opts.freq_alpha(2);
        
        trials = [];
        trials.pwrspctrm    = spectra(:,:,foi);
        trials.f            = f(foi);
        trials.events       = events;
        trials.channels     = channels;
        
        %%-- set trial name 'PRF' & set specs
        trials.events.trial_name(stmlist==1) = {'PRF'};
        eventList   = table(trials.events.trial_name,stmlist,'VariableNames',{'trial_name','stmlist'});
        eventList   = sortrows(eventList,'stmlist');
        eventList   = unique(eventList.trial_name,'stable');
        
        specs = [];
        specs.plot      = opts.plot;
        specs.plot.XLim = [min(f(foi)) max(f(foi))];
        
        if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
        else,         	whichElectrodes = trials.channels.name;
        end
        
        figureName = sprintf('spectra-prf%s_%s_alpha', postfix,subject);
        
        ecog_plotGridSpectra(trials, whichElectrodes, eventList,[], specs);
        hgexport(gcf, fullfile(opts.plotsavedir, figureName), hgexport('factorystyle'), 'Format', 'png'); close;
          if hasHDgrid
          [~,p]=ecog_plotGridSpectra(trials, 'GB', eventList, [], specs);
          for pp=1:length(p)
            if length(p)==1,  hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB',figureName)), hgexport('factorystyle'), 'Format', 'png');
            else,             hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB-%02d',figureName,pp)), hgexport('factorystyle'), 'Format', 'png');
            end
          end
          close(p);
          end
    end
    if opts.doplots_broadband
        foi  = f>=opts.freq_broadband(1) & f<=opts.freq_broadband(2);
        
        trials = [];
        trials.pwrspctrm    = spectra(:,:,foi);
        trials.f            = f(foi);
        trials.events       = events;
        trials.channels     = channels;
        
        %%-- set trial name 'PRF' & set specs
        trials.events.trial_name(stmlist==1) = {'PRF'};
        eventList   = table(trials.events.trial_name,stmlist,'VariableNames',{'trial_name','stmlist'});
        eventList   = sortrows(eventList,'stmlist');
        eventList   = unique(eventList.trial_name,'stable');
        
        specs = [];
        specs.plot          = opts.plot;
        specs.plot.XScale   = 'log';
        specs.plot.XLim     = [min(f(foi)) max(f(foi))];
        
        if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
        else,         	whichElectrodes = trials.channels.name;
        end
        
        figureName = sprintf('spectra-prf%s_%s_broadband', postfix,subject);
        
        ecog_plotGridSpectra(trials, whichElectrodes, eventList,[], specs);
        hgexport(gcf, fullfile(opts.plotsavedir, figureName), hgexport('factorystyle'), 'Format', 'png'); close;
          if hasHDgrid
          [~,p]=ecog_plotGridSpectra(trials, 'GB', eventList, [], specs);
          for pp=1:length(p)
            if length(p)==1,  hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB',figureName)), hgexport('factorystyle'), 'Format', 'png');
            else,             hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB-%02d',figureName,pp)), hgexport('factorystyle'), 'Format', 'png');
            end
          end
          close(p);
          end
    end
    if opts.doplots_whole
        if isempty(opts.plot.XLim)
            foi  = true(size(f));
        else
            foi  = f>=opts.plot.XLim(1) & f<=opts.plot.XLim(2);
        end
        
        trials = [];
        trials.pwrspctrm    = spectra(:,:,foi);
        trials.f            = f(foi);
        trials.events       = events;
        trials.channels     = channels;
        
        %%-- set trial name 'PRF' & set specs
        trials.events.trial_name(stmlist==1) = {'PRF'};
        eventList   = table(trials.events.trial_name,stmlist,'VariableNames',{'trial_name','stmlist'});
        eventList   = sortrows(eventList,'stmlist');
        eventList   = unique(eventList.trial_name,'stable');
        
        specs = [];
        specs.plot      = opts.plot;
        
        if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
        else,         	whichElectrodes = trials.channels.name;
        end
        
        figureName = sprintf('spectra-prf%s_%s_whole', postfix,subject);
        
        ecog_plotGridSpectra(trials, whichElectrodes, eventList,[], specs);
        hgexport(gcf, fullfile(opts.plotsavedir, figureName), hgexport('factorystyle'), 'Format', 'png'); close;
          if hasHDgrid
          [~,p]=ecog_plotGridSpectra(trials, 'GB', eventList, [], specs);
          for pp=1:length(p)
            if length(p)==1,  hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB',figureName)), hgexport('factorystyle'), 'Format', 'png');
            else,             hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB-%02d',figureName,pp)), hgexport('factorystyle'), 'Format', 'png');
            end
          end
          close(p);
          end
    end
    if opts.doplots_wholelog
        if isempty(opts.plot.XLim)
            foi = true(size(f));
        else
            foi     = f>=opts.plot.XLim(1) & f<=opts.plot.XLim(2);
        end
        
        trials = [];
        trials.pwrspctrm    = spectra(:,:,foi);
        trials.f            = f(foi);
        trials.events       = events;
        trials.channels     = channels;
        
        %%-- set trial name 'PRF' & set specs
        trials.events.trial_name(stmlist==1) = {'PRF'};
        eventList   = table(trials.events.trial_name,stmlist,'VariableNames',{'trial_name','stmlist'});
        eventList   = sortrows(eventList,'stmlist');
        eventList   = unique(eventList.trial_name,'stable');
        
        specs = [];
        specs.plot      = opts.plot;
        specs.plot.XScale   = 'log';
        
        if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
        else,         	whichElectrodes = trials.channels.name;
        end
        
        figureName = sprintf('spectra-prf_%s_wholelog', subject);
        
        ecog_plotGridSpectra(trials, whichElectrodes, eventList,[], specs);
        hgexport(gcf, fullfile(opts.plotsavedir, figureName), hgexport('factorystyle'), 'Format', 'png'); close;
          if hasHDgrid
          [~,p]=ecog_plotGridSpectra(trials, 'GB', eventList, [], specs);
          for pp=1:length(p)
            if length(p)==1,  hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB',figureName)), hgexport('factorystyle'), 'Format', 'png');
            else,             hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB-%02d',figureName,pp)), hgexport('factorystyle'), 'Format', 'png');
            end
          end
          close(p);
          end
    end
    
end
if ~cellinput,  freq = freq{1};  end
fprintf('[%s] Done! \n',mfilename);
end

%% Sub function
function postfix = cnstpostfix(fileid,allowlag)
postfix = '';
%%-- lag
addlag = allowlag && ~contains(fileid,'regresslag');
if addlag,   postfix = sprintf('%s_regresslag',postfix);	end
end