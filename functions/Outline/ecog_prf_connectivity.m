function [conn] = ecog_prf_connectivity(xfreq, freq, opts)

% Description: 
%
% [conn] = ecog_prf_connectivity(xfreq,[opts])
% [conn] = ecog_prf_connectivity(xfreq,freq,[opts])
%
% Input
% - xfreq           = Nx1 cell-array of data structure with the following fields:
%   - subject
%   - crossspectra  = channels x channels x events x f
%  (- spectra       = channels x events x f, refer when freq is not specified)
%   - f
%   - events
%   - channels
% - freq            = Nx1 cell-array of data structure with the following fields:
%   - subject
%   - spectra       = channels x events x f
%   - f
%   - events
%   - channels
% - opts
%   - method        = strings, can be
%       - 'coh'     = coherence (default)
%       - 'mscoh'   = magnitude-squared coherence
%       - 'imcoh'   = imaginary part of coherence
%       - 'plv'     = phase locking value
%       - 'pli'     = phase lag index
%       - 'wpli'    = weighted phase lag index
%   - average       = 'none','sessions','runs','stimuli'(default) or 'trials'
%   - isavgfreq     = true (default) or false, if ture, take avarage across frequency
%   - issave
%   - outputDir
%   ---------------------
%   - fileid
%
% Output
% - conn            = Nx1 cell-array of data structure with the following fields:
%   - subject
%   - connectivity  = channels x channels x events x f
%                     apply 'abs' to get magnitude
%   - f
%   - events
%   - channels
%   - average
%   - n_avg

% Hidden usage
% - opts.compute = false;        % load or bypass data
%   - opts.allowlag = false;
% 
% [conn] = ecog_prf_connectivity(conn, [opts])
% [conn] = ecog_prf_connectivity(subjectList, [opts])

% Dependency: <ECoG_utils>, ecog_regressout, SetDefault, cellstrfind, saveauto

% 20201116 - Yuasa
% 20201201 - Yuasa: bug fix - avoid to exclude channels which include NaN
%                             in any table components
% 20210513 - Yuasa: minor update

%% Set inputs
if nargin==2 && isstruct(freq) && ~(isfield(freq,'spectra') && isfield(freq,'f'))
    opts    = freq;
    freq    = {};
end
SetDefault('freq',{},'cell');
hasfreq   = ~isempty(freq);

%% Set options
%--Define inputs 
% <opts>
SetDefaultAnalysisPath('DATA','xSpectrum','opts.outputDir');
SetDefault('opts.average','stimuli');
SetDefault('opts.isavgfreq',false);
SetDefault('opts.issave',false);
% <method>
SetDefault('opts.method','coh');
% <hidden opts>
SetDefault('opts.targetBAND','');
SetDefault('opts.target_freq',[]);
SetDefault('opts.gammafit',false);      % just use to set default value in target_freq
SetDefault('opts.stimNames',{},'cell');
SetDefault('opts.compute',[]);
SetDefault('opts.skipexist',[]);
SetDefault('opts.fileid','freq_connectivity');
SetDefault('opts.allowlag',false);                      % just for loading files

% <StimNames>
if isempty(opts.stimNames) || isempty(opts.stimNames{1}) || strcmpi(opts.stimNames{1},'all')
    opts.stimNames = {'*'};
end

%-- check inputs and outputs
assert(~isempty(xfreq), 'Please provide the data struct');
if isempty(opts.compute)
    SetDefault('opts.skipexist',true);
    opts.compute = true;
elseif ~opts.compute
    opts.skipexist = false;
else
    SetDefault('opts.skipexist',false);
end
if opts.issave && ~exist(opts.outputDir, 'dir'),     mkdir(opts.outputDir); end

%% Loop across subjects
cellinput = iscell(xfreq);
if ~cellinput,     xfreq = {xfreq};  end
if hasfreq
    if ~iscell(freq),  freq  = {freq};   end
    assert(length(xfreq)==length(freq),'xfreq and freq must have include the same subjects');
end

conn = cell(size(xfreq));
for ii = 1 : numel(xfreq)
    iconn = [];
    %-- Set data
    if ~opts.compute && ~isstruct(xfreq{ii}) && ischar(xfreq{ii})
        subject     = xfreq{ii};
        isloadfile  = true;
        compute     = false;
        %-- Takeover parameters
        allowlag        = opts.allowlag;
        targetBAND      = opts.targetBAND;
        target_freq     = opts.target_freq;
        stimNames       = opts.stimNames;
    else
        ifrq        = xfreq{ii};
        subject     = ifrq.subject;
        isSave      = opts.issave;
        isloadfile  = opts.skipexist;
        compute     = opts.compute;
        %-- Takeover parameters
        SetDefault('ifrq.allowlag',opts.allowlag);
        if isempty(opts.targetBAND)
            SetDefault('ifrq.targetBAND','other');
        else
            SetDefault('ifrq.targetBAND',opts.targetBAND);
        end
        SetDefault('ifrq.target_freq',opts.target_freq);
        SetDefault('ifrq.stimNames',opts.stimNames);
        allowlag        = ifrq.allowlag;
        targetBAND      = ifrq.targetBAND;
        target_freq     = ifrq.target_freq;
        stimNames       = ifrq.stimNames;
    end
    [targetBAND,target_freq] = updatefreq(targetBAND,target_freq,opts);
    method          = opts.method;
    average         = opts.average;
    isavgfreq       = opts.isavgfreq;
    %-- Try to load files
    if isloadfile
        postfix = cnstpostfix(opts.fileid,targetBAND,target_freq,stimNames,allowlag,method,average,isavgfreq);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        %-- load files from directory
        if exist(filename,'file') || ~compute
            fprintf('[%s] Loading connectivity for subject %s from %s ',mfilename, subject, opts.outputDir);
            ifrq        = load(filename);
            isSave      = false;
            fprintf('\n');
            compute = false;
        end
    end
    %-- Main
    if compute
        fprintf('[%s] Computing connectivity for subject %s \n',mfilename, subject);
        
        %-- check input type
        hasspectra   = isfield(ifrq,'spectra');

        %-- Prepare arguments
        crossspectra    = ifrq.crossspectra;
        f               = ifrq.f;
        channels        = ifrq.channels;
        events          = ifrq.events;
        
        isfullx = ndims(crossspectra) > 3;
        if ~isfullx
          chancmb       = ifrq.chancmb;
          chanidx       = (chancmb-[1,0])*[height(channels);1];
          chanidxT      = (chancmb-[0,1])*[1;height(channels)];
          diagidx       = chancmb(:,1)==chancmb(:,2);
          hasdiag       = all(ismember(1:height(channels),chancmb(diagidx,1)));
        else
          hasdiag       = false;
        end
        
        %-- get power spectra
        if hasfreq               % if freq is specified
            iifreq = freq{ii};
            [events,trlidx1,trlidx2]    = intersect(events,iifreq.events,'stable');
            chktbl = intersect(intersect({'name','type','sampling_frequency','reference','group','subject_name'},...
                                    channels.Properties.VariableNames),iifreq.channels.Properties.VariableNames);
            [~,chnidx1,chnidx2]  = intersect(channels(:,chktbl),iifreq.channels(:,chktbl),'stable');
            channels             = channels(chnidx1,:);
            [f,fidx1,fidx2]  = intersect(f,iifreq.f,'stable');
            if isfullx
                crossspectra = crossspectra(chnidx1,chnidx1,trlidx1,fidx1);
            else
                chnidx3 = all(ismember(chancmb,chnidx1),2);
                chancmb      = chancmb(chnidx3,:);
                crossspectra = crossspectra(chnidx3,trlidx1,fidx1);
            end
            spectra    	 = iifreq.spectra(chnidx2,trlidx2,fidx2);
        elseif hasspectra        % if xfreq has field of spectra
            spectra    	 = ifrq.spectra;
        elseif isfullx
            spectra      = zeros(size(crossspectra,2:ndims(crossspectra)));
            for ich = 1:height(channels)
                spectra(ich,:,:) = crossspectra(ich,ich,:,:);
            end
        elseif hasdiag
            spectra      = crossspectra(diagidx,:,:);
        else
            error('Power spectra is required');
        end
        haspow   = any(~isnan(spectra(:)));
        
        %-- Check power spectra
        switch method
            case {'coh','abscoh','imcoh'}
                assert(hasfreq|hasspectra|haspow,'''%s'' requires power spectra',opts.method);
        end 
        
        %-- Target frequency
        if ~isempty(target_freq)
            fidx = f >= target_freq(1) & f <= target_freq(2);
            f(~fidx) = [];
            spectra(:,:,~fidx) = [];
            if isfullx,     crossspectra(:,:,:,~fidx) = [];
            else,           crossspectra(:,:,~fidx) = [];
            end
        end
        
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
        
        %-- Connectivity analysis
        avgidx = 2;      nf = length(f);
%         if isavgfreq,       avgidx = [2 3];  nf = 1;
%         end
        connectivity = nan(height(channels),height(channels),length(n_avg),nf);
        powspectra   = nan(height(channels),length(n_avg),nf);
        events_idx   = zeros(length(n_avg),1);
        switch method
            case 'coh'
                for iavg = 1:length(n_avg)
                    if isfullx
                      tmpspectra  = mean(crossspectra(:,:,avg_group==iavg,:),avgidx+1,'omitnan');
                    else
                      tmpspectra  = zeros(height(channels).^2,1,nf);
                      tmpspectra(chanidx,:,:)   = mean(crossspectra(:,avg_group==iavg,:),avgidx,'omitnan');
                      tmpspectra(chanidxT,:,:)  = mean(crossspectra(:,avg_group==iavg,:),avgidx,'omitnan');
                      tmpspectra = reshape(tmpspectra,height(channels),height(channels),1,nf);
                    end
                    connectivity(:,:,iavg,:)  = tmpspectra;
                    powspectra(:,iavg,:)      = sqrt(mean(spectra(:,avg_group==iavg,:),avgidx,'omitnan'));
                    events_idx(iavg)          = find(avg_group==iavg,1);
                end
                for ich = 1:height(channels)
                    connectivity(ich,:,:,:)  = permute(connectivity(ich,:,:,:),[2,3,4,1]) ./ powspectra;
                    connectivity(:,ich,:,:)  = permute(connectivity(:,ich,:,:),[1,3,4,2]) ./ powspectra;
                end
            case 'mscoh'
                for iavg = 1:length(n_avg)
                    if isfullx
                      tmpspectra  = crossspectra(:,:,avg_group==iavg,:).*conj(crossspectra(:,:,avg_group==iavg,:));
                    else
                      tmpspectra  = zeros(height(channels).^2,sum(avg_group==iavg),nf);
                      tmpspectra(chanidx,:,:)   = crossspectra(:,avg_group==iavg,:).*conj(crossspectra(:,avg_group==iavg,:));
                      tmpspectra(chanidxT,:,:)  = crossspectra(:,avg_group==iavg,:).*conj(crossspectra(:,avg_group==iavg,:));
                      tmpspectra = reshape(tmpspectra,height(channels),height(channels),sum(avg_group==iavg),nf);
                    end
                    powspectra                = spectra(:,avg_group==iavg,:);
                    for ich = 1:height(channels)
                        tmpspectra(ich,:,:,:)  = permute(tmpspectra(ich,:,:,:),[2,3,4,1]) ./ powspectra;
                        tmpspectra(:,ich,:,:)  = permute(tmpspectra(:,ich,:,:),[1,3,4,2]) ./ powspectra;
                    end
                    connectivity(:,:,iavg,:)  = mean(tmpspectra,avgidx+1,'omitnan');
                    events_idx(iavg)          = find(avg_group==iavg,1);
                end
            case 'imcoh'
                for iavg = 1:length(n_avg)
                    if isfullx
                      tmpspectra  = mean(imag(crossspectra(:,:,avg_group==iavg,:)),avgidx+1,'omitnan');
                    else
                      tmpspectra  = zeros(height(channels).^2,1,nf);
                      tmpspectra(chanidx,:,:)   = mean(imag(crossspectra(:,avg_group==iavg,:)),avgidx,'omitnan');
                      tmpspectra(chanidxT,:,:)  = mean(imag(crossspectra(:,avg_group==iavg,:)),avgidx,'omitnan');
                      tmpspectra = reshape(tmpspectra,height(channels),height(channels),1,nf);
                    end
                    connectivity(:,:,iavg,:)  = tmpspectra;
                    powspectra(:,iavg,:)      = sqrt(mean(spectra(:,avg_group==iavg,:),avgidx),'omitnan');
                    events_idx(iavg)          = find(avg_group==iavg,1);
                end
                for ich = 1:height(channels)
                    connectivity(ich,:,:,:)  = permute(connectivity(ich,:,:,:),[2,3,4,1]) ./ powspectra;
                    connectivity(:,ich,:,:)  = permute(connectivity(:,ich,:,:),[1,3,4,2]) ./ powspectra;
                end
            case 'plv'
                for iavg = 1:length(n_avg)
                    if isfullx
                      tmpspectra  = mean(crossspectra(:,:,avg_group==iavg,:)...
                                    ./abs(crossspectra(:,:,avg_group==iavg,:)),avgidx+1,'omitnan');
                    else
                      tmpspectra  = zeros(height(channels).^2,1,nf);
                      tmpspectra(chanidx,:,:)   = mean(crossspectra(:,avg_group==iavg,:)...
                                                  ./abs(crossspectra(:,avg_group==iavg,:)),avgidx,'omitnan');
                      tmpspectra(chanidxT,:,:)  = mean(crossspectra(:,avg_group==iavg,:)...
                                                  ./abs(crossspectra(:,avg_group==iavg,:)),avgidx,'omitnan');
                      tmpspectra = reshape(tmpspectra,height(channels),height(channels),1,nf);
                    end
                    connectivity(:,:,iavg,:)  = tmpspectra;
                    events_idx(iavg)          = find(avg_group==iavg,1);
                end
            case 'pli'
                for iavg = 1:length(n_avg)
                    if isfullx
                      tmpspectra  = mean(sign(imag(crossspectra(:,:,avg_group==iavg,:))),avgidx+1,'omitnan');
                    else
                      tmpspectra  = zeros(height(channels).^2,1,nf);
                      tmpspectra(chanidx,:,:)   = mean(sign(imag(crossspectra(:,avg_group==iavg,:))),avgidx,'omitnan');
                      tmpspectra(chanidxT,:,:)  = mean(sign(imag(crossspectra(:,avg_group==iavg,:))),avgidx,'omitnan');
                      tmpspectra = reshape(tmpspectra,height(channels),height(channels),1,nf);
                    end
                    connectivity(:,:,iavg,:)  = tmpspectra;
                    events_idx(iavg)          = find(avg_group==iavg,1);
                end
            case 'wpli'
                for iavg = 1:length(n_avg)
                    if isfullx
                      tmpspectra  = mean(imag(crossspectra(:,:,avg_group==iavg,:)),avgidx+1,'omitnan')...
                                    ./mean(abs(imag(crossspectra(:,:,avg_group==iavg,:))),avgidx+1,'omitnan');
                    else
                      tmpspectra  = zeros(height(channels).^2,1,nf);
                      tmpspectra(chanidx,:,:)   = mean(imag(crossspectra(:,avg_group==iavg,:)),avgidx,'omitnan')...
                                                  ./mean(abs(imag(crossspectra(:,avg_group==iavg,:))),avgidx,'omitnan');
                      tmpspectra(chanidxT,:,:)  = mean(imag(crossspectra(:,avg_group==iavg,:)),avgidx,'omitnan')...
                                                  ./mean(abs(imag(crossspectra(:,avg_group==iavg,:))),avgidx,'omitnan');
                      tmpspectra = reshape(tmpspectra,height(channels),height(channels),1,nf);
                    end
                    connectivity(:,:,iavg,:)  = tmpspectra;
                    events_idx(iavg)          = find(avg_group==iavg,1);
                end
            otherwise,          error('''%s'' is unknown connectivity method',method);     
        end
        if isavgfreq
            connectivity = mean(connectivity,4,'omitnan');
            f = round(mean(f),3,'significant');
        end
        events = events(events_idx,:);
    else
        %-- Load parameters
        connectivity    = ifrq.connectivity;
        f               = ifrq.f;
        channels        = ifrq.channels;
%         events          = ifrq.events(events_idx,:);
        events          = ifrq.events;
        method          = ifrq.method;
        average         = ifrq.average;
        n_avg           = ifrq.n_avg;
        if isfield(ifrq,'isavgfreq')
          isavgfreq     = ifrq.isavgfreq;
        else
          isavgfreq     = length(f) == 1;
        end
    end
    
    %-- Collect into an output struct
    iconn.subject        = subject;
    iconn.connectivity   = connectivity;
    iconn.f              = f;
    iconn.events         = events;
    iconn.channels       = channels;
    iconn.method         = method;
    iconn.average        = average;
    iconn.n_avg          = n_avg;
    iconn.isavgfreq      = isavgfreq;
    iconn.targetBAND     = targetBAND;
    iconn.target_freq    = target_freq;
    iconn.stimNames      = stimNames;
    iconn.allowlag       = allowlag;
    conn{ii}  = iconn;
    
    %-- Save out the data
    postfix = cnstpostfix(opts.fileid,targetBAND,target_freq,stimNames,allowlag,method,average,isavgfreq);
    if isSave
        fprintf('[%s] Saving connectivity for subject %s to %s \n',mfilename, subject, opts.outputDir);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        saveauto(filename,'-struct','iconn');
    end
    
end
if ~cellinput,  conn = conn{1};  end
fprintf('[%s] Done! \n',mfilename);
end

%% Sub function
function postfix = cnstpostfix(fileid,targetBAND,target_freq,stimNames,allowlag,method,average,isavgfreq)
postfix = '';
%%-- lag
addlag = allowlag && ~contains(fileid,'regresslag');
if addlag,   postfix = sprintf('%s_regresslag',postfix);	end
%-- average
postfix         = sprintf('%s_%s_avg-%s',postfix,method,average);
%%-- freq
if ~isempty(targetBAND) && ~strcmpi(targetBAND,'other')
    postfix         = sprintf('%s-%s',postfix,targetBAND);
elseif length(target_freq)==2
    postfix         = sprintf('%s-%g-%gHz',postfix,target_freq(1),target_freq(2));
end
%%-- freq average
if isavgfreq
    postfix         = sprintf('%s-averaged',postfix);
end
%-- stim
if ~ismember('*',stimNames)
    if length(stimNames)==1
        postfix     = [postfix '_' strrep(stimNames{1},'*','')];
    elseif length(stimNames)>1
        postfix     = [postfix '_' strrep(stimNames{1},'*','') '-' strrep(stimNames{end},'*','')];
    end
end
end

function [targetBAND,target_freq] = updatefreq(targetBAND,target_freq,opts)
if ~isempty(opts.target_freq) && length(opts.target_freq)<=2
    if length(opts.target_freq) == 1
        target_freq     = [opts.target_freq opts.target_freq];
    else
        target_freq     = opts.target_freq;
    end
elseif ~isempty(opts.targetBAND)
    targetBAND      = opts.targetBAND;
    switch targetBAND
        case {'alpha'}
            target_freq = [6 17];
        case {'broadband'}
            if opts.gammafit
                target_freq = [30 194];
            else
                target_freq = [70 180];
            end
    end
end
end