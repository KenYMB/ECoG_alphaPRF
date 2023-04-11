function [modeldata] = ecog_prf_constructTimeSeries(params, opts)

% Description: 
%
% [modeldata] = ecog_prf_constructTimeSeries(powspctrm, [opts])
%
% Input
% - powspctrm           = Nx1 cell-array of power-spectrum structure with the following fields:
%   - subject
%   - resamp_parms      = {channels} events x parameters x iteration
%   - events
%   - channels
%   - average
%   - n_avg
% - opts
%   - issave
%   - outputDir
%   - stimulus
%     - apertures
%     - res             = [x,y]: pixels to desampling
%     - size            = (degree): visual angle of stimulus
%   - targetBAND        = 'bb','bbS','a','aC','FaCLb',...
%                         broadband: 'bb','bbS'
%                           'bb'   = broadband estimated with linear regression
%                                    (gammafit = true in ecog_prf_fitalpha)
%                           'bbS'  = broadband estimated as mean power
%                                    (gammafit = false in ecog_prf_fitalpha)
%                         alpha: 'a','aC','aCR','aCL',... See also, ALPHACOMPUTATIONTYPES
%                           'a'    = spectral power change at alpha frequnecy
%                           'aC'   = estimated alpha power change with gaussian fitting
%                                    (estimateIAF = false, allownegfit = false in ecog_prf_fitalpha)
%                           'aCR'  = inverse value of 'aC'
%                           'aCL'  = 'aC' in log scale
%   - smoothingMode     = 'none', 'smooth' or 'decimate'(default)
%   - smoothingN        = scalar (default = 3)
%   ---------------------
%   - fileid
%
% Output
% - modeldata          	= Nx1 cell-array of time-series structure with the following fields:
%   - subject
%   - channels
%   - datats            = IxM cell-array of time-series data (I=iterations, M=runs)
%   - stimulus          = 1xM cell-array of stimulus data

% Hidden options
% - opts
%   - compute       = [];
%   - skipexist     = [];
%   - fileid        = 'freq_spectra-timeseries';

% Hidden usage
% - opts.compute    = false;        % load or bypass data
%   - opts.average  = 'runs';
%   - opts.allowlag = false;
% 
% [modeldata] = ecog_prf_constructTimeSeries(modeldata, [opts])
% [modeldata] = ecog_prf_constructTimeSeries(subjectList, [opts])

% Dependency: <ECoG_utils>, <analyzePRF>, <analyzePRFdog>,
%             SetDefault, saveauto, decimate_col

% 20200304 - Yuasa: separate ecog_prf_analyzePRF into ecog_prf_constructTimeSeries and ecog_prf_analyzePRF
% 20200921 - Yuasa: update for data with iterations
% 20210812 - Yuasa: fix uncorrelated alpha compuation
% 20210908 - Yuasa: update uncorrelated alpha compuation
% 20210914 - Yuasa: enable multiple targetBAND
% 20211117 - Yuasa: enable to estimate low_broadband

%% Set options
%--Define inputs 
% <opts>
SetDefaultAnalysisPath('DATA','pRFmodel','opts.outputDir');
SetDefault('opts.issave',false);
SetDefault('opts.stimulus.apertures','bar_apertures.mat');
SetDefault('opts.stimulus.res',[100 100]);
SetDefault('opts.stimulus.size',16.6);
SetDefault('opts.targetBAND',{},'cell');
SetDefault('opts.smoothingMode','decimate');    % 'none','smooth','decimate'(default)
SetDefault('opts.smoothingN',3);                % integar: width of smoothing window (default = 3)
% <Check Validity>
assert(~isempty(params), 'Please provide spectra parameters');
assert(~isempty(opts.targetBAND), 'opts.targetBAND is required');
if all(strcmpi(opts.targetBAND,'ERP'))||all(strcmpi(opts.targetBAND,'ERPlag'))
    opts.targetBAND = opts.targetBAND{1};
    warning('targetBAND is specified as ''%s''\nCalling %s...',...
                opts.targetBAND,'ecog_prf_constructTimeSeriesERP');
    modeldata = ecog_prf_constructTimeSeriesERP(params, opts);
    return;
end
% <hidden opts>
SetDefault('opts.compute',[]);
SetDefault('opts.skipexist',[]);
SetDefault('opts.fileid','freq_spectra-timeseries');
SetDefault('opts.average','runs');                      % just for saving and loading files
SetDefault('opts.allowlag',false);                      % just for saving and loading files
SetDefault('opts.maxlag',nan);                          % just for saving and loading files
        
%-- check inputs and outputs
if isempty(opts.compute)
    SetDefault('opts.skipexist',true);
    opts.compute = true;
elseif ~opts.compute
    opts.skipexist = false;
else
    SetDefault('opts.skipexist',false);
end
if opts.issave && ~exist(opts.outputDir, 'dir'),     mkdir(opts.outputDir); end
    
%% load stimulus set of pRF experiment
if opts.compute
    %-- load stimulus apertures
    res     = opts.stimulus.res;
    if ischar(opts.stimulus.apertures)
        apertures = load(opts.stimulus.apertures);
        tmp = fieldnames(apertures);
        apertures = apertures.(tmp{1});
    else
        apertures = opts.stimulus.apertures;
    end
    apertures = imresize(apertures, res, 'nearest');
end

%% Loop across subjects
cellinput = iscell(params);
if ~cellinput,  params = {params};  end

modeldata = cell(size(params));
for ii = 1 : numel(params)
    imodeldata = [];
    %-- Set data
    if ~opts.compute && ~isstruct(params{ii}) && ischar(params{ii})
        subject      = params{ii};
        isloadfile   = true;
        compute      = false;
        %-- Takeover parameters
        allowlag     = opts.allowlag;
        maxlag       = opts.maxlag;
        average      = opts.average;
    else
        iparams      = params{ii};
        subject      = iparams.subject;
        isSave       = opts.issave;
        isloadfile   = opts.skipexist;
        compute      = opts.compute;
        %-- Takeover parameters
        SetDefault('iparams.allowlag',opts.allowlag);
        SetDefault('iparams.maxlag',opts.maxlag);
        allowlag     = iparams.allowlag;
        maxlag       = iparams.maxlag;
        average      = iparams.average;
    end
    targetBAND      = getBANDname(opts.targetBAND,allowlag,ii,numel(params));
    smoothingMode   = opts.smoothingMode;
    smoothingN      = opts.smoothingN;
    if smoothingN == 1,  smoothingMode = 'none';    end
    %-- Try to load files
    if isloadfile
        postfix = cnstpostfix(average,targetBAND,smoothingMode,smoothingN);
        filename     = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        %-- load files from directory
        if exist(filename,'file') || ~compute
            fprintf('[%s] Loading time-series data for subject %s from %s <%s> ',mfilename, subject, opts.outputDir,postfix(2:end));
            iparams      = load(filename);
            isSave       = false;
            fprintf('\n');
            compute = false;
        end
    end
    n_avg        = iparams.n_avg;
    %-- Main
    if compute
        fprintf('[%s] Constructing time-series data for subject %s \n',mfilename, subject);
            
        %-- Prepare arguments
        resamp_parms = iparams.resamp_parms;
        channels     = iparams.channels;
        events       = iparams.events;
            %% Extract broadband response to each PRF stimulus
            %-- Compute average power spectra change in time window
            %--  {chs}(stims x params x boots) -> (chs x stims x boots x params)

            assert(numel(resamp_parms) == height(channels), 'Data structure looks wrong');
            bandparms      =permute(cat(4,resamp_parms{:}),[4,1,3,2]);

            %-- Actual alpha amplitude is amplitude./width bacause of the formula of fitting function
            if size(bandparms,4)<14     % when not have raw alpha
                m_base  = bandparms(:,:,:,8);      % baseline for each trial (bb_amp_low)
                bandparms(:,:,:,14) = bandparms(:,:,:,9)./bandparms(:,:,:,11) + m_base; % modeled alpha hight
            end
            switch targetBAND
                case broadbandComputationTypes('L')
            %-- low broadband
            PRFband = 100*(10.^bandparms(:,:,:,8)-1);
                case broadbandComputationTypes
            %-- broadband
            m_base  = bandparms(:,1,1,6);      % same value for all stimulus & all bootstrap trials
            PRFband = 100*(10.^bsxfun(@minus,bandparms(:,:,:,2),m_base)-1);
                case alphaComputationTypes('')
            %-- alpha
            PRFband = -100*(10.^bandparms(:,:,:,14)-1);
                case alphaComputationTypes('R')
            %-- reversed alpha
            PRFband = 100*(10.^-bandparms(:,:,:,14)-1);
                case alphaComputationTypes('L')
            %-- alpha in log
            PRFband = -bandparms(:,:,:,14);
                case alphaComputationTypes('C')
            %-- corrected alpha
            PRFband = -100*(10.^(bandparms(:,:,:,9)./bandparms(:,:,:,11))-1);
                case alphaComputationTypes('CR')
            %-- reversed corrected alpha
            PRFband = 100*(10.^-(bandparms(:,:,:,9)./bandparms(:,:,:,11))-1);
                case alphaComputationTypes('CL')
            %-- corrected alpha in log
            PRFband = -bandparms(:,:,:,9)./bandparms(:,:,:,11);
                otherwise
            error('Specified targetBAND is unknown');
            end

            %% seperate data by sessions and runs (chs x stims x boots) -> (chs x stims/run x runs x boots)
            n_sets       = numel(unique(events.session_name)) * numel(unique(events.run_name));
            n_stimsprun  = height(events)./n_sets;
            assert(isint(n_stimsprun),'[%s] %s do not consist with the same time series',mfilename,'runs');
            PRFband = reshape(PRFband,height(channels),n_stimsprun,n_sets,[]);
            events       = events(1:n_stimsprun,:);

            %% select valid trials of stimuli (stimulus = {1 x runs})
            stimulus = repmat({apertures(:,:,events.stim_file_index)},1,n_sets);

            %% apply smoothing
            stimtypelist = grp2idx(cellfun(@(x) strtok(x,'-'),events.trial_name,'UniformOutput',false));
            % boundaries = [0 28 40 56 84 96 112 140 152 168 196 208 224];
            boundaries = [0 reshape(find(diff(stimtypelist)~=0),1,[]) height(events)];
                switch smoothingMode
                    case {'smooth'}
                      for jj=1:length(boundaries)-1
                        PRFband(:,(boundaries(jj)+1):boundaries(jj+1),:,:) ...
                            = smoothdata(PRFband(:,(boundaries(jj)+1):boundaries(jj+1),:,:),2,'movmean',smoothingN);
                      end
                    case {'decimate'}
                      %--- decimate data
                      dimts = 2;              
                      PRFband = decimate_col(PRFband,smoothingN,[],dimts,'omitnan');

                      %--- decimate stimulus
                      stimulusPP = cell(size(stimulus));
                      for p=1:numel(stimulus)
                          stimulusPP{p} = squish(stimulus{p},2);
                      end
                      stimulusPP = decimate_col(double(catcell(3,stimulusPP)),smoothingN,[],2);
            %           stimulusPP(stimulusPP<0) = 0;
            %           stimulusPP = logical(round(stimulusPP));

                      stimulus = reshape(...
                                      mat2cell(reshape(stimulusPP,[res size(stimulusPP,2:ndims(stimulusPP))]),...
                                          res(1),res(2),...
                                          numel(stimulusPP)./prod(res)./numel(stimulus),...
                                          ones(1,numel(stimulus))),...
                                      size(stimulus));
                    otherwise
                        smoothingMode   = 'none';
                        smoothingN      = 1;
                end

            %% reshape datats (PRFband_mean = chs x stims x runs x boots -> datats = {boots x runs}(chs x stims))
            n_boot    = size(PRFband,4);
            datats   = reshape(...
                            mat2cell(PRFband,size(PRFband,1),size(PRFband,2),ones(n_sets,1),ones(n_boot,1)),...
                            n_sets,n_boot)';
            fprintf('[%s] Complete constructing time-series data for subject %s <avg-%s_%s-%s%d>\n',mfilename,subject,average,targetBAND,smoothingMode,smoothingN);

    else
        %-- Load parameters
        subject      = iparams.subject;
        datats       = iparams.datats;
        stimulus     = iparams.stimulus;
        channels     = iparams.channels;
        events       = iparams.events;
        targetBAND   = iparams.targetBAND;
            SetDefault('iparams.smoothingMode',opts.smoothingMode);
            SetDefault('iparams.smoothingN',opts.smoothingN);
        smoothingMode   = iparams.smoothingMode;
        smoothingN      = iparams.smoothingN;
    end
    
    %-- Collect into an output struct
    imodeldata.subject       = subject;
    imodeldata.datats        = datats;
    imodeldata.stimulus      = stimulus;
    imodeldata.events        = events;
    imodeldata.channels      = channels;
    imodeldata.average       = average;
    imodeldata.n_avg         = n_avg;
    imodeldata.allowlag      = allowlag;
    imodeldata.maxlag        = maxlag;
    imodeldata.targetBAND    = targetBAND;
    imodeldata.smoothingMode = smoothingMode;
    imodeldata.smoothingN    = smoothingN;
    modeldata{ii} = imodeldata;
    
    %-- Save out the powspctrm
    postfix = cnstpostfix(average,targetBAND,smoothingMode,smoothingN);
    if isSave
        fprintf('[%s] Saving time-series data for subject %s to %s <%s> \n',mfilename, subject, opts.outputDir,postfix(2:end));
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        saveauto(filename,'-struct','imodeldata');
    end
    
end
if ~cellinput,  modeldata = modeldata{1};  end
fprintf('[%s] Done! \n',mfilename);
end

%% Sub function
function postfix = cnstpostfix(average,targetBAND,smoothingMode,smoothingN)
postfix = '';
postfix      = sprintf('%s_avg-%s_%s',postfix,average,targetBAND);
issmooth     = ~ismember(smoothingMode,{'none'});
if issmooth,       postfix = sprintf('%s-%s%d',postfix,smoothingMode,smoothingN);  end
end

function targetBAND = getBANDname(targetBANDs,allowlag,idx,nidx)
if length(targetBANDs)==1
    targetBAND = targetBANDs{1};
elseif length(targetBANDs)==nidx
    targetBAND = targetBANDs{idx};
else
    error('%s must be a charactor-array or cell-array of the same length as the 1st argument','opts.targetBAND');
end
if allowlag && ~endsWith(targetBAND,'lag')
    targetBAND = [targetBAND 'lag'];
end
end