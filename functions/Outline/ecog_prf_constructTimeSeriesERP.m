function [modeldata] = ecog_prf_constructTimeSeriesERP(regressor, opts)

% Description: 
%
% [modeldata] = ecog_prf_constructTimeSeriesERP(regressor, [opts])
%
% Input
% - regressor           = Nx1 cell-array of regressor structure with the following fields:
%   - subject
%   - coef          = channels x events
%   - predictor     = t x 1 x channels
%   - t
%   - events
%   - channels
%   - fsample
% - opts
%   - issave
%   - outputDir
%   - stimulus
%     - apertures
%     - res             = [x,y]: pixels to desampling
%     - size            = (degree): visual angle of stimulus
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
%   - fileid        = 'data_erp-timeseries';
%   - targetBAND    = 'ERP';

% Hidden usage
% - opts.compute    = false;        % load or bypass data
%   - opts.average  = 'runs';
%   - opts.allowlag = false;
% 
% [modeldata] = ecog_prf_constructTimeSeriesERP(modeldata, [opts])
% [modeldata] = ecog_prf_constructTimeSeriesERP(subjectList, [opts])

% Dependency: <ECoG_utils>, <analyzePRF>, <analyzePRFdog>,
%             SetDefault, saveauto, decimate_col

% 20210414 - Yuasa
% 20220803 - Yuasa: delete unused options


%% Set options
%--Define inputs 
% <opts>
SetDefaultAnalysisPath('DATA','pRFmodel','opts.outputDir');
SetDefault('opts.issave',false);
SetDefault('opts.stimulus.apertures','bar_apertures.mat');
SetDefault('opts.stimulus.res',[100 100]);
SetDefault('opts.stimulus.size',16.6);
SetDefault('opts.smoothingMode','decimate');    % 'none','smooth','decimate'(default)
SetDefault('opts.smoothingN',3);                % integar: width of smoothing window (default = 3)
SetDefault('opts.targetBAND','ERP');
% <Check Validity>
assert(~isempty(regressor), 'Please provide regressor data');
if ~startsWith(opts.targetBAND,'ERP','IgnoreCase',true)
    warning('targetBAND is not ''EPR''\nCalling %s...','ecog_prf_constructTimeSeries');
    modeldata = ecog_prf_constructTimeSeries(regressor, opts);
    return;
end
% <hidden opts>
SetDefault('opts.compute',[]);
SetDefault('opts.skipexist',[]);
SetDefault('opts.fileid','data_erp-timeseries');
SetDefault('opts.average','runs');                      % just for saving and loading files
SetDefault('opts.allowlag',false);                      % just for saving and loading files
SetDefault('opts.maxlag',nan);                          % just for saving and loading files
        
%-- check inputs and outputs
if endsWith(opts.targetBAND,'lag')
    opts.allowlag = true;
end
if ~opts.allowlag
    opts.maxlag = 0;
end
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
cellinput = iscell(regressor);
if ~cellinput,  regressor = {regressor};  end

modeldata = cell(size(regressor));

for ii = 1 : numel(regressor)
    imodeldata = [];
    %-- Set data
    if ~opts.compute && ~isstruct(regressor{ii}) && ischar(regressor{ii})
        subject      = regressor{ii};
        isloadfile   = true;
        compute      = false;
        %-- Takeover parameters
        average      = opts.average;
        allowlag     = opts.allowlag;
        maxlag       = opts.maxlag;
    else
        idat         = regressor{ii};
        subject      = idat.subject;
        isSave       = opts.issave;
        isloadfile   = opts.skipexist;
        compute      = opts.compute;
        %-- Takeover parameters
        average      = idat.average;
        SetDefault('idat.allowlag',opts.allowlag);
        SetDefault('idat.maxlag',opts.maxlag);
        allowlag     = idat.allowlag;
        maxlag       = idat.maxlag;
    end
    targetBAND = fixBANDname(opts.targetBAND,allowlag);
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
            idat      = load(filename);
            isSave    = false;
            fprintf('\n');
            compute = false;
        end
    end
    n_avg        = idat.n_avg;
    %-- Main
    if compute
        fprintf('[%s] Constructing time-series data for subject %s \n',mfilename, subject);
            
        %-- Prepare arguments
        coef         = idat.coef;       % channels x events x  boots
        channels     = idat.channels;
        events       = idat.events;
            %% Extract ERP coefficients to each PRF stimulus
            %--  (chs x stims x boots)
            PRFband  = coef;

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
        subject      = idat.subject;
        datats       = idat.datats;
        stimulus     = idat.stimulus;
        channels     = idat.channels;
        events       = idat.events;
        targetBAND   = idat.targetBAND;
            SetDefault('idat.smoothingMode',opts.smoothingMode);
            SetDefault('idat.smoothingN',opts.smoothingN);
        smoothingMode   = idat.smoothingMode;
        smoothingN      = idat.smoothingN;
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

function targetBAND = fixBANDname(targetBAND,allowlag)
if allowlag && ~endsWith(targetBAND,'lag')
    targetBAND = [targetBAND 'lag'];
end
end