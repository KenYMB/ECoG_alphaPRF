function [prf] = ecog_prf_analyzePRF(modeldata, opts)

% Description: 
%
% [prf] = ecog_prf_analyzePRF(modeldata, [opts])
%
% Input
% - modeldata          	= Nx1 cell-array of time-series structure with the following fields:
%   - subject
%   - channels
%   - datats            = IxM cell-array of time-series data (I=iterations, M=runs)
%   - stimulus          = 1xM cell-array of stimulus data
% - opts
%   - issave
%   - outputDir
%   - hrf
%   - prfmodel      = 'css','fixexpt' or 'linear' (default)
%   - gaussianmode  = 'dog','og' or 'gs' (default)
%   ---------------------
%   - fileid
%
% Output
% - prf                = Nx1 cell-array of pRF information structure with the following fields:
%   - ecc
%   - ang
%   - expt
%   - rfsize
%   - R2
%   - params
%   - subject
%   - channels

% Hidden options
% - opts
%   - compute       = [];
%   - skipexist     = [];
%   - fileid        = 'prf';

% Hidden usage
% - opts.compute         = false;   % load or bypass data
%   - opts.average       = 'runs';
%   - opts.allowlag      = false;
%   - opts.targetBAND    = 'FaCLb';
%   - opts.smoothingMode = 'decimate';
%   - opts.smoothingN    = 3;
% 
% [prf] = ecog_prf_analyzePRF(prf, [opts])
% [prf] = ecog_prf_analyzePRF(subjectList, [opts])

% Dependency: <ECoG_utils>, <analyzePRF>, <analyzePRFdog>,
%             SetDefault, saveauto, decimate_col

% 20200224 - Yuasa
% 20200304 - Yuasa: separate ecog_prf_analyzePRF into ecog_prf_constructTimeSeries and ecog_prf_analyzePRF
% 20200921 - Yuasa: update for data with iterations
% 20210914 - Yuasa: enable multiple targetBAND

%% Set options
%--Define inputs 
% <opts>
SetDefaultAnalysisPath('DATA','pRFmodel','opts.outputDir');
SetDefault('opts.issave',false);
% <options for analyzePRF>
SetDefault('opts.hrf',1);
SetDefault('opts.maxpolydeg',0);
SetDefault('opts.prfmodel','linear');       %'css','fixexpt','linear'
SetDefault('opts.gaussianmode','gs');       %'dog','og','gs'
SetDefault('opts.isdouble',true);
SetDefault('opts.forcebounds',2);
SetDefault('opts.typicalgain',[]);
SetDefault('opts.typicalexpt',[]);
% <hidden opts>
SetDefault('opts.compute',[]);
SetDefault('opts.skipexist',[]);
SetDefault('opts.fileid','prf');
SetDefault('opts.average','runs');                      % just for loading files
SetDefault('opts.allowlag',false);                      % just for loading files
SetDefault('opts.targetBAND','FaCLb','cell');           % just for loading files
SetDefault('opts.smoothingMode','decimate');            % just for loading files
SetDefault('opts.smoothingN',3);                        % just for loading files
        
%-- check inputs and outputs
assert(~isempty(modeldata), 'Please provide time-series data');
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
cellinput = iscell(modeldata);
if ~cellinput,  modeldata = {modeldata};  end

prf = cell(size(modeldata));

for ii = 1 : length(modeldata)
    %-- Set data
    if ~opts.compute && ~isstruct(modeldata{ii}) && ischar(modeldata{ii})
        subject      = modeldata{ii};
        isloadfile   = true;
        compute      = false;
        %-- Takeover parameters
        average         = opts.average;
        targetBAND      = getBANDname(opts.targetBAND,opts.allowlag,ii,length(modeldata));
        smoothingMode   = opts.smoothingMode;
        smoothingN      = opts.smoothingN;
    else
        imodeldata   = modeldata{ii};
        subject      = imodeldata.subject;
        isSave       = opts.issave;
        isloadfile   = opts.skipexist;
        compute      = opts.compute;
        %-- Takeover parameters
        average         = imodeldata.average;
        targetBAND      = imodeldata.targetBAND;
        SetDefault('imodeldata.smoothingMode',opts.smoothingMode);
        SetDefault('imodeldata.smoothingN',opts.smoothingN);
        smoothingMode   = imodeldata.smoothingMode;
        smoothingN      = imodeldata.smoothingN;
    end
    if smoothingN == 1,  smoothingMode = 'none';    end
    if opts.isdouble,       compacc = @double;
    else,                   compacc = @single;
    end
    prfmodel      	= opts.prfmodel;
    gaussianmode   	= opts.gaussianmode;
    %-- Try to load files
    if isloadfile
        postfix = cnstpostfix(average,targetBAND,smoothingMode,smoothingN,prfmodel,gaussianmode,compacc);
        filename     = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        %-- load files from directory
        if exist(filename,'file') || ~compute
            fprintf('[%s] Loading pRF parameters for subject %s from %s <%s> ',mfilename, subject, opts.outputDir,postfix(2:end));
            imodeldata      = load(filename);
            isSave          = false;
            fprintf('\n');
            compute = false;
        end
    end
    n_avg        = imodeldata.n_avg;
    
    %-- Main
    if compute
        fprintf('[%s] Applying pRF analysis for subject %s \n',mfilename, subject);
            
        %-- Prepare arguments
        datats          = imodeldata.datats;
        stimulus        = imodeldata.stimulus;
        channels        = imodeldata.channels;

            %% typicalgain (datats =  {boots x runs}(chs x stims) -> chs x stims x iterations)
            if ~isempty(opts.typicalgain)
                datatsmean = nanmean(cat(3,datats{:}),3);     % average across all datasets
                datatsmean(datatsmean<=prctile(datatsmean,95,2)) = nan; % 95% percentile in each ch
                opts.typicalgain = nanmean(nanrms(datatsmean,2));    % rms across stimuli, mean across chs
            end

            %% Proceed to fitting....
            tr = 1;
            optsflds = fieldnames(opts);
            opts2 = rmfield(opts,optsflds(ismember(optsflds,...
                                {'stimulus','compute','skipexist','fileid','issave','outputDir','doplots','plotsavedir','plot',...
                                 'average','targetBAND','smoothingMode','smoothingN','isdouble'})));

                %-- Temporal file (back up)
                postfix = cnstpostfix(average,targetBAND,smoothingMode,smoothingN,prfmodel,gaussianmode,compacc);
                dd = 0; fileconflict = true;
                while (fileconflict)
                    dd= dd+1;
                    filename0      = fullfile(opts.outputDir, sprintf('%s_%s%s_temp%d.mat', subject,opts.fileid,postfix,dd));
                    fileconflict  = exist(filename0,'file');
                end

                %-- Run analyzePRF
                nboot = size(datats,1);
                resultsN      = cell(nboot,1);
                resultsN_xval = cell(nboot,1);
                for iboot=1:nboot

                      opts2.xvalmode = 0; % no cross-validation
                      resultsN{iboot} = analyzePRFdog(cellfun(compacc,stimulus,'UniformOutput',false),datats(iboot,:),tr,opts2);
                      fprintf('[%s] Complete pRF analysis <%s-%s@%s> N=%d \n',mfilename, prfmodel, gaussianmode,class(compacc(1)),iboot);

                      opts2.xvalmode = 2; % cross-validation within runs
                      resultsN_xval{iboot} = analyzePRFdog(cellfun(compacc,stimulus,'UniformOutput',false),datats(iboot,:),tr,opts2);
                      fprintf('[%s] Complete pRF analysis w/ cross-validation <%s-%s@%s> N=%d \n',mfilename, prfmodel, gaussianmode,class(compacc(1)),iboot);

                      resultsN{iboot}.xval = resultsN_xval{iboot}.xval;
                      
                      %-- Temporal file (back up)
                      if ~mod(iboot,100)
                          saveauto(filename0,'resultsN','resultsN_xval','iboot');
                      end
                end

                %-- Extract cell-array
                results = resultsN{1};  results_xval = resultsN_xval{1};
                prfflds = fieldnames(results);
                for ifld = prfflds'
                    if ~isempty(results.(ifld{:}))&&~isstruct(results.(ifld{:}))
                        iterdim = ndims(results.(ifld{:}))+1;
                        catfld = cellfun(@(C) getfield(C,ifld{:}),resultsN,'UniformOutput',false);
                        results.(ifld{:}) = cat(iterdim,catfld{:});
                        catfld = cellfun(@(C) getfield(C,ifld{:}),resultsN_xval,'UniformOutput',false);
                        results_xval.(ifld{:}) = cat(iterdim,catfld{:});
                    end
                end
                fprintf('[%s] Complete whole pRF analysis for subject %s <%s-%s_avg-%s_%s-%s%d@%s>\n',mfilename,subject,prfmodel,gaussianmode,average,targetBAND,smoothingMode,smoothingN,class(compacc(1)));
        
    else
        %-- Load parameters
        results         = rmfield(imodeldata,{'subject','channels','results_xval'});
        results_xval    = imodeldata.results_xval;
        subject         = imodeldata.subject;
        channels        = imodeldata.channels;
        compacc         = imodeldata.compacc; 
        prfmodel        = imodeldata.options.prfmodel;
        gaussianmode    = imodeldata.options.gaussianmode;
    end
    
    %-- Collect into an output struct
    results.results_xval    = results_xval;
    results.subject         = subject;
    results.channels        = channels;
    results.targetBAND      = targetBAND;
    results.smoothingMode   = smoothingMode;
    results.smoothingN      = smoothingN;
    results.compacc         = compacc;
    results.average         = average;
    results.n_avg           = n_avg;
    prf{ii}   = results;
    
    %-- Save out the powspctrm
    postfix = cnstpostfix(average,targetBAND,smoothingMode,smoothingN,prfmodel,gaussianmode,compacc);
    if isSave
        fprintf('[%s] Saving pRF parameters for subject %s to %s <%s> \n',mfilename, subject, opts.outputDir,postfix(2:end));
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        saveauto(filename,'-struct','results');
    end
    
    %-- Temporal file (back up)
    if exist('filename0','var') && exist(filename0,'file')
        delete(filename0);
    end
        
    %-- Plot figures (-> Outsourced to ecog_prf_plotPRFs)
    
end
if ~cellinput,  prf = prf{1};  end
fprintf('[%s] Done! \n',mfilename);
end

%% Sub function
function postfix = cnstpostfix(average,targetBAND,smoothingMode,smoothingN,prfmodel,gaussianmode,compacc)
postfix = '';
postfix       = sprintf('%s-%s-%s_avg-%s_%s',postfix,prfmodel,gaussianmode,average,targetBAND);
issmooth     = ~ismember(smoothingMode,{'none'});
if issmooth,       postfix = sprintf('%s-%s%d',postfix,smoothingMode,smoothingN);  end
isdouble     = isa(compacc(1),'double');
if ~isdouble, postfix = sprintf('%s_%s',postfix,class(compacc(1)));  end
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