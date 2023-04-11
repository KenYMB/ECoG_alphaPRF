function ecog_prf_plotPRFs(varargin)

% Description: Plots pRF timecourse and pRF map across subjects
%
% ecog_prf_plotPRFs(modeldata1, prf1, modeldata2, prf2, ..., [opts])
%
% Input
% - modeldata          	= Nx1 cell-array of time-series structure with the following fields:
%   - subject
%   - channels
%   - datats            = Mx1 cell-array of time-series data
%   - stimulus          = Mx1 cell-array of stimulus data
% - prf                 = Nx1 cell-array of pRF information structure with the following fields:
%   - ecc
%   - ang
%   - expt
%   - rfsize
%   - R2
%   - params
%   - subject
%   - channels
% - opts                = structure with the following fields:
%   - visualfield       = size in visual angle (default = 16.6)
%   - doplots
%    - doplots_prfmap
%    - doplots_prfts
%   - skipprojection
%   - plotbenson        = 'yes' or 'no'(default)  % plot results estimated from benson atlas
%   - closefig          = 'yes' or 'no'(default)  % automatically close figures
%   - plotsavedir
%   - plot
%   ---------------------
%   - 
%

% Dependency: <ECoG_utils>, <analyzePRF>, <analyzePRFdog>,
%             SetDefault, saveauto, decimate_col

% 20200304 - Yuasa
% 20220222 - Yuasa: load HDgrid information from subjectlist.tsv

%% Set options
%--Define inputs
narginchk(2,inf);
if isodd(nargin),   opts = varargin{end}; end
ndatas = floor(nargin./2);
modeldata = {}; results = {};
for jj=1:ndatas
    assert(length(varargin{jj*2-1})==length(varargin{1}),'inputs are inconsistent');
    assert(length(varargin{jj*2})==length(varargin{1}),  'inputs are inconsistent');
    for ii=1:length(varargin{jj*2})
        if iscell(varargin{jj*2-1}(ii))
            modeldata{ii,jj} = varargin{jj*2-1}{ii};
        else
            modeldata{ii,jj} = varargin{jj*2-1}(ii);
        end
        if iscell(varargin{jj*2}(ii))
            results{ii,jj} = varargin{jj*2}{ii};
        else
            results{ii,jj} = varargin{jj*2}(ii);
        end
    end
end

% <opts>
SetDefaultAnalysisPath('FIGURE','pRF','opts.plotsavedir');
SetDefault('opts.visualfield',16.6);
SetDefault('opts.doplots',true);
SetDefault('opts.doplots_prfmap',opts.doplots);
SetDefault('opts.doplots_prfts',opts.doplots);
SetDefault('opts.skipprojection',true);
SetDefault('opts.plotbenson','no');
SetDefault('opts.closefig','no');
SetDefault('opts.plot.addR2ToTitle','yes');
SetDefault('opts.plot.RotGrid',true,1);
SetDefault('opts.plot.fontSize',14);

%-- check inputs and outputs
if (opts.doplots_prfmap || opts.doplots_prfts) && ~exist(opts.plotsavedir, 'dir'), mkdir(opts.plotsavedir); end

%-- Collect Subject Information
subjectList_fname = 'subjectlist.tsv';
SbjInfo    = loadSbjInfo(subjectList_fname,'all');
hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');

%% Loop across subjects
clsfig = strcmp(opts.closefig,'yes');

for ii = 1 : size(modeldata,1)
    
    imodeldata   = modeldata(ii,:);
    iresults     = results(ii,:);
    
        subject         = imodeldata{1}.subject;
        channels        = imodeldata{1}.channels;
        
        average         = imodeldata{1}.average;
%         n_avg           = imodeldata{1}.n_avg;
        if all(isfield(imodeldata{1},{'smoothingMode','smoothingN'}))
            smoothingMode   = imodeldata{1}.smoothingMode;
            smoothingN      = imodeldata{1}.smoothingN;
        else
            smoothingMode   = 'none';
            smoothingN      = 1;
        end
        
        if isfield(iresults{1}.options,'compacc')
            compacc         = iresults{1}.options.compacc;
        else
            compacc         = iresults{1}.compacc;
        end
        if isfield(iresults{1}.options,'prfmodel')
            prfmodel         = iresults{1}.options.prfmodel;
        else
            prfmodel         = iresults{1}.prfmodel;
        end
        if isfield(iresults{1}.options,'gaussianmode')
            gaussianmode         = iresults{1}.options.gaussianmode;
        else
            gaussianmode         = iresults{1}.gaussianmode;
        end
        
        stim_tmp = imodeldata{1}.stimulus;
        while iscell(stim_tmp)
            stim_tmp = stim_tmp{1};
        end
        res = size(stim_tmp,[1,2]);
        
    %-- plot Benson reference
    plotbenson   = false;
    if strcmpi(opts.plotbenson,'yes')
        if all(ismember({'bensoneccen','bensonangle','bensonsigma'},channels.Properties.VariableNames))
            plotbenson   = true;
        else
            warning('Benson atlas information is not found for %s',subject);
        end
    end
            
    %-- HD grid flag
    if hasSbjInfo && istablefield(SbjInfo,'hasHDgrid')
        hasHDgrid = any(strcmpi(SbjInfo.hasHDgrid(ismember(SbjInfo.participant_id,subject)),'yes'));
    else
        hasHDgrid = false;
    end
    
    %-- Plot figures
    if hasHDgrid,  	whichElectrodes = find(~contains(channels.name,'GB'));
    else,         	whichElectrodes = 1:height(channels);
    end

    %%-- Plot model time-series
    if opts.doplots_prfts
    for jj=1:ndatas
        targetBAND      = imodeldata{jj}.targetBAND;
        switch targetBAND
            case alphaComputationTypes
                label = 'alpha';
            case broadbandComputationTypes
                label = 'broadband';
            otherwise
                label = targetBAND;
        end
    
        postfix = '';        % prepare for filenames
        postfix       = sprintf('%s-%s-%s_avg-%s_%s',postfix,prfmodel,gaussianmode,average,label);
        issmooth     = ~ismember(smoothingMode,{'none'});
        if issmooth,       postfix = sprintf('%s-%s%d',postfix,smoothingMode,smoothingN);  end
        isdouble     = isa(compacc(1),'double');
        if ~isdouble, postfix = sprintf('%s_%s',postfix,class(compacc(1)));  end
        
        figureName = sprintf('pRFmodelts-prf_%s%s', subject,postfix);

        opt2plot = [];
        opt2plot.plot        = opts.plot;
        opt2plot.plot.pix2deg   = opts.visualfield./max(res);
        if plotbenson
            opt2plot.plotbenson     = 'yes';
            opt2plot.plot.legend    = {label,'estimation','benson'};
        else
            opt2plot.plotbenson     = 'no';
            opt2plot.plot.legend    = label;
        end
        opt2plot.skipprojection = opts.skipprojection;
        ecog_plotGridPRFts(imodeldata{jj}.datats, imodeldata{jj}.stimulus, iresults{jj}, whichElectrodes,opt2plot);
        hgexport(gcf, fullfile(opts.plotsavedir, figureName), hgexport('factorystyle'), 'Format', 'png'); if clsfig, close; end
          if hasHDgrid
          [p]=ecog_plotGridPRFts(imodeldata{jj}.datats, imodeldata{jj}.stimulus, iresults{jj}, 'GB',opt2plot);
          for pp=1:length(p)
            if length(p)==1,  hgexport(p{pp}, fullfile(opts.plotsavedir, sprintf('%s-GB',figureName)), hgexport('factorystyle'), 'Format', 'png');
            else,             hgexport(p{pp}, fullfile(opts.plotsavedir, sprintf('%s-GB-%02d',figureName,pp)), hgexport('factorystyle'), 'Format', 'png');
            end
          end
          if clsfig, close(p{:}); end
          end
    end
    end

    %%-- Plot pRF
    if opts.doplots_prfmap
    label = {};
    for jj=1:ndatas
        targetBAND      = imodeldata{jj}.targetBAND;
        switch targetBAND
            case alphaComputationTypes
                label{jj} = 'alpha';
            case broadbandComputationTypes
                label{jj} = 'broadband';
            otherwise
                label{jj} = targetBAND;
        end
    end
    
    postfix = '';        % prepare for filenames
    postfix       = sprintf('%s-%s-%s_avg-%s_%s',postfix,prfmodel,gaussianmode,average,strjoin(label,'-'));
    issmooth     = ~ismember(smoothingMode,{'none'});
    if issmooth,       postfix = sprintf('%s-%s%d',postfix,smoothingMode,smoothingN);  end
    isdouble     = isa(compacc(1),'double');
    if ~isdouble, postfix = sprintf('%s_%s',postfix,class(compacc(1)));  end
        
    figureName = sprintf('pRFposition-prf_%s%s', subject,postfix);
    
    opt2plot = [];
    opt2plot.plot        = opts.plot;
    opt2plot.plot.pix2deg   = opts.visualfield./max(res);
    opt2plot.plot.XLim      = [-1 1].*opts.visualfield;
    opt2plot.plot.YLim      = [-1 1].*opts.visualfield;
    if plotbenson
        opt2plot.plotbenson     = 'yes';
        opt2plot.plot.legend    = [label,{'benson'}];
    else
        opt2plot.plotbenson     = 'no';
        opt2plot.plot.legend    = label;
    end
    ecog_plotGridPRF(whichElectrodes, opt2plot, iresults{:});
    hgexport(gcf, fullfile(opts.plotsavedir, figureName), hgexport('factorystyle'), 'Format', 'png'); if clsfig, close; end
      if hasHDgrid
      [p]=ecog_plotGridPRF('GB', opt2plot, iresults{:});
      for pp=1:length(p)
        if length(p)==1,  hgexport(p{pp}, fullfile(opts.plotsavedir, sprintf('%s-GB',figureName)), hgexport('factorystyle'), 'Format', 'png');
        else,             hgexport(p{pp}, fullfile(opts.plotsavedir, sprintf('%s-GB-%02d',figureName,pp)), hgexport('factorystyle'), 'Format', 'png');
        end
      end
      if clsfig, close(p{:}); end
      end
    end
    
end
fprintf('[%s] Done! \n',mfilename);
end
