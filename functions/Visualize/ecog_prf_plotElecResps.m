function ecog_prf_plotElecResps(varargin)

% Description: 
%
% ecog_prf_plotElecResps(modeldata1,[modeldata2],[opts])
%
% Input
% - modeldata          	= Nx1 cell-array of time-series structure with the following fields:
%   - subject
%   - channels
%   - datats            = Mx1 cell-array of time-series data
%   - stimulus          = Mx1 cell-array of stimulus data
% - opts
%   - alpha             = criterion for confidence interval (default = 0.05)
%   - plotsavedir
%   - plot
%   ---------------------
%   - 
%
% 
% If modeldata is NxM cell-array, M datasets were averaged for plot

% Dependency: ECoG_utils, SetDefault

% 20200321 - Yuasa
% 20200407 - Yuasa: add errorbar for bootstrap data
% 20220222 - Yuasa: load HDgrid information from subjectlist.tsv

%% Set options
%--Define inputs
narginchk(1,3);
modeldata1 = varargin{1};
if nargin>2
    modeldata2 = varargin{2};
    opts = varargin{3};
elseif iscell(varargin{2})
    modeldata2 = varargin{2};
    assert(length(modeldata1)==length(modeldata2),'Modeldatas must have same length');
else
    opts = varargin{2};
end

% <opts>
SetDefault('opts.alpha',0.05);
SetDefault('opts.plotsavedir',fullfile(analysisRootPath, 'Figures', 'wholeVF'));
SetDefault('opts.plot.fontSize',14);

SetDefault('opts.plot.RotGrid', true);

%-- check inputs and outputs
if  ~exist(opts.plotsavedir, 'dir'), mkdir(opts.plotsavedir); end

%-- check plot type
dotwoside = exist('modeldata2','var');
flglgnd   = isfield(opts.plot,'legend');

%-- Collect Subject Information
subjectList_fname = 'subjectlist.tsv';
SbjInfo    = loadSbjInfo(subjectList_fname,'all');
hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');

%% Main
ndatspsbj = size(modeldata1,2);
for ii = 1 : size(modeldata1,1)

    %--
    imodeldata1 = modeldata1{ii,1};
    subject  = imodeldata1.subject;
    channels = imodeldata1.channels;
    datats1     = zeros(height(channels),ndatspsbj);
    for jj = 1:ndatspsbj,   datats1(:,jj) = [modeldata1{ii,jj}.datats{:}];   end
    isalpha1 = ismember(imodeldata1.targetBAND,alphaComputationTypes);
    datats1  = datats1.*(-1).^isalpha1;
    if dotwoside
      imodeldata2 = modeldata2{ii,1};
      datats2     = zeros(height(channels),ndatspsbj);
      for jj = 1:ndatspsbj,   datats2(:,jj) = [modeldata2{ii,jj}.datats{:}];   end
      isalpha2 = ismember(imodeldata2.targetBAND,alphaComputationTypes);
      datats2  = datats2.*(-1).^isalpha2;
    end
    

    
    
    %-- Not show area label
    if any(ismember(channels.Properties.VariableNames,'bensonarea'))
        channels.bensonarea = [];
    end
    if any(ismember(channels.Properties.VariableNames,'wangarea'))
        channels.wangarea = [];
    end
    
    %-- HD grid flag
    if hasSbjInfo && istablefield(SbjInfo,'hasHDgrid')
        hasHDgrid = any(strcmpi(SbjInfo.hasHDgrid(ismember(SbjInfo.participant_id,subject)),'yes'));
    else
        hasHDgrid = false;
    end
    
    if hasHDgrid,  	selElectrodes = ~contains(channels.name,'GB');
    else,         	selElectrodes = true(height(channels),1);
    end
    
    %-- Channel #
    nchs = numel(find(selElectrodes));
    
    %-- Plot
    p=figure('Position',[150 600 500+nchs*10 420]);
    hl=[];
    if dotwoside,   yyaxis left;    end
    resps1   = nanmean(datats1,2);
    lerr1    = resps1-prctile(datats1,opts.alpha/2*100,2);    uerr1    = prctile(datats1,(1-opts.alpha/2)*100,2)-resps1;
%     hl(1)=plot(resps(selElectrodes),'o:','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','auto');     % line plot 
    if dotwoside
      hl(1)=bar((1:nchs)-0.2,resps1(selElectrodes),0.4);       % bar plot
      if ndatspsbj>1
          hold on;
          eb = errorbar((1:nchs)-0.2,resps1(selElectrodes),lerr1(selElectrodes),uerr1(selElectrodes));
          eb.Color = [0 0 0];
          eb.LineStyle = 'none';
          hold off;
      end
    else
      hl(1)=bar(resps1(selElectrodes),0.8);       % bar plot
    end
    set(gca,'XTick',1:nchs,'XTickLabel',channels.name(selElectrodes));
    set(gca,'FontSize',opts.plot.fontSize,'XTickLabelRotation',90,'XLim',[0,nchs]+0.5);
    yyl = ylim;
    if isalpha1
        axis ij;
        yyl = fliplr(-yyl);
    end
    if dotwoside
      yyaxis right
      resps2   = nanmean(datats2,2);
      lerr2    = resps2-prctile(datats2,2.5,2);    uerr2    = prctile(datats2,97.5,2)-resps2;
%       hl(2)=plot(resps(selElectrodes),'o:','LineWidth',2,'MarkerSize',8,'MarkerFaceColor','auto');     % line plot 
      hl(2)=bar((1:nchs)+0.2,resps2(selElectrodes),0.4);       % bar plot
      if ndatspsbj>1
          hold on;
          eb = errorbar((1:nchs)+0.2,resps2(selElectrodes),lerr2(selElectrodes),uerr2(selElectrodes));
          eb.Color = [0 0 0];
          eb.LineStyle = 'none';
          hold off;
      end
      straightline(0,'h','k--');
      yyrstep = mean(diff(yticks));
      yyr = [nanmin(resps2(selElectrodes)), nanmax(resps2(selElectrodes))];
      if isalpha2
          axis ij;
          yyr = fliplr(-yyr);
      end
      yyr = [floor(yyr(1)./yyrstep), ceil(yyr(2)./yyrstep)].*yyrstep;
      lrratio = yyr./yyl;
      if lrratio(1)>lrratio(2)
          yyr(2) = yyl(2)./yyl(1).*yyr(1);
      else
          yyr(1) = yyl(1)./yyl(2).*yyr(2);
      end
      if isalpha2
          yyr = fliplr(-yyr);
      end
      ylim(yyr);
      if flglgnd
      hlg=legend(hl,opts.plot.legend,'Location','northwest','NumColumns',2,'FontSize',opts.plot.fontSize);
      end
    elseif flglgnd
      hlg=legend(hl,opts.plot.legend,'Location','northwest','FontSize',opts.plot.fontSize);
    end
    hA=gca;
    setsubplotaxes(hA);
    hA.Position = hA.Position + [0.02,0.01,-0.04,hA.TightInset(4)*1/5];
    if flglgnd,  hlg.Position(2) = 0.5 - hlg.Position(4)/2 + hA.Position(2)/2 + hA.Position(4)/2;   end
    
    figureName = sprintf('ElecResponses_%s', subject);
    hgexport(p, fullfile(opts.plotsavedir, figureName), hgexport('factorystyle'), 'Format', 'png'); close(p);
    

    %% plot Grid
      if hasHDgrid
          whichHDgrid = 'GB';
          
%           isalpha1 = ismember(imodeldata1.targetBAND,alphaComputationTypes);
%           resps1   = [imodeldata1.datats{:}].*(-1).^isalpha1;
          
          if dotwoside
%             isalpha2 = ismember(imodeldata2.targetBAND,alphaComputationTypes);
%             resps2   = [imodeldata2.datats{:}].*(-1).^isalpha2;
            
            opts2 = [];
            opts2.channels    = channels;
            opts2.plot        = opts.plot;
            opts2.plot.axis   = {'xy','xy'};
            if isalpha1,    opts2.plot.axis{1} = 'ij';  end
            if isalpha2,    opts2.plot.axis{2} = 'ij';  end
            [p] = ecog_plotGridBarTwoSide(resps1, resps2, whichHDgrid, opts2);
          else
            opts2 = [];
            opts2.channels    = channels;
            opts2.plot        = opts.plot;
            opts2.plot.colors  = 1;
            if flglgnd,     opts2.timelabel = opts.plot.legend; end
            [p] = ecog_plotGridBarwitherr(resps1, [], [], whichHDgrid, opts2);
            
            if isalpha1
              for jj=1:length(p)
                hpanel = get(p{jj},'Children');
                hpanel(~ismember(get(hpanel,'Type'),{'axes'})) = [];
                axis(hpanel,'ij');
              end
            end
          end
          
          for pp=1:length(p)
            if length(p)==1,  hgexport(p{pp}, fullfile(opts.plotsavedir, sprintf('%s-GB',figureName)), hgexport('factorystyle'), 'Format', 'png');
            else,             hgexport(p{pp}, fullfile(opts.plotsavedir, sprintf('%s-GB-%02d',figureName,pp)), hgexport('factorystyle'), 'Format', 'png');
            end
          end
          close(p{:});
      end
    
end
