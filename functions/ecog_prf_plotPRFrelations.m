function   varargout = ecog_prf_plotPRFrelations(prfs,prfs_corr,reltarget,varargin)
% h = ecog_prf_plotPRFrelations(prfs,prfs_corr,Target,[parms,conds,selrois,roiList,opts])
%   plots pRF relations between two pRF parameters or pRF conditions.
% Output figure has NxM subplots, where N is the numbers of pRF parameters
% or pRF conditions and M is the number of ROIs.
% 
% Input
%   prfs:       structure consists of pRF parameters, 
%               the field names are assumed as "ParameterName_ConditionName", 
%               each filed is cell-array consists of ROIs.
%   prfs_corr:  structure consists of correlations of pRF parameters,
%               the field names are assumed as "ParameterName1_ParameterName2_fit", 
%               or as "ParameterName_ConditionName1_ConditionName2_fit",
%               required when plotting linear regressed line.
%               (if opts.method = 'line' or opts.showfitted = 'fitlm')
%   reltarget:  'Prop1_Prop2', characters consist of two properties to show relations
% 
%   parms:      pRF parameters list to plot (cf. 'ang','ecc')
%   conds:      pRF model conditions list to plot (cf. 'bb','a')
%   selrois:    ROI names to plot
%   roiList:    full ROI names to identify which index of prfs for plotting
% 
%   opts:
%     method:         'histogram','scatter','line', or 'none'
%     condNames:      cell-array of characters, names of conds
%     pix2deg:        coefficience to update units
%     overplot:       [](default),true,false, if overplotting multiple properties in the same panel
%     showdiag:       [](default),true,false, if showing diagonal line
%     scatt_jitter:   scalar, rate of jitter size to axis in scatter plot, default = 0.01
%     ang_limpad:     scalar, padding axis for plotting 'ang' to show circular correlation
%     ang_circshift:  false(default),true, if shifting centers of axes for plotting 'ang'
% 
%     showfitted:     'ellipse','fitlm','tls', or 'none'
%     fit_sd:         scalar, standard deviation to show covariance ellipses
%                     or confidence intervals of fitting
%     fit_alpha:      scalar value, indicates as alpha instead of sd
%     fit_excludeang: true(default),false, if avoiding to show fitting for 'ang'
%     fitlm_domean:   true(default),false, take mean if true, take median if false
%     fitlm_showref:  true(default),false, if showing fMRI reference for 'fitlm'
% 
%     FontSize:       font size in the figure
%     Color:          plot color
%     colormap:       colormap
%     showcolorbar:   false(default),true
% 
%     (histogram)
%     NBin:           (default = 25)
%     (scatter)
%     MarkerSize:     (default = [])
%     (line)
%     LineWidth:      (default = 1.5)
% 
%         
% Example 1:
%   Target  = {'bb_a'};
%   parms   = {'ang','ecc','rfsize'};
%   conds   = {'bb','a'};
%   selrois = {'V1-V3','Dorsolateral'};
%   opts = [];
%   opts.mothod      = 'hist';
%   opts.showfitted  = 'ellipse';
%   opts.fit_sd      = 1;        % 1sd ellipse
%   opts.ang_limpad  = 90;       % degree for angle
%   opts.pix2deg     = 16.6/100;
%   ecog_prf_plotPRFrelations(prfs,prfs_corr,Target,parms,conds,selrois,roiList,opts)
% 
% Example 2:
%   Target  = {'ecc_rfsize'};
%   parms   = {'ang','ecc','rfsize'};
%   conds   = {'bb','a'};
%   selrois = {'V1','V2','V3','Dorsolateral'};
%   opts = [];
%   opts.mothod      = 'line';
%   opts.fit_sd      = 1;        % 1sd CI
%   opts.ang_limpad  = 90;       % degree for angle
%   opts.pix2deg     = 16.6/100;
%   ecog_prf_plotPRFrelations(prfs,prfs_corr,Target,parms,conds,selrois,roiList,opts)
% 

% Dependency: keepfields(fieldTrip), SetDefault, realignYLabels, error_ellipse, fitline2derror

% 20220607 Yuasa
% 20221026 Yuasa - add robust option
% 20221122 Yuasa - enable to plot other parameters

%-- Set inputs
narginchk(2,inf);
if ~isempty(prfs_corr)&&~isstruct(prfs_corr)
    varargin  = [{reltarget} varargin];
    reltarget = prfs_corr;
    prfs_corr = [];
end
if ~exist('reltarget','var') || isempty(reltarget)
    error('Target is required');
elseif iscell(reltarget)
    reltarget = reshape(reltarget,1,[]);
else
    reltarget = {reltarget};
end

%-- Get parameters and conditions lists
prffields = fieldnames(prfs);
[allparms,allconds] = strtok(prffields,'_');
[allconds,remain]   = strtok(allconds,'_');
extrafiled = cellfun(@isempty,allconds)|~cellfun(@isempty,remain)|...
             endsWith(prffields,{'_dist','_fit','_corr','_R2'});
allparms = unique(allparms(~extrafiled),'stable');
allconds = unique(allconds(~extrafiled),'stable');

%-- Solve varargin
parms = [];
conds = [];
selrois = [];
roiList = [];
opts = [];

for ii = 1:numel(varargin)
    if ii == numel(varargin) && isstruct(varargin{end})
        opts = varargin{end};
    elseif isSTRs(varargin{ii})
        if isempty(parms) && any(ismember(varargin{ii},allparms))
            parms = varargin{ii};
        elseif isempty(conds) && any(ismember(varargin{ii},allconds))
            conds = varargin{ii};
        elseif isempty(selrois) && ...
                  (ii+1)<=numel(varargin) && isSTRs(varargin{ii+1}) && ...
                  any(ismember(varargin{ii},varargin{ii+1}))
            selrois = varargin{ii};
            roiList = varargin{ii+1};
            varargin{ii+1} = [];
        elseif isempty(roiList)
            roiList = varargin{ii};
        end
    elseif isnumeric(varargin{ii}) || islogical(varargin{ii})
        if isempty(selrois)
            selrois = varargin{ii};
        end
    end
end

%-- Check specified parameters and conditions
if isempty(parms)
    parms = reshape(allparms,1,[]);
else
    if ~iscell(parms),  parms = {parms};    end
    chkinp = ismember(parms,allparms);
    if ~all(chkinp)
        warning('Ignore unknown parameter:%s',sprintf(' ''%s''',parms{~chkinp}));
    end
    parms = reshape(parms(chkinp),1,[]);
end
if isempty(conds)
    conds = reshape(allconds,1,[]);
else
    if ~iscell(conds),  conds = {conds};    end
    chkinp = ismember(conds,allconds);
    if ~all(chkinp)
        warning('Ignore unknown condition:%s',sprintf(' ''%s''',conds{~chkinp}));
    end
    conds = reshape(conds(chkinp),1,[]);
end

%-- Check ROIs
if isempty(roiList)
    nroi = size(prfs.(prffields{find(~extrafiled,1)}),1);
    roiList = arrayfun(@(n)sprintf('%d',n),1:nroi,'UniformOutput',false);
else
%     nroi = numel(roiList);
    roiList = reshape(roiList,1,[]);
end
if isempty(selrois)
    selrois = roiList;
elseif isnumeric(varargin{ii}) || islogical(varargin{ii})
    selrois = reshape(roiList(selrois),1,[]);
else
    selrois = reshape(selrois,1,[]);
end

%-- Options
SetDefault('opts.mothod','hist');
SetDefault('opts.condNames',[]);
SetDefault('opts.overplot',[]);
SetDefault('opts.showdiag',[]);
SetDefault('opts.scatt_jitter',0.01);
SetDefault('opts.ang_limpad',90);
SetDefault('opts.ang_circshift',false);
SetDefault('opts.showfitted','none');
SetDefault('opts.fit_sd',1);
SetDefault('opts.fit_alpha',[]);        % 0.32 for 1sd, 0.05 for 2sd, 0.003 for 3sd
SetDefault('opts.fit_excludeang',true);
SetDefault('opts.fitlm_domean',true);
SetDefault('opts.fitlm_showref',[]);
SetDefault('opts.pix2deg',1);
SetDefault('opts.Color',[]);
SetDefault('opts.FontSize',18);
SetDefault('opts.colormap',[[0,0,0];hot].^0.3);
SetDefault('opts.showcolorbar',false);
SetDefault('opts.NBin',24);
SetDefault('opts.MarkerSize',[]);
SetDefault('opts.LineWidth',1.5);

%-- Robustness
isrobustfit = contains(lower(opts.showfitted),'robust');
opts.showfitted = strrep(lower(opts.showfitted),'robust','');

%-- Set main options
switch lower(opts.mothod)
    case {'line','linfit','fitlm'}
        opts.mothod     = 'none';
        opts.showfitted = 'fitlm';
        opts.fit_excludeang = false;
    case {'hist','histogram'}
        opts.mothod     = 'histogram';
        opts.overplot   = false;
end
condNames   = opts.condNames;
  if isempty(condNames)
      condNames = conds;
      condNames(ismember(condNames,'bb')) = {'Broadband'};
      condNames(ismember(condNames,'a'))  = {'Alpha'};
      condNames(ismember(condNames,'e'))  = {'ERP'};
  else
      assert(numel(condNames)==numel(conds),...
          'opts.condNames and conds must have the same size');
  end

%-- Solve options
plotmethod  = opts.mothod;
if isempty(opts.overplot)
    isoverplt = strcmpi(plotmethod,'none');
else
    isoverplt = opts.overplot;
end
jitterrate  = opts.scatt_jitter;
iscircshift = opts.ang_circshift;
anglimpad   = opts.ang_limpad;
cfactor     = opts.pix2deg;
FontSize    = opts.FontSize;
ColMap      = opts.colormap;
showclbar   = opts.showcolorbar;
plt_NBin        = opts.NBin;
plt_SizeData    = opts.MarkerSize;
plt_LineWidth   = opts.LineWidth;

%-- Solve fit options
switch opts.showfitted
    case {'tls'},               showfitted = 1;
    case {'ellipse'},           showfitted = 2;
    case {'fitlm','linfit'},    showfitted = 3;
    otherwise,                  showfitted = 0;
end
if isempty(opts.fit_alpha)
    fit_alpha = (1-(1+erf(opts.fit_sd*2.^(-1/2)))/2)*2;
else
    fit_alpha  = opts.fit_alpha;
end
fit_exang  = opts.fit_excludeang;
fit_domean = opts.fitlm_domean;

%%

%-- 
selroi = reshape(find(ismember(roiList,selrois)),1,[]);

if ~exist('prfs_corr','var')||isempty(exist('prfs_corr','var'))
    prfs_corr = keepfields(prfs, prffields(endsWith(prffields,{'_fit','_corr','_R2'})));
end


%% Plot figure
%-- Loop for Targets
hF=gobjects(0);
for trgt = reltarget

%-- Solve 'VS'
sbtrgt    =  strsplit(trgt{:},'_');
rel4conds = all(ismember(sbtrgt,allconds));
rel4parms = ~rel4conds && all(ismember(sbtrgt,allparms));
assert(rel4conds || rel4parms, 'Target ''%s'' is invalid',trgt{:});

%-- Check direction
if rel4conds,   lpprops = parms;    % rel4conds
else,           lpprops = conds;    % rel4parms
end

%-- Check specific options
showref     = opts.fitlm_showref;
  if isempty(showref)
      showref = showfitted == 3;  % default = true if showfitted = 'fitlm'
  end
showdiag    = opts.showdiag;
  if isempty(showdiag)
      showdiag = ~showref;
  end

%-- Figure shape
noclmn = length(lpprops)==1 || isoverplt;
if noclmn
    ix = ceil(sqrt(length(selroi)));    iy = ceil(length(selroi)./ix);
else
    ix = length(selroi);  iy = length(lpprops);
end

switch lower(plotmethod)    % if ~noclmn,     figSiz = [ix*300+100 iy*300+100];
    case {'histogram'},     figSiz = [ix*300+100 iy*330+10];
    otherwise,              figSiz = [ix*450+40 iy*330+40];
end
hF(end+1) = figure('MenuBar','none');
hF(end).Position = [max([1,1],hF(end).Position(1:2)-[figSiz - hF(end).Position(3:4)]./2), figSiz];
hT=tiledlayout(iy,ix,'Padding','compact','TileSpacing','compact');

%-- Line properties
if isempty(opts.Color)
    plcol = get(hF(end),'defaultAxesColorOrder');
    if strcmpi(plotmethod,'histogram')
        plcol = circshift(plcol,-2,1);
    end
else
    plcol = matlab.graphics.internal.convertToRGB(opts.Color);
end
ncol = size(plcol,1);
set(hF(end),'defaultAxesColorOrder',plcol);
if strcmpi(plotmethod,'histogram')
    lineSpec = 'w-';  sublineSpec = 'w--';  reflineCol = [1 1 1].*0.7;
else
    lineSpec = 'k--'; sublineSpec = 'k:';   reflineCol = [1 1 1].*0.3;
end

%%-- Loop in Figure
ii = 1; newax = true;
yL=gobjects(0);
for iroi = selroi
    hR=gobjects(0); hl=gobjects(0); legtitle = {};
    for sbprop = lpprops
        hax = nexttile(ii);
        if rel4conds
            parmX = sbprop{:};  parmY = sbprop{:};
            condX = sbtrgt{1};  condY = sbtrgt{2};
            sbfldF  = [sbprop{:} '_fit'];
            sbfldR  = [sbprop{:} '_corr'];
            sbfldR2 = [sbprop{:} '_R2'];
        else %rel4parms
            parmX = sbtrgt{1};  parmY = sbtrgt{2};
            condX = sbprop{:};  condY = sbprop{:};
            sbfldF  = [trgt{:} '_' sbprop{:} '_fit'];
            sbfldR  = [trgt{:} '_' sbprop{:} '_corr'];
            sbfldR2 = [trgt{:} '_' sbprop{:} '_R2'];
        end
        sbfldX  = [parmX '_' condX];
        sbfldY  = [parmY '_' condY];
        %%% roi %%%
        %-- set data
        if iscell(prfs.(sbfldX))
            datX = prfs.(sbfldX){iroi};
            datY = prfs.(sbfldY){iroi};
        else
            datX = prfs.(sbfldX)(iroi,:);
            datY = prfs.(sbfldY)(iroi,:);
        end
        if isfield(prfs_corr,sbfldF)
            if iscell(prfs_corr.(sbfldF))
                datF = prfs_corr.(sbfldF){iroi};
            else
                datF = prfs_corr.(sbfldF)(iroi,:);
            end
            if showfitted==3 && size(datF,1)<2     % error if datF has invalid data
                close(hF(end));
                error('prfs_corr must have the %s field with at least two rows',sbfldF);
            end
        else
            if showfitted==3     % warn if showfitted = 'fitlm'
                warning('prfs_corr does not have the %s field',sbfldF);
            end
            datF = [nan;nan];
        end
        
        %-- unit change
        switch parmX
            case {'ecc', 'rfsize'}
            %-- (pixel -> degree)
                datX = datX .* cfactor;
                if size(datF,1)>=2
                datF(2,:) = datF(2,:) ./ cfactor;
                end
        end
        switch parmY
            case {'ecc', 'rfsize'}
            %-- (pixel -> degree)
                datY = datY .* cfactor;
                if size(datF,1)>=2
                datF(:,:) = datF(:,:) .* cfactor;
                end
        end
        %-- for circular relations
        datShiftX = 0;    datShiftY = 0;
        if rel4conds
            showlim = anglimpad;
            switch parmX
              case {'ang'}
                %-- circular shift
                if iscircshift
%                     datShift = round(nanmedian([datX,datY]))-180;
                    datShiftX = round(nanmedian(datX))-180;
                    datShiftY = round(nanmedian(datY))-180;
                    datX = mod(datX-datShiftX,360)+datShiftX;
                    datY = mod(datY-datShiftY,360)+datShiftY;
                else
                    datShiftX = 0;     datShiftY = 0;
                end
                %-- replicate around 90 deg
                if isscalar(showlim), showlim = [-1 1].*showlim;    end
                datX = [datX, datX,     datX,     datX+360, datX+360, datX+360, datX-360, datX-360, datX-360];
                datY = [datY, datY+360, datY-360, datY,     datY+360, datY-360, datY,     datY+360, datY-360];
                showch = datX<(360+showlim(2)) & datX>(showlim(1)) & datY<(360+showlim(2)) & datY>(showlim(1));
                datX = datX(showch);
                datY = datY(showch);
            end
        else
            showlim = [0 0];
        end
        
        %-- set axis
        switch parmX
            case {'R2','xR2','xval'}
                xmin = min(0,floor(min(datX,[],'all','omitnan')/10)*10);
                if     xmin<-50,  xgaps = 50;
                elseif xmin<0,    xgaps = 25;
                else,             xgaps = 20;
                end
                xlim([xmin 100]);  xticks(fliplr(100:-xgaps:xmin));
                xNbin = plt_NBin;
                xptitle = 'R^2 (%)';
            case {'ecc'}
                xlim([0 8.5]);  xticks(0:2:8);
                xNbin = plt_NBin;
                xptitle = 'Eccentricity (deg)';
            case {'rfsize'}
                xlim([0 10]);  xticks(0:2:10);
                xNbin = plt_NBin;
                xptitle = 'Size (deg)';
            case {'ang'}
%                 xlim([0 360]); xticks([0 90 180 270 360]);
                xlim([0 360]+[showlim(1) showlim(2)]+datShiftX); xticks([0 90 180 270 360]+datShiftX);
%                 xlim([0 720]+datShift); xticks([-180 -90 0 90 180 270 360 450 540]);
%                                     xticklabels([180 270 0 90 180 270 360 90 180]);
                xNbin = round(plt_NBin*(diff(showlim(1:2))./360+1));
                xptitle = 'Angle (deg)';
            otherwise
                xlim(getdatarange(datX));
                xNbin = plt_NBin;
                xptitle = parmX;
        end
        switch parmY
            case {'R2','xR2','xval'}
                ymin = min(0,floor(min(datY,[],'all','omitnan')/10)*10);
                if     ymin<-50,  ygaps = 50;
                elseif ymin<0,    ygaps = 25;
                else,             ygaps = 20;
                end
                ylim([ymin 100]);  yticks(fliplr(100:-ygaps:ymin));
                yNbin = plt_NBin;
                yptitle = 'R^2 (%)';
            case {'ecc'}
                ylim([0 8.5]);  yticks(0:2:8);
                yNbin = plt_NBin;
                yptitle = 'Eccentricity (deg)';
            case {'rfsize'}
                ylim([0 10]);  yticks(0:2:10);
                yNbin = plt_NBin;
                yptitle = 'Size (deg)';
            case {'ang'}
%                 ylim([0 360]); yticks([0 90 180 270 360]);
                ylim([0 360]+[showlim(1) showlim(2)]+datShiftY); yticks([0 90 180 270 360]+datShiftY);
%                 ylim([0 720]+datShift); yticks([-180 -90 0 90 180 270 360 450 540]);
%                                     yticklabels([180 270 0 90 180 270 360 90 180]);
                yNbin = round(plt_NBin*(diff(showlim(1:2))./360+1));
                yptitle = 'Angle (deg)';
            otherwise
                ylim(getdatarange(datY));
                yNbin = plt_NBin;
                yptitle = parmY;
        end
        xctitle = condNames{find(ismember(conds,condX),1)};
        yctitle = condNames{find(ismember(conds,condY),1)};
        if rel4conds
            clmntitle = xptitle;
            xglabel = xctitle;  yglabel = yctitle;
        else % rel4parms
            clmntitle = xctitle;
            xglabel = xptitle;  yglabel = yptitle;
        end
        legtitle{end+1} = clmntitle;
                
        %-- plot main
        hold on;
        pli = mod(max(numel(hR),numel(hl)),ncol)+1;
        switch lower(plotmethod)
            case {'histogram'}
                histopts = {'XBinLimits',xlim,'YBinLimits',ylim,'NumBins',[xNbin yNbin]};
                %-- plot histogram & set axis
                hR(end+1) = histogram2(datX,datY,histopts{:},'DisplayStyle','tile','ShowEmptyBins','on');
                hR(end).EdgeColor = 'none';
                view([0 90]); axis square;
            case {'scatter'}
                histopts = {'XBinLimits',xlim,'YBinLimits',ylim,'NumBins',[xNbin yNbin].*2};
                %-- plot scatter
                roughN = histcounts2(datX,datY,histopts{:});  roughN(roughN == 0) = [];
                roughN = prctile(roughN,0.85);
                MFalpha = 0.5*10.^-sqrt(log10(max(roughN-20,1)).*1.2).*.1; MFalpha = 1;
                datXpl = datX + randn(size(datX)).*diff(xlim)*jitterrate./5;
                datYpl = datY + randn(size(datY)).*diff(ylim)*jitterrate./5;
                hR(end+1) = scatter(datXpl,datYpl,[],plcol(pli,:),'o',...
                    'MarkerFaceColor','flat','MarkerFaceAlpha',MFalpha,'MarkerEdgeAlpha',MFalpha);
                if ~isempty(plt_SizeData)
                    hR(end).SizeData   = plt_SizeData;
                end
        end
        
        %-- plot fitted line or ellipse
        fitname = '';
        if showfitted && ~(fit_exang && ismember(sbprop,{'ang'}))
            t = xlim;
            t = linspace(t(1),t(2),100);
            
            %-- exclude circular data for anlge & exclude outliers if robust fit is required
            switch sbprop{:}
                case {'ang'}
                    ngch = datY < (datX-180) | datY > (datX+180) |...
                           (datX<(0+datShiftX) & datY<(0+datShiftX)) |...
                           (datX>(360+datShiftX) & datY>(360+datShiftX));
                    datX = datX(~ngch);
                    datY = datY(~ngch);
                otherwise
                    if isrobustfit
                        robustalpha = 0.0455;   % 2sigma
                        datDist = sqrt((datX - median(datX)).^2 + (datY - median(datY)).^2);
                        datDist(datDist<=0) = max(min(datDist).*1e-3,eps);
                        pd   = fitdist(datDist','Lognormal');
                        ngch = datDist > pd.icdf(1-robustalpha);
                        datX = datX(~ngch);
                        datY = datY(~ngch);
                    end
            end
                    
            if showfitted == 1      % Total Least Square
                fitname = '-tls';
%                 [eVc,eVl]  = eig(cov(dat1,dat2));
%                 [~,mainax] = diag(eVl);
                
                x = fitline2derror(datX,datY);
                y = x(1)*t+x(2);
                
                hl(end+1) = plot(t,y,'Color',plcol(pli,:),'LineWidth',plt_LineWidth);
                
            elseif showfitted == 2  % Covariance Ellipse
                fitname = '-ellipse';
            
                hl(end+1) = error_ellipse(cov(datX,datY),[median(datX),median(datY)],'conf',1-fit_alpha);
                hl(end).LineWidth = plt_LineWidth;      hl(end).Color=plcol(pli,:);
                hl(end).ZData = [];
            
            elseif showfitted == 3  % Linear Regression (bootstrap)
                fitname = '-fits';
                plfun = @(t,x) x(2,:)'*t + x(1,:)';
                if size(datF,1)>=2
                ys = plfun(t,datF);
                
                if fit_domean,  y1 = mean(ys,1);
                else,       y1 = median(ys,1);
                end
                y2 = prctile(ys,(fit_alpha/2)*100,1);
                y3 = prctile(ys,(1-fit_alpha/2)*100,1);
                
                hl(end+1) = plot(t,y1,'Color',plcol(pli,:),'LineWidth',plt_LineWidth);
                fill([t, fliplr(t)], [y2, fliplr(y3)],plcol(pli,:),'EdgeColor','none','FaceAlpha',0.3);
                end
            end
            
        end
        
        if  newax
          %-- set column title
          if ii<=length(selroi)
              title(roiList(iroi));
          end
          if ~noclmn && (mod(ii,length(selroi))==1 || length(selroi)==1)
              yL(end+1) = ylabel(clmntitle,'VerticalAlignment','bottom','FontWeight','bold');
          end
          set(hax,'FontSize',FontSize)
          
          %-- set diagonal
          if showdiag
              axlim = [min([xlim,ylim]), max([xlim,ylim])];
              plot(axlim,axlim,lineSpec,'LineWidth',0.8);
          end
        
          %-- plot separator
          switch sbprop{:}
              case {'ang'}
                  if showdiag
                  plot(xlim,ylim+360,sublineSpec,'LineWidth',0.8);
                  plot(xlim+360,ylim,sublineSpec,'LineWidth',0.8);
                  end
                  plot(xlim,[1 1].*(0+datShiftX),lineSpec,'LineWidth',1.0);
                  plot([1 1].*(0+datShiftX),ylim,lineSpec,'LineWidth',1.0);
                  plot(xlim,[1 1].*(360+datShiftX),lineSpec,'LineWidth',1.0);
                  plot([1 1].*(360+datShiftX),ylim,lineSpec,'LineWidth',1.0);
              otherwise
                if prod(xlim)<0
                  plot([0 0],ylim,sublineSpec,'LineWidth',0.8);
                end
                if prod(ylim)<0
                  plot(xlim,[0 0],sublineSpec,'LineWidth',0.8);
                end
          end
          
        end
        
        %-- color set
        if strcmpi(plotmethod,'histogram')
            colormap(ColMap);
            if showclbar
                hc=colorbar;
                hc.Ticks = hc.Limits;
                hc.TickLabels = {'0'; 'Max'};
            end
        end
        
        %-- set next axis
        newax = false;
        if ~noclmn,  ii = ii + ix; newax = true;  end
        
        %-- show reference line [Himmelberg 2021] (if next plots in the new axis)
        hasref = false;
        if newax
          if showref
              [hasref,href,reflebel] = plot_ref(trgt{:},iroi,roiList,reflineCol,plt_LineWidth);
              if hasref
                  hl(end+1)         = href;
                  legtitle{end+1}   = reflebel;
              end
          end
        end
        
        %-- set fitting at the top
        uistack(hl,'top');
        %-- reset fit lines
        if newax,  hR=gobjects(0); hl=gobjects(0); end
    end
    
    %-- show reference line [Himmelberg 2021] (if not has reference yet)
    if ~hasref
        if showref
            [hasref,href,reflebel] = plot_ref(trgt{:},iroi,roiList,reflineCol,plt_LineWidth);
            if hasref
                hl(end+1)         = href;
                legtitle{end+1}   = reflebel;
            end
        end
    end
    
    %-- set Legend
    if isoverplt && ii==1
        if isempty(hl) || numel(hR) > numel(hl)
            legend(hR,legtitle,'Location','best','AutoUpdate','off');
        else
            legend(hl,legtitle,'Location','best','AutoUpdate','off');
        end
    end
    
    %-- set next axis
    ii = mod(ii,ix.*iy)+1; newax = true;
end
% if rel4conds && ~isoverplt
    realignYLabels(yL);
    xlabel(hT,xglabel,'HorizontalAlignment','center','FontSize',FontSize);
    if noclmn
    ylabel(hT,yglabel,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',FontSize);
    else
    axPos = hT.Children(1).Position;
    txpos = [1.01,(0.5-axPos(2))/axPos(4)];
    text(hax,txpos(1),txpos(2),yglabel,'Unit','normalized','HorizontalAlignment','center','VerticalAlignment','bottom','Rotation',-90,'FontSize',FontSize);
    end
% end
%     
    set(hF(end),'Name','pRF relations');
    
    
end

if nargout > 0
    varargout{1} = hF;
end

end

function retval = isSTRs(inp)
% Check if characters-like variable or not

retval = iscellstr(inp) || isstring(inp) || ischar(inp);
end

function [hasref,h,label] = plot_ref(hax,Target,iroi,roiList,reflineCol,LineWidth)
% Plot reference  

narginchk(3,5);

%-- get axis
if ~exist('reflineCol','var')
    reflineCol = [];
    if ~ishandle(hax)
        if exist('roiList','var')
            reflineCol = roiList;
        end
        roiList    = iroi;
        iroi       = Target;
        Target     = hax;
        hax        = gca;
    end
end
          
%-- set ref
hasref  = true;
h       = [];
label   = [];
if strcmpi(Target,'ecc_rfsize')
    label = 'fMRI';
    switch roiList{iroi}
        case 'V1', Sl = 0.17; Int = 0.57;
        case 'V2', Sl = 0.25; Int = 0.65;
        case 'V3', Sl = 0.30; Int = 1.24;
        otherwise, hasref = false; label = [];
    end
else
    hasref = false;
end

%-- plot ref
if hasref
    t = xlim(hax);   t = linspace(t(1),t(2),100);
    h = plot(hax, t,t*Sl + Int,'-.','Color',reflineCol,'LineWidth',LineWidth);
end

end

function datrange = getdatarange(data,scaletic)
% Set data range for plot

dmax = max(data,[],'all','omitnan');
dmin = min(data,[],'all','omitnan');
if ~exist('scaletic','var')||isempty(scaletic)
scaletic = round((dmax-dmin)./10,1,'significant');
end
dmax = ceil(dmax/scaletic)*scaletic;
dmin = floor(dmin/scaletic)*scaletic;
if dmax*dmin > 0 && min(abs(dmax),abs(dmin))/(dmax-dmin)<0.1
    if dmin < 0, dmax = 0;
    else,        dmin = 0;
    end
end

datrange = [dmin dmax];

end

