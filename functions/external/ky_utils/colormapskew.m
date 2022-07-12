function    cm  = colormapskew(cm_name,cm_skew,varargin)

% COLORMAPSKEW(colormap_name,skewness,[range,M])
%   returns an Mx3 skewed colormap array of colormap_name.
%   If range is specified as [min max] in [0 1], it returns the part of
%   colormap array.
%   The colomap is skewed to maximum value if skewness is larger than 1,
%   while it is skewed to minimum value if skewness is smaller than 1.
%   
% COLORMAPSKEW(colormap_name,skewness,[range,M],'inverse')
%   returns a flipped colormap array.
% 
% 20201117 Yuasa

%% Check inputs
narginchk(2,5);
assert(cm_skew>0, 'cm_skew must be larger than 0.');

narginplus = nargin - 2;
if narginplus > 0 && ischar(varargin{narginplus})
  if strcmp(varargin{narginplus},'inverse')
    isinv       = true;
    narginplus  = narginplus - 1;
  else
    isinv       = false;
    warning('''%s'' is unknown parameter');
  end
else
    isinv       = false;
end

if narginplus > 0
    cm_range    = varargin{1};
else
    cm_range    = [0 1];
end
assert(isvector(cm_range) & length(cm_range)==2 & diff(cm_range)>0,...
    'cm_range must be a 2-element vector of increasing numeric values.');
assert(cm_range(1)>=0&cm_range(2)<=1, 'cm_range must be in [0,1].');

if narginplus > 1
    cm_res      = varargin{2};
else
   hf = get(groot,'CurrentFigure');
   if isempty(hf)
      cm_res = size(get(groot,'DefaultFigureColormap'),1);
   else
      cm_res = size(hf.Colormap,1);
   end
end

res_ratio   = 4./diff(cm_range);

%% Make colormap
%-- get colormap
cm_size     = round(cm_res*res_ratio);
cm          = eval(sprintf('%s(%d)',cm_name,cm_size));

%-- get range
cmrange_idx = round(cm_range.*(cm_size-1))+1;

%-- make skewed colormap
cmskew_idx  = round((linspace(0,1,cm_res).^(cm_skew.^((-1).^~isinv)))...
                    .*diff(cmrange_idx)+cmrange_idx(1));
cm          = cm(cmskew_idx,:);

%-- inverse
if isinv
    cm      = flipud(cm);
end

