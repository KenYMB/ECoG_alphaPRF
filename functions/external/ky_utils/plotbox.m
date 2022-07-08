function    p = plotbox(h,xbounds,ybounds,boxcol,varargin)

%   PLOTBOX([xmin xmax],[ymin ymax],C) creates box (rectangle) specified with
%   [xmin xmax] and [ymin ymax]. If either of [xmin xmax] or [ymin ymax] is
%   empty, the width or the highet of the rectangle is set to the full range
%   of the current axis. C determines the box colors.
% 
%   PLOTBOX({[x1min x1max],[x2min x2max],...},[ymin ymax],C) 
%   PLOTBOX([xmin xmax],{[y1min y1max},[y2min y2max],...},C) 
%   PLOTBOX({[x1min x1max],[x2min x2max],...},{[y1min y1max},[y2min y2max],...},C) 
%   creats multiple boxes (rectangles) specified with each pair of [x1min x1max]
%   and [y1min y1max], [x2min x2max] and [y2min y2max], ...
% 
%   PLOTBOX(...,Name,Value) specifies patch properties using one or more
%   Name,Value pair arguments.
% 
%   PLOTBOX(container,...) creates box in the axes, group, or transform
%   specified by container, instead of in the current axes.
%
%   P = PLOTBOX(...) returns the patch object that contains the data for all
%   the boxes.
% 
%   Exaple:
%       t = 0:1/250:10;
%       x = 2*sin(2*pi*20*t) + 0.5*sin(2*pi*80*t) + 20*sin(2*pi*60*t) + (rand(1,length(t))-0.5);
%       pwelch(x,[],[],[],250);
%       plotbox([55 65],[],'k');
% 
% See also PATCH, RECTANGLE.

% 20200203 Yuasa
% 20220106 Yuasa - minor bug fix

%-- check arguments
narginchk(2,inf);

if ishandle(h(1))
    narginchk(3,inf);
    h = h(1);
    if nargin < 4,      boxcol   = {'none'};     end
else
    if nargin > 4,      varargin = [{boxcol} varargin];     end
    if nargin > 3,      boxcol   = ybounds;
    else,               boxcol   = {'none'};
    end
    ybounds  = xbounds;
    xbounds  = h;
    h        = gca;
end

if ~iscell(xbounds),    xbounds = {xbounds};    end
if ~iscell(ybounds),    ybounds = {ybounds};    end
if ~iscell(boxcol),     boxcol  = {boxcol};     end
assert(numel(xbounds)==1 || numel(ybounds)==1 || numel(xbounds)==numel(ybounds),...
    'the bounds for x and y must have the same size');

nboxes = max(numel(xbounds),numel(ybounds));
if numel(xbounds)==1,   xbounds = repmat(xbounds,1,nboxes);  end
if numel(ybounds)==1,   ybounds = repmat(ybounds,1,nboxes);  end
if numel(boxcol)==1,    boxcol  = repmat(boxcol,1,nboxes);   end

%-- plot boxes
pp = [];

for iplt = 1:nboxes
    if isempty(xbounds{iplt}), xbounds{iplt} = xlim;    end
    if isempty(ybounds{iplt}), ybounds{iplt} = ylim;    end
    if isempty(boxcol{iplt}),  boxcol{iplt}  = 'none';  end
    xpos = xbounds{iplt}([1 2 2 1]);
    ypos = ybounds{iplt}([1 1 2 2]);
    fcol = boxcol{iplt};
    if ischar(fcol)
      switch fcol
        case {'none'}
            fcol = 'white';
            options = {'FaceColor','none','EdgeColor','k'};
        case {'gray'}
            fcol = 'white';
            options = {'FaceColor',[0.5 0.5 0.5],'EdgeColor','none'};
        otherwise
            %-- not plot edge as default
            options = {'EdgeColor','none'};
      end
    else
            %-- not plot edge as default
            options = {'EdgeColor','none'};
    end
    
    pp(iplt) = patch(h,xpos,ypos,fcol,options{:},varargin{:});
end

if nargout > 0
    p = handle(pp);
end
