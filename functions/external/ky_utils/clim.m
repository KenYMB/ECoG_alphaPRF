function varargout = clim(varargin)
%CLIM COLOR limits.
%   CL = CLIM             gets the color limits of the current axes.
%   CLIM([CMIN CMAX])     sets the color limits.
%   CLMODE = CLIM('mode') gets the color limits mode.
%   CLIM(mode)            sets the color limits mode.
%                            (mode can be 'auto' or 'manual')
%   CLIM(AX,...)          uses axes AX instead of current axes.
%
%   CLIM sets or gets the CLim or CLimMode property of an axes.
%
%   See also XLIM, YLIM, ZLIM.

%   20170205 Yuasa: modify XLIM @matlab2016a
%   20200311 Yuasa: update for matlab2016b or later

if verLessThan('matlab','9.1')     % matlab2016a or earlier
    narginchk(0,2);
    nargoutchk(0,1);
    if nargin == 0
        a = get(gca,'CLim');
    else
        arg1 = varargin{1};
        if nargin==2, arg2 = varargin{2}; end
        if isscalar(arg1) && ishghandle(arg1) && isprop(arg1,'CLim')
            ax = arg1;
            if nargin==2
                val = arg2;
            else
                a = get(ax,'CLim');
                varargout{1} = a;
                return
            end
        else
            if nargin==2
                error(message('MATLAB:xlim:InvalidNumberArguments'))
            else
                ax = gca;
                val = arg1;
            end
        end

        matlab.graphics.internal.markFigure(ax);
        if ischar(val)
            if(strcmp(val,'mode'))
                a = get(ax,'CLimMode');
            else
                set(ax,'CLimMode',val);
            end
        else
            set(ax,'CLim',val);
        end
    end
    if exist('a','var'), varargout{1} = a; end
else
    varargout = matlab.graphics.internal.ruler.rulerFunctions(mfilename, nargout, varargin);
end