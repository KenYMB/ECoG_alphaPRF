function savefigauto(fig, filename, varargin)

% SAVEFIGAUTO save figure with the format of hgexport('factorystyle')
% using the Painters renderer.
%   SAVEFIGAUTO(H,'FILENAME')
%     saves figure with the formats of 'fig', 'png', and 'eps'.
%   SAVEFIGAUTO(H,'FILENAME',{'FORMAT'})
%     saves figure with the formats specified in 'FORMAT'.
%   any options
% 
% see also, hgexport

% 20210623 Yuasa

%-- parameter
narginchk(2,inf);

if nargin < 3 || ~mod(nargin,2)
    formats = {'fig','png','eps'};
else
    formats = varargin{1};
    varargin(1) = [];
end
if ~iscell(formats)
    formats = {formats};
end

%-- save
for fmt = reshape(formats,1,[])
    switch fmt{:}
        case {'fig','m','mfig'}
            saveas(fig, filename, fmt{:});
        case {'pdf','eps','epsc','eps2','epsc2','meta','svg','ps','psc','ps2','psc2'}
            hgexport(fig, filename, hgexport('factorystyle'), ...
                'Renderer', 'painters', 'Format', fmt{:}, varargin{:});
        otherwise
            hgexport(fig, filename, hgexport('factorystyle'), ...
                'Format', fmt{:}, varargin{:});
    end
end
