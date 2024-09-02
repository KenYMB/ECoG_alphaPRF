function savefigauto(fig, filename, varargin)

% SAVEFIGAUTO save figure with the format of hgexport('factorystyle')
% using the Painters renderer.
%   SAVEFIGAUTO(H,'FILENAME')
%     saves figures with the formats of 'fig', 'png', and 'eps'.
%   SAVEFIGAUTO(H,'FILENAME',{'FORMAT'})
%     saves figures with the formats specified in 'FORMAT'.
%   SAVEFIGAUTO(H,'FILENAME',{'VECTOR_FORMAT'},'-vector')
%     saves figures using vector rendering even for complex figures.
%     This option is same as {'Render','painters'}
%   any options compatible with SAVEAS or HGEXPORT can be used.
% 
% see also, saveas, hgexport, exportgraphics

% 20210623 Yuasa
% 20220621 Yuasa - add '-vector' option
% 20240527 Yuasa - update to use exportgraphics in case hgexport is out-of-date
% 20240901 Yuasa - suppress warning of hgexport

%-- parameter
narginchk(2,inf);

optidx = find(cellfun(@ischar,varargin));
isvec = ismember(varargin(optidx),'-vector');
varargin(optidx(isvec)) = [];
nopts = nargin - 2 - sum(isvec);
isvec = any(isvec);

if nargin < 3 || ~mod(nopts,2)  % nopts is even number
    formats = {'fig','png','eps'};
else
    formats = varargin{1};
    varargin(1) = [];
end
if ~iscell(formats)
    formats = {formats};
end

%-- suppress warning of hgexport
warnid = 'MATLAB:hgexport:StyleSheetNotSupportedInFutureRelease';
warnst = warning('query',warnid);

%-- save
try
 for fmt = reshape(formats,1,[])
    switch fmt{:}
        case {'fig','m','mfig'}
            saveas(fig, filename, fmt{:});
        case {'pdf','eps','epsc','eps2','epsc2','meta','svg','ps','psc','ps2','psc2'}
          if exist('hgexport','file')
            if isvec
              hgexport(fig, filename, hgexport('factorystyle'), ...
                  'Renderer', 'painters', 'Format', fmt{:}, varargin{:});
            else
              hgexport(fig, filename, hgexport('factorystyle'), ...
                  'Format', fmt{:}, varargin{:});
            end
          else
            if isvec
              exportgraphics(fig, sprintf('%s.%s',filename, fmt{:}),...
                  'ContentType', 'vector', varargin{:});
            else
              exportgraphics(fig, sprintf('%s.%s',filename, fmt{:}), varargin{:});
            end
          end
        otherwise
          if exist('hgexport','file')
            hgexport(fig, filename, hgexport('factorystyle'), ...
                'Format', fmt{:}, varargin{:});
          else
            exportgraphics(fig, sprintf('%s.%s',filename, fmt{:}), varargin{:});
          end
    end
 end
 warning(warnst);
catch ME
 warning(warnst);
 rethrow(ME);
end
