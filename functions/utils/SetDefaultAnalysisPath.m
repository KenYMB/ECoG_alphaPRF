function varargout = SetDefaultAnalysisPath(PthName, defval, varname)
% SETDEFAULTANALYSISPATH 
%   sets path variables in your workspace directly. 
%   Variables which are already exist are not updated.
% 
% path = SETDEFAULTANALYSISPATH('BIDS, DAT, or FIG')
%   reads and return the path of the specified pass type from your workspace.
%   or returns the default path if it does not exist in your workspace.
% 
% path = SETDEFAULTANALYSISPATH('BIDS, DAT, or FIG', DirectoryName, [VarName])
%   returns the path with DirectoryName if VarName is not exist or epmty,
%   returns the path with VarName if VarName is a directory name,
%   or passthrough VarName if VarName is a path to a directory.
% 
% SETDEFAULTANALYSISPATH('BIDS, DAT, or FIG', DirectoryName, VarName)
%   sets the path directory to VarName instead of returning the path.
% 
% SETDEFAULTANALYSISPATH('SHOW')
%   prints current path to 'BIDS', 'DAT', and 'FIG'

% 20220224 Yuasa
% 20220412 Yuasa - Updates to return path
% 20220428 Yuasa - refer bidsPth/datPth/figPth in the workspace
% 20220518 Yuasa - add 'SHOW' mode

narginchk(0,3);

% Set default pass
bidsPth = fullfile(bidsRootPath,'derivatives',filesep);
datPth  = fullfile(analysisRootPath,'Data',filesep);
figPth  = fullfile(analysisRootPath,'Figures',filesep);

% If argument is exist
if nargin
    
  if strcmpi(PthName,'SHOW')
    
    bidsPth = SetDefaultAnalysisPath('BIDS');
    datPth  = SetDefaultAnalysisPath('DAT');
    figPth  = SetDefaultAnalysisPath('FIG');
    
    fprintf(['--------------------------------------------------------------------------\n'...
             'bidsPth = %s\n'...
             'datPth  = %s\n'...
             'figPth  = %s\n'...
             '--------------------------------------------------------------------------\n'],...
             bidsPth,datPth,figPth);
      
  else
      
    switch upper(PthName)
        case {'BIDS'}
            targetname = 'bidsPth';
        case {'DAT','DATA'}
            targetname = 'datPth';
        case {'FIG','FIGURE'}
            targetname = 'figPth';
    end
    if evalin('base',sprintf('exist(''%s'',''var'')&&~isempty(%s)&&ischar(%s)',...
                                targetname,targetname,targetname))
        target = evalin('base',targetname);
    else
        target = eval(targetname);
    end
    
    % If defval is specified
    if nargin > 1
        %-- Get variable value if exist
        existval = false;
        if exist('varname','var') && ~isempty(varname)
            assert(ischar(varname) || (isstring(varname) && isscalar(varname)), ...
                'VarName must be either a character vector or a string scalar.');
            
            varnameparts = strsplit(varname,'.');
            existval = evalin('caller',sprintf('exist(''%s'',''var'')',varnameparts{1}));
            %--- loop for struct
            for ifld = 2:length(varnameparts)
                if existval
                    existval = evalin('caller',sprintf('isfield(%s,''%s'')', ...
                        strjoin(varnameparts(1:(ifld-1)),'.'), varnameparts{ifld}));
                end
            end
        end
        if existval
            inpval = evalin('caller',varname);
        else
            inpval = '';
        end
        
        %-- If inpval is empty, then output defval
        if isempty(inpval)
            target = fullfile(target,defval);
            
        else
            assert(ischar(inpval) || (isstring(inpval) && isscalar(inpval)), ...
                '%s must be either a character vector or a string scalar.',varname);
            
        %-- If inpval is a directory name, then concatenate with path and directory
            if numel(strsplit(inpval,filesep)) == 1
                target = fullfile(target,inpval);
                
        %-- If inpval is a path to a directory name, then passthrough inpval
            else
                target = inpval;
            end
        end
    end

    %-- Set variable directory ui variable name is specified & no return value is required
    if ~nargout && exist('varname','var')
        evalin('caller',sprintf('%s = ''%s'';',varname,target));
    else
        varargout = {target};
    end
    
  end
  
% If argument is not exist
else
    
    numsetpath = 0;

    %-- Input path
    numsetpath = numsetpath + ...
        SetDefault('bidsPth', bidsPth,'base');

    %-- Output path 
    numsetpath = numsetpath + ...
        SetDefault('datPth', datPth,'base');
    numsetpath = numsetpath + ...
        SetDefault('figPth', figPth,'base');

    %-- Return numbers of path set in this function
    if nargout~=0
        varargout = {numsetpath};
    end
    
end
end
