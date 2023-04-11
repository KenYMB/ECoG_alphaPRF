function    tbpass = lookaroundToolbox(funcname,issilent)
% ToolboxPath = lookaroundToolbox(FunctionName)
% ToolboxPath = lookaroundToolbox(FunctionName,'silent')
% Look around for function name, and return the path

% 20220523 Yuasa

narginchk(1,2);
if isempty(funcname)
    error('FunctionName is required');
end
issilent = exist('issilent','var') && ~isempty(issilent) && ...
            ischar(issilent) && strcmpi(issilent,'silent');

%-- Set search locations
defPath = userpath;
if isempty(defPath)
    if ispc,    defPath = getenv('userprofile');
    else,       defPath = '~';
    end
    tmp = dir(fullfile(defPath,'matlab'));
    if isempty(tmp)
        tmp = dir(fullfile(defPath,'matlab'));
    end
    if isempty(tmp)
        tmp = dir(defPath);
    end
    defPath = tmp(1).folder;
end
finddepth   = {'', '*', fullfile('*','*'), fullfile('*','*','*')};
ndepth = length(finddepth);

%-- Check path
tbpass = dir(which(funcname));
isSetPath = ~isempty(tbpass);

%-- Search file
defPath     = {fullfile(defPath,'toolboxes'),fullfile(defPath,'TOOLBOXES'),...
                fullfile(defPath,'toolbox'),fullfile(defPath,'TOOLBOX')};
iloc = 1;
while ~isSetPath && iloc <= length(defPath)*ndepth
    tbpass = dir(fullfile(defPath{ceil(iloc./ndepth)},finddepth{mod(iloc-1,ndepth)+1},funcname));
    isSetPath = ~isempty(tbpass);
    iloc = iloc + 1;
end
iloc = 1;
while ~isSetPath && iloc <= length(defPath)*ndepth
    tbpass = dir(fullfile(defPath{ceil(iloc./ndepth)},finddepth{mod(iloc-1,ndepth)+1},[funcname '.m']));
    isSetPath = ~isempty(tbpass);
    iloc = iloc + 1;
end

%-- Additional search file
defPath     = {'.','..',fullfile('..','..')};
iloc = 1;
while ~isSetPath && iloc <= length(defPath)*ndepth
    if mod(iloc-1,ndepth)+1 < (6 - ceil(iloc./ndepth))
    tbpass = dir(fullfile(defPath{ceil(iloc./ndepth)},finddepth{mod(iloc-1,ndepth)+1},funcname));
    isSetPath = ~isempty(tbpass);
    end
    iloc = iloc + 1;
end
iloc = 1;
while ~isSetPath && iloc <= length(defPath)*ndepth
    if mod(iloc-1,ndepth)+1 < (6 - ceil(iloc./ndepth))
    tbpass = dir(fullfile(defPath{ceil(iloc./ndepth)},finddepth{mod(iloc-1,ndepth)+1},[funcname '.m']));
    isSetPath = ~isempty(tbpass);
    end
    iloc = iloc + 1;
end

%-- Output
if isSetPath
    tbpass = tbpass(1).folder;
else
    tbpass = '';
    if ~issilent,  warning('%s is not found',funcname);  end
end
