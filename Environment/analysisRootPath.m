function rootPath = analysisRootPath()

% Users of the ECoG_alpha toolbox should update this path to point to
% where they want the processed data and figures to be written to.

%-- Set your root path for analysis
rootPath = '';

%-- Set current directory as the root path if rootPath does not exist and
%-- if this function is placed on the working directory,  
%-- or set the toolbox root directory as the root path (default)
if ~exist(rootPath,'dir')
    rootPath = fileparts(mfilename('fullpath'));
    if ~isequal(rootPath,pwd)
    rootPath = fullfile(rootPath,'..');
    end
    [~,fileinfo] = fileattrib(rootPath);
    rootPath = fileinfo.Name;
end