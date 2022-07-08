function rootPath = analysisRootPath()

% Users of the ECoG_alpha toolbox should update this path to point to
% where they want the processed data and figures to be written to.

%-- Set your root path for analysis
rootPath = '';

%-- When you copy analysisRootPath to working space, then please comment out next line.
rootPath = regexprep(mfilename('fullpath'),mfilename,'');

%-- Set toolbox root directory if rootPath does not exist (default)
if ~exist(rootPath,'dir')
    rootPath = fullfile(fileparts(mfilename('fullpath')),'..','..');
    [st,fileinfo] = fileattrib(rootPath);
    rootPath = fileinfo.Name;
end

