function rootPath = analysisRootPath()

% Users of the ECoG_alphaPRF toolbox should update this path to point to
% where they want the processed data and figures to be written to.

rootPath = '';
if ~exist(rootPath,'dir')
    rootPath = fullfile(fileparts(mfilename('fullpath')),'..','..');
    [st,fileinfo] = fileattrib(rootPath);
    rootPath = fileinfo.Name;
end
