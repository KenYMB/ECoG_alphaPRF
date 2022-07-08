function    sbjectInfo = loadSbjInfo(filename,opt)

% sbjectInfo = LOADSBJINFO(filename)
%   loads subject information from 'filename'.
%   'filename' should be strings or cell array of strings including file
%   extensions.
% 
% sbjectInfo = LOADSBJINFO(filename,'all')
%   loads subject information from 'name*.ext' files.
% 
% See also, SetSubjectsList

% 20220222 Yuasa

%-- Option
narginchk(1,2);
if nargin < 2,  opt = '';  end

%-- Get filenames
if ~iscell(filename)
    filename = {filename};
end
if strcmpi(opt,'all')
    basefilename = filename;
    filename = [];
    for ifile=1:numel(basefilename)
        filepath = which(basefilename{ifile});
        if ~isempty(filepath)
            [fpath,fname,ext] = fileparts(filepath);
            filelist = dir(fullfile(fpath,[fname '*' ext]));
            filename = vertcat(filename,arrayfun(@(S) S.name, filelist, 'UniformOutput', false));
        end
    end
end

%-- Load files
sbjectInfo = table;
for ifile=1:numel(filename)
    sbjectInfo = vertcat(sbjectInfo, readtable(filename{ifile}, 'FileType', 'text'));
end
sbjectInfo = unique(sbjectInfo,'stable');