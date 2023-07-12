function rootPath = bidsRootPath()

% Users of the ECoG_alpha toolbox should update this path to point to
% where they downloaded the BAIR dataset.

%-- Set your root path to download dataset
rootPath = '';

%-- Set the root path as 'analysisRootPath/BIDS' if rootPath does not exist (default)
if ~exist(rootPath,'dir')
    rootPath = fullfile(analysisRootPath, 'BIDS'); 
end

end