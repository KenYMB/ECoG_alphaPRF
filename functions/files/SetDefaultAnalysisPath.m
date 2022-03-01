function varargout = SetDefaultAnalysisPath()
% SETDEFAULTANALYSISPATH directly returns some variables in your workspace
% for the functions in ECoG_alphaPRF toolbox.
% Variables which are already exist are not updated.

% 20220224 Yuasa

numsetpath = 0;

%-- Input path
numsetpath = numsetpath + ...
    SetDefault('dataPth', fullfile(bidsRootPath,'derivatives',filesep),'base');

%-- Output path 
numsetpath = numsetpath + ...
    SetDefault('savePth', fullfile(analysisRootPath,'Data',filesep),'base');
numsetpath = numsetpath + ...
    SetDefault('figPth', fullfile(analysisRootPath,'Figures',filesep),'base');

%-- Return numbers of path set in this function
if nargout~=0
    varargout = {numsetpath};
end