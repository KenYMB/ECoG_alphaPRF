function    checkPath()
% Check path setting

% 20220224 Yuasa

%-- Check availability of ToolboxToolbox
if exist('tbUse','file')
    % check some functions
    if ~exist('analysisRootPath','file') || ~exist('ecog_plotGrid','file') || ~exist('analyzePRF','file')
        tbUse('ECoG_alphaPRF');
    end
else
    % check toolboxes
    if ~exist('analysisRootPath','file')
        error('Please set path to ECoG_alphaPRF')
    elseif ~exist('ecog_plotGrid','file')
        error('Please set path to ECoG_utils')
    elseif ~exist('analyzePRF','file')
        error('Please set path to analyzePRF')
    end
end

