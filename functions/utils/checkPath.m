function    checkPath()
% Check path setting

% Dependency: lookaroundToolbox, ProjectName

% 20220224 Yuasa

%-- Set ProjectName
if isempty(which('ProjectName'))
    envPath = dir(fullfile(fileparts(mfilename('fullpath')),'..','Environment','ProjectName.m'));
    addpath(envPath(1).folder);
end

%-- Check availability of ToolboxToolbox
if ~isempty(which('tbUse'))
    % check some functions
    if isempty(which('ecog_prf_analyzePRF')) || isempty(which('bidsEcogGetMatchedAtlas')) || isempty(which('ft_defaults')) || isempty(which('analyzePRF'))
        tbUse(ProjectName);
    end
else
    % check toolboxes
    if isempty(which('ecog_prf_analyzePRF'))
        tbPath = fileparts(mfilename('fullpath'));
        [~,tbPath] = fileattrib(fullfile(tbPath,'..','..'));
        tbPath = tbPath.Name;
        warning('Setting path to %s',ProjectName);
        addpath(genpath(tbPath));
    end
    if isempty(which('bidsEcogGetMatchedAtlas'))
        tbPath = lookaroundToolbox('bidsEcogGetMatchedAtlas.m','silent');
        if ~isempty(tbPath)
            [~,tbPath] = fileattrib(fullfile(tbPath,'..','..'));
            tbPath = tbPath.Name;
            warning('Setting path to ECoG_utils');
            addpath(genpath(tbPath));
        else
            error('Please set path to ECoG_utils');
        end
    end
    if isempty(which('ft_defaults'))
        tbPath = lookaroundToolbox('ft_defaults.m','silent');
        if ~isempty(tbPath)
            [~,tbPath] = fileattrib(fullfile(tbPath,'.'));
            tbPath = tbPath.Name;
            warning('Setting path to fieldTrip');
            addpath(tbPath);
            ft_defaults;
        else
            error('Please set path to fieldTrip');
        end
    end
    if isempty(which('analyzePRF'))
        tbPath = lookaroundToolbox('analyzePRF.m','silent');
        if ~isempty(tbPath)
            [~,tbPath] = fileattrib(fullfile(tbPath,'.'));
            tbPath = tbPath.Name;
            warning('Setting path to analyzePRF');
            addpath(genpath(tbPath));
        else
            error('Please set path to analyzePRF');
        end
    end
end

