function    checkPath(Project)
% Check path setting

% Dependency: lookaroundToolbox, ProjectName

% 20220224 Yuasa
% 20220608 Yuasa - update

%-- Set ProjectName
if ~exist('Project','var')||isempty(Project)
    if isempty(which('ProjectName'))
        envPath = dir(fullfile(fileparts(mfilename('fullpath')),'..','..','Environment','ProjectName.m'));
        assert(~isempty(envPath),'ProjectName not found');
        addpath(envPath(1).folder);
    end
    Project = ProjectName;
end

%-- Check availability of ToolboxToolbox
if ~isempty(which('tbUse'))
    % check some functions
    if isempty(which('ecog_prf_analyzePRF')) || isempty(which('bidsEcogGetMatchedAtlas')) || isempty(which('ft_defaults')) || isempty(which('analyzePRF'))
        tbUse(Project);
    end
else
    % check toolboxes
    if isempty(which('ecog_prf_analyzePRF'))
        tbPath = fileparts(mfilename('fullpath'));
        [~,tbPath] = fileattrib(fullfile(tbPath,'..','..'));
        tbPath = tbPath.Name;
        warning('Setting path to %s',Project);
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

