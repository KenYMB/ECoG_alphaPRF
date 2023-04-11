function run_checkPath()

% Users of the toolbox should update this path to point to where
% they copy the toolbox.

%-- Toolbox Root Path
Project    = 'ECoG_alphaPRF';
tbrootPath = fullfile(userpath,'toolboxes',Project);

%-- Call
if isempty(which('checkPath')) || isempty(which('ProjectName')) || ~strcmp(ProjectName,Project)
    run(fullfile(tbrootPath,'functions','utils','checkPath'));
else
    checkPath;
end