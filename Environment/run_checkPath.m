function run_checkPath()

% Users of the toolbox should update this path to point to where
% they copy the toolbox.

%-- Toolbox Root Path
tbrootPath = fullfile(userpath,'toolboxes','ECoG_alphaPRF');

%-- Call
if isempty(which('checkPath'))
    run(fullfile(tbrootPath,'functions','utils','checkPath'));
else
    checkPath;
end