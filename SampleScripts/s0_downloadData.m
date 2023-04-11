%% s0_downloadData
% Downloading ECoG data from fileserver

% 20220223 Yuasa

%% run in system
assert(~ispc, 'Data downloading is currently not available for Windows.');
run_checkPath;
system([which('ecog_APRFF_00_downloadData.sh') ' ' bidsRootPath]);
