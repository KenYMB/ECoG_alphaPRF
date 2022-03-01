%% s0_downloadData
% Downloading ECoG data from fileserver

% 20220223 Yuasa

%% run in system
assert(~ispc, 'Data downloading is currently not available for Windows.');
system([which('ecog_APRF_00_downloadData.sh') ' ' analysisRootPath]);
