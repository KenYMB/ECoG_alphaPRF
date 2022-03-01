%% run_tbUse
% Initialize for users using ToolboxToolbox

% 20220225 Yuasa

if isempty(which('checkPath'))
    if exist('tbUse','file')
        tbUse ECoG_alphaPRF;
    else
        error('Please set path to ECoG_alphaPRF')
    end
end

checkPath;
