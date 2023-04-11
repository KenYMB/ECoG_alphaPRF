% ecog_APRFF_INITc_postfix
% General script to set postfix about pRF models for file or figure name

% 20201014 Yuasa
% 20210825 Yuasa - add options

%% Parameters
% Inputs: 
%   subjectList (if allowmixbeta is true)
% 
%   average        ='runs';
%   smoothingMode  ='decimate';
%   smoothingN     = 3;
%   prfmodel       ='linear';
%   gaussianmode   ='gs';
% 
% %   allowbeta      = false;     % if true, load FaCLbB
% %   allowwide      = false;     % if true, load FaCLbW, FaCLbBW
% %     OR
% %   allowmixbeta   = false;     % if true, set beta options automatically
% % 
% %   allowlag       = false;     % if true, load 'lag' data
% %   linalpha       = false;     % if true, load FaCRb
% %   noACorrect     = false;     % if true, load FaLb
% 
% Outputs: 
%   postfix
%   [modeldataID]
%   [prfID]
%   [alphaType]


%% load analyzePRF
switch smoothingMode
    case {'none'}
        smoothingN = 1;
end

postfix      = '';
postfix      = sprintf('%s-%s-%s_avg-%s',postfix,prfmodel,gaussianmode,average);
if smoothingN~=1,  postfix = sprintf('%s_%s%d',postfix,smoothingMode,smoothingN);  end

%% Compute options
ecog_APRFF_INITsub_loadoptions;

%-- process options
if noACorrect
    postfix     = sprintf('%s_noACorr',postfix);
end
if linalpha 
    postfix     = sprintf('%s_linAlpha',postfix);
end
if allowmixbeta
    postfix     = sprintf('%s_mixbeta',postfix);
elseif allowbeta && allowwide
    postfix     = sprintf('%s_betawide',postfix);
elseif allowbeta
    postfix     = sprintf('%s_beta',postfix);
elseif allowwide
    postfix     = sprintf('%s_wide',postfix);
end
if allowlag
    postfix     = sprintf('%s_lag',postfix);
end
