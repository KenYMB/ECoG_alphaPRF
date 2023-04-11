% ecog_APRFF_INITd_threshold
% General script to set postfix about pRF models for file or figure name

% 20210903 Yuasa
% 20210916 Yuasa - update to load cod-perm files
% 20211124 Yuasa - outsource getthreshold

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
% %   usefulltsxR2   = false;
% % 
% %   broadbandType  = 'bbS';
% %   alphaType      = 'FaCLb';
% % 
% %   allowlag       = false;
% %   recmpt         = false;
% 
% Outputs:
%   threshold_bb   	: threshold for broadband [%]
%   threshold_a     : threshold for alpha [%]
%   eclimit         : threshold in eccentricity [pixel]
% 
% Example:
% elec_ok = ~(prf_all_bb.xval<=threshold_bb | prf_all_a.xval<=threshold_a ...
%          | prf_all_bb.ecc >= thresh_ecc | prf_all_a.ecc >= thresh_ecc); 

%% Set threshold
%%%% LIST (broadband alpha)
% [linear gs wangprobchs]
%   no-lag mix       29 20
%   lag    mix       29 24
% [linear gs visualchs]
%   no-lag mix       29 19
%   lag    mix       29 23
% [linear gs allchs]
%   no-lag mix       21 19
%   lag    mix       24 23
% [linear og wangprobchs]
%   no-lag mix       30 21
%   lag    mix       27 24

%-- prepare postfix
if ~exist('postfix','var')
    ecog_APRFF_INITc_postfix;
end

%-- set default path
SetDefault('prfstatPth','pRFanalysis');

%-- Load & set threshold
alpha_cod = 0.95;
[threshold_bb,threshold_a,threshold_erp] = ecog_prf_getthreshold(prfstatPth,R2mode,selectchs,postfix,alpha_cod);
eclimit   = 50;
