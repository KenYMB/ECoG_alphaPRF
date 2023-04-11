% ecog_APRFF_INITa_loaddata
% General script to load model time-series and pRF data
%   Need to specify parameters listed in "Parameters"

% 20201008 Yuasa
% 20210421 Yuasa - use ecog_updatePRFdata
% 20210914 Yuasa - use ecog_APRFF_INITsub_loadoptions
% 20211117 Yuasa - separate main processing into ecog_prf_loadprfs
% 20221203 Yuasa - enable to load alpha/broabdan separately

%% Parameters
% Inputs: 
%   subjectList
%   prfPth
% 
%   average        ='runs';
%   smoothingMode  ='decimate';
%   smoothingN     = 3;
%   prfmodel       ='linear';
%   gaussianmode   ='gs';
% 
% %   selectchs      ='allchs';
% % % selectch_exFEF   = true;
% % % selectch_thresh  = 0.05;
% % 
% %   usefulltsR2    = false;
% %   usefulltsxR2   = false;
% % 
% % modeldataID    = [];
% % prfID          = [];
% % 
% % % alphaType      = 'FaCLb';
% % % broadbandType  = 'bbS';
% % % 
% % % allowbeta    = false;     % if true, load FaCLbB
% % % allowwide    = false;     % if true, load FaCLbW, FaCLbBW
% % %   OR
% % % allowmixbeta = false;     % if true, set beta options automatically
% % % 
% % % allowlag     = false;     % if true, load 'lag' data
% % % linalpha     = false;     % if true, load FaCRb
% % % noACorrect   = false;     % if true, load FaLb
% % 
% % skipsummarizeROIs = false;
% 
% Outputs:
%   modeldata_a
%   modeldata_bb
%   prf_params_a
%   prf_params_bb
% 
%   R2mode

%% Prepare parameters
SetDefault('selectchs','allchs');
SetDefault('selectch_exFEF',true);
SetDefault('selectch_thresh',0.05);

SetDefault('usefulltsR2',false);
SetDefault('usefulltsxR2',false);

SetDefault('modeldataID',[]);
SetDefault('prfID',[]);

SetDefault('skipsummarizeROIs',false);

if strcmp(smoothingMode,'none')||smoothingN==1
    if usefulltsR2 || usefulltsxR2
        warning('Variance explained has been computed from full time-series');
    end
    usefulltsR2 = false;
    usefulltsxR2 = false;
end

%-- prepare other parameters
ecog_APRFF_INITsub_loadoptions;

%-- set default path
SetDefault('prfPth','pRFmodel');

%% main
if ~exist('INITload','var')||isempty(INITload)||any(ismember({'a','alpha','all'},lower(INITload)))
[modeldata_a, prf_params_a, xR2fldname, usefulltsxR2]   = ecog_prf_loadprfs(subjectList,alphaType,prfPth,modeldataID,prfID,average,smoothingMode,smoothingN,prfmodel,gaussianmode,selectchs,selectch_exFEF,selectch_thresh,usefulltsR2,usefulltsxR2,skipsummarizeROIs);
end
if ~exist('INITload','var')||isempty(INITload)||any(ismember({'bb','broadband','all'},lower(INITload)))
[modeldata_bb, prf_params_bb, xR2fldname, usefulltsxR2] = ecog_prf_loadprfs(subjectList,broadbandType,prfPth,modeldataID,prfID,average,smoothingMode,smoothingN,prfmodel,gaussianmode,selectchs,selectch_exFEF,selectch_thresh,usefulltsR2,usefulltsxR2,skipsummarizeROIs);
end
