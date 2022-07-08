function [th_bb,th_a,th_erp] = ecog_prf_getthreshold(prfstatPth,R2mode,selectchs,postfix,alpha_cod)
% [th_bb,th_a,th_erp] = ECOG_PRF_GETTHRESHOLD(prfstatPth,R2mode,selectchs,postfix,alpha_cod)
%    loads permutation data and get threshold at alpha_cod
% 
% Outputs:
%   th_bb       = threshold R2 value for broadband
%   th_a        = threshold R2 value for alpha
%   th_erp      = threshold R2 value for erp
% 
% Inputs:
%   [fine name related]
%   prfstatPth      = Directory where mat files are saved 
%   R2mode          = 
%   selectchs       = 
%   postfix         = 
% 
%   [threshold parameter]
%   alpha_cod       = threshold criterion (default: 0.95)

% 20211124 Yuasa

% Dependency: SetDefault (ky_utils)

%% Load permutation data and get threshold
%-- Set default value of alpha_cod
SetDefault('alpha_cod',0.95);

%-- Get 1st value if selectchs is cell
if iscell(selectchs),   selectchs = selectchs{1};   end

%-- Default path of all_codperm data
SetDefaultAnalysisPath('DATA','pRFanalysis','prfstatPth');
%-- broadband & alpha
filename = sprintf('all_cod-permhalves%s-%s%s.mat',R2mode,selectchs,postfix);
filepath = fullfile(prfstatPth,filename);
if exist(filepath,'file')
    load(filepath,'cod_bb','cod_a');
    th_bb = round(prctile(cod_bb(:),alpha_cod*100));
    th_a  = round(prctile(cod_a(:),alpha_cod*100));
else
    th_bb = nan;
    th_a  = nan;
end
%-- erp
filename = sprintf('all_cod-permhalves%s-%s%s-%s',R2mode,selectchs,postfix,'ERP');
filepath = fullfile(prfstatPth,filename);
if exist(filepath,'file')
    load(filepath,'cod_erp');
    th_erp = round(prctile(cod_erp(:),alpha_cod*100));
else
    th_erp = nan;
end
end