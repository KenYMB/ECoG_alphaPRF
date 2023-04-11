% ecog_APRFF_INITb_mergedata
% General script to merge pRF data across subjects
%   Need to specify parameters listed in "Parameters"

% Dependency: copyfields (fieldTrip)

% 20201008 Yuasa
% 20211117 Yuasa - separate main processing into ecog_prf_mergeprfs
% 20221203 Yuasa - enable to load alpha/broabdan separately

%% Parameters
% va_area = 'wangarea';
% 
% % usexvalparams = false;

% Inputs: 
%   modeldata_a
%   modeldata_bb
%   prf_params_a
%   prf_params_bb
% 
%   va_area = 'wangarea';
% 
% %   usexvalparams = false;
% 
% Outputs:
%   model_all_a
%   model_all_bb
%   prf_all_a
%   prf_all_bb
% 
%   channels
%   nchan

SetDefault('usexvalparams',false);

%% estimate va_area
if exist('selectchs','var')
    if ismember(selectchs,{'bensonchs'})
         SetDefault('va_area','bensonarea');
    elseif ismember(selectchs,{'wangchs','wangprobchs'})
         SetDefault('va_area','wangarea');
    elseif ismember(selectchs,{'hcpchs','hcpprobchs'})
         SetDefault('va_area','hcparea');
    else
         SetDefault('va_area','wangarea');
    end
    if ( (ismember({'wangarea'},va_area) && ismember({'wangprobchs'},selectchs)) || ...
         (ismember({'hcparea'},va_area) && ismember({'wangprobchs'},selectchs)) ) && ...
            exist('selectch_thresh','var') && ~isempty(selectch_thresh)
        rearropt = {'norm',selectch_thresh};     % reassign visual area labels
    else
        rearropt = {};
    end
else
    assert(exist('va_area','var')&&~isempty(va_area),'Variable ''va_area'' is required.');
    SetDefault('rearropt',{},'cell');
end

%% main
%-- merge data
if ~exist('INITload','var')||isempty(INITload)||any(ismember({'a','alpha','all'},lower(INITload)))
[model_all_a, prf_all_a, model_va_a, prf_va_a, usexvalparams]     = ecog_prf_mergeprfs(modeldata_a,prf_params_a,va_area,usexvalparams,rearropt);
end
if ~exist('INITload','var')||isempty(INITload)||any(ismember({'bb','broadband','all'},lower(INITload)))
[model_all_bb, prf_all_bb, model_va_bb, prf_va_bb, usexvalparams] = ecog_prf_mergeprfs(modeldata_bb,prf_params_bb,va_area,usexvalparams,rearropt);
end

%-- set useful variables
if ~exist('prf_all_bb','var')
channels = prf_all_a.channels;
else
channels = prf_all_bb.channels;
end
nchan    = height(channels);

if exist('prf_all_bb','var')&&exist('prf_all_a','var')
prf_table = ecog_prf_prftable(channels,va_area,prf_all_bb,'bb',prf_all_a,'a');
end
