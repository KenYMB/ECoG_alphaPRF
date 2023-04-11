function [model_all, prf_all, model_va, prf_va, usexvalparams] = ecog_prf_mergeprfs(modeldata,prf_params,va_area,usexvalparams,rearropt)
% [model_all, prf_all, model_va, prf_va, usexvalparams] = ECOG_PRF_MERGEPRFS(modeldata,prf_params,va_area,usexvalparams,rearrange_opt)
%    merge pRF results separated in each subject
% 
% Outputs:
%   model_all       = merged model timecourse
%   prf_all         = merged pRF parameters
%   model_va        = model timecourses in each visual area
%   prf_va          = pRF parameterss in each visual area
%   usexvalparams   = update 'usexvalparams' if need
% 
% Inputs:
%   [pRF results]
%   modeldata       = model timecourse
%   prf_params      = pRF parameters
% 
%   [merge options]
%   va_area         = target visual area name
%   usexvalparams   = if mearge prf_params.results_xval instead of prf_params
%   rearrange_opt   = options
%   See also: ecog_rearrangePRF

% 20211117 Yuasa

%% merge data
SetDefault('rearropt',{},'cell');

%-- merge modeldata
[model_va, model_all]    = ecog_rearrangePRF(ecog_realignTimeSeries(modeldata,1),va_area);

%-- merge prf results
if usexvalparams && ~isfield(prf_params{1},'results_xval')
  warning('pRF results do not include cross-validated results');
  usexvalparams = false;
end
if usexvalparams
  %-- copy structure information to results_xval
  addinfo = setdiff(fieldnames(prf_params{1}),[fieldnames(prf_params{1}.results_xval);{'results_xval'}]);
  prf_params_xval  = cellfun(@(x) copyfields(x,x.results_xval,addinfo),prf_params);
  %-- merge results.results_xval
  [prf_va, prf_all]      = ecog_rearrangePRF(prf_params_xval,va_area,rearropt{:});
else
  %-- merge results
  [prf_va, prf_all]      = ecog_rearrangePRF(prf_params,va_area,rearropt{:});
end
