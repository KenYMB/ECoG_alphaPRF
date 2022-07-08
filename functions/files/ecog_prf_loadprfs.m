function [modeldata, prf_params, xR2fldname, usefulltsxR2] = ecog_prf_loadprfs(subjectList,targetBAND,prfPth,modeldataID,prfID,average,smoothingMode,smoothingN,prfmodel,gaussianmode,selectchs,selectch_exFEF,selectch_thresh,usefulltsR2,usefulltsxR2,skipsummarizeROIs)
% [modeldata, prf_params, xR2fldname, usefulltsxR2] = ECOG_PRF_LOADPRFS(subjectList,targetBAND,prfPth,modeldataID,prfID,average,smoothingMode,smoothingN,prfmodel,gaussianmode,selectchs,selectch_exFEF,selectch_thresh,usefulltsR2,usefulltsxR2,skipsummarizeROIs)
%    loads pRF results and rearrange data for further analysis
% 
% Outputs:
%   modeldata       = model timecourse
%   prf_params      = pRF parameters
%   xR2fldname      = fieldname of cross-validated R2 in prf_params
%   usefulltsxR2    = update 'usefulltsxR2' if need
% 
% Inputs:
%   [fine name related]
%   subjectList     = Subjects name list
%   prfPth          = Directory where mat files are saved 
%   modeldataID     = file name id for model timecourse
%   prfID           = file name id for pRF parameters
% 
%   [pRF parameters]
%   targetBAND      =
%   average         =
%   smoothingMode   =
%   smoothingN      =
%   prfmodel        =
%   gaussianmode    =
%   See also: ecog_prf_constructTimeSeries, ecog_prf_analyzePRF
% 
%   [rearrange options]
%   selectchs       =
%   selectch_exFEF  =
%   selectch_thresh =
%   usefulltsR2     =
%   usefulltsxR2    = 
%   See also: ecog_updatePRFdata, ecog_channelSelection
% 
%   [other options]
%   skipsummarizeROIs = true/false, if apply ecog_summarizeROIs or not

% 20211117 Yuasa

%% load analyzePRF
usefllts = usefulltsR2 | usefulltsxR2;
%-- load model time series
opts = [];
opts.fileid         = modeldataID;
opts.average        = average;
if usefllts                     % load full time-series for modeldata
opts.smoothingMode  = 'none';
opts.smoothingN     = [];
else
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;
end
opts.issave         = false;
opts.compute        = false;
opts.outputDir      = prfPth;

opts.targetBAND     = targetBAND;
modeldata   = ecog_prf_constructTimeSeries(subjectList, opts);

%-- load pRF results
opts.fileid         = prfID;
opts.prfmodel       = prfmodel;
opts.gaussianmode   = gaussianmode;
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;

opts.targetBAND     = targetBAND;
prf_params  = ecog_prf_analyzePRF(subjectList, opts);

%% check cross-validated data
if usefulltsxR2 && ~isfield(prf_params{1},'results_xval')
    warning('pRF results do not include cross-validated results');
    usefulltsxR2 = false;
end

%% load new channels for the data computed without wangprob
%-- update channel, R2 & xR2
if min([cellfun(@(x) size(x.channels,2),modeldata);cellfun(@(x) size(x.channels,2),prf_params)])<41
            updatechan = true;
else,       updatechan = false;
end
[modeldata, prf_params, xR2fldname] = ecog_updatePRFdata(modeldata, prf_params, updatechan, usefulltsR2, usefulltsxR2);

%-- summarize ROIs
if ~skipsummarizeROIs
  modeldata   = ecog_summarizeROIs(modeldata);
  prf_params  = ecog_summarizeROIs(prf_params);
end

%% apply channel selection
[modeldata]   = ecog_channelSelection(modeldata,selectchs,selectch_exFEF,selectch_thresh);
[prf_params]  = ecog_channelSelection(prf_params,selectchs,selectch_exFEF,selectch_thresh);

