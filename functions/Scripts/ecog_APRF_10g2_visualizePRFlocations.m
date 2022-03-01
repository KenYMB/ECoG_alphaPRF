% visualize pRF correlations computed in 06c (equivalant 03c2 + 03d2)
%   Update xR2 with full time-series (no decimation)
%   Folk from ecog_APRF_06d_visualizePRFparamsFPM_bootall_fullts.m
%         and ecog_APRF_06e_visualizePRFsFPM_bootall_fullts.m

% 20210107 Yuasa - Update from 06d,06e
% 20210809 Yuasa - modify for paper
% 20210916 Yuasa - for APRFF

%% Define paths and dataset
% close all; clear all;
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
SetDefault('issaveplot',true);
if issaveplot
    plotsavePth    = fullfile(figPth, 'pRFrelations-representative');
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%% load analyzePRF & recompute full-ts R2
average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';
selectchs      = 'wangprobchs';
    allowlag       = false;
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;

va_area = 'wangarea';
usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

clear alphaType broadbandType
ecog_APRF_INITa_loaddata;
ecog_APRF_INITc_postfix;
ecog_APRF_INITd_threshold;

%% Load correlations

% nboot     = 5000;
boottype  = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling

filename = sprintf('all_prfparams-%sall%s-%s%s_thresh%02d',boottype,R2mode,selectchs,postfix,threshold_bb);
if exist('threshold_a','var')&&~isnan(threshold_a)
    filename = sprintf('%s-%02d',filename,threshold_a); end
if exist('eclimit','var')&&~isnan(eclimit)
    filename = sprintf('%s-ecc%02d',filename,eclimit);     end
filepath = fullfile(savePth,'pRFanalysis',filename);

if exist([filepath '.mat'],'file')
 fprintf('loading %s...',filename);
 load(filepath);
 fprintf('\n');
 nroi = length(rois);
 nboot = size(prfs.num,2);
 SetDefault('threshold_a',nan); SetDefault('eclimit',nan);
else
 error('Failed to find %s',filename);
end

%-- update ROI labels
rois(ismember(rois,'low'))={'V1-V3'};
rois(ismember(rois,'high'))={'Dorsolateral'};
rois(ismember(rois,'dorsolateral'))={'Dorsolateral'};

%-- reject inaccurate ROIs
rejthr = nboot .* 2;
for ifld = fieldnames(prfs)'
    for iroi = 1:nroi
        if ~ismember(ifld{:},{'num','chanidx'}) && length(prfs.(ifld{:}){iroi})<rejthr
        prfs.(ifld{:}){iroi} = nan;
        end
    end
end

%%% prepare visualize
nsubjects = length(prf_params_bb);
nroi  = length(rois);
selroi = (nroi-1):nroi;
% nvisarea  = length(prf_va_bb);
% areaList  = unique(prf_all_bb.channels.(va_area));

res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');

% plotsavedir    = fullfile(plotsavePth, postfix(2:end),'arrange');
% if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end


%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations
% close all;
fontSiz = 18;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot PRF locations
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ecog_APRF_INITb_mergedata;
%% Wang prob
wangprob = channels(:,startsWith(channels.Properties.VariableNames,'wangprob')&endsWith(channels.Properties.VariableNames,rois));
%-- low
catrois = {'V1';'V2';'V3'};
catidx  = endsWith(wangprob.Properties.VariableNames,catrois);
wangprob.wangprob_low = sum(wangprob{:,catidx},2);
%-- high
catrois = {'V3a';'V3b';'LO1';'LO2';'TO';'IPS'};
catidx  = endsWith(wangprob.Properties.VariableNames,catrois);
wangprob.wangprob_high = sum(wangprob{:,catidx},2);

%-- update Subjects
subjects = prf_all_bb.subjects;
nsbj = length(subjects);
subjectsnum = cellfun(@(s) ['S' s],cellstr(num2str([1:nsbj]')),'UniformOutput',false);
for isbj=1:nsbj
    prf_all_bb.channels.subject_name(ismember(prf_all_bb.channels.subject_name,subjects{isbj})) = subjectsnum(isbj);
    prf_all_a.channels.subject_name(ismember(prf_all_a.channels.subject_name,subjects{isbj})) = subjectsnum(isbj);
end
prf_all_bb.subjects = subjectsnum;
prf_all_a.subjects  = subjectsnum;


%% prf in visual areas (and in low visual area)
% close all

hF = {};
for ii=1:1
if ii==1
plroi = selroi;
else
plroi = 1:3;
end

% overlap = 'both';       % include overlapped channels in both ROIs
% overlap = 'none';       % exclude overlapped channels
overlap = 'prob';       % assign overlapped channels based on Wang probability

chidx = [];
for iroi = plroi
    prf_loc_bb = [prfs.ecc_bb{iroi}', prfs.ang_bb{iroi}', prfs.rfsize_bb{iroi}'];
    prf_loc_a  = [prfs.ecc_a{iroi}',  prfs.ang_a{iroi}',  prfs.rfsize_a{iroi}'];

    prf_loc_bb = unique(prf_loc_bb,'rows','stable');
    prf_loc_a  = unique(prf_loc_a,'rows','stable');
    
    chidx{iroi} = nearlyeq(prf_all_bb.ecc,prf_loc_bb(:,1));
    assert(isequal(chidx{iroi}, nearlyeq(prf_all_bb.ang,prf_loc_bb(:,2))),'need to check channles');
    assert(isequal(chidx{iroi}, nearlyeq(prf_all_a.ecc,prf_loc_a(:,1))),'need to check channles');
    chidx{iroi} = unique(chidx{iroi});
    
    switch overlap
        case {'none'}
            chidx{iroi}(wangprob{chidx{iroi},iroi} < sum(wangprob{chidx{iroi},plroi},2)) = [];
        case {'prob'}
            chidx{iroi}(wangprob{chidx{iroi},iroi} < max(wangprob{chidx{iroi},plroi},[],2)) = [];
    end
end
for iroi = plroi
    if     iroi == 12, iwdth = 7;
    elseif iroi == 13, iwdth = 9;
    else,              iwdth = 5;
    end
if ~isempty(chidx{iroi})
    
    opts = [];
    opts.plot.pix2deg = cfactor;
    opts.plot.XLim    = [-1 1].*12;
    opts.plot.YLim    = [-1 1].*12;
    opts.plot.addChsToTitle     = 'yes';
    opts.plot.addSbjToTitle     = 'yes';
    opts.plot.addBensonToTitle  = 'no';
    opts.plot.addWangToTitle    = 'no';
    opts.plot.fontSize          = 14;
    opts.plot.nSubPlots = [0 min(iwdth,numel(chidx{iroi}))];
%     ecog_plotGridPRF(chidx{iroi}, opts, prf_all_bb,prf_all_a);
%     
%     set(gcf,'MenuBar','none');
%     
% figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_location_size-%s_overlap-%s',threshold_bb,threshold_a,eclimit,rois{iroi},overlap);
% savefigauto(gcf,fullfile(plotsavedir, figname));

    opts.plot.showaxis = 0;         % if show X & Y axis
    hF = horzcat(hF,ecog_plotGridPRF(chidx{iroi}, opts, prf_all_bb,prf_all_a));
    
    set(gcf,'MenuBar','none');
    
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_location_size-%s_overlap-%s_noaxis',threshold_bb,threshold_a,eclimit,rois{iroi},overlap);
% savefigauto(gcf,fullfile(plotsavedir,sprintf('width%d',iwdth), figname));
end
end
end
