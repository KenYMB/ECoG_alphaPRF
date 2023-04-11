% visualize pRF locations based on computation in 06c (equivalant 03c2 + 03d2)
%   Update xR2 with full time-series (no decimation)
%   Folk from ecog_APRF_06d_visualizePRFparamsFPM_bootall_fullts.m
%         and ecog_APRF_06e_visualizePRFsFPM_bootall_fullts.m

% 20210107 Yuasa - Update from 06d,06e
% 20210809 Yuasa - modify for paper
% 20210916 Yuasa - for APRFF

%% prefix
% close all; clear all;
% startupToolboxToolbox;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'pRFrelations-representative');
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
else
plotsavePth    = 'pRFrelations-representative';
prfPth         = 'pRFmodel';
prfstatPth     = 'pRFanalysis';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth));
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%-- Plotting Setting
FntSiz  = 14;

%%
%%% load analyzePRF & recompute full-ts R2
clear alphaType broadbandType
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

ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

%% Load correlations

% nboot     = 5000;
boottype  = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling
robustfit = true;

[prfs,prfs_corr,prfs_corrx,rois] = ...
    ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit);

%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations

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

hF = gobjects(0);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot PRF locations
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ecog_APRFF_INITb_mergedata;
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
if ~issaveplot && all(startsWith(subjects,'p'))
subjectsnum = cellfun(@(C) sprintf('P%d',str2double(C{1})),regexp(subjects,'\d+','match'),'UniformOutput',false);
else
subjectsnum = cellfun(@(s) ['S' s],cellstr(num2str([1:nsbj]')),'UniformOutput',false);
end
for isbj=1:nsbj
    prf_all_bb.channels.subject_name(ismember(prf_all_bb.channels.subject_name,subjects{isbj})) = subjectsnum(isbj);
    prf_all_a.channels.subject_name(ismember(prf_all_a.channels.subject_name,subjects{isbj})) = subjectsnum(isbj);
end
prf_all_bb.subjects = subjectsnum;
prf_all_a.subjects  = subjectsnum;


%% prf in visual areas (and in low visual area)

plroi = selroi;

% overlap = 'both';       % include overlapped channels in both ROIs
% overlap = 'none';       % exclude overlapped channels
overlap = 'prob';       % assign overlapped channels based on Wang probability

chidx = [];
for iroi = plroi
    if iscell(prfs.chanidx)
        chidx{iroi} = [prfs.chanidx{iroi,:}];
    else
        prf_loc_bb = [prfs.ecc_bb{iroi}', prfs.ang_bb{iroi}', prfs.rfsize_bb{iroi}'];
        prf_loc_a  = [prfs.ecc_a{iroi}',  prfs.ang_a{iroi}',  prfs.rfsize_a{iroi}'];

        prf_loc_bb = unique(prf_loc_bb,'rows','stable');
        prf_loc_a  = unique(prf_loc_a,'rows','stable');

        chidx{iroi} = nearlyeq(prf_all_bb.ecc,prf_loc_bb(:,1));
        assert(isequal(chidx{iroi}, nearlyeq(prf_all_bb.ang,prf_loc_bb(:,2))),'need to check channles');
        assert(isequal(chidx{iroi}, nearlyeq(prf_all_a.ecc,prf_loc_a(:,1))),'need to check channles');
    end
    chidx{iroi} = unique(chidx{iroi});
    
    switch overlap
        case {'none'}
            chidx{iroi}(wangprob{chidx{iroi},iroi} < sum(wangprob{chidx{iroi},plroi},2)) = [];
        case {'prob'}
            chidx{iroi}(wangprob{chidx{iroi},iroi} < max(wangprob{chidx{iroi},plroi},[],2)) = [];
    end
end
for iroi = plroi
if ~isempty(chidx{iroi})
    
    nplch = numel(chidx{iroi});
    if issaveplot
      if nplch > 21,     iwdth = 9;
      elseif nplch > 18, iwdth = 7;
      elseif nplch > 10, iwdth = 6;
      else,    iwdth = min(5,nplch);
      end
    else
      if nplch > 21,     iwdth = 9;
      elseif nplch > 10, iwdth = 7;
      else,    iwdth = min(5,nplch);
      end
    end
    
    opts = [];
    opts.plot.pix2deg = cfactor;
    opts.plot.XLim    = [-1 1].*12;
    opts.plot.YLim    = [-1 1].*12;
    opts.plot.addChsToTitle     = 'yes';
    opts.plot.addSbjToTitle     = 'yes';
    opts.plot.addBensonToTitle  = 'no';
    opts.plot.addWangToTitle    = 'no';
    opts.plot.fontSize          = FntSiz;
    opts.plot.nSubPlots = [0 iwdth];
%     p = ecog_plotGridPRF(chidx{iroi}, opts, prf_all_bb,prf_all_a);
%     hF = [hF p{:}];
%     
%     set(gcf,'MenuBar','none');
%     
% figureName = sprintf('prf-%02d%%-%02d%%-ecc%02d_location_size-%s_overlap-%s',threshold_bb,threshold_a,eclimit,rois{iroi},overlap);
% savefigauto(gcf,fullfile(plotsavedir, figureName));

    opts.plot.showaxis = 0;         % if show X & Y axis
    p = ecog_plotGridPRF(chidx{iroi}, opts, prf_all_bb,prf_all_a);
    hF = [hF p{:}];
    
    set(gcf,'MenuBar','none');
    
figureName = sprintf('prf-%02d%%-%02d%%-ecc%02d_location_size-%s_overlap-%s_noaxis',threshold_bb,threshold_a,eclimit,rois{iroi},overlap);
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName));    end
end
end
