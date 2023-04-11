%%% explorer channels 
% close all; clear all;

%% Define dataset
close all; clearvars;
%-- Set path
checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('chanPth',       'Channels');
else
chanPth        = 'Channels';
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

%% load full-channels data

channelsC = cell(0);
fldnames = [];
for subject= reshape(string(subjectList),1,[])
    load(fullfile(SetDefaultAnalysisPath('DAT',chanPth),sprintf('%s-channels',subject)),'channels');
    channelsC{end+1,1} = channels;
end
channelsAll = vertcat(channelsC{:});

%-- channel selection
wangprobflds = startsWith(fieldnames(summary(channelsAll)),'wangprob');
wangnormflds = wangprobflds & ~endsWith(fieldnames(summary(channelsAll)),'FEF');
normthresh   = 0.05;
bensonchs    = ~ismember(channelsAll.bensonarea,'none');
wangchs      = ~ismember(channelsAll.wangarea,'none');
wangprobchs  = sum(channelsAll{:,wangprobflds},2)>0;
vischs       = bensonchs | wangchs | wangprobchs;
wangnormchs  = sum(channelsAll{:,wangnormflds},2)>normthresh;

%% Show numbers of channels
%-- All
channels = channelsAll;
fprintf('--------------------\n');
fprintf('%9s:\t%5d\t(good)%5d\n','Total chs',height(channels),sum(ismember(channels.status,'good')));
GPN = 'strip';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
GPN = 'grid';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
GPN = 'HDgrid';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
GPN = 'depth';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
fprintf('%9s:\t%5d\t(w/o %s)%5d\n','bad chs',sum(ismember(channels.status,'bad')),GPN,sum(ismember(channels.status,'bad')&~ismember(channels.group,GPN)));
fprintf('--------------------\n');

%-- Visual
channels = channelsAll(vischs,:);
fprintf('--------------------\n');
fprintf('%9s:\t%5d\t(good)%5d\n','Total chs',height(channels),sum(ismember(channels.status,'good')));
GPN = 'strip';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
GPN = 'grid';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
GPN = 'HDgrid';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
GPN = 'depth';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
fprintf('%9s:\t%5d\t(w/o %s)%5d\n','bad chs',sum(ismember(channels.status,'bad')),GPN,sum(ismember(channels.status,'bad')&~ismember(channels.group,GPN)));
fprintf('--------------------\n');

%-- Wang norm
channels = channelsAll(wangnormchs,:);
fprintf('--------------------\n');
fprintf('%9s:\t%5d\t(good)%5d\n','Total chs',height(channels),sum(ismember(channels.status,'good')));
GPN = 'strip';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
GPN = 'grid';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
GPN = 'HDgrid';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
GPN = 'depth';
fprintf('%9s:\t%5d\t(good)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&ismember(channels.status,'good')));
fprintf('%9s:\t%5d\t(w/o %s)%5d\n','bad chs',sum(ismember(channels.status,'bad')),GPN,sum(ismember(channels.status,'bad')&~ismember(channels.group,GPN)));
fprintf('--------------------\n');


%% %%%%%%%%%%%%%%%%%%
%% Analyzed electrode
%% %%%%%%%%%%%%%%%%%%
%% Load data
%%% Raw data
outputDir  = 'Raw';
[data_pre] = ecog_prf_getData(false, subjectList, [], [], [], [], [], [], outputDir);

%%% selected data
opts = [];
opts.doplots        = false;
opts.outputDir      = 'Preprocessed';
opts.issave         = false;
opts.compute        = false;
[data_post] = ecog_prf_selectData(subjectList,[],opts);

%%% channels
    channels  = cellfun(@(C) C.channels,data_pre,'UniformOutput',false);
    fldnames = [];
    for ii = 1:length(channels)
        %-- Get common field names
        if isempty(fldnames),   fldnames = channels{ii}.Properties.VariableNames;
        else,                   fldnames = intersect(fldnames, channels{ii}.Properties.VariableNames,'stable');
        end
    end
    for ii = 1:length(channels)
        %-- Exclude subject specific field
        channels{ii} = channels{ii}(:,fldnames);
        %-- Rename group
        channels{ii}.group(~ismember(lower(channels{ii}.group),{'strip','grid','hdgrid','depth'}) & ismember(lower(channels{ii}.type),{'ecog'})) = {'strip'};
        channels{ii}.group(~ismember(lower(channels{ii}.group),{'strip','grid','hdgrid','depth'}) & ismember(lower(channels{ii}.type),{'seeg'})) = {'depth'};
    end
channels_pre  = vertcat(channels{:});

    channels  = cellfun(@(C) C.channels,data_post,'UniformOutput',false);
    fldnames = [];
    for ii = 1:length(channels)
        %-- Get common field names
        if isempty(fldnames),   fldnames = channels{ii}.Properties.VariableNames;
        else,                   fldnames = intersect(fldnames, channels{ii}.Properties.VariableNames,'stable');
        end
    end
    for ii = 1:length(channels)
        %-- Exclude subject specific field
        channels{ii} = channels{ii}(:,fldnames);
        %-- Rename group
        channels{ii}.group(~ismember(lower(channels{ii}.group),{'strip','grid','hdgrid','depth'}) & ismember(lower(channels{ii}.type),{'ecog'})) = {'strip'};
        channels{ii}.group(~ismember(lower(channels{ii}.group),{'strip','grid','hdgrid','depth'}) & ismember(lower(channels{ii}.type),{'seeg'})) = {'depth'};
    end
channels_post = vertcat(channels{:});

%% Show numbers of channels
normthresh   = 0.05;

%-- Analyzed electrodes
channels = channels_pre;
wangnormflds = startsWith(fieldnames(summary(channels)),'wangprob') & ~endsWith(fieldnames(summary(channels)),'FEF');
wangnormchs  = sum(channels{:,wangnormflds},2)>normthresh;
fprintf('--------------------\n');
fprintf('%9s:\t%5d\t(Wang)%5d\n','Total chs',height(channels),sum(wangnormchs));
GPN = 'strip';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
GPN = 'grid';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
GPN = 'HDgrid';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
GPN = 'depth';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
fprintf('--------------------\n');

%-- Selected electrodes
channels = channels_post;
wangnormflds = startsWith(fieldnames(summary(channels)),'wangprob') & ~endsWith(fieldnames(summary(channels)),'FEF');
wangnormchs  = sum(channels{:,wangnormflds},2)>normthresh;
fprintf('--------------------\n');
fprintf('%9s:\t%5d\t(Wang)%5d\n','Total chs',height(channels),sum(wangnormchs));
GPN = 'strip';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
GPN = 'grid';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
GPN = 'HDgrid';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
GPN = 'depth';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
fprintf('--------------------\n');


%-- Wang
GPN = {'strip','grid','HDgrid'};
fprintf('--------------------\n');
channels = channelsAll;
wangnormflds = startsWith(fieldnames(summary(channels)),'wangprob') & ~endsWith(fieldnames(summary(channels)),'FEF');
wangnormchs  = sum(channels{:,wangnormflds},2)>normthresh;
fprintf('%9s:\t%5d\t(Strip/Grid)%5d\t(HDgrid)%5d\n','Electrodes',sum((ismember(channels.group,GPN)&wangnormchs)|ismember(channels.group,'HDgrid'))...
                                                       ,sum(ismember(channels.group,GPN)&~ismember(channels.group,'HDgrid')&wangnormchs)...
                                                       ,sum(ismember(channels.group,'HDgrid')));
fprintf('%9s:\t%5d\t(Strip/Grid)%5d\t(HDgrid)%5d\n','(Visual)',sum(ismember(channels.group,GPN)&wangnormchs)...
                                                       ,sum(ismember(channels.group,GPN)&~ismember(channels.group,'HDgrid')&wangnormchs)...
                                                       ,sum(ismember(channels.group,'HDgrid')&wangnormchs));
fprintf('%9s:\t%5d\t(Strip/Grid)%5d\t(HDgrid)%5d\n','Good',sum(((ismember(channels.group,GPN)&wangnormchs)|ismember(channels.group,'HDgrid'))&ismember(channels.status,'good'))...
                                                       ,sum(ismember(channels.group,GPN)&~ismember(channels.group,'HDgrid')&wangnormchs&ismember(channels.status,'good'))...
                                                       ,sum(ismember(channels.group,'HDgrid')&ismember(channels.status,'good')));
fprintf('%9s:\t%5d\t(Strip/Grid)%5d\t(HDgrid)%5d\n','(Visual)',sum(ismember(channels.group,GPN)&wangnormchs&ismember(channels.status,'good'))...
                                                       ,sum(ismember(channels.group,GPN)&~ismember(channels.group,'HDgrid')&wangnormchs&ismember(channels.status,'good'))...
                                                       ,sum(ismember(channels.group,'HDgrid')&wangnormchs&ismember(channels.status,'good')));
                                                   
channels = channels_post;
wangnormflds = startsWith(fieldnames(summary(channels)),'wangprob') & ~endsWith(fieldnames(summary(channels)),'FEF');
wangnormchs  = sum(channels{:,wangnormflds},2)>normthresh;
fprintf('%9s:\t%5d\t(Strip/Grid)%5d\t(HDgrid)%5d\n','Selected',sum((ismember(channels.group,GPN)&wangnormchs)|ismember(channels.group,'HDgrid'))...
                                                       ,sum(ismember(channels.group,GPN)&~ismember(channels.group,'HDgrid')&wangnormchs)...
                                                       ,sum(ismember(channels.group,'HDgrid')));
fprintf('%9s:\t%5d\t(Strip/Grid)%5d\t(HDgrid)%5d\n','(Visual)',sum(ismember(channels.group,GPN)&wangnormchs)...
                                                       ,sum(ismember(channels.group,GPN)&~ismember(channels.group,'HDgrid')&wangnormchs)...
                                                       ,sum(ismember(channels.group,'HDgrid')&wangnormchs));
fprintf('--------------------\n');

%-- HDgrid
GPN = 'HDgrid';
fprintf('--------------------\n');
channels = channelsAll;
wangnormflds = startsWith(fieldnames(summary(channels)),'wangprob') & ~endsWith(fieldnames(summary(channels)),'FEF');
wangnormchs  = sum(channels{:,wangnormflds},2)>normthresh;
fprintf('%9s:\t%5d\t(Wang)%5d\t(Other)%5d\n','HDgrid',sum(ismember(channels.group,GPN))...
                                          ,sum(ismember(channels.group,GPN)&wangnormchs)...
                                          ,sum(ismember(channels.group,GPN)&~wangnormchs));
fprintf('%9s:\t%5d\t(Wang)%5d\t(Other)%5d\n','Good',sum(ismember(channels.group,GPN)&ismember(channels.status,'good'))...
                                       ,sum(ismember(channels.group,GPN)&ismember(channels.status,'good')&wangnormchs)...
                                       ,sum(ismember(channels.group,GPN)&ismember(channels.status,'good')&~wangnormchs));
channels = channels_post;
wangnormflds = startsWith(fieldnames(summary(channels)),'wangprob') & ~endsWith(fieldnames(summary(channels)),'FEF');
wangnormchs  = sum(channels{:,wangnormflds},2)>normthresh;
fprintf('%9s:\t%5d\t(Wang)%5d\t(Other)%5d\n','Selected',sum(ismember(channels.group,GPN))...
                                           ,sum(ismember(channels.group,GPN)&wangnormchs)...
                                           ,sum(ismember(channels.group,GPN)&~wangnormchs));
fprintf('--------------------\n');


%% %%%%%%%%%%%%%%%%%%
%% pRF electrode
%% %%%%%%%%%%%%%%%%%%
%% Load data
%%% pRF data
opts = [];
opts.outputDir      = 'pRFmodel';
opts.issave         = false;
opts.compute        = false;
opts.prfmodel       = 'linear';
opts.gaussianmode   = 'gs';
opts.targetBAND     = 'bbS';
opts.smoothingMode  = 'decimate';
opts.smoothingN     = 3;
opts.average        = 'runs';
[prf] = ecog_prf_analyzePRF(subjectList, opts);

%%% channels
    channels  = cellfun(@(C) C.channels,prf,'UniformOutput',false);
    fldnames = [];
    for ii = 1:length(channels)
        %-- Get common field names
        if isempty(fldnames),   fldnames = channels{ii}.Properties.VariableNames;
        else,                   fldnames = intersect(fldnames, channels{ii}.Properties.VariableNames,'stable');
        end
    end
    for ii = 1:length(channels)
        %-- Exclude subject specific field
        channels{ii} = channels{ii}(:,fldnames);
        %-- Rename group
        channels{ii}.group(~ismember(lower(channels{ii}.group),{'strip','grid','hdgrid','depth'}) & ismember(lower(channels{ii}.type),{'ecog'})) = {'strip'};
        channels{ii}.group(~ismember(lower(channels{ii}.group),{'strip','grid','hdgrid','depth'}) & ismember(lower(channels{ii}.type),{'seeg'})) = {'depth'};
    end
channels_prf  = vertcat(channels{:});

%% Show numbers of channels
normthresh   = 0.05;

%-- pRF electrodes
channels = channels_prf;
wangnormflds = startsWith(fieldnames(summary(channels)),'wangprob') & ~endsWith(fieldnames(summary(channels)),'FEF');
wangnormchs  = sum(channels{:,wangnormflds},2)>normthresh;
fprintf('--------------------\n');
fprintf('%9s:\t%5d\t(Wang)%5d\n','Total chs',height(channels),sum(wangnormchs));
GPN = 'strip';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
GPN = 'grid';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
GPN = 'HDgrid';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
GPN = 'depth';
fprintf('%9s:\t%5d\t(Wang)%5d\n',GPN,sum(ismember(channels.group,GPN)),sum(ismember(channels.group,GPN)&wangnormchs));
fprintf('--------------------\n');


%% %%%%%%%%%%%%%%%%%%
%% pRF electrode with threshold
%% %%%%%%%%%%%%%%%%%%
clear alphaType broadbandType

average        ='runs';
prfmodel       ='linear';
smoothingMode  ='decimate';
smoothingN     = 3;
gaussianmode   ='gs';
selectchs      = 'wangprobchs';
% selectchs      = 'HDgridchs';
% selectchs      = {'wangprobchs','HDgridchs'};
    allowlag       = false;
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;

va_area = 'wangarea';
usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;

%%% Load correlations
boottype  = 'boot';         % 'boot': w/ resampling; 'iter': w/o resampling
robustfit = true;

[prfs,prfs_corr,prfs_corrx,rois,nroi] = ...
    ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit);


%%
elec_ok = ~(prf_all_bb.xval<=threshold_bb | prf_all_a.xval<=threshold_a ...
     | prf_all_bb.ecc >= eclimit | prf_all_a.ecc >= eclimit);  % use tilde for nan
elec_bb = ~(prf_all_bb.xval<=threshold_bb | prf_all_bb.ecc >= eclimit);  % use tilde for nan

wang_prob    = startsWith(fieldnames(summary(channels)),'wangprob_');
wang_early   = wang_prob & endsWith(fieldnames(summary(channels)),{'V1','V2','V3'});
wang_dorsal  = wang_prob & endsWith(fieldnames(summary(channels)),{'V3a','V3b','LO1','LO2','TO','IPS'});
wang_ventral = wang_prob & endsWith(fieldnames(summary(channels)),{'hV4','VO','PHC'});
wang_other   = wang_prob & ~(wang_early | wang_dorsal | wang_ventral);

%-- Subject Table
fprintf('\n=============================================================\n');
for subject = string(subjectList')

    elec_sbj = ismember(channels.subject_name,subject);
    fprintf('%s\t# Vis:%3d\tArea:%-32s\t# Sel:%3d\tArea:%s\n',subject,...
        sum(elec_sbj),sprintf(' %s',unique(channels.wangarea(elec_sbj))),...
        sum(elec_sbj&elec_ok),sprintf(' %s',unique(channels.wangarea(elec_sbj&elec_ok))));
    
end


%-- Analyzed channels (all pRF)
fprintf('\n========================= All pRFs =========================\n');
%%-- Subject
fprintf('%3d electrodes from %d subjects:%s\n',length(elec_ok),numel(unique(channels.subject_name(:))),...
                                              sprintf(' %s',unique(channels.subject_name(:))));
for subject = unique(channels.subject_name(:))'
    fprintf('\t%3d electrodes from %s\n',sum(ismember(channels.subject_name,subject)),subject);
end
%%-- Probabilistic ROIs
Label = 'V1-V3'; chanidx = sum(channels{:,wang_early},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));
Label = 'Dorsolateral'; chanidx = sum(channels{:,wang_dorsal},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));
Label = 'Either'; chanidx = sum(channels{:,wang_early},2)>0|sum(channels{:,wang_dorsal},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));
Label = 'Both'; chanidx = sum(channels{:,wang_early},2)>0&sum(channels{:,wang_dorsal},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
Label = 'Ventral'; chanidx = sum(channels{:,wang_ventral},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
Label = 'SPL'; chanidx = sum(channels{:,wang_other},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
Label = 'Either others'; chanidx = sum(channels{:,wang_ventral|wang_other},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));
Label = 'Total'; chanidx = sum(channels{:,wang_prob},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));


%-- Analyzed channels (broadband Threshold)
fprintf('\n======================= BB Threshold =======================\n');
%%-- Subject
fprintf('%3d electrodes from %d subjects:%s\n',sum(elec_bb),numel(unique(channels.subject_name(elec_bb))),...
                                              sprintf(' %s',unique(channels.subject_name(elec_bb))));
for subject = unique(channels.subject_name(elec_bb))'
    fprintf('\t%3d electrodes from %s\n',sum(elec_bb&ismember(channels.subject_name,subject)),subject);
end
%%-- Probabilistic ROIs
Label = 'V1-V3'; chanidx = elec_bb&sum(channels{:,wang_early},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));
Label = 'Dorsolateral'; chanidx = elec_bb&sum(channels{:,wang_dorsal},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));
Label = 'Either'; chanidx = elec_bb&(sum(channels{:,wang_early},2)>0|sum(channels{:,wang_dorsal},2)>0);
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));
Label = 'Both'; chanidx = elec_bb&sum(channels{:,wang_early},2)>0&sum(channels{:,wang_dorsal},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
Label = 'Ventral'; chanidx = elec_bb&sum(channels{:,wang_ventral},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));


%-- Analyzed channels (pRF Threshold)
fprintf('\n======================= pRF Threshold =======================\n');
%%-- Subject
fprintf('%3d electrodes from %d subjects:%s\n',sum(elec_ok),numel(unique(channels.subject_name(elec_ok))),...
                                              sprintf(' %s',unique(channels.subject_name(elec_ok))));
for subject = unique(channels.subject_name(elec_ok))'
    fprintf('\t%3d electrodes from %s\n',sum(elec_ok&ismember(channels.subject_name,subject)),subject);
end
%%-- Probabilistic ROIs
Label = 'V1-V3'; chanidx = elec_ok&sum(channels{:,wang_early},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));
Label = 'Dorsolateral'; chanidx = elec_ok&sum(channels{:,wang_dorsal},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));
Label = 'Either'; chanidx = elec_ok&(sum(channels{:,wang_early},2)>0|sum(channels{:,wang_dorsal},2)>0);
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
    % summary(channels.subject_name(chanidx));
Label = 'Both'; chanidx = elec_ok&sum(channels{:,wang_early},2)>0&sum(channels{:,wang_dorsal},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
Label = 'Ventral'; chanidx = elec_ok&sum(channels{:,wang_ventral},2)>0;
fprintf('%3d electrodes in %-12s from%s\n',sum(chanidx),Label,sprintf(' %s',unique(channels.subject_name(chanidx))));
fprintf('Averaged %.2f electrodes in %s\n',nanmean(prfs.num(end-1,:)),'V1-V3');
fprintf('Averaged %.2f electrodes in %s\n',nanmean(prfs.num(end,:)),'Dorsolateral');

