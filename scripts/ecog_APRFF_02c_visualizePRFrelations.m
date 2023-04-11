% visualize pRF correlations computed in 06c (equivalant 03c2 + 03d2)
%   Update xR2 with full time-series (no decimation)
%   Folk from ecog_APRF_06d_visualizePRFparamsFPM_bootall_fullts.m
%         and ecog_APRF_06e_visualizePRFsFPM_bootall_fullts.m

% 20210107 Yuasa - Update from 06d,06e
% 20210809 Yuasa - modify for paper
% 20210916 Yuasa - Update for new environment

%% prefix
% close all; clear all;
% startupToolboxToolbox;
checkPath;
%-- Input & Output path
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'pRFrelations');
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
else
plotsavePth    = 'pRFrelations';
prfPth         = 'pRFmodel';
prfstatPth     = 'pRFanalysis';
end
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

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

[prfs,prfs_corr,prfs_corrx,rois,nroi,nboot] = ...
    ecog_prf_loadprfrelations(prfstatPth,R2mode,selectchs,postfix,boottype,threshold_bb,threshold_a,eclimit,robustfit);


%%% prepare visualize
% nsubjects = length(prf_params_bb);
% nroi  = length(rois);
% selroi = (nroi-1):nroi;
% nvisarea  = length(prf_va_bb);
% areaList  = unique(prf_all_bb.channels.(va_area));

res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

alpha = 0.32;   % 0.32 for 1sd, 0.05 for 2sd, 0.003 for 3sd

plcol = get(groot,'defaultAxesColorOrder');

plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth), postfix(2:end));
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end


%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations
close all;
FntSiz   = 18;
LFntSiz  = 22;

fitting = '-ellipse';

showref = true;
domean  = true;

%% plot 2D histogram (replication around 90 deg for angle)
%% %%%%%%%%%%%%%%%%%%%%
%%% arrange in one figure
relt    = 'bb_a';
flds    = {'ang','ecc','rfsize'};
selroi  = {'V1-V3','Dorsolateral'};
opts = [];
opts.mothod      = 'hist';
opts.ang_limpad  = 90;       % degree for angle
opts.pix2deg     = cfactor;
opts.FontSize    = FntSiz;
opts.showfitted  = strtok(fitting,'-');
opts.fit_alpha   = alpha;
ecog_prf_plotPRFrelations(prfs,prfs_corr,relt,flds,selroi,rois,opts);

figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_2Dhist%s_%s',threshold_bb,threshold_a,eclimit,fitting,'all');
savefigauto(gcf,fullfile(plotsavedir, figname));

%% %%%%%%%%%%%%%%%%%%%%
%%% cross-correlations

%% plot fitted line in ROIs
relt    = 'ecc_rfsize';
flds    = {'bb','a'};
selroi  = {'V1','V2','V3','V1-V3','Dorsolateral'};
opts = [];
opts.mothod        = 'line';
opts.pix2deg       = cfactor;
opts.FontSize      = LFntSiz;
opts.fit_alpha     = alpha;
opts.fitlm_domean  = domean;
opts.fitlm_showref = showref;
ecog_prf_plotPRFrelations(prfs,prfs_corrx,relt,flds,selroi,rois,opts);

if showref,     isrefname = '_withref';
else,           isrefname = '';
end
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_line-roi_%s%s-largeFont',threshold_bb,threshold_a,eclimit,relt,isrefname);
if domean,  figname = sprintf('%s-mean',figname);   end
savefigauto(gcf,fullfile(plotsavedir, figname));

%% plot xcorr in 2D histogram
relt    = 'ecc_rfsize';
flds    = {'bb','a'};
selroi  = {'V1','V2','V3','V1-V3','Dorsolateral'};
opts = [];
opts.mothod        = 'histogram';
opts.pix2deg       = cfactor;
opts.FontSize      = FntSiz;
opts.showdiag      = true;
opts.showfitted    = strtok(fitting,'-');
opts.fit_alpha     = alpha;
opts.fitlm_showref = showref;
ecog_prf_plotPRFrelations(prfs,prfs_corrx,relt,flds,selroi,rois,opts);

if showref,     isrefname = '_withref';
else,           isrefname = '';
end
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_2Dhist%s_%s%s',threshold_bb,threshold_a,eclimit,fitting,relt,isrefname);
savefigauto(gcf,fullfile(plotsavedir, figname));

%%
close all;


%%
%% ColorBar for 2Dhist
meshres = 300;
maxrad  = meshres/2;

figure;

[X,Y] = meshgrid(0:meshres,0:meshres/10);
Z = X;
surf(X,Y,Z,'EdgeColor','none');
axis equal;
colormap([[0,0,0];hot].^0.3);
set(gcf, 'InvertHardcopy', 'off', 'Color', [1 1 1], 'Position', [200 100 50 250]);
view([90 0.1]); box on;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);

figname = 'colorbar';
savefigauto(gcf, fullfile(plotsavedir,figname),{'png','eps'});




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
subjects = subjectList;
nsbj = length(subjects);
subjectsnum = cellfun(@(s) ['S' s],cellstr(num2str([1:nsbj]')),'UniformOutput',false);
for isbj=1:nsbj
    prf_all_bb.channels.subject_name(ismember(prf_all_bb.channels.subject_name,subjects{isbj})) = subjectsnum(isbj);
    prf_all_a.channels.subject_name(ismember(prf_all_a.channels.subject_name,subjects{isbj})) = subjectsnum(isbj);
end
prf_all_bb.subjects = subjectsnum;
prf_all_a.subjects  = subjectsnum;

selroi = [12:13];       % low, high
% selroi = [1:3];         % V1,V2,V3

%% prf center location (or polar angle relation)
close all

flds = {'cntr','ang'};

for sbfld = flds
[ix,iy] = arrangeinrect(length(selroi));
figure('Position',[1000/2 900/2 900 440],'MenuBar','none');
ii = 1;
for iroi = selroi
    prf_loc_bb = [prfs.ecc_bb{iroi}', prfs.ang_bb{iroi}', prfs.rfsize_bb{iroi}'];
    prf_loc_a  = [prfs.ecc_a{iroi}',  prfs.ang_a{iroi}',  prfs.rfsize_a{iroi}'];

    prf_loc_bb = unique(prf_loc_bb,'rows','stable');
    prf_loc_a  = unique(prf_loc_a,'rows','stable');
    
    switch sbfld{:}
        case 'ang'
            prf_loc_bb(:,1) = 30;
            prf_loc_a(:,1)  = 50;
    end

    subplot_er(iy,ix,ii);
    switch sbfld{:}
        case 'cntr'
            rectangle('Position',[-0.5 -0.5 1 1].*cfactor*100,'Curvature',[1 1],'LineWidth',1,'LineStyle',':');
        case 'ang'
            rectangle('Position',[-0.5 -0.5 1 1].*cfactor*prf_loc_bb(1)*2,'Curvature',[1 1],'LineWidth',1,'LineStyle',':');
            hold on;
            rectangle('Position',[-0.5 -0.5 1 1].*cfactor*prf_loc_a(1)*2,'Curvature',[1 1],'LineWidth',1,'LineStyle',':');
            xticks([]); yticks([]);
    end
    axis([-1 1 -1 1 -1 1].*10, 'equal');
    hold on;
    plot([0 0],ylim,'k--','LineWidth',1);
    plot(xlim,[0 0],'k--','LineWidth',1);
    
    prf_loc_bb = prf_loc_bb(:,1)*cfactor.*[cos(deg2rad(prf_loc_bb(:,2))),...
                                    sin(deg2rad(prf_loc_bb(:,2)))];
    prf_loc_a  = prf_loc_a(:,1)*cfactor.*[cos(deg2rad(prf_loc_a(:,2))),...
                                    sin(deg2rad(prf_loc_a(:,2)))];
                                
    for jj = 1:size(prf_loc_bb,1)
        plot([prf_loc_bb(jj,1),prf_loc_a(jj,1)],[prf_loc_bb(jj,2),prf_loc_a(jj,2)],'k-','LineWidth',1);
    end
        h = [];
        h(1)=scatter(prf_loc_bb(:,1),prf_loc_bb(:,2),60,'o','LineWidth',1.0,'MarkerEdgeColor',plcol(1,:),'MarkerFaceColor','w','MarkerFaceAlpha',1);
        h(2)=scatter(prf_loc_a(:,1),prf_loc_a(:,2),60,'s','LineWidth',1.0,'MarkerEdgeColor',plcol(2,:),'MarkerFaceColor','w','MarkerFaceAlpha',1);
    
    title(rois{iroi});
%     legend(h,{'broadband','alpha'});
    set(gca,'FontSize',FntSiz);
    
    ii = ii+1;
end
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_location_%s',threshold_bb,threshold_a,eclimit,sbfld{:});
savefigauto(gcf,fullfile(plotsavedir, figname));
end

%% prf center location (or polar angle relation) in low visual area
close all

flds = {'cntr','ang'};

plroi = 1:3;

for sbfld = flds
[ix,iy] = arrangeinrect(length(plroi),3);
figure('Position',[1000/4 900/2 1400 440],'MenuBar','none');
ii = 1;
for iroi = plroi
    prf_loc_bb = [prfs.ecc_bb{iroi}', prfs.ang_bb{iroi}', prfs.rfsize_bb{iroi}'];
    prf_loc_a  = [prfs.ecc_a{iroi}',  prfs.ang_a{iroi}',  prfs.rfsize_a{iroi}'];

    prf_loc_bb = unique(prf_loc_bb,'rows','stable');
    prf_loc_a  = unique(prf_loc_a,'rows','stable');
    
    switch sbfld{:}
        case 'ang'
            prf_loc_bb(:,1) = 30;
            prf_loc_a(:,1)  = 50;
    end

    subplot_er(iy,ix,ii);
    switch sbfld{:}
        case 'cntr'
            rectangle('Position',[-0.5 -0.5 1 1].*cfactor*100,'Curvature',[1 1],'LineWidth',1,'LineStyle',':');
        case 'ang'
            rectangle('Position',[-0.5 -0.5 1 1].*cfactor*prf_loc_bb(1)*2,'Curvature',[1 1],'LineWidth',1,'LineStyle',':');
            hold on;
            rectangle('Position',[-0.5 -0.5 1 1].*cfactor*prf_loc_a(1)*2,'Curvature',[1 1],'LineWidth',1,'LineStyle',':');
            xticks([]); yticks([]);
    end
    axis([-1 1 -1 1 -1 1].*10, 'equal');
    hold on;
    plot([0 0],ylim,'k--','LineWidth',1);
    plot(xlim,[0 0],'k--','LineWidth',1);
    
    prf_loc_bb = prf_loc_bb(:,1)*cfactor.*[cos(deg2rad(prf_loc_bb(:,2))),...
                                    sin(deg2rad(prf_loc_bb(:,2)))];
    prf_loc_a  = prf_loc_a(:,1)*cfactor.*[cos(deg2rad(prf_loc_a(:,2))),...
                                    sin(deg2rad(prf_loc_a(:,2)))];
                                
    for jj = 1:size(prf_loc_bb,1)
        plot([prf_loc_bb(jj,1),prf_loc_a(jj,1)],[prf_loc_bb(jj,2),prf_loc_a(jj,2)],'k-','LineWidth',1);
    end
        h = [];
        h(1)=scatter(prf_loc_bb(:,1),prf_loc_bb(:,2),60,'o','LineWidth',1.0,'MarkerEdgeColor',plcol(1,:),'MarkerFaceColor','w','MarkerFaceAlpha',1);
        h(2)=scatter(prf_loc_a(:,1),prf_loc_a(:,2),60,'s','LineWidth',1.0,'MarkerEdgeColor',plcol(2,:),'MarkerFaceColor','w','MarkerFaceAlpha',1);
    
    title(rois{iroi});
%     legend(h,{'broadband','alpha'});
    set(gca,'FontSize',FntSiz);
    
    ii = ii+1;
end
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_location_%s-low',threshold_bb,threshold_a,eclimit,sbfld{:});
savefigauto(gcf,fullfile(plotsavedir, figname));
end


%% prf in visual areas (and in low visual area)
close all

for ii=1:2
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
    opts.plot.nSubPlots = [0 min(9,numel(chidx{iroi}))];
%     ecog_plotGridPRF(chidx{iroi}, opts, prf_all_bb,prf_all_a);
%     
%     set(gcf,'MenuBar','none');
%     
% figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_location_size-%s_overlap-%s',threshold_bb,threshold_a,eclimit,rois{iroi},overlap);
% savefigauto(gcf,fullfile(plotsavedir, figname));

    opts.plot.showaxis = 0;         % if show X & Y axis
    ecog_plotGridPRF(chidx{iroi}, opts, prf_all_bb,prf_all_a);
    
    set(gcf,'MenuBar','none');
    
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_location_size-%s_overlap-%s_noaxis',threshold_bb,threshold_a,eclimit,rois{iroi},overlap);
savefigauto(gcf,fullfile(plotsavedir, figname));
end
end
end

%% prf time series in visual areas
close all

stimulus = model_all_bb.stimulus;

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
if ~isempty(chidx{iroi})
    
    opts = [];
    opts.plot.addChsToTitle     = 'yes';
    opts.plot.addSbjToTitle     = 'yes';
    opts.plot.addBensonToTitle  = 'no';
    opts.plot.addWangToTitle    = 'no';
    opts.plot.fontSize          = 14;
    opts.plot.nSubPlots = [0 min(5,numel(chidx{iroi}))];
    
    %-- broadband
    opts.plot.model_options    = {'Color',plcol(1,:)};
%     opts.skipprojection        = 'no';
%     ecog_plotGridPRFts(model_all_bb.datats,stimulus,prf_all_bb,chidx{iroi}, opts);
%     set(gcf,'MenuBar','none');
%     
% figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_timeseries-%s-broadband_overlap-%s',threshold_bb,threshold_a,eclimit,rois{iroi},overlap);
% savefigauto(gcf,fullfile(plotsavedir, figname));

    %---- no polynomial projection
    opts.skipprojection        = 'yes';
    ecog_plotGridPRFts(model_all_bb.datats,stimulus,prf_all_bb,chidx{iroi}, opts);
    set(gcf,'MenuBar','none');
    
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_timeseries-noproj-%s-broadband_overlap-%s',threshold_bb,threshold_a,eclimit,rois{iroi},overlap);
savefigauto(gcf,fullfile(plotsavedir, figname));


    %-- alpha
    opts.plot.model_options    = {'Color',plcol(2,:)};
%     opts.skipprojection        = 'no';
%     ecog_plotGridPRFts(model_all_a.datats,stimulus,prf_all_a,chidx{iroi}, opts);
%     set(gcf,'MenuBar','none');
%     
% figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_timeseries-%s-alpha_overlap-%s',threshold_bb,threshold_a,eclimit,rois{iroi},overlap);
% savefigauto(gcf,fullfile(plotsavedir, figname));

    %---- no polynomial projection
    opts.skipprojection        = 'yes';
    ecog_plotGridPRFts(model_all_a.datats,stimulus,prf_all_a,chidx{iroi}, opts);
    set(gcf,'MenuBar','none');
    
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_timeseries-noproj-%s-alpha_overlap-%s',threshold_bb,threshold_a,eclimit,rois{iroi},overlap);
savefigauto(gcf,fullfile(plotsavedir, figname));
end
end
end

%%
close all;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Percent of inclusion
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-- pickup each pRF
plroi = selroi;
% overlap = 'both';       % include overlapped channels in both ROIs
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
%%
%-- Percent of inclusion (average across electrode with probability asignment)
bbINa = [];
for iroi = selroi
    [~,bbINa{iroi}] = circintarea(...
                          reshape(prf_all_bb.params(1,1,chidx{iroi}),[],1),...
                          reshape(prf_all_bb.params(1,2,chidx{iroi}),[],1),...
                          reshape(prf_all_bb.params(1,3,chidx{iroi}),[],1),...
                          reshape(prf_all_a.params(1,1,chidx{iroi}),[],1),...
                          reshape(prf_all_a.params(1,2,chidx{iroi}),[],1),...
                          reshape(prf_all_a.params(1,3,chidx{iroi}),[],1));
end

% cellfun(@mean,bbINa)
% cellfun(@median,bbINa)
                     
%-- Percent of inclusion (average across ramdomize assigned electrode)
bbINa_boot = [];
for iroi = selroi
    [~,bbINa_boot{iroi}] = circintarea(...
                          reshape(prfs.ecc_bb{iroi} .* cos(prfs.ang_bb{iroi}*pi/180),[],1),...
                          reshape(prfs.ecc_bb{iroi} .* sin(prfs.ang_bb{iroi}*pi/180),[],1),...
                          reshape(prfs.rfsize_bb{iroi},[],1),...
                          reshape(prfs.ecc_a{iroi} .* cos(prfs.ang_a{iroi}*pi/180),[],1),...
                          reshape(prfs.ecc_a{iroi} .* sin(prfs.ang_a{iroi}*pi/180),[],1),...
                          reshape(prfs.rfsize_a{iroi},[],1));
end

% cellfun(@mean,bbINa_boot)
% cellfun(@median,bbINa_boot)

%% Permutation

%-- Percent of inclusion (average across electrode with probability asignment)
figure('MenuBar','none');
set(gcf,'Position',get(gcf,'Position').*[.8 .8 1.5 1]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for iroi = selroi
perm_bbINa = [];
for iboot = 1:nboot
permidx = Shuffle(chidx{iroi});
    [~,bbINai] = circintarea(...
                          reshape(prf_all_bb.params(1,1,permidx),[],1),...
                          reshape(prf_all_bb.params(1,2,permidx),[],1),...
                          reshape(prf_all_bb.params(1,3,permidx),[],1),...
                          reshape(prf_all_a.params(1,1,chidx{iroi}),[],1),...
                          reshape(prf_all_a.params(1,2,chidx{iroi}),[],1),...
                          reshape(prf_all_a.params(1,3,chidx{iroi}),[],1));
    perm_bbINa(iboot) = mean(bbINai);
end 
mean_bbINa = mean(bbINa{iroi},'omitnan');

nexttile;
histogram(perm_bbINa.*100,'BinWidth',5);
set(gca,'FontSize',18);
hold on;
plot(median(perm_bbINa,'omitnan').*[1 1].*100,ylim,'k-.','LineWidth',1.3);
plot(prctile(perm_bbINa,95).*[1 1].*100,ylim,'k:','LineWidth',1.6);
plot(mean_bbINa.*[1 1].*100,ylim,'-','LineWidth',1.6,'Color',plcol(1,:));
hold off;
title(rois{iroi});
xlim([0 100]);
if mean_bbINa>0.8
    text(mean_bbINa.*100-diff(xlim)*.03,diff(ylim)*.8,sprintf('%.1f%%',mean_bbINa.*100),'FontSize',18,'HorizontalAlignment','right');
else
    text(mean_bbINa.*100+diff(xlim)*.03,diff(ylim)*.8,sprintf('%.1f%%',mean_bbINa.*100),'FontSize',18);
end
end

figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_inclusion-bbINa_overlap-%s_noaxis',threshold_bb,threshold_a,eclimit,'prob');
set(gcf,'Name',figname);
savefigauto(gcf,fullfile(plotsavedir, figname));


%-- Percent of inclusion (average across ramdomize assigned electrode)
figure('MenuBar','none');
set(gcf,'Position',get(gcf,'Position').*[.8 .8 1.5 1]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for iroi = selroi
perm_bbINa_boot = [];
for iboot = 1:nboot
bootidx = sum(prfs.num(iroi,1:(iboot-1)))+[1:prfs.num(iroi,iboot)];
permidx = randperm(sum(prfs.num(iroi,:)),prfs.num(iroi,iboot));
 
    [~,bbINai_boot] = circintarea(...
                          reshape(prfs.ecc_bb{iroi}(permidx) .* cos(prfs.ang_bb{iroi}(permidx)*pi/180),[],1),...
                          reshape(prfs.ecc_bb{iroi}(permidx) .* sin(prfs.ang_bb{iroi}(permidx)*pi/180),[],1),...
                          reshape(prfs.rfsize_bb{iroi}(permidx),[],1),...
                          reshape(prfs.ecc_a{iroi}(bootidx) .* cos(prfs.ang_a{iroi}(bootidx)*pi/180),[],1),...
                          reshape(prfs.ecc_a{iroi}(bootidx) .* sin(prfs.ang_a{iroi}(bootidx)*pi/180),[],1),...
                          reshape(prfs.rfsize_a{iroi}(bootidx),[],1));
    perm_bbINa_boot(iboot) = mean(bbINai_boot,'omitnan');
end 
mean_bbINa_boot = mean(bbINa_boot{iroi},'omitnan');

nexttile;
histogram(perm_bbINa_boot.*100,'BinWidth',5);
set(gca,'FontSize',18);
hold on;
plot(median(perm_bbINa_boot,'omitnan').*[1 1].*100,ylim,'k-.','LineWidth',1.3);
plot(prctile(perm_bbINa_boot,95).*[1 1].*100,ylim,'k:','LineWidth',1.6);
plot(mean_bbINa_boot.*[1 1].*100,ylim,'-','LineWidth',1.6,'Color',plcol(1,:));
hold off;
title(rois{iroi});
xlim([0 100]);
if mean_bbINa_boot>0.8
    text(mean_bbINa_boot.*100-diff(xlim)*.03,diff(ylim)*.8,sprintf('%.1f%%',mean_bbINa_boot.*100),'FontSize',18,'HorizontalAlignment','right');
else
    text(mean_bbINa_boot.*100+diff(xlim)*.03,diff(ylim)*.8,sprintf('%.1f%%',mean_bbINa_boot.*100),'FontSize',18);
end
end

figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_inclusion-bbINa_overlap-%s_noaxis',threshold_bb,threshold_a,eclimit,'boot');
set(gcf,'Name',figname);
savefigauto(gcf,fullfile(plotsavedir, figname));

%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Percent of inclusion (Reverse)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%-- Percent of inclusion (average across electrode with probability asignment)
aINbb = [];
for iroi = selroi
    [~,aINbb{iroi}] = circintarea(...
                          reshape(prf_all_a.params(1,1,chidx{iroi}),[],1),...
                          reshape(prf_all_a.params(1,2,chidx{iroi}),[],1),...
                          reshape(prf_all_a.params(1,3,chidx{iroi}),[],1),...
                          reshape(prf_all_bb.params(1,1,chidx{iroi}),[],1),...
                          reshape(prf_all_bb.params(1,2,chidx{iroi}),[],1),...
                          reshape(prf_all_bb.params(1,3,chidx{iroi}),[],1));
end

% cellfun(@mean,aINbb)
% cellfun(@median,aINbb)
                     
%-- Percent of inclusion (average across ramdomize assigned electrode)
aINbb_boot = [];
for iroi = selroi
    [~,aINbb_boot{iroi}] = circintarea(...
                          reshape(prfs.ecc_a{iroi} .* cos(prfs.ang_a{iroi}*pi/180),[],1),...
                          reshape(prfs.ecc_a{iroi} .* sin(prfs.ang_a{iroi}*pi/180),[],1),...
                          reshape(prfs.rfsize_a{iroi},[],1),...
                          reshape(prfs.ecc_bb{iroi} .* cos(prfs.ang_bb{iroi}*pi/180),[],1),...
                          reshape(prfs.ecc_bb{iroi} .* sin(prfs.ang_bb{iroi}*pi/180),[],1),...
                          reshape(prfs.rfsize_bb{iroi},[],1));
end

% cellfun(@mean,aINbb_boot)
% cellfun(@median,aINbb_boot)

%% Permutation

%-- Percent of inclusion (average across electrode with probability asignment)
figure('MenuBar','none');
set(gcf,'Position',get(gcf,'Position').*[.8 .8 1.5 1]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for iroi = selroi
perm_aINbb = [];
for iboot = 1:nboot
permidx = Shuffle(chidx{iroi});
    [~,aINbbi] = circintarea(...
                          reshape(prf_all_a.params(1,1,chidx{iroi}),[],1),...
                          reshape(prf_all_a.params(1,2,chidx{iroi}),[],1),...
                          reshape(prf_all_a.params(1,3,chidx{iroi}),[],1),...
                          reshape(prf_all_bb.params(1,1,permidx),[],1),...
                          reshape(prf_all_bb.params(1,2,permidx),[],1),...
                          reshape(prf_all_bb.params(1,3,permidx),[],1));
    perm_aINbb(iboot) = mean(aINbbi);
end 
mean_aINbb = mean(aINbb{iroi},'omitnan');

nexttile;
histogram(perm_aINbb.*100,'BinWidth',5);
set(gca,'FontSize',18);
hold on;
plot(median(perm_aINbb,'omitnan').*[1 1].*100,ylim,'k-.','LineWidth',1.3);
plot(prctile(perm_aINbb,95).*[1 1].*100,ylim,'k:','LineWidth',1.6);
plot(mean_aINbb.*[1 1].*100,ylim,'-','LineWidth',1.6,'Color',plcol(1,:));
hold off;
title(rois{iroi});
xlim([0 100]);
if mean_aINbb>0.8
    text(mean_aINbb.*100-diff(xlim)*.03,diff(ylim)*.8,sprintf('%.1f%%',mean_aINbb.*100),'FontSize',18,'HorizontalAlignment','right');
else
    text(mean_aINbb.*100+diff(xlim)*.03,diff(ylim)*.8,sprintf('%.1f%%',mean_aINbb.*100),'FontSize',18);
end
end

figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_inclusion-aINbb_overlap-%s_noaxis',threshold_bb,threshold_a,eclimit,'prob');
set(gcf,'Name',figname);
savefigauto(gcf,fullfile(plotsavedir, figname));


%-- Percent of inclusion (average across ramdomize assigned electrode)
figure('MenuBar','none');
set(gcf,'Position',get(gcf,'Position').*[.8 .8 1.5 1]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
for iroi = selroi
perm_aINbb_boot = [];
for iboot = 1:nboot
bootidx = sum(prfs.num(iroi,1:(iboot-1)))+[1:prfs.num(iroi,iboot)];
permidx = randperm(sum(prfs.num(iroi,:)),prfs.num(iroi,iboot));
 
    [~,aINbbi_boot] = circintarea(...
                          reshape(prfs.ecc_a{iroi}(bootidx) .* cos(prfs.ang_a{iroi}(bootidx)*pi/180),[],1),...
                          reshape(prfs.ecc_a{iroi}(bootidx) .* sin(prfs.ang_a{iroi}(bootidx)*pi/180),[],1),...
                          reshape(prfs.rfsize_a{iroi}(bootidx),[],1),...
                          reshape(prfs.ecc_bb{iroi}(permidx) .* cos(prfs.ang_bb{iroi}(permidx)*pi/180),[],1),...
                          reshape(prfs.ecc_bb{iroi}(permidx) .* sin(prfs.ang_bb{iroi}(permidx)*pi/180),[],1),...
                          reshape(prfs.rfsize_bb{iroi}(permidx),[],1));
    perm_aINbb_boot(iboot) = mean(aINbbi_boot,'omitnan');
end 
mean_aINbb_boot = mean(aINbb_boot{iroi},'omitnan');

nexttile;
histogram(perm_aINbb_boot.*100,'BinWidth',5);
set(gca,'FontSize',18);
hold on;
plot(median(perm_aINbb_boot,'omitnan').*[1 1].*100,ylim,'k-.','LineWidth',1.3);
plot(prctile(perm_aINbb_boot,95).*[1 1].*100,ylim,'k:','LineWidth',1.6);
plot(mean_aINbb_boot.*[1 1].*100,ylim,'-','LineWidth',1.6,'Color',plcol(1,:));
hold off;
title(rois{iroi});
xlim([0 100]);
if mean_aINbb_boot>0.8
    text(mean_aINbb_boot.*100-diff(xlim)*.03,diff(ylim)*.8,sprintf('%.1f%%',mean_aINbb_boot.*100),'FontSize',18,'HorizontalAlignment','right');
else
    text(mean_aINbb_boot.*100+diff(xlim)*.03,diff(ylim)*.8,sprintf('%.1f%%',mean_aINbb_boot.*100),'FontSize',18);
end
end

figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_inclusion-aINbb_overlap-%s_noaxis',threshold_bb,threshold_a,eclimit,'boot');
set(gcf,'Name',figname);
savefigauto(gcf,fullfile(plotsavedir, figname));

%%



