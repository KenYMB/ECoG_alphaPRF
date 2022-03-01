% visualize pRF correlations computed in 06c (equivalant 03c2 + 03d2)
%   Update xR2 with full time-series (no decimation)
%   Folk from ecog_APRF_06d_visualizePRFparamsFPM_bootall_fullts.m
%         and ecog_APRF_06e_visualizePRFsFPM_bootall_fullts.m

% 20210107 Yuasa - Update from 06d,06e
% 20210809 Yuasa - modify for paper
% 20210916 Yuasa - Update for new environment

%% Define paths and dataset
% close all; clear all;
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
SetDefault('issaveplot',true);
if issaveplot
    plotsavePth    = fullfile(figPth, 'pRFrelations-representative');
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

alpha = 0.32;   % 0.32 for 1sd, 0.05 for 2sd, 0.003 for 3sd

plcol = get(groot,'defaultAxesColorOrder');

if issaveplot
    plotsavedir    = fullfile(plotsavePth, postfix(2:end));
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%% %%%%%%%%%%%%%%%%%%%%
%%% Visualize correlations
% close all;
fontSiz = 18;

%% plot 2D histogram (replication around 90 deg for angle)
%% %%%%%%%%%%%%%%%%%%%%
%%% arrange in one figure

flds = {'ang','ecc','rfsize'};
iscircshift = false;
showfitted  = 3;        % 0=none, 1=total least square, 2=covariance ellipse, 3=2(except for ang->0)
showlim     = 90;       % degree for angle

ix = length(selroi);  iy = length(flds);
hF = figure('Position',[1000/2 900/2 700 1000],'MenuBar','none');
tiledlayout(iy,ix,'TileSpacing','compact','Padding','compact');
ii = 1;
for sbfld = flds
    for iroi = selroi
        %-- set data
        if iscell(prfs.([sbfld{:} '_bb']))
            dat1 = prfs.([sbfld{:} '_bb']){iroi};
            dat2 = prfs.([sbfld{:} '_a']){iroi};
        else
            dat1 = prfs.([sbfld{:} '_bb'])(iroi,:);
            dat2 = prfs.([sbfld{:} '_a'])(iroi,:);
        end
        if exist('prfs_corr','var') && iscell(prfs_corr.([sbfld{:} '_fit']))
            dat3 = prfs_corr.([sbfld{:} '_fit']){iroi};
        else
            dat3 = nan;
        end
        switch sbfld{:}
            case {'ecc', 'rfsize'}
                %-- unit change (pixel -> degree)
                dat1 = dat1 .* cfactor;
                dat2 = dat2 .* cfactor;
                dat3(1,:) = dat3(1,:) .* cfactor;
            case {'ang'}
                %-- circular shift
                if iscircshift
                    datShift = round(nanmedian([dat1,dat2]))-180;
                    dat1 = mod(dat1-datShift,360)+datShift;
                    dat2 = mod(dat2-datShift,360)+datShift;
                else
                    datShift = 0;
                end
                %-- replicate around 90 deg
                if isscalar(showlim), showlim = [-1 1].*showlim;    end
                dat1 = [dat1, dat1,     dat1,     dat1+360, dat1+360, dat1+360, dat1-360, dat1-360, dat1-360];
                dat2 = [dat2, dat2+360, dat2-360, dat2,     dat2+360, dat2-360, dat2,     dat2+360, dat2-360];
                showch = dat1<(360+showlim(2)) & dat1>(showlim(1)) & dat2<(360+showlim(2)) & dat2>(showlim(1));
                dat1 = dat1(showch);
                dat2 = dat2(showch);
        end
        
        nexttile;
        %-- plot histogram & set axis
        h1 = histogram2(dat1,dat2,'DisplayStyle','tile','ShowEmptyBins','on');
        switch sbfld{:}
            case {'R2'}
                xlim([0 100]);  xticks(0:20:100);
                ylim([0 100]);  yticks(0:20:100);
                Nbin = 25;
                tname = 'R^2 (%)';
            case {'ecc'}
                xlim([0 8.5]);  xticks(0:2:8);
                ylim([0 8.5]);  yticks(0:2:8);
                Nbin = 25;
                tname = 'Eccentricity (deg)';
            case {'rfsize'}
                xlim([0 10]);  xticks(0:2:10);
                ylim([0 10]);  yticks(0:2:10);
                Nbin = 25;
                tname = 'Size (deg)';
            case {'ang'}
                xlim([0 360]+[showlim(1) showlim(2)]+datShift); xticks([0 90 180 270 360]+datShift);
                ylim([0 360]+[showlim(1) showlim(2)]+datShift); yticks([0 90 180 270 360]+datShift);
%                 xlim([0 720]+datShift); xticks([-180 -90 0 90 180 270 360 450 540]);
%                                     xticklabels([180 270 0 90 180 270 360 90 180]);
%                 ylim([0 720]+datShift); yticks([-180 -90 0 90 180 270 360 450 540]);
%                                     yticklabels([180 270 0 90 180 270 360 90 180]);
                Nbin = 25*2;
                tname = 'Angle (deg)';
        end
        h1.XBinLimits = xlim;
        h1.YBinLimits = ylim;
        h1.NumBins = [Nbin Nbin];
        h1.EdgeColor = 'none';
        view([0 90]); axis square;
      if ii<=length(selroi)
        title(rois(iroi));
      end
      if mod(ii,length(selroi))==1
        txpos = [xlim*[1;0] - diff(xlim)*0.22, mean(ylim)];
        text(txpos(1),txpos(2),tname,'HorizontalAlignment','center','Rotation',90,'FontSize',fontSiz*1.1,'FontWeight','bold');
      end
      if ii == (length(flds)-1)*length(selroi) + 1
        txpos = [mean(xlim), ylim*[1;0] - diff(ylim)*0.14];
        text(txpos(1),txpos(2),'Broadband','HorizontalAlignment','center','FontSize',fontSiz);
        txpos = [xlim*[1;0] - diff(xlim)*0.12, mean(ylim)];
        text(txpos(1),txpos(2),'Alpha','HorizontalAlignment','center','Rotation',90,'FontSize',fontSiz);
      end
        set(gca,'FontSize',fontSiz)
        
        %-- set diagonal
        hold on
        axlim = [min([xlim,ylim]), max([xlim,ylim])];
        plot(axlim,axlim,'w-','LineWidth',0.8);
        
        %-- plot fitted line or ellipse
        fitname = '';
        if showfitted && ~(showfitted == 3 && ismember(sbfld,{'ang'}))
            t = xlim;
            t = linspace(t(1),t(2),100);
            
            %-- exclude circular data for anlge
            switch sbfld{:}
                case {'ang'}
                    ngch = dat2 < (dat1-180) | dat2 > (dat1+180) |...
                           (dat1<(0+datShift) & dat2<(0+datShift)) |...
                           (dat1>(360+datShift) & dat2>(360+datShift));
                    dat1 = dat1(~ngch);
                    dat2 = dat2(~ngch);
            end
                    
            if showfitted == 1  % Total Least Square
                fitname = '-tls';
%                 [eVc,eVl]  = eig(cov(dat1,dat2));
%                 [~,mainax] = diag(eVl);
                if isempty(which('fitline2derror')), tbUse('knkUtils'); end
                
                x = fitline2derror(dat1,dat2);
                y = x(2)*t+x(1);
                
                plot(t,y,'Color',plcol(1,:),'LineWidth',1.6);
                
            elseif showfitted == 2 || showfitted == 3 % Covariance Ellipse
                fitname = '-ellipse';
            
                h=error_ellipse(cov(dat1,dat2),[median(dat1),median(dat2)],'conf',1-alpha);
                h.LineWidth = 1.5;
            
            else    % Least Square Method
                fitname = '-fits';
                plfun = @(t,x) x(2,:)'*t + x(1,:)';
                ys = plfun(t,dat3);
                
                y1 = median(ys,1);
                y2 = prctile(ys,(alpha/2)*100,1);
                y3 = prctile(ys,(1-alpha/2)*100,1);
                
                fill([t, fliplr(t)], [y2, fliplr(y3)],plcol(1,:),'EdgeColor','none','FaceAlpha',0.3);
                plot(t,y1,'Color',plcol(1,:),'LineWidth',1.5);
            end
            
        end
        
        %-- plot separator
        switch sbfld{:}
            case {'ang'}
                plot(xlim,ylim+360,'w--','LineWidth',0.8);
                plot(xlim+360,ylim,'w--','LineWidth',0.8);
                plot(xlim,[1 1].*(0+datShift),'w-','LineWidth',1.0);
                plot([1 1].*(0+datShift),ylim,'w-','LineWidth',1.0);
                plot(xlim,[1 1].*(360+datShift),'w-','LineWidth',1.0);
                plot([1 1].*(360+datShift),ylim,'w-','LineWidth',1.0);
        end
        
        %-- color set
        colormap([[0,0,0];hot].^0.3);
%         hc=colorbar;
%         hc.Ticks = hc.Limits;
%         hc.TickLabels = {'0'; 'Max'};
        
        ii = ii +1;
    end
end
    set(gcf,'Name','pRF relations');
    
    figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_2Dhist%s_%s',threshold_bb,threshold_a,eclimit,fitname,'all');
    if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));  end

%% %%%%%%%%%%%%%%%%%%%%
%%% cross-correlations
fontSiz = 22;
% close all;
%% plot fitted line in ROIs

flds = {'ecc_rfsize'};
showref = true;
domean  = true;

plroi = [1 2 3 nroi];   % V1,V2,V3,high

for sbfld = flds

    sbflds = strsplit(sbfld{:},'_');
    
    [ix,iy] = arrangeinrect(length(plroi));
    hF(2) = figure('Position',[1000/2 900/2 940 700],'MenuBar','none');
    tiledlayout(iy,ix,'Padding','compact','TileSpacing','compact');
    
    ii = 1;
for iroi = plroi
    nexttile;
    jj=1;
    hl = [];
    for frqbnd = {'bb','a'}
        %-- set data
        dat1 = prfs.([sbflds{1} '_' frqbnd{:}]){iroi};
        dat2 = prfs.([sbflds{2} '_' frqbnd{:}]){iroi};
        dat3 = prfs_corrx.([sbfld{:} '_' frqbnd{:} '_fit']){iroi};
        datR = prfs_corrx.([sbfld{:} '_' frqbnd{:} '_corr']){iroi};
        datR2 = prfs_corrx.([sbfld{:} '_' frqbnd{:} '_R2']){iroi};

        %-- unit change (pixel -> degree)
        switch sbflds{1}
            case {'ecc', 'rfsize'}
                dat1 = dat1 .* cfactor;
                dat3(2,:) = dat3(2,:) ./ cfactor;
        end
        switch sbflds{2}
            case {'ecc', 'rfsize'}
                dat2 = dat2 .* cfactor;
                dat3(:,:) = dat3(:,:) .* cfactor;
        end

        %-- set axis
        switch sbflds{1}
            case {'R2'}
                xlim([0 100]);
                xlabel('R^2 (%)');
            case {'ecc'}
                xlim([0 8.5]);
                xlabel('Eccentricity (deg)');
            case {'rfsize'}
                xlim([0 10]);
                xlabel('Size (deg)');
            case {'ang'}
                xlim([0 360]); xticks([0 90 180 270 360]);
                xlabel('Angle (deg)');
        end
        switch sbflds{2}
            case {'R2'}
                ylim([0 100]);
                ylabel('R^2 (%)');
            case {'ecc'}
                ylim([0 8.5]);
                ylabel('Eccentricity (deg)');
            case {'rfsize'}
                ylim([0 10]);
                ylabel('Size (deg)');
            case {'ang'}
                ylim([0 360]); yticks([0 90 180 270 360]);
                ylabel('Angle (deg)');
        end
    
        %-- plot fitted line
        t = xlim;
        t = linspace(t(1),t(2),100);
        plfun = @(t,x) x(2,:)'*t + x(1,:)';
        ys = plfun(t,dat3);

        if domean,  y1 = mean(ys,1);
        else,       y1 = median(ys,1);
        end
        y2 = prctile(ys,(alpha/2)*100,1);
        y3 = prctile(ys,(1-alpha/2)*100,1);

        hl(jj) = plot(t,y1,'Color',plcol(jj,:),'LineWidth',1.5);
        hold on;
        fill([t, fliplr(t)], [y2, fliplr(y3)],plcol(jj,:),'EdgeColor','none','FaceAlpha',0.3);
        
%         fprintf('slope = %.2f,%.2f\n',median(dat3(2,:),2),median(datx3(2,:),2));

        jj=jj+1;
    end
    uistack(hl,'top');
    
    %-- plot reference line (Marc)
    hasref = false;
    if showref
        hasref = true;
        switch rois{iroi}
            case 'V1', Sl = 0.17; Int = 0.57; 
            case 'V2', Sl = 0.25; Int = 0.65; 
            case 'V3', Sl = 0.30; Int = 1.24; 
            otherwise, hasref = false;
        end
        if hasref
            hl(jj) = plot(t,t*Sl + Int,'--','Color',[1 1 1].*0.3,'LineWidth',1.5);
        end
    end
    
%     xlabel(sbflds{1}); ylabel(sbflds{2});
    title(rois(iroi))
    set(gca,'FontSize',fontSiz)
    box off;

%     %-- set diagonal
%     hold on
%     axlim = [min([xlim,ylim]), max([xlim,ylim])];
%     plot(axlim,axlim,'k-','LineWidth',0.6);
        
    if hasref
      legend(hl,{'Broadband','Alpha','fMRI'},'Location','northwest');
    else
      legend(hl,{'Broadband','Alpha'},'Location','northwest');
    end
    ii=ii+1;
end
set(gcf,'Name',[sbfld{:}]);

if showref,     isrefname = '_withref';
else,           isrefname = '';
end
    
figname = sprintf('prf-%02d%%-%02d%%-ecc%02d_line-roi_%s-%s%s-largeFont',threshold_bb,threshold_a,eclimit,sbflds{1},sbflds{2},isrefname);
if domean,  figname = sprintf('%s-mean',figname);   end
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname));  end
end
