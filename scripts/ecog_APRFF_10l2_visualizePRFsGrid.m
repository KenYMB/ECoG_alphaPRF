% visualize pRF parameters on HD-Grid computed in 06c (equivalant 03c2 + 03d2)
%   Applicable xR2 with full time-series (no decimation)

% 202011111 Yuasa - update from 06f


%% prefix
% close all; clearvars;
% startupToolboxToolbox;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'pRF_Grid-representative');
SetDefault('prfPth',        'pRFmodel');
else
plotsavePth    = 'pRF_Grid-representative';
prfPth         = 'pRFmodel';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'pRF');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'hasHDgrid', 'yes');

%-- Plotting Setting
FntSiz  = 16;

%% Set pRF parameters
%-- Channels and Grid Setting
whichHDgrid     = 'GB';

%%% Set parameters for load analyzePRF
clear alphaType broadbandType
average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';
selectchs      = [whichHDgrid '*'];
    allowlag       = false;
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;

va_area = 'wangarea';

usefulltsR2   = false;
usefulltsxR2  = false;

ecog_APRFF_INITa_loaddata;

%% Plot parameters
%%% Set threshold
allowmixbeta   = true;
selectchs      = 'wangprobchs';
ecog_APRFF_INITc_postfix;
ecog_APRFF_INITd_threshold;
selectchs      = 'gridchs';

res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;

alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');

%% %%%%%%%%
%% Grid plot

%% plot pRF parameters in Grid
% close all;

if issaveplot
flds    = {'xval','ang','ecc','rfsize'};
else
flds    = {'xval','ang','ecc'};
end
bands   = {'broadband','alpha'};
setmask = 1;     % 0: no threshold, 1: threshold for each, 2: both threshold
eccmask = 0;     % 0: no threshold, 1: set eccentricity threshold

%-- Plot
hF = gobjects(0);
for isbj = 1:length(subjectList)
% isbj = 1;

for tarBAND = bands

switch tarBAND{:}
    case {'bb','broadband'}
        iprf    = prf_params_bb{isbj};
        thresh  = threshold_bb;
    case {'a','alpha'}
        iprf    = prf_params_a{isbj};
        thresh  = threshold_a;
end
subject = iprf.subject;

if eccmask,     thresh_ecc = eclimit;
else,           thresh_ecc = nan;
end
switch setmask
    case 1
    elec_ok = ~(iprf.xval<=thresh | iprf.ecc >= thresh_ecc); 
    case 2
    elec_ok = ~(prf_params_bb{isbj}.xval<=threshold_bb | prf_params_a{isbj}.xval<=threshold_a ...
         | prf_params_bb{isbj}.ecc >= thresh_ecc | prf_params_a{isbj}.ecc >= thresh_ecc); 
end

for sbfld = flds

specs = [];
specs.channels          = iprf.channels;
specs.plot.nSubPlots    = [];
specs.plot.RotGrid      = true;
specs.plot.fontSize     = FntSiz;
specs.plot.FigName      = sprintf('pRF property');

% specs.plot.ShowLabel    = true;
specs.plot.ShowLabel    = false;

switch sbfld{:}
    case {'xval','xR2','R2'}
        %-- Variance Explained
        specs.plot.CLim         = [0 100];
        specs.plot.colorMap     = colormapskew('hot',1.4',[0 0.95]);
        specs.timelabel         = {'Variance Explained'};
        pltdat      = iprf.(sbfld{:});
        if setmask && ~ismember(sbfld,{'xval','xR2','R2'})
            threshidx = fix((thresh - specs.plot.CLim(1)) ./ diff(specs.plot.CLim)...
                                .*size(specs.plot.colorMap,1));
            specs.plot.colorMap(1:threshidx,:)     = ones(threshidx,3);
        end
    case {'ecc'}
        %-- Eccentricity
        if eccmask
          specs.plot.CLim         = [0 eclimit].*cfactor;
        else
          specs.plot.CLim         = round([0 resmx/2].*cfactor,-1);
        end
%         specs.plot.colorMap     = colormapskew('hot',1.4',[0 0.95]);
        load colormap_benson14_eccen
        specs.plot.colorMap     = cmap(1:450,2:4);
        specs.timelabel         = {'Eccentricity'};
        pltdat      = iprf.(sbfld{:}).*cfactor;
    case {'rfsize'}
        %-- pRF size
        if eccmask
          specs.plot.CLim         = round([0 eclimit].*cfactor,-1);
        else
          specs.plot.CLim         = round([0 resmx/2].*cfactor,-1);
        end
%         specs.plot.colorMap     = colormapskew('hot',1.4',[0 0.95]);
        load colormap_benson14_sigma
        specs.plot.colorMap     = cmap(1:450,2:4);
        specs.timelabel         = {'Size'};
        pltdat      = iprf.(sbfld{:}).*cfactor;
    case {'ang'}
        %-- Polar angle
        specs.plot.CLim         = [0 360];
        specs.plot.CTick        = [0:45:360];
%         specs.plot.colorMap     = circshift([jet(128); flipud(jet(128))],64,1);
%         specs.plot.colorMap     = circshift(hsv(256),0,1);
        load colormap_benson14_angle
        specs.plot.colorMap     = cmap(:,2:4);
        specs.timelabel         = {'Polar angle'};
        pltdat      = iprf.(sbfld{:});
end
if setmask && ~ismember(sbfld,{'xval','xR2','R2'})
    specs.plot.AlphaData = elec_ok;
end

specs.plot.colorBar = false;
hFc = ecog_plotGridSC(pltdat, whichHDgrid, specs);
hF = [hF hFc{:}];

%-- rotation
iflefthemi = mean(ismember(prf_params_bb{isbj}.channels.hemisphere,'L')) > 0.7;
if iflefthemi                   % if left hemisphere
    for ii=1:numel(hFc)
        set(hFc{ii}.Children,'View',[180 90]);
    end
end

if eccmask,     eccfix = sprintf('-ecc%02d',eclimit);
else,           eccfix = '';
end
switch setmask .* ~ismember(sbfld,{'xval','xR2','R2'})
    case 1
    figname = sprintf('prfGrid_%s-%02d%%%s_%s-%s',subject,thresh,eccfix,tarBAND{:},sbfld{:});
    case 2
    figname = sprintf('prfGrid_%s-%02d%%-%02d%%%s_%s-%s',subject,threshold_bb,threshold_a,eccfix,tarBAND{:},sbfld{:});
    otherwise
    figname = sprintf('prfGrid_%s_%s-%s',subject,tarBAND{:},sbfld{:});
end
if ~specs.plot.ShowLabel
    figname = sprintf('%s_nolabel',figname);
end
if issaveplot
if numel(hFc)==1
    savefigauto(hFc{1},fullfile(plotsavedir,figname),'-vector');
else
    for ii=1:numel(hFc)
    savefigauto(hFc{ii},fullfile(plotsavedir,sprintf('%s-%d',figname,ii)),'-vector');
    end
end
end

end
end
end

%% %%%%%%%%%%%%%%%%%%%%%%
%% color label
%% %%%%%%%%%%%%%%%%%%%%%%
% close all;
if issaveplot

meshres = 300;
maxrad  = meshres/2;

%% Variance Explained
figure;

[X,Y] = meshgrid(0:meshres,0:meshres/10);
Z = X;
surf(X,Y,Z,'EdgeColor','none');
axis equal;
colormap(colormapskew('hot',1.4',[0 0.95]));
set(gcf, 'InvertHardcopy', 'off', 'Color', [1 1 1], 'Position', [200 100 50 250]);
view([90 0.1]); box on;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);

figname = 'Colorbar-R2';
savefigauto(gcf, fullfile(plotsavedir,figname),'-vector');
    
%% Eccentricity
figure;

[X,Y] = meshgrid(0:meshres,0:meshres);
Z = sqrt((X-meshres/2).^2+(Y-meshres/2).^2);
Z(Z>maxrad) = nan;
surf(X,Y,Z,'EdgeColor','none');
axis equal;
load colormap_benson14_eccen
colormap(cmap(1:450,2:4));
set(gcf, 'InvertHardcopy', 'off', 'Color', [1 1 1], 'Position', [200 100 250 250]);
view([0 90]); axis off;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);

figname = 'Colorbar-ecc';
savefigauto(gcf, fullfile(plotsavedir,figname),'-vector');


%% pRF Size
figure;

[X,Y] = meshgrid(0:meshres,0:meshres/10);
Z = X;
surf(X,Y,Z,'EdgeColor','none');
axis equal;
load colormap_benson14_sigma
colormap(cmap(1:450,2:4));
set(gcf, 'InvertHardcopy', 'off', 'Color', [1 1 1], 'Position', [200 100 50 250]);
view([90 0.1]); box on;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);

figname = 'Colorbar-rfsize';
savefigauto(gcf, fullfile(plotsavedir,figname),'-vector');


%% Polar angle
figure;

[X,Y] = meshgrid(0:meshres,0:meshres);
Z = angle((X-meshres/2)+1i*(Y-meshres/2));
Z(Z<0) = Z(Z<0) + 2*pi;
Z2 = sqrt((X-meshres/2).^2+(Y-meshres/2).^2);
Z(Z2>maxrad) = nan;
surf(X,Y,Z,'EdgeColor','none');
axis equal;
load colormap_benson14_angle
colormap(cmap(:,2:4));
set(gcf, 'InvertHardcopy', 'off', 'Color', [1 1 1], 'Position', [200 100 250 250]);
view([0 90]); axis off;
set(gca,'XTick',[],'YTick',[],'ZTick',[]);

figname = 'Colorbar-ang';
savefigauto(gcf, fullfile(plotsavedir,figname),'-vector');

end
