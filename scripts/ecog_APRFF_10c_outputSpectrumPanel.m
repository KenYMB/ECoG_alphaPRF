% Show spectrum shapes to demonstrate computation of alpha suppression
% Show spectrum shapes for pRF stimuli
% for figure 3

% 20210107 Yuasa - update from test_Output_Spectrum2
% 20210507 Yuasa - 3 rows

%%
% close all; clear all;
% startupToolboxToolbox;
run_checkPath;
if contains(which('resample'),'Fieldtrip')
    rmpath(fileparts(which('resample')));
end

%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'spectra-representative');
SetDefault('spctrmPth',     'Spectrum');
else
plotsavePth    = 'spectra-representative';
spctrmPth      = 'Spectrum';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'Bar');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Plotting Setting
FontSize = 24;

%-- Load Bar Spectra
ecog_APRFF_10b0_prepSpectrumBar;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot RAW spectrum w/ and W/o baseline adjustment
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xplshow = [3 25];
xplshow = [1 30];

hF = gobjects(0);
if issaveplot,  corrslope = [true, false];
else,           corrslope = true;
end
for corrslope = corrslope

f_alpha     = ismember(f,f_alpha4fit);
hF(end+1) = figure('Position',[500 500 910 1050],'MenuBar','none');
ht=tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

enames =  {'Oc18', 'Oc17', 'Oc25'};
elist = zeros(size(enames));
for ee = 1:numel(elist)
    elist(ee) = find(ismember({pltdata.elec},enames(ee)),1);
end
for ee = elist  % Oc18, Oc17, Oc25

        data_fit    = pltdata(ee).fit;
        data_base   = pltdata(ee).base;
        iaf         = pltdata(ee).iaf;
        f           = pltdata(ee).f;
        elec_sel    = pltdata(ee).elec;
        bb_amp      = pltdata(ee).bb_amp;
        slope       = pltdata(ee).slope;
        iaw         = pltdata(ee).iaw;
        
%% get amplitude at alpha
[aa_base, f_detail] = resample(data_base(f_alpha),f_alpha4fit,50);
[aa_base] = interp1(f_alpha4fit,data_base(f_alpha),f_detail,'spline','extrap');
[aa_fit]  = interp1(f_alpha4fit,data_fit(f_alpha),f_detail,'spline','extrap');

iafidx = nearest(f_detail,iaf);   % f_detail(iaidx)

% figure; semilogx(f_alpha4fit,data_base(f_alpha)); hold on; plot(f_detail,aa_base,'--');
% figure; semilogx(f_alpha4fit,data_fit(f_alpha)); hold on; plot(f_detail,aa_fit,'--');

[~,iafidx_ord] = findpeaks(log10(aa_base)-log10(aa_fit));
iafidx_ord = iafidx_ord((iafidx_ord-iafidx) == min(iafidx_ord-iafidx));
iaf_ord = f_detail(iafidx_ord);
if length(iaf_ord)>1,   [~,midx] = min(abs(iaf_ord-10)); iafidx_ord = iafidx_ord(midx); end

%-- plot option
plotiaf = true;
plotarrow = false;
% plotarrow = true;
useord  = false;
% useord  = true;

%% plot in linear-log
% close all

if useord,  iaidx = iafidx_ord;
else,       iaidx = iafidx;
end

func_model = @(X,F,C) -(X(1) - C*F - X(2)*sqrt(2*pi)*normpdf(F,X(3),X(4)));

%-- raw
nexttile;
h=semilogy(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=semilogy(f,data_fit,'LineWidth',4,'Color','#A06440');
xlim(xplshow)
y_lin = ylim;
if plotiaf, plot(f_detail([iaidx iaidx]),y_lin,'k:','LineWidth',2); end
if plotarrow, plot_arrow(f_detail(iaidx),aa_base(iaidx),f_detail(iaidx), aa_fit(iaidx),'LineWidth',1.6,'absheadheight',abs(diff([aa_base(iaidx), aa_fit(iaidx)]))*0.2,'absheadwidth',1.0); end
% xticks(unique([min(xlim) xticks]));
ylim(y_lin);

set(gca,'FontSize',FontSize)
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
if ee==elist(1), title(sprintf('Power Spectrum\nwithout Correction')); legend(h,{'BLANK','STIMULUS'});end
if ee==elist(end), xlabel('Frequency (Hz)'); end
if ee~=elist(end), xticklabels([]); end
% ylabel('Power (\muV^2/Hz)');

%-- shift
nexttile;
h=semilogy(f,data_base,'LineWidth',4,'Color','k');
hold on;
if corrslope
    bb_comp = func_model([-bb_amp,0,log10(iaf),iaw],log10(f)-log10(iaf),slope);
else
    bb_comp = bb_amp;
end
h(2)=semilogy(f,10.^(log10(data_fit)-bb_comp),'LineWidth',4,'Color','#A06440');
xlim(xplshow)
ylim(y_lin);
if plotiaf, plot(f_detail([iaidx iaidx]),y_lin,'k:','LineWidth',2); end
if plotarrow, plot_arrow(f_detail(iaidx),aa_base(iaidx),f_detail(iaidx), 10.^(log10(aa_fit(iaidx))-bb_amp),'LineWidth',1.6,'absheadheight',abs(diff([aa_base(iaidx), 10.^(log10(aa_fit(iaidx))-bb_amp)]))*0.2,'absheadwidth',1.0); end
% xticks(unique([min(xlim) xticks]));
ylim(y_lin);

set(gca,'FontSize',FontSize)
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
if ee==elist(1), title(sprintf('Baseline Adjusted\nPower Spectrum')); legend(h,{'BLANK','STIMULUS'});end
if ee==elist(end), xlabel('Frequency (Hz)'); end
if ee~=elist(end), xticklabels([]); end
% ylabel('Power (\muV^2/Hz)');
yticklabels([]);

end

% xlabel(ht,'Frequency (Hz)','FontSize',FontSize);
ylabel(ht,'Power (\muV^2/Hz)','FontSize',FontSize);

if corrslope,  subname = '_corrslope';
else,          subname = '';
end
figureName = sprintf('alpha_computation%s_%s_%s',subname,subject,'3rows');
if plotiaf,   figureName = sprintf('%s_iaf',figureName);
elseif plotarrow,   figureName = sprintf('%s_arrow',figureName);
end
if useord&&(plotiaf||plotarrow),   figureName = sprintf('%s-ord',figureName); end
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end

end
