% Show spectrum shapes to demonstrate computation of alpha suppression
% Show spectrum shapes for pRF stimuli
% for figure 1

% 20210107 Yuasa - update from test_Output_Spectrum2

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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot difference of spectra w/o and w/ fitting model in log-log
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plotcore = false;
adjustbb = false;

xplshow = [3 30];
xpltick = [3 10 20 30];
ypltick = [.5 1 2 3];
yplminr = [0.5:0.1:2 3];

f_alpha     = ismember(f,f_alpha4fit);
hF = gobjects(0);
for ee = find(ismember({pltdata.elec},'Oc17'))

        data_fit    = pltdata(ee).fit;
        data_base   = pltdata(ee).base;
        iaf         = pltdata(ee).iaf;
        f           = pltdata(ee).f;
        elec_sel    = pltdata(ee).elec;
        bb_amp      = pltdata(ee).bb_amp;
        bb0         = pltdata(ee).bb0;
        slope       = pltdata(ee).slope;
        iap         = pltdata(ee).iap;
        iaw         = pltdata(ee).iaw;
        
%% get amplitude at alpha
[~, f_detail] = resample(data_base(f_alpha),f_alpha4fit,50);
[aa_base] = interp1(f_alpha4fit,data_base(f_alpha),f_detail,'spline','extrap');
[aa_fit]  = interp1(f_alpha4fit,data_fit(f_alpha),f_detail,'spline','extrap');

iafidx = nearest(f_detail,iaf);   % f_detail(iaidx)

% figure; semilogx(f_alpha4fit,data_base(f_alpha)); hold on; plot(f_detail,aa_base,'--');
% figure; semilogx(f_alpha4fit,data_fit(f_alpha)); hold on; plot(f_detail,aa_fit,'--');

[~,iafidx_ord] = findpeaks(log10(aa_base)-log10(aa_fit));
iafidx_ord = iafidx_ord((iafidx_ord-iafidx) == min(iafidx_ord-iafidx));
iaf_ord = f_detail(iafidx_ord);

%-- plot option
plotiaf = false;
useord  = false;

%% plot in log-log
% close all

if adjustbb, bb_shift = bb_amp;
else,        bb_shift = 0;
end

if useord,  iaidx = iafidx_ord;
else,       iaidx = iafidx;
end

func_model = @(X,F,C) -(X(1) - C*F - X(2)*sqrt(2*pi)*normpdf(F,X(3),X(4)));

%-- raw spectrum
hF(end+1) = figure('Position',[500 500 500 420]);
h=loglog(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=loglog(f,10.^(log10(data_fit)-bb_shift),'LineWidth',4,'Color','#A06440');
if plotiaf, plot(f_detail([iaidx iaidx]),[aa_base(iaidx), 10.^(log10(aa_fit(iaidx))-bb_shift)],'k:','LineWidth',2); end
xlim(xplshow); xticks(xpltick);
y_lin = ylim;
legend(h,{'BLANK','STIMULUS'});

set(gca,'FontSize',FontSize)
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
% title('Correct Broadband Shift');
title('Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

%%%
figureName = sprintf('alpha_model-%s-%s_%s_%s','log','Raw',subject,elec_sel);
if plotiaf,   figureName = sprintf('%s_iaf',figureName);
end
if useord,   figureName = sprintf('%s-ord',figureName); end
if plotcore,   figureName = sprintf('%s-core',figureName); end
if ~adjustbb,  figureName = sprintf('%s-noshift',figureName); end
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end


%-- ratio in log scale
if issaveplot
for plotarrow = [true false]
hF(end+1) = figure('Position',[500 500 500 420]);
semilogx([1 150],[0 0],'k--','LineWidth',2);
hold on;
h=semilogx(f,-(log10(data_base)-log10(data_fit)+bb_shift),'-','LineWidth',4,'Color','#77AC30');
if plotiaf, plot(f_detail([iaidx iaidx]),[bb_amp-bb_shift -(log10(aa_base(iaidx))-log10(aa_fit(iaidx)))-bb_shift],'k:','LineWidth',2); end
if ~adjustbb && plotarrow
    plot(f,func_model([bb0,0,log10(iaf),iaw],log10(f),slope)-bb_shift,':','LineWidth',2,'Color','#A2142F');
    plot_arrow(f_detail(iaidx),0,f_detail(iaidx), bb_amp-bb_shift,'LineWidth',1.6,'absheadheight',0.022,'absheadwidth',1.0);
end
xlim(xplshow); xticks(xpltick);
ylim([-1 1].*0.33+round(bb_amp-bb_shift,2));
yticks(log10(ypltick));
legend(h,{'STIMULUS / BLANK'});

% set(gca,'FontSize',18,'YTick',[0.1 0.5 1 2 4]);
set(gca,'FontSize',FontSize);
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(round(10.^get(gca,'YTick'),2)));
set(get(gca,'YAxis'),'MinorTick','on','MinorTickValues',log10(yplminr));
if ~adjustbb && plotarrow,  title('Broadband Elevation');
else,                       title('Ratio of Power Spectra');
end
xlabel('Frequency (Hz)');
ylabel('Ratio');

%%%
figureName = sprintf('alpha_model-%s-%s_%s_%s','log','Ratio',subject,elec_sel);
if plotiaf,   figureName = sprintf('%s_iaf',figureName);
elseif plotarrow,   figureName = sprintf('%s_arrow',figureName); end
if useord,   figureName = sprintf('%s-ord',figureName); end
if plotcore,   figureName = sprintf('%s-core',figureName); end
if ~adjustbb,  figureName = sprintf('%s-noshift',figureName); end
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end
end


%-- ratio in log scale with model
if issaveplot,  plotarrow = [true false];
else,           plotarrow = false;
end
for plotarrow = plotarrow
hF(end+1) = figure('Position',[500 500 500 420]);
semilogx([1 150],[0 0],'k--','LineWidth',2);
hold on;
h=semilogx(f,-(log10(data_base)-log10(data_fit)+bb_shift),'-','LineWidth',4,'Color','#77AC30');
h2=semilogx(f,func_model([bb0,iap,log10(iaf),iaw],log10(f),slope)-bb_shift,'--','LineWidth',4,'Color','#A2142F');
%   plot(f_detail([iaidx iaidx]),[0 func_model([bb0,iap,log10(iaf),iaw],log10(iaf),slope)-bb_shift],'k:','LineWidth',2);
if plotarrow
    plot(f,func_model([bb0,0,log10(iaf),iaw],log10(f),slope)-bb_shift,':','LineWidth',2,'Color','#A2142F');
    plot_arrow(f_detail(iaidx),bb_amp-bb_shift,f_detail(iaidx), func_model([bb0,iap,log10(iaf),iaw],log10(iaf),slope)-bb_shift,'LineWidth',1.6,'absheadheight',0.022,'absheadwidth',1.0);
end
xlim(xplshow); xticks(xpltick);
ylim([-1 1].*0.33+round(bb_amp-bb_shift,2));
yticks(log10(ypltick));
legend([h,h2],{'STIMULUS / BLANK','MODEL'});

% set(gca,'FontSize',18,'YTick',[0.1 0.5 1 2 4]);
set(gca,'FontSize',FontSize);
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(round(10.^get(gca,'YTick'),2)));
set(get(gca,'YAxis'),'MinorTick','on','MinorTickValues',log10(yplminr));
if ~adjustbb && plotarrow,  title('Alpha Suppression');
else,                       title('Ratio of Power Spectra');
end
xlabel('Frequency (Hz)');
ylabel('Ratio');

%%%
figureName = sprintf('alpha_model-%s-%s_%s_%s','log','Model',subject,elec_sel);
if plotiaf,   figureName = sprintf('%s_iaf',figureName);
elseif plotarrow,   figureName = sprintf('%s_arrow',figureName); end
if useord,   figureName = sprintf('%s-ord',figureName); end
if plotcore,   figureName = sprintf('%s-core',figureName); end
if ~adjustbb,  figureName = sprintf('%s-noshift',figureName); end
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end


%-- ratio in log scale with model w/ shift
if issaveplot
for plotarrow = [true false]
hF(end+1) = figure('Position',[500 500 500 420]);
semilogx([1 150],[0 0],'k--','LineWidth',2);
hold on;
h=semilogx(f,-(log10(data_base)-log10(data_fit)+bb_amp),'-','LineWidth',4,'Color','#77AC30');
h2=semilogx(f,func_model([bb0,iap,log10(iaf),iaw],log10(f),slope)-bb_amp,'--','LineWidth',4,'Color','#A2142F');
plot(f,func_model([bb0,0,log10(iaf),iaw],log10(f),slope)-bb_amp,':','LineWidth',2,'Color','#A2142F');
%   plot(f_detail([iaidx iaidx]),[0 func_model([bb0,iap,log10(iaf),iaw],log10(iaf),slope)-bb_shift],'k:','LineWidth',2);
if plotarrow, plot_arrow(f_detail(iaidx),bb_amp-bb_amp,f_detail(iaidx), func_model([bb0,iap,log10(iaf),iaw],log10(iaf),slope)-bb_amp,'LineWidth',1.6,'absheadheight',0.022,'absheadwidth',1.0); end
xlim(xplshow); xticks(xpltick);
ylim([-1 1].*0.33+round(bb_amp-bb_amp,2));
yticks(log10(ypltick)); 
legend([h,h2],{'STIMULUS / BLANK','MODEL'});

% set(gca,'FontSize',18,'YTick',[0.1 0.5 1 2 4]);
set(gca,'FontSize',FontSize);
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(round(10.^get(gca,'YTick'),2)));
set(get(gca,'YAxis'),'MinorTick','on','MinorTickValues',log10(yplminr));
if ~adjustbb && plotarrow,  title('Alpha Suppression');
else,                       title('Model Fit');
end
xlabel('Frequency (Hz)');
ylabel('Ratio');

%%%
figureName = sprintf('alpha_model-%s-%s_%s_%s','log','ModelShift',subject,elec_sel);
if plotiaf,   figureName = sprintf('%s_iaf',figureName);
elseif plotarrow,   figureName = sprintf('%s_arrow',figureName); end
if useord,   figureName = sprintf('%s-ord',figureName); end
if plotcore,   figureName = sprintf('%s-core',figureName); end
if ~adjustbb,  figureName = sprintf('%s-noshift',figureName); end
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end
end


%-- broadband elevation
plotarrow = true;
hF(end+1) = figure('Position',[500 500 500 420]);
semilogx([1 150],[0 0],'k--','LineWidth',2);
hold on;
plot(f,func_model([bb0,0,log10(iaf),iaw],log10(f),slope)-bb_shift,':','LineWidth',2.5,'Color','#A2142F');
%   plot(f_detail([iaidx iaidx]),[0 func_model([bb0,iap,log10(iaf),iaw],log10(iaf),slope)-bb_shift],'k:','LineWidth',2);
if plotarrow, plot_arrow(f_detail(iaidx),0,f_detail(iaidx), bb_amp-bb_shift,'LineWidth',1.6,'absheadheight',0.022,'absheadwidth',1.0); end
xlim(xplshow); xticks(xpltick);
ylim([-1 1].*0.33+round(bb_amp-bb_amp,2));
yticks(log10(ypltick));
legend([h,h2],{'STIMULUS / BLANK','MODEL'});

% set(gca,'FontSize',18,'YTick',[0.1 0.5 1 2 4]);
set(gca,'FontSize',FontSize);
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(round(10.^get(gca,'YTick'),2)));
set(get(gca,'YAxis'),'MinorTick','on','MinorTickValues',log10(yplminr));
title('Broadband Elevation');
xlabel('Frequency (Hz)');
ylabel('Ratio');

%%%
figureName = sprintf('alpha_model-%s-%s_%s_%s','log','bbElevation',subject,elec_sel);
if plotiaf,   figureName = sprintf('%s_iaf',figureName);
elseif plotarrow,   figureName = sprintf('%s_arrow',figureName); end
if useord,   figureName = sprintf('%s-ord',figureName); end
if plotcore,   figureName = sprintf('%s-core',figureName); end
if ~adjustbb,  figureName = sprintf('%s-noshift',figureName); end
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end


%-- alpha suppression
plotarrow = true;
hF(end+1) = figure('Position',[500 500 500 420]);
semilogx([1 150],[0 0],'k--','LineWidth',2);
hold on;
h2=semilogx(f,func_model([bb0,iap,log10(iaf),iaw],log10(f),slope)-func_model([bb0,0,log10(iaf),iaw],log10(f),slope),'--','LineWidth',2.5,'Color','#A2142F');
%   plot(f_detail([iaidx iaidx]),[0 func_model([bb0,iap,log10(iaf),iaw],log10(iaf),slope)-bb_shift],'k:','LineWidth',2);
if plotarrow, plot_arrow(f_detail(iaidx),bb_amp-bb_amp,f_detail(iaidx), func_model([bb0,iap,log10(iaf),iaw],log10(iaf),slope)-bb_amp,'LineWidth',1.6,'absheadheight',0.022,'absheadwidth',1.0); end
xlim(xplshow); xticks(xpltick);
ylim([-1 1].*0.33+round(bb_amp-bb_amp,2));
yticks(log10(ypltick));

% set(gca,'FontSize',18,'YTick',[0.1 0.5 1 2 4]);
set(gca,'FontSize',FontSize);
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(round(10.^get(gca,'YTick'),2)));
set(get(gca,'YAxis'),'MinorTick','on','MinorTickValues',log10(yplminr));
title('Alpha Suppression');
xlabel('Frequency (Hz)');
ylabel('Ratio');

%%%
figureName = sprintf('alpha_model-%s-%s_%s_%s','log','aSuppression',subject,elec_sel);
if plotiaf,   figureName = sprintf('%s_iaf',figureName);
elseif plotarrow,   figureName = sprintf('%s_arrow',figureName); end
if useord,   figureName = sprintf('%s-ord',figureName); end
if plotcore,   figureName = sprintf('%s-core',figureName); end
if ~adjustbb,  figureName = sprintf('%s-noshift',figureName); end
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end

end

