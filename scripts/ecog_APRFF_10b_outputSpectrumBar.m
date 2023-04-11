% Show spectrum shapes for large Visul Filed stimuli
% for figure 2

% 20210119 Yuasa - update from ecog_APRF_10a_outputSpectrum
% 20210430 Yuasa - same figure from bar stimuli

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
FontSize = 18;

alpha_lim = [1 30];
% alpha_lim = [3 25];

%-- Load Bar Spectra
ecog_APRFF_10b0_prepSpectrumBar;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot RAW spectrum in broadband and alpha
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_alpha     = ismember(f,f_alpha4fit);
if issaveplot
    pltel = 1:length(pltdata);
else
    pltel = find(ismember(channels.name(el),{'Oc17'}));
end
hF = gobjects(0);
for ee = pltel

        data_fit    = pltdata(ee).fit;
        data_base   = pltdata(ee).base;
        iaf         = pltdata(ee).iaf;
        f           = pltdata(ee).f;
        elec_sel    = pltdata(ee).elec;
        bb_amp      = pltdata(ee).bb_amp;
        
%% get amplitude at alpha
[aa_base, f_detail] = resample(data_base(f_alpha),f_alpha4fit,50);
[aa_base] = interp1(f_alpha4fit,data_base(f_alpha),f_detail,'spline','extrap');
[aa_fit]  = interp1(f_alpha4fit,data_fit(f_alpha),f_detail,'spline','extrap');

iafidx = nearest(f_detail,iaf);   % f_detail(iaidx)

% figure; semilogx(f_alpha4fit,data_base(f_alpha)); hold on; plot(f_detail,aa_base,'--');
% figure; semilogx(f_alpha4fit,data_fit(f_alpha)); hold on; plot(f_detail,aa_fit,'--');

[~,iafidx_ord] = findpeaks(log10(aa_base)-log10(aa_fit));
iaf_ord = f_detail(iafidx_ord);
if length(iaf_ord)>1,   [~,midx] = min(abs(iaf_ord-10)); iafidx_ord = iafidx_ord(midx); end

%-- plot option
% subplot = @(m,n,p) subaxis(m,n,p,'MR',0.03,'ML',0.05, 'SH',0.04);
subplot = @(m,n,p) subaxis(m,n,p,'MR',0.02,'ML',0.08, 'SH',0.08,'MB',0.13); % for 2 H plots
plotiaf = false;        % show IAF by dotted line
plotarrow = false;      % show IAF by solid arrow 
useord  = false;        % show IAF recomputed from plotting instead of which by gaussian fitting
showbox = true;
% showbox = false;

%% plot in linear-log
% close all

if useord,  iaidx = iafidx_ord;
else,       iaidx = iafidx;
end

hF(end+1) = figure('Position',[500 500 930 420]);
%-- broadband
subplot(1,2,1);
h=semilogy(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=semilogy(f,data_fit,'LineWidth',4,'Color','#A06440');
if plotiaf, plot(f_detail([iaidx iaidx]),[aa_base(iaidx), aa_fit(iaidx)],'k:','LineWidth',2); end
if plotarrow, plot_arrow(f_detail(iaidx),aa_base(iaidx),f_detail(iaidx), aa_fit(iaidx),'LineWidth',1.6,'absheadheight',abs(diff([aa_base(iaidx), aa_fit(iaidx)]))*0.2,'absheadwidth',1.0); end
xlim([1 150])
y_lin = ylim;
hbb = gca;
legend(h,{'BLANK','STIMULUS'},'AutoUpdate','off');

set(gca,'FontSize',FontSize)
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Broadband');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

%-- alpha
subplot(1,2,2);
h=semilogy(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=semilogy(f,data_fit,'LineWidth',4,'Color','#A06440');
if plotiaf, plot(f_detail([iaidx iaidx]),[aa_base(iaidx), aa_fit(iaidx)],'k:','LineWidth',2); end
if plotarrow, plot_arrow(f_detail(iaidx),aa_base(iaidx),f_detail(iaidx), aa_fit(iaidx),'LineWidth',1.6,'absheadheight',abs(diff([aa_base(iaidx), aa_fit(iaidx)]))*0.2,'absheadwidth',1.0); end
xlim(alpha_lim)
y_lin = ylim;
y_linbb = [max(y_lin, ylim(hbb))*[1 0]', min(y_lin, ylim(hbb))*[0 1]'];
if showbox
  plot(hbb,alpha_lim * [1 0;1 0;0 1;0 1;1 0]',y_linbb * [1 0;0 1;0 1;1 0;1 0]','k:','LineWidth',4);
  plot(f_detail([iaidx iaidx]),y_lin,'k:','LineWidth',2.5);
end
legend(h,{'BLANK','STIMULUS'},'AutoUpdate','off');

set(gca,'FontSize',FontSize)
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Alpha');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

figureName = sprintf('spectra-log_Bar_%s_%s',subject,elec_sel);
if plotiaf,   figureName = sprintf('%s_iaf',figureName);
elseif plotarrow,   figureName = sprintf('%s_arrow',figureName);
elseif showbox,   figureName = sprintf('%s_box',figureName); end
if useord&&(plotiaf||plotarrow),   figureName = sprintf('%s-ord',figureName); end
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end

%% plot in linear-linear
% close all
if issaveplot
if useord,  iaidx = iafidx_ord;
else,       iaidx = iafidx;
end

hF(end+1) = figure('Position',[500 500 930 420]);
%-- broadband
subplot(1,2,1);
h=semilogy(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=semilogy(f,data_fit,'LineWidth',4,'Color','#A06440');
if plotiaf, plot(f_detail([iaidx iaidx]),[aa_base(iaidx), aa_fit(iaidx)],'k:','LineWidth',2); end
if plotarrow, plot_arrow(f_detail(iaidx),aa_base(iaidx),f_detail(iaidx), aa_fit(iaidx),'LineWidth',1.6,'absheadheight',abs(diff([aa_base(iaidx), aa_fit(iaidx)]))*0.2,'absheadwidth',1.0); end
xlim([1 150])
y_lin = ylim;
hbb = gca;
legend(h,{'BLANK','STIMULUS'},'AutoUpdate','off');

set(gca,'FontSize',FontSize)
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Broadband');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

%-- alpha
subplot(1,2,2);
h=plot(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=plot(f,data_fit,'LineWidth',4,'Color','#A06440');
if plotiaf, plot(f_detail([iaidx iaidx]),[aa_base(iaidx), aa_fit(iaidx)],'k:','LineWidth',2); end
if plotarrow, plot_arrow(f_detail(iaidx),aa_base(iaidx),f_detail(iaidx), aa_fit(iaidx),'LineWidth',1.6,'absheadheight',abs(diff([aa_base(iaidx), aa_fit(iaidx)]))*0.2,'absheadwidth',1.0); end
xlim(alpha_lim)
y_lin = ylim;
if showbox
  plot(hbb,alpha_lim * [1 0;1 0;0 1;0 1;1 0]',y_linbb * [1 0;0 1;0 1;1 0;1 0]','k:','LineWidth',4);
  plot(f_detail([iaidx iaidx]),y_lin,'k:','LineWidth',2.5);
end
legend(h,{'BLANK','STIMULUS'},'AutoUpdate','off');

set(gca,'FontSize',FontSize)
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Alpha');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

figureName = sprintf('spectra-linear_Bar_%s_%s',subject,elec_sel);
if plotiaf,   figureName = sprintf('%s_iaf',figureName);
elseif plotarrow,   figureName = sprintf('%s_arrow',figureName);
elseif showbox,   figureName = sprintf('%s_box',figureName); end
if useord&&(plotiaf||plotarrow),   figureName = sprintf('%s-ord',figureName); end
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end

end
%%
% close all;
