% Show spectrum shapes at representative time point
% for figure 4

% 20210301 Yuasa

%%
% close all; clear all;
% startupToolboxToolbox;
run_checkPath;

%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'spectra-representative');
SetDefault('spctrmPth',     'Spectrum');
SetDefault('prfPth',        'pRFmodel');
else
plotsavePth    = 'spectra-representative';
spctrmPth      = 'Spectrum';
prfPth         = 'pRFmodel';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'timepoint');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
sbjectInfo  = loadSbjInfo(subjectList_fname);
subjectList = sbjectInfo.participant_id;
subjectSite = sbjectInfo.site;
isbetawide  = ismember(alphaFitTypes(subjectList,'name'),'betawide');

%-- Plotting Setting
FntSiz = 24;

%% Representative channels
if issaveplot
repelec = {'name',{'Oc18','GB103','Oc17','GB102'},...
            'subject_name',{'p02','p10','p02','p10'}};
ipk = 1;    % pick up peak #
else
% repelec = {'name',{'Oc17'},'subject_name',{'p02'}}; ipk = 5;
repelec = {'name',{'GB103'},'subject_name',{'p10'}}; ipk = 4;
end

  average        = 'runs';
  gammafit       = false;
  
  allowlag       = false;
  allowbetafit   = false;
  allowwidefit   = false;
  
%% Load spctrum & parameter & model
selsbj      = ismember(subjectList,repelec{4});
subjectList = subjectList(selsbj);
subjectSite = subjectSite(selsbj);
isbetawide  = isbetawide(selsbj);

selset      = ismember(repelec{4},subjectList);
repelec{2}  = repelec{2}(selset);
repelec{4}  = repelec{4}(selset);

opts = [];
opts.outputDir      = spctrmPth;
opts.compute        = false;
opts.allowlag       = allowlag;
[freq] = ecog_prf_spectra(subjectList, opts);

opts = [];
opts.outputDir      = spctrmPth;
opts.compute        = false;
opts.average        = average;
opts.gammafit       = false;
opts.estimateIAF    = true;
opts.allownegfit    = true;
opts.allowbetafit   = isbetawide;
opts.allowwidefit   = isbetawide;
spcrm_params = ecog_prf_fitalpha(subjectList, opts);

opts = [];
opts.outputDir      = prfPth;
opts.compute        = false;
opts.targetBAND     ='bbS';
opts.smoothingMode  ='none';
modeldata_bb = ecog_prf_constructTimeSeries(subjectList, opts);

    
%%
% isbj = 2;
hF = gobjects(0);
for isbj = 1:length(subjectList)
%% Take average
  %-- set parameters
    f = freq{isbj}.f;
    channels = freq{isbj}.channels;
    events   = freq{isbj}.events;
    spectra      = freq{isbj}.spectra;
    spectra_off  = freq{isbj}.spectra_off;
    resamp_parms   =  spcrm_params{isbj}.resamp_parms;
    datats  = cat(2,modeldata_bb{isbj}.datats{:});
    
    %-- take average across repeats & change dims (chan x events x f -> f x averaged events x channels)
    [data_spctr, events] = ecog_averageEvents(spectra,events,average,@geomean);
    data_spctr      = permute(data_spctr,[3,2,1]);
    data_spctr_off  = permute(spectra_off,[3,2,1]);
    
%% plot
el = findtable(channels,repelec{:});

ypltick = [.1 .3 1 3 10 30 100];
ypltick2 = [.25 .5 1 2 4 8];
yplminr = [.01:.01:.09 0.1:0.1:.9 1:1:9 10:10:100];

for elec = el'
    
%-- get peak
[~,peaks] = findpeaks(datats(elec,:),'MinPeakHeight',prctile(datats(elec,:),85),'MinPeakDistance',10,'NPeaks',6,'SortStr','descend');
peaks = sort(peaks);

data_fits   = data_spctr(:,peaks(ipk),elec);
data_base   = data_spctr_off(:,1,elec);

%-- get model parameters

func_model = @(X,F,C) -(X(1) - C*F - X(2)*sqrt(2*pi)*normpdf(F,X(3),X(4)));
func_model2 = @(X,F,C,X2) -(X(1) - C*F - X(2)*sqrt(2*pi)*normpdf(F,X(3),X(4)) + X2(1)*sqrt(2*pi)*normpdf(F,X2(2),X2(3)));    
    
    tmp_params = resamp_parms{elec}(peaks(ipk),:);
    
    bb_amp_low  = tmp_params(8);
    alpha_amp   = tmp_params(9);
    alpha_freq  = tmp_params(10);
    alpha_width = tmp_params(11);
    slope       = tmp_params(12);
    bb0         = tmp_params(13);
    beta_params = tmp_params(15:end)';
    
    
% bb_amp = resamp_parms{elec}(peaks(ipk),8);
bb_amp = bb_amp_low;

xplshow_bb = [1 150];
xplshow_a  = [3 30];

%%% broadband
if issaveplot
xplshow = xplshow_bb;

hF(end+1) = figure('Position',[500 500 500 420]);
%-- raw
h=semilogy(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=semilogy(f,data_fits,'LineWidth',4,'Color','#A06440');
xlim(xplshow);
y_lin = ylim;
legend(h,{'BLANK','STIMULUS'});

set(gca,'FontSize',FntSiz)
yTic = yticks; yticks(yTic(yTic<1e3));
xticklabels(string(xticks));    yticklabels(string(yticks));
title('Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

figureName = sprintf('%s_%s_broadband-%s',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end

%%% broadband w/o power-line
if true
xplshow = xplshow_bb;

hF(end+1) = figure('Position',[500 500 500 420]);
%--  set power-line
if strcmpi(subjectSite{isbj},'NYU')    % NYU subjects: 60Hz line noise
%     f_sel = [1 52; 70 115; 126 175; 186 230];
    f_sel = [1 52; 70 max(f)];
else                                   % Utrecht subjects: 50Hz line noise
%     f_sel = [1 42; 60 95; 106 145; 156 190; 210 240];
    f_sel = [1 42; 60 max(f)];
end
%-- raw
for ifreq = 1:size(f_sel,1)
  f_plot = f>=f_sel(ifreq,1) & f<=f_sel(ifreq,2);
  h=semilogy(f(f_plot),data_base(f_plot),'LineWidth',4,'Color','k');
  hold on;
  h(2)=semilogy(f(f_plot),data_fits(f_plot),'LineWidth',4,'Color','#A06440');
end
xlim(xplshow);
y_lin = ylim;
legend(h,{'BLANK','STIMULUS'});

set(gca,'FontSize',FntSiz)
yTic = yticks; yticks(yTic(yTic<1e3));
xticklabels(string(xticks));    yticklabels(string(yticks));
title('Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

figureName = sprintf('%s_%s_broadband_nopower-%s',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end

%%% alpha
xplshow = xplshow_a;

%-- raw
if issaveplot
hF(end+1) = figure('Position',[500 500 500 420]);
h=semilogy(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=semilogy(f,data_fits,'LineWidth',4,'Color','#A06440');
xlim(xplshow);
y_lin = ylim;
legend(h,{'BLANK','STIMULUS'});

set(gca,'FontSize',FntSiz)
yTicL = [string(yticks)'; string(yticklabels)]; yTicL([yticks>=1e3, yticks<1e3]) = [];
xticklabels(string(xticks));    yticklabels(yTicL);
title('Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

figureName = sprintf('%s_%s_alpha-%s',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end

%-- shift
if issaveplot
hF(end+1) = figure('Position',[500 500 500 420]);
h=semilogy(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=semilogy(f,data_fits./10.^func_model([bb0,0,alpha_freq,alpha_width],log10(f),slope),'LineWidth',4,'Color','#A06440');
xlim(xplshow);
y_lin = ylim;
legend(h,{'BLANK','STIMULUS'});

set(gca,'FontSize',FntSiz)
yTicL = [string(yticks)'; string(yticklabels)]; yTicL([yticks>=1e3, yticks<1e3]) = [];
xticklabels(string(xticks));    yticklabels(yTicL);
title('Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

figureName = sprintf('%s_%s_alpha-shift-%s',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end

%-- raw + shift
if true
hF(end+1) = figure('Position',[500 500 500 420]);
h=semilogy(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=semilogy(f,data_fits,'LineWidth',4,'Color','#A06440');
semilogy(f,data_fits./10.^func_model([bb0,0,alpha_freq,alpha_width],log10(f),slope),'--','LineWidth',4,'Color','#A06440');
xlim(xplshow);
legend(h,{'BLANK','STIMULUS'});

set(gca,'FontSize',FntSiz)
title('Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

%%%%% log-log
set(get(gca,'XAxis'),'Scale','log','TickValue',[min(xplshow) 10:10:max(xplshow)]);xlim(xplshow);
yTicL = [string(yticks)'; string(yticklabels)]; yTicL([yticks>=1e3, yticks<1e3]) = [];
xticklabels(string(xticks));    yticklabels(yTicL);

figureName = sprintf('%s_%s_alpha-both-%s',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end

%%% alpha difference
xplshow = xplshow_a;

%-- difference
if issaveplot
hF(end+1) = figure('Position',[500 500 500 420]);
plot([1 150],[0 0],'k--','LineWidth',2);
hold on;
h=plot(f,-(log10(data_base)-log10(data_fits)),'LineWidth',4,'Color','#77AC30');
xlim(xplshow);
legend(h,{'STIMULUS / BLANK'},'Location','northwest');

set(gca,'FontSize',FntSiz)
% set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Ratio of Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Log Ratio');

%%%% log
ylabel('Ratio');
if min(ylim)>log10(.1) && max(ylim)<log10(3)
    yticks(log10(ypltick2)); yticklabels(num2str(ypltick2'));
else
    yticks(log10(ypltick)); yticklabels(num2str(ypltick'));
end
set(get(gca,'YAxis'),'MinorTick','on','MinorTickValues',log10(yplminr));

%%%% w/ baseline
plot(f,func_model([bb0,0,alpha_freq,alpha_width],log10(f),slope),':','LineWidth',2,'Color','#A2142F');
legend(h,{'STIMULUS / BLANK'},'Location','northwest');

%%%%% log-log
set(get(gca,'XAxis'),'Scale','log','TickValue',[min(xplshow) 10:10:max(xplshow)]);

figureName = sprintf('%s_%s_alpha_diff-loglog-base-%s',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end

%%% alpha model
xplshow = xplshow_a;

%-- difference
if issaveplot
hF(end+1) = figure('Position',[500 500 500 420]);
plot([1 150],[0 0],'k--','LineWidth',2);
hold on;
h=plot(f,-(log10(data_base)-log10(data_fits)+bb_amp),'LineWidth',4,'Color','#77AC30');
xlim(xplshow);
y_lin = ylim;
legend(h,{'STIMULUS / BLANK'},'Location','northwest');

set(gca,'FontSize',FntSiz)
% set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Ratio of Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Log Ratio');

%%%% log
ylabel('Ratio');
if min(ylim)>log10(.1) && max(ylim)<log10(3)
    yticks(log10(ypltick2)); yticklabels(num2str(ypltick2'));
else
    yticks(log10(ypltick)); yticklabels(num2str(ypltick'));
end
set(get(gca,'YAxis'),'MinorTick','on','MinorTickValues',log10(yplminr));

%%%% w/ baseline
plot(f,func_model([bb0,0,alpha_freq,alpha_width],log10(f),slope)-bb_amp,':','LineWidth',2,'Color','#A2142F');
legend(h,{'STIMULUS / BLANK'},'Location','northwest');

%%%% log-log
set(get(gca,'XAxis'),'Scale','log','TickValue',[min(xplshow) 10:10:max(xplshow)]);

figureName = sprintf('%s_%s_alpha_model-loglog-base-%s',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end


%%% alpha model check
xplshow = xplshow_a;

%-- difference
if issaveplot
hF(end+1) = figure('Position',[500 500 500 420]);
plot([1 150],[0 0],'k--','LineWidth',2);
hold on;
h=plot(f,-(log10(data_base)-log10(data_fits)+bb_amp),'LineWidth',4,'Color','#77AC30');
h2=plot(f,func_model([bb0,alpha_amp,alpha_freq,alpha_width],log10(f),slope)-bb_amp,'--','LineWidth',4,'Color','#A2142F');
plot(f,func_model([bb0,0,alpha_freq,alpha_width],log10(f),slope)-bb_amp,':','LineWidth',2,'Color','#A2142F');
xlim(xplshow);
y_lin = ylim;
legend(h,{'STIMULUS / BLANK'},'Location','northwest');

set(gca,'FontSize',FntSiz)
% set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Ratio of Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Log Ratio');

figureName = sprintf('%s_%s_alpha_modelcheck-%s',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end


%%% alpha model check with beta
xplshow = xplshow_a;

%-- difference
if issaveplot
hF(end+1) = figure('Position',[500 500 500 420]);
plot([1 150],[0 0],'k--','LineWidth',2);
hold on;
h=plot(f,-(log10(data_base)-log10(data_fits)+bb_amp),'LineWidth',4,'Color','#77AC30');
h2=plot(f,func_model2([bb0,alpha_amp,alpha_freq,alpha_width],log10(f),slope,beta_params)-bb_amp,'--','LineWidth',4,'Color','#A2142F');
plot(f,func_model2([bb0,0,alpha_freq,alpha_width],log10(f),slope,beta_params.*[0,1,1])-bb_amp,':','LineWidth',2,'Color','#A2142F');
xlim(xplshow);
y_lin = ylim;
legend(h,{'STIMULUS / BLANK'},'Location','northwest');

set(gca,'FontSize',FntSiz)
% set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Ratio of Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Log Ratio');

figureName = sprintf('%s_%s_alpha_modelcheck-beta-%s',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end

if issaveplot,  close(hF);  end
end
end

