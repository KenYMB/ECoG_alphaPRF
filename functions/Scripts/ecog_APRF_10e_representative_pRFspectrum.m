% Show spectrum shapes at representative time point
% for figure 4

% 20210301 Yuasa


%% Define paths and dataset
% close all; clear all;
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
SetDefault('issaveplot',true);
if issaveplot
    plotsavedir    = fullfile(figPth, 'spectra-representative','timepoint');
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

FontSize = 24;

%% Set representative subjects
repelec = {'name',{'GB103'},'subject_name',{'som726'}};

  average        = 'runs';
  gammafit       = false;
  
  allowlag       = false;
  allowbetafit   = false;
  allowwidefit   = false;
  
%% Load spctrum & parameter & model
subjectList = unique(repelec{4});
isbetawide = ismember(alphaFitTypes(subjectList,'name'),'betawide');

opts = [];
opts.compute        = false;
opts.issave         = false;
opts.doplots        = false;
opts.outputDir      = fullfile(savePth, 'Spectrum');
[freq] = ecog_prf_spectra(subjectList, opts);

opts = [];
opts.compute        = false;
opts.average        = average;
opts.gammafit       = false;
opts.estimateIAF    = true;
opts.allownegfit    = true;
  spcrm_params = cell(size(subjectList));
  if any(~isbetawide)
      opts.allowbetafit   = false;
      opts.allowwidefit   = false;
      tmp = ecog_prf_fitalpha(subjectList(~isbetawide), opts);
      spcrm_params(~isbetawide) = tmp;
  end
  if any(isbetawide)
      opts.allowbetafit   = true;
      opts.allowwidefit   = true;
      tmp = ecog_prf_fitalpha(subjectList(isbetawide), opts);
      spcrm_params(isbetawide) = tmp;
  end

opts = [];
opts.targetBAND     ='bbS';
opts.smoothingMode  ='none';
opts.compute        = false;
opts.issave         = false;
opts.doplots        = false;
[modeldata_bb] = ecog_prf_constructTimeSeries(subjectList, opts);
    
%%
isbj = 1;
%% Take average
  %-- set parameters
    f = freq{isbj}.f;
    channels = freq{isbj}.channels;
    events   = freq{isbj}.events;
    spectra      = freq{isbj}.spectra;
    spectra_off  = freq{isbj}.spectra_off;
    resamp_parms   =  spcrm_params{isbj}.resamp_parms;
    datats  = cat(2,modeldata_bb{isbj}.datats{:});
    
    switch average
        case {'none'},      avg_group = 1:height(events);
        case {'sessions'},  avg_group = findgroups(events.task_name,events.run_name,events.stim_file_index);
        case {'runs'},      avg_group = findgroups(events.task_name,events.stim_file_index);
        case {'stimuli'},   avg_group = findgroups(events.task_name,events.trial_type);
        case {'trials'},    avg_group = ones(height(events),1);
                            bslIndex  = contains(events.trial_name, 'BLANK');
                            avg_group(bslIndex) = 2;
        otherwise,          error('''%s'' is unknown average type',average);
    end
    n_avg     = groupcounts(avg_group);

    %-- change dims & take average across repeats (chan x events x f -> f x averaged events x channels)
    data_spctr      = zeros(length(f),length(n_avg),height(channels));
    events_idx      = zeros(length(n_avg),1);
    for iavg = 1:length(n_avg)
        data_spctr(:,iavg,:) = geomean(permute(spectra(:,avg_group==iavg,:),[3,2,1]),2,'omitnan');
        events_idx(iavg)     = find(avg_group==iavg,1);
    end
    data_spctr_off  = permute(spectra_off,[3,2,1]);
    
%% plot
el = findtable(channels,repelec{:});

ypltick = [.1 .3 1 3 10 30 100];
ypltick2 = [.25 .5 1 2 4 8];
yplminr = [.01:.01:.09 0.1:0.1:.9 1:1:9 10:10:100];

for elec = el'
    
%-- get peak
[~,peaks] = findpeaks(datats(elec,:),'MinPeakHeight',prctile(datats(elec,:),85),'MinPeakDistance',3);
ipk = 1;
    
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

%%% broadband w/o power-line
xplshow = xplshow_bb;

hF = figure('Position',[500 500 500 420]);
%--  set power-line
if startsWith(channels.subject_name{elec},'som')    % NYU subjects: 60Hz line noise
%     f_sel = [1 52; 70 115; 126 175; 186 230];
    f_sel = [1 52; 70 max(f)];
else                                                % Utrecht subjects: 50Hz line noise
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

set(gca,'FontSize',FontSize)
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

figureName = sprintf('%s_%s_broadband_nopower-1st',channels.subject_name{elec},channels.name{elec});
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName));  end


%%% alpha difference
xplshow = xplshow_a;

hF(2) = figure('Position',[500 500 500 420]);
%-- difference
plot([1 150],[0 0],'k--','LineWidth',2);
hold on;
h=plot(f,-(log10(data_base)-log10(data_fits)),'LineWidth',4,'Color','#77AC30');
xlim(xplshow);
legend(h,{'STIMULUS / BLANK'},'Location','northwest');

set(gca,'FontSize',FontSize)
% set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Spectra Ratio');
xlabel('Frequency (Hz)');
ylabel('Log Ratio');

%%% log
ylabel('Ratio');
if min(ylim)>log10(.1) && max(ylim)<log10(3)
    yticks(log10(ypltick2)); yticklabels(num2str(ypltick2'));
else
    yticks(log10(ypltick)); yticklabels(num2str(ypltick'));
end
set(get(gca,'YAxis'),'MinorTick','on','MinorTickValues',log10(yplminr));

%%% w/ baseline
plot(f,func_model([bb0,0,alpha_freq,alpha_width],log10(f),slope),':','LineWidth',2,'Color','#A2142F');
legend(h,{'STIMULUS / BLANK'},'Location','northwest');

%%% log-log
set(get(gca,'XAxis'),'Scale','log','TickValue',[min(xplshow) 10:10:max(xplshow)]);

figureName = sprintf('%s_%s_alpha_diff-loglog-base-1st',channels.subject_name{elec},channels.name{elec});
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName));  end

end
