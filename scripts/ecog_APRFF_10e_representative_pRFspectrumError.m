% Show spectrum shapes at representative time point
% for figure 4

% 20210301 Yuasa
% 20231031 Yuasa - Add error shading
% 20241115 Yuasa - Add option for error

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
SetDefault('isbooterror',0);    % if N>0, use bootstrapping for error with N iteration

%% Representative channels
if issaveplot
% repelec = {'name',{'Oc18','GB103','Oc17','GB102'},...
            % 'subject_name',{'p02','p10','p02','p10'}};
repelec = {'name',{'GB103'},'subject_name',{'p10'}};
ipks = 1:6;    % pick up peak #
else
% repelec = {'name',{'Oc17'},'subject_name',{'p02'}}; ipks = 5;
repelec = {'name',{'GB103'},'subject_name',{'p10'}}; ipks = 4;
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
%%%%%% Loop for subjects %%%%%%
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
    
%% plot
el = findtable(channels,repelec{:});

%-- plot parameters
xplshow_bb = [1 150];
xplshow_a  = [3 30];
ypltick = [.1 .3 1 3 10 30 100];
ypltick2 = [.25 .5 1 2 4 8];
yplminr = [.01:.01:.09 0.1:0.1:.9 1:1:9 10:10:100];

%%%%%% Loop for electrodes %%%%%%
for elec = el'
    
%-- get peak
[~,peaks] = findpeaks(datats(elec,:),'MinPeakHeight',prctile(datats(elec,:),85),'MinPeakDistance',10,'NPeaks',6,'SortStr','descend');
peaks = sort(peaks);
%   data_fits   = data_spctr(:,peaks(ipk),elec);
%   data_base   = data_spctr_off(:,1,elec);

%%%%%% Loop for peaks %%%%%%
for ipk = ipks

%-- get error and average
ipeaks = events.stim_file_index == peaks(ipk);
iblnks = ismember(events.trial_name,'BLANK');

epochs_fits = permute(spectra(elec,ipeaks,:),[3,2,1]);
epochs_base = permute(spectra(elec,iblnks,:),[3,2,1]);
data_fits   = geomean(epochs_fits,2,'omitnan');
data_base   = geomean(epochs_base,2,'omitnan');
if isbooterror
    boot_fits = zeros(length(f),isbooterror);
    boot_base = zeros(length(f),isbooterror);
    stimbaserat = round(sum(iblnks)/sum(ipeaks));
    for iboot = 1:isbooterror
        bootsmpl = randsample(1:sum(ipeaks),sum(ipeaks),true);
        boot_fits(:,iboot) = geomean(epochs_fits(:,bootsmpl,:),2,'omitnan');
        boot_base(:,iboot) = geomean(epochs_base(:,reshape((1:stimbaserat)+stimbaserat*(bootsmpl-1)',1,[]),:),2,'omitnan');
    end
    ci_alpha = 0.6827;
    [se_fits]     = quantile(boot_fits,[0.5-ci_alpha/2, 0.5+ci_alpha/2],2);
    [se_base]     = quantile(boot_base,[0.5-ci_alpha/2, 0.5+ci_alpha/2],2);
    se_l_fits = data_fits-se_fits(:,1);   se_u_fits = se_fits(:,2)-data_fits; 
    se_l_base = data_base-se_base(:,1);   se_u_base = se_base(:,2)-data_base; 
else
    [se_l_fits,se_u_fits]     = geosem(epochs_fits,0,2,'omitnan');
    [se_l_base,se_u_base]     = geosem(epochs_base,0,2,'omitnan');
end


%--  exclude power-line
if strcmpi(subjectSite{isbj},'NYU')    % NYU subjects: 60Hz line noise
%     f_sel = [1 52; 70 115; 126 175; 186 230];
%     f_sel = [1 52; 70 max(f)];
    f_nan = f<1 | (f>52 & f<70);
else                                   % Utrecht subjects: 50Hz line noise
%     f_sel = [1 42; 60 95; 106 145; 156 190; 210 240];
%     f_sel = [1 42; 60 max(f)];
    f_nan = f<1 | (f>42 & f<60);
end
data_fits(f_nan) = nan;
data_base(f_nan) = nan;

%-- get model parameters
func_model = @(X,F,C) -(X(1) - C*F - X(2)*sqrt(2*pi)*normpdf(F,X(3),X(4)));
% func_model2 = @(X,F,C,X2) -(X(1) - C*F - X(2)*sqrt(2*pi)*normpdf(F,X(3),X(4)) + X2(1)*sqrt(2*pi)*normpdf(F,X2(2),X2(3)));    
    tmp_params = resamp_parms{elec}(peaks(ipk),:);
%     bb_amp_low  = tmp_params(8);
%     alpha_amp   = tmp_params(9);
    alpha_freq  = tmp_params(10);
    alpha_width = tmp_params(11);
    slope       = tmp_params(12);
    bb0         = tmp_params(13);
%     beta_params = tmp_params(15:end)';
    

%%% broadband w/o power-line
if true
xplshow = xplshow_bb;

hF(end+1) = figure('Position',[500 500 500 420]);
%-- raw
f_plot = true(size(f));
h=errorshade([f(f_plot),f(f_plot)],[data_base(f_plot),data_fits(f_plot)],...
      [se_l_base(f_plot),se_l_fits(f_plot)],[se_u_base(f_plot),se_u_fits(f_plot)],...
      ["k","#A06440"],'LineWidth',4);
h(2,:) = [];
set(gca,'YScale','log');
xlim(xplshow);
y_lin = ylim;
legend(h,{'BLANK','STIMULUS'});

set(gca,'FontSize',FntSiz)
yTic = yticks; yticks(yTic(yTic<1e3));
xticklabels(string(xticks));    yticklabels(string(yticks));
title('Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

figureName = sprintf('%s_%s_broadband_nopower-%s-error',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end


%%% alpha
xplshow = xplshow_a;

%-- raw + shift
if true
hF(end+1) = figure('Position',[500 500 500 420]);
f_plot = f>1&f<45;
bb_shift = 10.^func_model([bb0,0,alpha_freq,alpha_width],log10(f),slope);
h=errorshade([f(f_plot),f(f_plot),f(f_plot)],...
      [data_base(f_plot),data_fits(f_plot),data_fits(f_plot)./bb_shift(f_plot)],...
      [se_l_base(f_plot),se_l_fits(f_plot),se_l_fits(f_plot)./bb_shift(f_plot)],...
      [se_u_base(f_plot),se_u_fits(f_plot),se_u_fits(f_plot)./bb_shift(f_plot)],...
      ["k","#A06440","#A06440"],'LineWidth',4);
h(2,:) = []; h(end).LineStyle = '--';
set(gca,'YScale','log');
xlim(xplshow);
f_plot = f>=(xplshow(1)) & f<(xplshow(2));
y_lin = [min([data_base(f_plot)-se_l_base(f_plot),data_fits(f_plot)-se_l_fits(f_plot),(data_fits(f_plot)-se_l_fits(f_plot))./bb_shift(f_plot)],[],'all','omitnan'),...
         max([data_base(f_plot)+se_u_base(f_plot),data_fits(f_plot)+se_u_fits(f_plot),(data_fits(f_plot)+se_u_fits(f_plot))./bb_shift(f_plot)],[],'all','omitnan')];
ylim(exp(log(y_lin) + [-1 1]*0.02*range(log(y_lin))));
% legend(h,{'BLANK','STIMULUS',sprintf('STIMULUS,\nBASELINE\nADJUSTED')});
legend(h(1:2),{'BLANK','STIMULUS'});

set(gca,'FontSize',FntSiz)
title('Power Spectra');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

%%%%% log-log
set(get(gca,'XAxis'),'Scale','log','TickValue',[min(xplshow) 10:10:max(xplshow)]);xlim(xplshow);
yTicL = [string(yticks)'; string(yticklabels)]; yTicL([yticks>=1e3, yticks<1e3]) = [];
xticklabels(string(xticks));    yticklabels(yTicL);

figureName = sprintf('%s_%s_alpha-both-%s-error',channels.subject_name{elec},channels.name{elec},int2ordinal(ipk));
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName));    end
end

end
if issaveplot,  close(hF);  end
end
end

