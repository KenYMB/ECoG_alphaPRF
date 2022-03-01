% Show spectrum shapes for large Visul Filed stimuli
% for figure 1

% 20210119 Yuasa - update from ecog_APRF_10a_outputSpectrum
% 20210430 Yuasa - same figure from bar stimuli

%% Define paths and dataset
% close all; clear all;
checkPath;
%-- Input & Output path
SetDefaultAnalysisPath;
SetDefault('issaveplot',true);
if issaveplot
    plotsavedir    = fullfile(figPth, 'spectra-representative','Bar');
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%%
filename = fullfile('Data','Spectrum','Methods_AlphaSuppression_Bar');

if exist([filename '.mat'],'file')
  load(filename);

else
  %% load spectrum & fitting
  
  subject   = 'chaam';
  outputDir = fullfile(savePth, 'Spectrum');
  
  %-- spectrum
  average        = 'trials';
  gammafit       = false;
  allownegfit    = true;
  isbetawide     = alphaFitTypes(subject,'name');
  allowbetafit   = contains(isbetawide,'beta');
  allowwidefit   = contains(isbetawide,'wide');
  
  fileid    = 'freq_spectra';
  dataname    = fullfile(outputDir, sprintf('%s_%s.mat', subject,fileid));
  idat        = load(dataname);
  
  %-- set parameters
    f      	 = idat.f;
    channels = idat.channels;
    events   = idat.events;
    spectra  = idat.spectra;
        
  %% prepare for IAF estimation
  %-- load stimulus apertures
    res     = [100 100];
    resmx   = 100;
    pix2deg = 16.6./resmx;
    hrf   = 1;
    
      apertures = load(fullfile(analysisRootPath, 'Data','stimuli','bar_apertures.mat'));
      tmp = fieldnames(apertures);
      apertures = apertures.(tmp{1});

  apertures = imresize(apertures, res, 'nearest');
  [d,xx,yy] = makegaussian2d(resmx,2,2,2,2);
  
  %-- Prepare the stimuli for use in the model
  stimulusPP = squish(apertures,2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP = [stimulusPP 1*ones(size(stimulusPP,1),1)];  % this adds a dummy column to indicate run breaks
  
  %-- set OG-CSS model
  modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));
  PRF2CSS  = @(pp,res,gain,expt) [(1+res)/2 - pp(1)*sin(pp(2)*pi/180),...
                 pp(1)*sin(pp(2)*pi/180) + (1+res)/2,...
                 pp(3)*expt, gain, expt];

  %% get selected channels
  el = find(ismember(channels.name,{'Oc17','Oc18','Oc25'}))';
  el = [el find(ismember(channels.name,{'Oc19','Oc26','Oc27'}))'];
  %   el = [1:height(channels)];
    
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
%         data_spctr_off  = permute(spectra_off,[3,2,1]);
        
        %-- set new params
        n_stims     = size(data_spctr,2);
        n_params    = 11;
        events      = events(events_idx,:);
        n_iter = 1;
        
        %%% loop across electrodes
        pltdata = struct([]);
        ee = 1;
        for elec = el
            elec_sel = channels.name{elec};
            fprintf('Processing electrode %s for subject %s across %s \n', elec_sel, subject, average);
            
            data_fits   = data_spctr(:,:,elec);
%             data_base   = data_spctr_off(:,:,elec);
            
            %--  set parameters
            if startsWith(subject,'som')    % NYU subjects: 60Hz line noise
              if gammafit
                f_use4fit   = f((f>=30 & f<=52) | (f>=70 & f<=115) | (f>=126 & f<=175) | (f>=186 & f<=194));
              else
                f_use4fit   = f((f>=70 & f<=115) | (f>=126 & f<=175));
              end
            else                            % Utrecht subjects: 50Hz line noise
              if gammafit
                f_use4fit   = f((f>=30 & f<=42) | (f>=60 & f<=95) | (f>=106 & f<=145) | (f>=156 & f<=194));
              else
                f_use4fit   = f((f>=70 & f<=95) | (f>=106 & f<=145) | (f>=156 & f<=180));
              end
            end
            if allowwidefit
                f_alpha4fit = f(f>=3 & f<=34);
            else
                f_alpha4fit = f(f>=3 & f<=26);
            end
                
            %-- fit alpha with PRF trials
            taskIndex = ~ismember(events.trial_name,'BLANK');
            data_fit  = geomean(data_fits(:,taskIndex),2,'omitnan'); % task
            data_base = geomean(data_fits(:,~taskIndex),2,'omitnan'); % baseline
            [bb_amp_low,alpha_amp,alpha_freq,alpha_width,fit_fd2,out_exp,beta_params] = ...
                ecog_fitalpha(f,f_alpha4fit,[7 14; 8 13],data_base',data_fit',allownegfit,allowbetafit);

            iaf = 10.^alpha_freq; % [Hz]
                
            pltdata(ee).fit    = data_fit;
            pltdata(ee).base   = data_base;
            pltdata(ee).bb_amp = bb_amp_low;
            pltdata(ee).iaf    = iaf;
            pltdata(ee).f      = f;
            pltdata(ee).elec   = elec_sel;
            pltdata(ee).slope  = out_exp(1);
            pltdata(ee).bb0    = out_exp(2);
            pltdata(ee).iap    = alpha_amp;
            pltdata(ee).iaw    = alpha_width;
            pltdata(ee).beta   = beta_params;
            
            ee = ee + 1;
        end
        
  save(filename,'pltdata','el','f','f_alpha4fit','taskIndex','events','channels','subject')
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all;

alpha_lim = [1 30];
% alpha_lim = [3 25];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot RAW spectrum in broadband and alpha
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_alpha     = ismember(f,f_alpha4fit);
ee = find(ismember({pltdata.elec},'Oc17'));

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

hF = figure('Position',[500 500 930 420]);
tiledlayout(hF,1,2,'TileSpacing','Compact','Padding','compact');
%-- broadband
nexttile;
h=semilogy(f,data_base,'LineWidth',4,'Color','k');
hold on;
h(2)=semilogy(f,data_fit,'LineWidth',4,'Color','#A06440');
if plotiaf, plot(f_detail([iaidx iaidx]),[aa_base(iaidx), aa_fit(iaidx)],'k:','LineWidth',2); end
if plotarrow, plot_arrow(f_detail(iaidx),aa_base(iaidx),f_detail(iaidx), aa_fit(iaidx),'LineWidth',1.6,'absheadheight',abs(diff([aa_base(iaidx), aa_fit(iaidx)]))*0.2,'absheadwidth',1.0); end
xlim([1 150])
y_lin = ylim;
hbb = gca;
legend(h,{'BLANK','STIMULUS'},'AutoUpdate','off');

set(gca,'FontSize',18)
set(gca,'XTickLabel',string(get(gca,'XTick')),'YTickLabel',string(get(gca,'YTick')));
title('Broadband');
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2/Hz)');

%-- alpha
nexttile;
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

set(gca,'FontSize',18)
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
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figureName));  end

%%
