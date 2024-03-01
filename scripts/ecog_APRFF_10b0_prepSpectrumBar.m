% Show spectrum shapes for large Visul Filed stimuli
% for figure 2

% 20220523 Yuasa - separate this part
% 20231031 Yuasa - add error information

%%
% close all; clear all;
% startupToolboxToolbox;
run_checkPath;

%-- Input & Output path
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('spctrmPth',        'Spectrum');
else
spctrmPth      = 'Spectrum';
end
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%%
filename = fullfile(SetDefaultAnalysisPath('DAT',spctrmPth),'Methods_AlphaSuppression_Bar');

if exist([filename '.mat'],'file')
  load(filename);

else
  %% load spectrum & fitting
  
  subject   = intersect(subjectList,{'p02'});  subject = subject{1};
  outputDir = SetDefaultAnalysisPath('DAT',spctrmPth);
  
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
        
  %% check spectra in all electrodes  ->  Oc25, Oc17, Oc18 is better for explanation
  %{
  f_plot = [1 150];
        foi     = f>=f_plot(1) & f<=f_plot(2);
        
        trials = [];
        trials.pwrspctrm    = spectra(:,:,foi);
        trials.f            = f(foi);
        trials.events       = events;
        trials.channels     = channels;
        
        %%-- set 0 for prf BLANK, 1 for prf stimuli, and 10+trial_type for other stimuli
        stmlist = double(ismember(trials.events.task_name,'prf')&~ismember(trials.events.trial_name,'BLANK')) + ...
                  ~ismember(trials.events.task_name,'prf').*(trials.events.trial_type+10); 

        %%-- set trial name 'PRF' & set specs
        trials.events.trial_name(stmlist==1) = {'STIMULUS'};
        eventList   = table(trials.events.trial_name,stmlist,'VariableNames',{'trial_name','stmlist'});
        eventList   = sortrows(eventList,'stmlist');
        eventList   = unique(eventList.trial_name,'stable');
                
  specs = [];
  specs.plot.XLim      = f_plot;
  specs.plot.fontSize  = 16;
  specs.plot.XScale   = 'linear';
  specs.plot.mean     = 'geometric';

  ecog_plotGridSpectra(trials, trials.channels.name, eventList,[],specs);

  %%
  close;
  %}
  %% prepare for IAF estimation
  %-- load stimulus apertures
    res     = [100 100];
    resmx   = 100;
    pix2deg = 16.6./resmx;
    hrf   = 1;
    
      apertures = load('bar_apertures.mat');
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

  %-- Collect Subject Information
  SbjInfo    = loadSbjInfo(subjectList_fname,'all');
  hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');
  
  %% get selected channels
  el = find(ismember(channels.name,{'Oc17','Oc18','Oc25'}))';
%   el = [el find(ismember(channels.name,{'Oc19','Oc26','Oc27'}))'];
%     el = [1:height(channels)];
    
        %-- take average across repeats & change dims (chan x averaged events x f)
        [data_spctr, events, ~, events_idx]  = ecog_averageEvents(spectra,events,average,@geomean);
        data_spctr      = permute(data_spctr,[3,2,1]);        
        
        %-- Recording site
        if hasSbjInfo && istablefield(SbjInfo,'site')
            RECsite = SbjInfo.site(ismember(SbjInfo.participant_id,subject));
            RECsite = RECsite{1};
        else
            RECsite = 'unknown';
        end
        
        %%% loop across electrodes
        pltdata = struct([]);
        ee = 1;
        for elec = el
            elec_sel = channels.name{elec};
            fprintf('Processing electrode %s for subject %s across %s \n', elec_sel, subject, average);
            
            %-- pickup electrode (chan x events x f -> f x events x 1)
            data_elec   = permute(data_spctr(elec,:,:),[3,2,1]);
            spctr_elec  = permute(spectra(elec,:,:),[3,2,1]);
            
            %--  set parameters
            switch upper(RECsite)
                case {'NYU'}    % NYU subjects: 60Hz line noise
                  if gammafit
                    f_use4fit   = f((f>=30 & f<=52) | (f>=70 & f<=115) | (f>=126 & f<=175) | (f>=186 & f<=194));
                  else
                    f_use4fit   = f((f>=70 & f<=115) | (f>=126 & f<=175));
                  end
                case {'UMCU'}   % UMCU subjects: 50Hz line noise
                  if gammafit
                    f_use4fit   = f((f>=30 & f<=42) | (f>=60 & f<=95) | (f>=106 & f<=145) | (f>=156 & f<=194));
                  else
                    f_use4fit   = f((f>=70 & f<=95) | (f>=106 & f<=145) | (f>=156 & f<=180));
                  end
                otherwise
                  if gammafit
                    f_use4fit   = f(f>=30 & f<=194);
                  else
                    f_use4fit   = f(f>=70 & f<=180);
                  end
            end
            if allowwidefit
                f_alpha4fit = f(f>=3 & f<=34);
            else
                f_alpha4fit = f(f>=3 & f<=26);
            end
                
            %-- fit alpha with PRF trials
            taskIndex = ismember(events.trial_name,'prf');
            blnkIndex = ismember(events.trial_name,'BLANK');
            data_fit  = data_elec(:,taskIndex); % task     (already averaged)
            data_base = data_elec(:,blnkIndex); % baseline (already averaged)
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

            %-- extra information
            [pltdata(ee).fit_se(:,1),pltdata(ee).fit_se(:,2)] ...
                               = geosem(spctr_elec(:,events_idx{taskIndex}),2,'omitnan');  % task
            [pltdata(ee).base_se(:,1),pltdata(ee).base_se(:,2)] ...
                               = geosem(spctr_elec(:,events_idx{blnkIndex}),2,'omitnan'); % baseline
            pltdata(ee).fit_n  = sum(~isnan(spctr_elec(:,events_idx{taskIndex})),2);  % task
            pltdata(ee).base_n = sum(~isnan(spctr_elec(:,events_idx{blnkIndex})),2); % baseline
            
            ee = ee + 1;
        end

        channels = channels(el,:);
        
  save(filename,'pltdata','f','f_alpha4fit','events','channels','subject')
end
