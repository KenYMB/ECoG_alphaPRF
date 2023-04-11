function [params] = ecog_prf_fitalpha(freq, opts)

% Description: 
%
% [powspctrm] = ecog_prf_fitalpha(freq, [opts])
%
% Input
% - freq            = Nx1 cell-array of freq structure with the following fields:
%   - subject
%   - spectra       = channels x events x f
%   - spectra_off   = channels x 1 x f
%   - f
%   - events
%   - channels
% - opts
%   - average       = 'none','sessions','runs'(default),'stimuli' or 'trials'
%   - issave
%   - outputDir
%   - gammafit      = true or false (default): if compute gamma power
%                     it also change the metrics to compute broadband power
%                     if true, return estimated broadband power at 0Hz,
%                        (intercept of linear regression in 30-194Hz)
%                     if false, return mean power in broadband (70-180Hz)
%   - estimateIAF   = true (default) or false
%                     estimate IAF and bounds alpha peaks around the IAF
%                     (±1Hz)
%   - allownegfit   = true (default) or false,
%                     allow negative coefficients in alpha fitting 
%   - allowbetafit  = true or false (default),
%                     allow to fit a beta bump in alpha fitting 
%   - allowwidefit  = true or false (default),
%                     allow to fit with wider frequency range in alpha fitting 
%   - stimulus      = structure  % required if estimateIAF=true
%     - apertures
%     - res         = [x,y]: pixels to desampling
%     - size        = (degree): visual angle of stimulus
%   ---------------------
%   - fileid
%   - hrf
%
% Output
% - powspctrm       = Nx1 cell-array of power-spectrum structure with the following fields:
%   - subject
%   - resamp_parms  = {channels} events x parameters x boots
%   - events
%   - channels
%   - average
%   - n_avg

% Hidden options
% - opts
%   - compute       = [];
%   - skipexist     = [];
%   - fileid        = 'freq_spectra-params';
%   - hrf           = 1;
%   - bootstrap     = 0;

% Hidden usage
% - opts.compute    = false;        % load or bypass data
% - opts.allowlag   = false;
% 
% - opts.bootstrap  = number (p) > 0: apply bootstrapping with p iterations
%                     0(default)    : not apply bootstrapping
% 
% [params] = ecog_prf_fitalpha(params, [opts])
% [params] = ecog_prf_fitalpha(subjectList, [opts])

% Dependency: <ECoG_utils>, <analyzePRF>, ecog_fitgamma, ecog_fitalpha, SetDefault, saveauto

% 20200224 - Yuasa
% 20200321 - Yuasa: add 'estimateIAF' option (to disable estimating IAF before spectral fitting)
% 20200414 - Yuasa: consider 50Hz line noise for Utrecht subjects
% 20200702 - Yuasa: add 'gammafit' option (combine folk function wrote in 20200414)
% 20201117 - Yuasa: update average computation
% 20201204 - Yuasa: update for bootstrapping
% 20201209 - Yuasa: add subjects information, change to use channles.hemisphere
% 20201211 - Yuasa: minor bug fix
% 20210420 - Yuasa: optimize
% 20210907 - Yuasa: ouput additional parameter for debug 
% 20220222 - Yuasa: load recording site information (& recorded hemisphere information if needed)
%                   from subjectlist.tsv
% 20220810 - Yuasa: update for dataset with trial-by-trial spectra-off

%% Set options
%--Define inputs 
% <opts>
SetDefaultAnalysisPath('DATA','Spectrum','opts.outputDir');
SetDefault('opts.average','runs');
SetDefault('opts.issave',false);
SetDefault('opts.gammafit',false);
SetDefault('opts.estimateIAF',true);
SetDefault('opts.allownegfit',true);
SetDefault('opts.allowbetafit',false);
SetDefault('opts.allowwidefit',false);
SetDefault('opts.stimulus.apertures','bar_apertures.mat');
SetDefault('opts.stimulus.res',[100 100]);
SetDefault('opts.stimulus.size',16.6);
% <hidden opts>
SetDefault('opts.hrf',1);               % used in model for IAF estimation
SetDefault('opts.bootstrap',0);
SetDefault('opts.compute',[]);
SetDefault('opts.skipexist',[]);
SetDefault('opts.fileid','freq_spectra-params');
SetDefault('opts.allowlag',false);                      % just for loading files
SetDefault('opts.maxlag',nan);                          % just for saving and loading files

%-- check inputs and outputs
assert(~isempty(freq), 'Please provide the freq struct');
if isempty(opts.compute)
    SetDefault('opts.skipexist',true);
    opts.compute = true;
elseif ~opts.compute
    opts.skipexist = false;
else
    SetDefault('opts.skipexist',false);
end
if opts.issave && ~exist(opts.outputDir, 'dir'),     mkdir(opts.outputDir); end

%-- Collect Subject Information
subjectList_fname = 'subjectlist.tsv';
SbjInfo    = loadSbjInfo(subjectList_fname,'all');
hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');

%% Prepare prf model to estimate pRF enter trials 
if opts.compute && opts.estimateIAF
  %-- load stimulus apertures
  res     = opts.stimulus.res;
  resmx   = max(res);
  pix2deg = opts.stimulus.size./resmx;
  if ischar(opts.stimulus.apertures)
      apertures = load(opts.stimulus.apertures);
      tmp = fieldnames(apertures);
      apertures = apertures.(tmp{1});
      clear tmp;
  else
      apertures = opts.stimulus.apertures;
  end
  apertures = imresize(apertures, res, 'nearest');
  [d,xx,yy] = makegaussian2d(resmx,2,2,2,2);
  
  %-- Prepare the stimuli for use in the model
  stimulusPP = squish(apertures,2)';  % this flattens the image so that the dimensionality is now frames x pixels
  stimulusPP = [stimulusPP 1*ones(size(stimulusPP,1),1)];  % this adds a dummy column to indicate run breaks
  
  %-- set OG-CSS model
  modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),opts.hrf,dd(:,prod(res)+1));
  PRF2CSS  = @(pp,res,gain,expt) [(1+res)/2 - pp(1)*sin(pp(2)*pi/180),...
                 pp(1)*sin(pp(2)*pi/180) + (1+res)/2,...
                 pp(3)*expt, gain, expt];
end
%% Loop across subjects
cellinput = iscell(freq);
if ~cellinput,  freq = {freq};  end

params = cell(size(freq));
for ii = 1 : length(freq)
    iparams = [];
    %-- Set data
    if ~opts.compute && ~isstruct(freq{ii}) && ischar(freq{ii})
        subject     = freq{ii};
        isloadfile  = true;
        compute     = false;
        %-- Takeover parameters
        allowlag    = opts.allowlag;
        maxlag      = opts.maxlag;
    else
        ifrq        = freq{ii};
        subject     = ifrq.subject;
        isSave      = opts.issave;
        isloadfile  = opts.skipexist;
        compute     = opts.compute;
        %-- Takeover parameters
        SetDefault('ifrq.allowlag',opts.allowlag);
        SetDefault('ifrq.maxlag',opts.maxlag);
        allowlag    = ifrq.allowlag;
        maxlag      = ifrq.maxlag;
    end
    average      = opts.average;
    gammafit     = opts.gammafit;
    allownegfit  = getparam(opts,'allownegfit',ii,length(freq));
    allowbetafit = getparam(opts,'allowbetafit',ii,length(freq));
    allowwidefit = getparam(opts,'allowwidefit',ii,length(freq));
    estimateIAF  = opts.estimateIAF;
    %-- Try to load files
    postfix = cnstpostfix(opts.fileid,allowlag,average,estimateIAF,allownegfit,gammafit,allowbetafit,allowwidefit);
    if isloadfile
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        %-- load files from directory
        if exist(filename,'file') || ~compute
            fprintf('[%s] Loading power spectra parameters for subject %s from %s <%s> ',mfilename, subject, opts.outputDir,postfix(2:end));
            ifrq        = load(filename);
            isSave      = false;
            fprintf('\n');
            compute = false;
        end
    end
    
    %-- Main
    if compute
        fprintf('[%s] Computing power spectra parameters for subject %s \n',mfilename, subject);
            
        %-- Prepare arguments
        spectra     = ifrq.spectra;
        spectra_off = ifrq.spectra_off;
        f           = ifrq.f;
        channels    = ifrq.channels;
        events      = ifrq.events;
        isboot      = logical(opts.bootstrap);
        
        isnonavgoff = size(spectra_off,2)>1;    % if spectra_off has individual data for each event
        
        %-- take average across repeats & change dims (chan x events x f -> f x averaged events x channels)
        if isnonavgoff
            [spectra_off]                    = ecog_averageEvents(spectra_off,events,average,@geomean);
        end
        [data_spctr, events, n_avg, stimIdx] = ecog_averageEvents(spectra,events,average,@geomean);
        data_spctr      = permute(data_spctr,[3,2,1]);
        data_spctr_off  = permute(spectra_off,[3,2,1]);
        
        %-- set new params
        n_stims     = size(data_spctr,2);
        n_params    = 17;   %11 + 3 + 3
        if isboot,  n_iter = opts.bootstrap;
        else,       n_iter = 1;
        end
        
        %-- Recording site
        if hasSbjInfo && istablefield(SbjInfo,'site')
            RECsite = SbjInfo.site(ismember(SbjInfo.participant_id,subject));
            RECsite = RECsite{1};
        else
            RECsite = 'unknown';
        end
        
        %-- ECoG hemisphere
        if hasSbjInfo && istablefield(SbjInfo,'hemi')
            ELEChemi = SbjInfo.hemi(ismember(SbjInfo.participant_id,subject));
            ELEChemi = ELEChemi{1};
        else
            ELEChemi = 'unknown';
        end
            
        %--  set frequency parameters
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
%           f_alpha4fit = f(f>=3 & f<=20);
          f_alpha4fit = f(f>=3 & f<=26);
        end
        f_alpha     = find(f>=8 & f<=13);
        
        %%% loop across electrodes
        %-- Temporal file (back up)
        filename_tmp  = fullfile(opts.outputDir, sprintf('%s_%s%s_temp.mat', subject,opts.fileid,postfix));
        if exist(filename_tmp,'file')
            fprintf('[%s] Temporary file is found. Loading %s\n',mfilename,filename_tmp);
            tmp = load(filename_tmp);
            resamp_parms = tmp.resamp_parms;
            iter         = tmp.iter + 1;
            clear tmp;
            fprintf('[%s] Continue from %s loop\n',mfilename,iptnum2ordinal(iter));
        else
            resamp_parms = repmat({NaN(n_stims,n_params,n_iter)},height(channels),1);
            iter         = 1;
        end
        for iter = iter:n_iter
          if isboot
            %-- recompute data_spctr
            parfor iavg = 1:length(n_avg)
                data_spctr(:,iavg,:) = geomean(permute(spectra(:,randsample(stimIdx{iavg},length(stimIdx{iavg}),true),:),[3,2,1]),2,'omitnan');
            end
          end
          for elec = 1:height(channels)
            elec_sel = channels.name{elec};
            fprintf('[%s] Processing electrode %s (%d) for subject %s across %s \n',mfilename, elec_sel, iter, subject, average);
            
            data_fits   = data_spctr(:,:,elec);
            data_bases  = data_spctr_off(:,:,elec);

            %-- estimate representative peak alpha frequency across stimuli
            if estimateIAF
                if all(ismember({'bensoneccen','bensonangle','bensonsigma'},channels.Properties.VariableNames))
                    if ismember('hemisphere',channels.Properties.VariableNames) && ...
                            ismember(channels.hemisphere(elec),{'L','R'})
                        islefths = ismember(channels.hemisphere(elec),'L');
                    else
                      switch ELEChemi
                        case 'R'  % right hemisphere
                            islefths = false;
                        case 'L'  % left hemisphere
                            islefths = true;
                        case 'LR' % bilateral
                            islefths = strcmpi(channels.name{elec}(1),'L');
                        otherwise
                            error('''%s'' is unknown subject', subject);
                      end
                    end
                    bensonPRF = [channels.bensoneccen(elec)./pix2deg, ((-1)^islefths).*channels.bensonangle(elec)+90, channels.bensonsigma(elec)./pix2deg];
                else
                    bensonPRF = nan(1,3);
                end
                if any(isnan(bensonPRF))
                    modelts = nan;
                else
                    %%-- convert benson pRF into CSS parameters with expt=0.05
                    bensonCSS = PRF2CSS(bensonPRF,resmx,1,0.05);

                    stimulusPPR = stimulusPP(events.stim_file_index,:);
                    modelts = modelfun(bensonCSS,stimulusPPR);
                end
                if all(isnan(modelts)) || var(modelts,'omitnan')==0  % flat signal because of large eccentricity
                    taskIndex = ~ismember(events.trial_name,'BLANK');
                else
                    taskIndex = modelts > prctile(modelts,85);
                end

                %-- fit alpha with selected trials to estimate peak alpha frequency
                data_fit  = geomean(data_fits(:,taskIndex),2,'omitnan');  % task
                if isnonavgoff, data_base  = geomean(data_bases(:,taskIndex),2,'omitnan');
                else,           data_base  = data_bases;
                end
                [~,~,alpha_freq] = ...
                    ecog_fitalpha(f,f_alpha4fit,[],data_base',data_fit',false,allowbetafit);

                iaf = 10.^alpha_freq + [-1 1]; % [Hz] allow +-1 range
            else
                iaf = []; % empty to estimate alpha peaks in whole alpha frequency range
            end

            %%% trials
            tmp_params   = zeros(n_stims,n_params);
            parfor istim = 1:n_stims       % 20.9s -> 12.8s (4 cores)
                data_fit = data_fits(:,istim);
                if isnonavgoff, data_base = data_bases(:,istim);
                else,           data_base = data_bases;
                end

                if any(~isnan(data_fit))
                % do the fitting
                if gammafit
                  [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width] = ...
                      ecog_fitgamma(f,f_use4fit,data_base',data_fit');
                else
                  f_sel       = ismember(f,f_use4fit);
                  out_exp     = [nan,nan];
                  out_exp(2)  = mean(log10(data_base(f_sel)'),'omitnan');
                  bb_amp      = mean(log10(data_fit(f_sel)'),'omitnan');
                  gamma_amp   = nan;
                  gamma_freq  = nan;
                  gamma_width = nan;
                end
                %-- fix F, both side
                [bb_amp_low,alpha_amp,alpha_freq,alpha_width,~,alpha_exp,beta_params] = ...
                    ecog_fitalpha(f,f_alpha4fit,iaf,data_base',data_fit',allownegfit,allowbetafit);
                [~,afidx] = min((f-10.^alpha_freq).^2);
                tmp_params(istim,:) = [out_exp(1); % this is the slope used in all cases
                                    bb_amp;
                                    gamma_amp;
                                    gamma_freq;
                                    gamma_width;
                                    out_exp(2); % this is the baseline intercept
                                    % calculate alpha change
                                    mean(log10(data_fit(f_alpha)) - log10(data_base(f_alpha)));
                                    bb_amp_low;
                                    alpha_amp;
                                    alpha_freq;
                                    alpha_width;
                                    alpha_exp(1);  % the baseline slope to reconstruct alpha fitting
                                    alpha_exp(2);  % the baseline intercept to reconstruct alpha fitting
                                    log10(data_fit(afidx)) - log10(data_base(afidx));  % alpha difference
                                    beta_params'];
                else
                tmp_params(istim,:)  = nan(n_params,1);
                end
            end
            resamp_parms{elec}(:,:,iter)  = tmp_params;
          end
          %-- Temporal file (back up)
          if ~mod(iter,10)
            saveauto(filename_tmp,'iter','resamp_parms');
          end
        end
    else
        %-- Load parameters
        channels     = ifrq.channels;
        events       = ifrq.events;
        resamp_parms = ifrq.resamp_parms;
        average      = ifrq.average;
        n_avg        = ifrq.n_avg;
            SetDefault('ifrq.gammafit',opts.gammafit);
            SetDefault('ifrq.allownegfit',opts.allownegfit);
            SetDefault('ifrq.allowbetafit',opts.allowbetafit);
            SetDefault('ifrq.allowwidefit',opts.allowwidefit);
            SetDefault('ifrq.estimateIAF',opts.estimateIAF);
        gammafit     = ifrq.gammafit;
        allownegfit  = ifrq.allownegfit;
        allowbetafit = ifrq.allowbetafit;
        allowwidefit = ifrq.allowwidefit;
        estimateIAF  = ifrq.estimateIAF;
    end
    
    %-- Collect into an output struct
    iparams.subject      = subject;
    iparams.resamp_parms = resamp_parms;
    iparams.events       = events;
    iparams.channels     = channels;
    iparams.average      = average;
    iparams.n_avg        = n_avg;
    iparams.gammafit     = gammafit;
    iparams.allownegfit  = allownegfit;
    iparams.allowbetafit = allowbetafit;
    iparams.allowwidefit = allowwidefit;
    iparams.estimateIAF  = estimateIAF;
    iparams.allowlag     = allowlag;
    iparams.maxlag       = maxlag;
    params{ii} = iparams;
    
    %-- Save out the freq
    postfix = cnstpostfix(opts.fileid,allowlag,average,estimateIAF,allownegfit,gammafit,allowbetafit,allowwidefit);
    if isSave
        fprintf('[%s] Saving power spectra parameters for subject %s to %s <%s> \n',mfilename, subject, opts.outputDir,postfix(2:end));
        filename    = fullfile(opts.outputDir, sprintf('%s_%s%s.mat', subject,opts.fileid,postfix));
        saveauto(filename,'-struct','iparams');
    end
    
    %-- Temporal file (back up)
    if exist('filename_tmp','var') && exist(filename_tmp,'file')
        delete(filename_tmp);
    end
    
end
if ~cellinput,  params = params{1};  end
fprintf('[%s] Done! \n',mfilename);
end

%% Sub function
function postfix = cnstpostfix(fileid,allowlag,average,estimateIAF,allownegfit,gammafit,allowbetafit,allowwidefit)
postfix = '';
%%-- lag
addlag = allowlag && ~contains(fileid,'regresslag');
if addlag,   postfix = sprintf('%s_regresslag',postfix);	end
%%-- parameters 
fitparams = '';
addbeta = allowbetafit && ~contains(fileid,{'_beta'});
addwide = allowwidefit && ~contains(fileid,{'wide_'}) && ~endsWith(fileid,{'wide'}) ;
if addbeta,         fitparams = sprintf('%sbeta',fitparams); end
if addwide,         fitparams = sprintf('%swide',fitparams); end
if ~estimateIAF, 	fitparams = sprintf('%s-freeAlpha',fitparams); end
if ~allownegfit,  	fitparams = sprintf('%s-oneside',fitparams); end
if ~gammafit,       fitparams = sprintf('%s-nogammafit',fitparams); end
fitparams    = regexprep(fitparams,'^-*','');
postfix      = sprintf('%s_%s_avg-%s',postfix,fitparams,average);
end

function param = getparam(opt,paramname,idx,nidx)

param = opt.(paramname);
if length(param)==nidx
    param = param(idx);
elseif length(param)~=1
    error('%s.%s must be a boolean or have the same length to the 1st argument',inputname(1),paramname);
end
end