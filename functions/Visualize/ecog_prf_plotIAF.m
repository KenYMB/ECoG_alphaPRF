function [h,fitness] = ecog_prf_plotIAF(freq, opts)

% Description: 
%
% [h,fitness] = ecog_prf_plotIAF(freq, [opts])
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
%   - allowbetafit  = true or false (default),
%                     allow to fit a beta bump in alpha fitting 
%   - allowwidefit  = true or false (default),
%                     allow to fit with wider frequency range in alpha fitting 
%   - consvfit      = true or false (default),
%                     if true, limit alpha frequency range strictly
%   - plotsavedir
%   - plot          % see ecog_plotGridSpectra for this option
%   ---------------------
%   - fileid
%

% Hidden options
% - opts
%   - fileid        = 'tsrfreq';


% Dependency: <ECoG_utils>, <analyzePRF>, ecog_fitgamma, ecog_fitalpha,
%             <fieldTrip> copyfields, removefields,
%             SetDefault, saveauto

% 202220823 - Yuasa

%% Set options
%--Define inputs 
% <opts>
SetDefaultAnalysisPath('FIGURE','spectrum','opts.plotsavedir');
SetDefault('opts.allowbetafit',false);
SetDefault('opts.allowwidefit',false);
SetDefault('opts.consvfit',false);
SetDefault('opts.doplots',true);
SetDefault('opts.plot.XScale','log');
SetDefault('opts.plot.fontSize',24);
SetDefault('opts.plot.showlegend',[],true);

% <hidden opts>
SetDefault('opts.fileid','tsrfreq');
SetDefault('opts.allowlag',false);                      % just for file names

%-- check inputs and outputs
assert(~isempty(freq), 'Please provide the freq struct');
if opts.doplots&&~exist(opts.plotsavedir, 'dir'), mkdir(opts.plotsavedir); end

%-- Collect Subject Information
subjectList_fname = 'subjectlist.tsv';
SbjInfo    = loadSbjInfo(subjectList_fname,'all');
hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');

%% Loop across subjects
cellinput = iscell(freq);
if ~cellinput,  freq = {freq};  end

h = gobjects(0);
fitness = cell(size(freq));
for ii = 1 : length(freq)
    if isempty(freq{ii}),   continue;   end
    %-- Set data
    assert(isstruct(freq{ii}), 'Please provide the freq struct');
    
    ifrq        = freq{ii};
    subject     = ifrq.subject;
    %-- Takeover parameters
    SetDefault('ifrq.allowlag',opts.allowlag);
    allowlag    = ifrq.allowlag;
                    
    allowbetafit  = getparam(opts,'allowbetafit',ii,length(freq));
    allowwidefit  = getparam(opts,'allowwidefit',ii,length(freq));
    consvfit = opts.consvfit;
    
    %-- Try to load files
    fileid = cnstpostfix(opts.fileid,allowlag,allowbetafit,allowwidefit);
    
    %-- HD grid flag
    if hasSbjInfo && istablefield(SbjInfo,'hasHDgrid')
        hasHDgrid = any(strcmpi(SbjInfo.hasHDgrid(ismember(SbjInfo.participant_id,subject)),'yes'));
    else
        hasHDgrid = false;
    end
    
    %%% Main
    %-- Prepare arguments
    subject     = ifrq.subject;
    spectra     = ifrq.spectra;
    spectra_off = ifrq.spectra_off;
    f           = ifrq.f;
    channels    = ifrq.channels;
    events      = ifrq.events;

    %-- use same value for all stimulus
    spectra_off = geomean(spectra_off,2,'omitnan');

    %-- take average across repeats & change dims (chan x events x f -> f x averaged events x channels)
    [data_spctr, events] = ecog_averageEvents(spectra,events,'trials',@geomean);
    isblank = ismember(events.trial_name,'BLANK');
    data_spctr(:,isblank,:) = [];    events(isblank,:) = [];
    if size(data_spctr,2) > 1
        [data_spctr, events] = ecog_averageEvents(data_spctr,events,'all',@geomean);
    end
    data_spctr      = permute(data_spctr,[3,2,1]);
    data_spctr_off  = permute(spectra_off,[3,2,1]);
        

    %--  set frequency parameters
    if allowwidefit
      f_alpha4fit = f(f>=3 & f<=34);
    else
      f_alpha4fit = f(f>=3 & f<=26);
    end
    if consvfit == -1
        f_alphalimit     = [6 17];       % for robust fitting
    elseif consvfit
        f_alphalimit     = [8 13];       % for conservative fitting
    else
        f_alphalimit     = [];           % for balance fitting
    end

    %-- get IAF
    model_spctr = zeros(size(data_spctr));
    model_gof   = zeros(height(channels),1);
    model_iaf   = zeros(height(channels),1);
    parfor elec = 1:height(channels)
        data_fit   = data_spctr(:,:,elec);
        data_base  = data_spctr_off(:,:,elec);
        
        %-- fit alpha with selected trials to estimate peak alpha frequency
        [~,~,alpha_freq,~,model,~,~,~,R2] = ...
            ecog_fitalpha(f,f_alpha4fit,f_alphalimit,data_base',data_fit',false,allowbetafit);
        
        model_iaf(elec)       = 10.^alpha_freq; % [Hz] allow +-1 range
        model_spctr(:,:,elec) = 10.^model(:,1);
        model_gof(elec)       = R2;
    end
    
    %-- Output
    channels.iaf          = model_iaf; % [Hz] allow +-1 range
    channels.gof_alpha    = model_gof; % max 100%
    fitness{ii}.subject   = subject;
    fitness{ii}.channels  = channels;
    fitness{ii}.events    = events;

    %%% Plot
    if opts.doplots
    %-- trials for plot
    trials = [];
    trials.pwrspctrm    = permute(cat(2,data_spctr,data_spctr_off,model_spctr),[3,2,1]);
    trials.pwrspctrm(:,4,:)    = trials.pwrspctrm(:,2,:) ./ trials.pwrspctrm(:,1,:);
    trials.pwrspctrm    = trials.pwrspctrm(:,:,f_alpha4fit);
    trials.f            = f(f_alpha4fit);
    trials.events       = repmat(events,4,1);
    trials.channels     = channels;

    %-- Event list
    trials.events.trial_type(2) = 0;
    trials.events.trial_name{2} = 'BLANK';
    trials.events.trial_type(4) = max(trials.events.trial_type)+1;
    trials.events.trial_name{4} = 'DIFFERENCE';
    trials.events.trial_type(3) = max(trials.events.trial_type)+1;
    trials.events.trial_name{3} = 'Model';
    trials.events.trial_name    = cellfun(@upper, trials.events.trial_name, 'UniformOutput', false);
    eventList   = trials.events(:,{'trial_type','trial_name'});
    eventList   = sortrows(eventList,'trial_type');
    eventList   = unique(eventList.trial_name,'stable');

    %%% plot IAF
    specs = [];
    specs.plot      = opts.plot;
    specs.plot.XLim = [min(f_alpha4fit) max(f_alpha4fit)];
    if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
    else,         	whichElectrodes = trials.channels.name;
    end
    maxplchan       = 24;
    numpltchan      = length(whichElectrodes);
    npltfig         = ceil(numpltchan./maxplchan);
    nSubPlots       = ceil(arrangeinrect(numpltchan,npltfig./1.2,npltfig./[1.0,1.8],'silent')./[npltfig,1]);
    specs.plot.nSubPlots  = nSubPlots;
    specs.plot.showlegend = setlegendpos(opts.plot.showlegend,numpltchan,nSubPlots);

    [~,hF] = ecog_plotGridSpectra(trials, whichElectrodes, eventList(1:2),[], specs);
    figureName = sprintf('%s_%s_IAF', fileid,subject);
    for ihF = hF
    %-- show IAF
    showIAF(ihF,channels,f_alpha4fit);
    %-- Save
    if length(hF) > 1, fignum = sprintf('-%d',find(ismember(hF,ihF)));  else, fignum = ''; end
    set(ihF,'Name',figureName);
    hgexport(ihF, fullfile(opts.plotsavedir, sprintf('%s%s',figureName,fignum)), hgexport('factorystyle'), 'Format', 'png'); % close;
    end
    h = [h hF];

    %-- HD grid
    if hasHDgrid
        specs.plot.nSubPlots  = [];
        specs.plot.showlegend = 'Outside';
        [~,hF] = ecog_plotGridSpectra(trials, 'GB', eventList(1:2),[], specs);
        figureName = sprintf('%s_%s_IAF-GB', fileid,subject);
        for ihF = hF
        %-- show IAF
        showIAF(ihF,channels,f_alpha4fit);
        %-- Save
        if length(hF) > 1, fignum = sprintf('-%d',find(ismember(hF,ihF)));  else, fignum = ''; end
        set(ihF,'Name',figureName);
        hgexport(ihF, fullfile(opts.plotsavedir, sprintf('%s%s',figureName,fignum)), hgexport('factorystyle'), 'Format', 'png'); % close;
        end
        h = [h hF];
    end
    
    %%% plot Difference
    if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
    else,         	whichElectrodes = trials.channels.name;
    end
    specs.plot.nSubPlots  = nSubPlots;
    specs.plot.showlegend = setlegendpos(opts.plot.showlegend,numpltchan,nSubPlots);
    [~,hF] = ecog_plotGridSpectra(trials, whichElectrodes, eventList(3:4),[], specs);
    figureName = sprintf('%s_%s_IAFdiff', fileid,subject);
    for ihF = hF
    %-- show IAF
    showIAF(ihF,channels,f_alpha4fit);
    %-- Save
    if length(hF) > 1, fignum = sprintf('-%d',find(ismember(hF,ihF)));  else, fignum = ''; end
    set(ihF,'Name',figureName);
    hgexport(ihF, fullfile(opts.plotsavedir, sprintf('%s%s',figureName,fignum)), hgexport('factorystyle'), 'Format', 'png'); % close;
    end
    h = [h hF];

    %-- HD grid
    if hasHDgrid
        specs.plot.nSubPlots  = [];
        specs.plot.showlegend = 'Outside';
        [~,hF] = ecog_plotGridSpectra(trials, 'GB', eventList(3:4),[], specs);
        figureName = sprintf('%s_%s_IAFdiff-GB', fileid,subject);
        for ihF = hF
        %-- show IAF
        showIAF(ihF,channels,f_alpha4fit);
        %-- Save
        if length(hF) > 1, fignum = sprintf('-%d',find(ismember(hF,ihF)));  else, fignum = ''; end
        set(ihF,'Name',figureName);
        hgexport(ihF, fullfile(opts.plotsavedir, sprintf('%s%s',figureName,fignum)), hgexport('factorystyle'), 'Format', 'png'); % close;
        end
        h = [h hF];
    end
    end
end
end

%% Sub function
function postfix = cnstpostfix(fileid,allowlag,allowbetafit,allowwidefit)
%%-- lag
if allowlag
    if allowlag && ~contains(fileid,'regresslag')
        fileid = sprintf('%s_regresslag',fileid);
    elseif ~contains(fileid,'regress')
        fileid = sprintf('%s_regress',fileid);
    end
end
%%-- parameters 
fitparams = '';
addbeta = allowbetafit && ~contains(fileid,{'_beta'});
addwide = allowwidefit && ~contains(fileid,{'wide_'}) && ~endsWith(fileid,{'wide'}) ;
if addbeta,         fitparams = sprintf('%sbeta',fitparams); end
if addwide,         fitparams = sprintf('%swide',fitparams); end
fitparams    = regexprep(fitparams,'^-*','');
if ~isempty(fitparams), fitparams = ['_' fitparams]; end
postfix      = sprintf('%s%s',fileid,fitparams);
end

function param = getparam(opt,paramname,idx,nidx)

param = opt.(paramname);
if length(param)==nidx
    param = param(idx);
elseif length(param)~=1
    error('%s.%s must be a boolean or have the same length to the 1st argument',inputname(1),paramname);
end
end

function showlegend = setlegendpos(showlegend,numpltchan,nSubPlots)

if isempty(showlegend)
    if numpltchan<=12,  showlegend = 'Inside';
    elseif numpltchan<=23 && numpltchan~=prod(nSubPlots)
                        showlegend = 'Last';
    else,               showlegend = 'best';
    end
end
end
    
function  showIAF(hf,channels,f_alpha4fit)
%-- Fix legend
set(findobj(hf,'Type','Legend'),'AutoUpdate','off');  % do not allow update legend
%-- Get Axis
hax = findobj(hf,'Type','Axes');
for iax = 1:length(hax)
    %-- Get channel number
    channelname = strtok(hax(iax).Title.String);
    channelidx  = find(ismember(channels.name,channelname),1);
    if ~isempty(channelname) && isempty(channelidx)
        [eleccat, elecnum] = strtok(channelname,int2str(0:9));
        channelname = {sprintf('%s%d',eleccat,str2double(elecnum)),...
                       sprintf('%s%02d',eleccat,str2double(elecnum)),...
                       sprintf('%s%03d',eleccat,str2double(elecnum)),...
                       sprintf('%s%04d',eleccat,str2double(elecnum))};
        channelidx  = find(ismember(channels.name,channelname),1);
    end
    if ~isempty(channelidx)
        %-- Plot IAF
        plot(hax(iax),[1 1].*channels.iaf(channelidx),ylim(hax(iax)),'k--','LineWidth',3);
    end
    %-- Simplify xtics
    xticks(hax(iax),unique([[5:5:20] minmax(reshape(f_alpha4fit,1,[]))]));
    xticklabels(hax(iax),arrayfun(@int2str,xticks(hax(iax)),'UniformOutput',false));    % update xticklabels
end
if length(hax)>40
    xlabel(hax,''); ylabel(hax,'');
end
end