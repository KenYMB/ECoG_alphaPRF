function [data] = ecog_prf_selectData(data, stimNames, opts)

% Description
%
% [data] = ecog_prf_selectData(data, [stimNames], [opts])
% 
% Outputs reduced version of data after following steps:
%
% Removes bad runs and channels with many bad epochs (epochOpts)
% Removes channels that do not match inclusion criteria (elecOpts)
% Converts to percent signal change (using baselineTime)
% Averages across trials (make optional?)
% Normalizes by max (optional)
% Averages across visual areas (optional)
% Makes plots (optional)
% 
% Input
% - data
%   - epochs     = t x events x channels
%   - events
%   - channels
% - stimNames     = cell-strings array    % trial names for analysis
% - opts
%   - baseline_time         % time period across which to compute normalization baseline
%   - stim_on_time          % time period across which stimulus evoked response is expected
%   - epoch_time            % time period across which are considered in data selection
%   - epoch_jump_thresh     % max jump in voltage allowed within stim_on_time period
%   - epoch_dist_thresh     % x-fold max magnitude above which epoch is ignored to estimate a distribution
%   - epoch_outlier_thresh  % percentile of a distribution in which epoch will be labeled as outlier
%   - elec_selection_method % variance(default), splithalf, meanpredict
%                             see also ecog_selectElectrodes except for 'variance'
%     - elec_thresh_var         % required for 'variance'
%                                 x-fold variance above which electrode will be labeled as outlier
%     - elec_splitmethod        % required for 'splithalf'
%                                 alternative(default), twohalves,
%                                 alternativeINsessions, twohalvesINsessions,
%                                 alternativeINruns, twohalvesINruns
%     - elec_splithalf_thresh   % required for 'splithalf'
%                                 minimum required R2 between split halves of data
%     - elec_meanpredict_thresh % required for 'meanpredict'
%                                 minimum required R2 for prediction by mean (1 - (SSEresidual/SSEtotal)
%   - elec_exclude_depth    % true(default) / false
%   - issave                % true / false(default)
%   - outputDir
%   - doplots
%   - plotsavedir
%   - plot
%   ---------------------
%   - fileid
% 
% Output
% - data

% Hidden options
% - opts
%   - compute       = [];
%   - skipexist     = [];
%   - fileid        = 'data_visualelecs';
%   - epoch_regress
%      - allowlag      = true or false(default), if allow lag in regression
%       - maxlag        = maximum lags (default = 0.01s, valid if allowlag is true)
%       - lag_maxabs    = true or false(default), if estimate time lag with maxabs or max

% To load selected data
%   use ecog_prf_getData as following:
% 
%   compute = false;
%   outputDir = 'Preprocessed';
%   [data] = ecog_prf_getData(compute, [], outputDir, subjectList);

% Dependency: ecog_regressout, SetDefault, cellstrfind, nanrms, copyfields

% 20200127 - Yuasa: modified from tde_selectData in temporalECoG
% 20220215 - Yuasa: enable to set time range to check epochs
% 20220222 - Yuasa: load HDgrid information from subjectlist.tsv

%% Set options
%--Define inputs 
% <opts>
SetDefaultAnalysisPath('DATA','Preprocessed','opts.outputDir');
SetDefaultAnalysisPath('FIGURE','dataselection','opts.plotsavedir');
SetDefault('opts.baseline_time',[-0.2 0]);
SetDefault('opts.stim_on_time',[0 0.5]);
SetDefault('opts.issave',false);
SetDefault('opts.epoch_time',[-0.2 0.8]);
SetDefault('opts.epoch_jump_thresh',250);
SetDefault('opts.epoch_dist_thresh',20);
SetDefault('opts.epoch_outlier_thresh',0.01*1e-2);
SetDefault('opts.elec_selection_method','variance');
SetDefault('opts.elec_thresh_var',3);
SetDefault('opts.elec_splitmethod','alternative');
SetDefault('opts.elec_splithalf_thresh',0.22);
SetDefault('opts.elec_meanpredict_thresh',0);
SetDefault('opts.elec_exclude_depth',true);
SetDefault('opts.doplots',false);
SetDefault('opts.plot.XLim',opts.epoch_time);
SetDefault('opts.plot.fontSize',16);
% <hidden opts>
SetDefault('opts.compute',[]);
SetDefault('opts.skipexist',[]);
SetDefault('opts.fileid','data_visualelecs');
SetDefault('opts.epoch_regress.allowlag',false);
SetDefault('opts.epoch_regress.maxlag',0.01);
SetDefault('opts.epoch_regress.lag_maxabs',false);
SetDefault('opts.plot.RotGrid',false);

% <stimNames>
SetDefault('stimNames',{},'cell');
if isempty(stimNames) || isempty(stimNames{1}) || strcmpi(stimNames{1},'all')
    stimNames = {'*'};
end

%% Initialize
%-- check inputs and outputs
assert(~isempty(data), 'Please provide the data struct');
if isempty(opts.compute)
    SetDefault('opts.skipexist',true);
    opts.compute = true;
elseif ~opts.compute
    opts.skipexist = false;
else
    SetDefault('opts.skipexist',false);
end
if opts.issave && ~exist(opts.outputDir, 'dir'),  	mkdir(opts.outputDir); end
if opts.doplots && ~exist(opts.plotsavedir, 'dir'), mkdir(opts.plotsavedir); end
if opts.doplots
  if ~exist(fullfile(opts.plotsavedir, 'preselection'), 'dir');       mkdir(fullfile(opts.plotsavedir, 'preselection'));       end
  if ~exist(fullfile(opts.plotsavedir, 'epochselection'), 'dir');     mkdir(fullfile(opts.plotsavedir, 'epochselection'));     end
  if ~exist(fullfile(opts.plotsavedir, 'electrodeselection'), 'dir'); mkdir(fullfile(opts.plotsavedir, 'electrodeselection')); end
  if ~exist(fullfile(opts.plotsavedir, 'runselection'), 'dir');       mkdir(fullfile(opts.plotsavedir, 'runselection'));       end
end

%-- check to save
isSave      = opts.issave;

%% Loop over subjects
for ii = 1:length(data) % Loop over subjects
    %-- Set data
    if ~opts.compute && ~isstruct(data{ii}) && ischar(data{ii})
        subject     = data{ii};
        isloadfile  = true;
        compute     = false;
    else
        idat        = data{ii};
        subject     = idat.subject;
        epochs      = idat.epochs;
        t           = idat.t;
        channels    = idat.channels;
        events      = idat.events;
        fsample     = idat.fsample;
        isloadfile  = opts.skipexist;
        compute     = opts.compute;
    end
    %-- Try to load files
    if isloadfile
        filename    = fullfile(opts.outputDir, sprintf('%s_%s.mat', subject, opts.fileid));
        %-- load files from directory
        if exist(filename,'file') || ~compute
            fprintf('[%s] Loading selected data for subject %s from %s ',mfilename, subject, opts.outputDir);
            idat        = load(filename);
            fprintf('\n');
            compute = false;
        end
    end
    
    %-- Main
    if compute
        
%% STEP 0-1 Restrict selection
    fprintf('[%s] Selecting data for subject %s \n',mfilename, subject);
    
    %-- HD grid flag
    hasHDgrid = istablefield(channels,'group') && any(ismember(channels.group,'HDgrid'));
    
    %-- Restrict selection to relevant stimuli only
    stimsForSelection = cellstrfind(events.trial_name, stimNames,1);
    epochs = epochs(:, stimsForSelection, :);
    events = events(stimsForSelection, :);
    
    %-- Restrict selection to included channels only
    %--- Exclude depth electrodes
    if opts.elec_exclude_depth 
        select_chs_pre = ~contains(lower(channels.type), 'seeg');
        channels = channels(select_chs_pre,:);
        epochs   = epochs(:,:,select_chs_pre);
    end  
    %--- Exclude runs based on experiment note
    runList = grp2idx(cellfun(@(x,y) cat(2,x,y),events.session_name,events.run_name,'UniformOutput',false));
    select_runs_pre = true(numel(unique(runList)),1);
    switch subject
        case {'p05'} 
            %---- exclude fiest 2 runs because of broke fixation
            select_runs_pre([1,2]) = false;
    end
    select_idx  = ismember(runList,find(select_runs_pre));
    epochs      = epochs(:,select_idx,:);
    events      = events(select_idx,:);
    runList     = runList(select_idx);
    
    %-- Get sessions x runs list
    runSet   = unique(runList);
    nruns    = numel(runSet);
    trlprun  = [];
    for kk = 1:length(runSet)
        trlprun(kk)  = numel(find(runList==runSet(kk)));
    end
    
    %-- Set 0 for prf BLANK, 1 for prf stimuli, and 10+trial_type for other stimuli
    stmlist = double(ismember(events.task_name,'prf')&~ismember(events.trial_name,'BLANK')) + ...
              ~ismember(events.task_name,'prf').*(events.trial_type+10); 
          
    %-- Copy events and change trial names for pRF
    events_prf  = events;
    events_prf.trial_name(stmlist==1) = {'PRF'};
    events_prf.isprf    = ismember(events_prf.task_name,'prf');
    eventList   = table(events_prf.trial_name,stmlist,'VariableNames',{'trial_name','stmlist'});
    eventList   = sortrows(eventList,'stmlist');
    eventList   = unique(eventList.trial_name,'stable');
          
%% STEP 0-2 Baseline correction
    fprintf('[%s] Apply baseline correction...\n',mfilename);
    epochs = ecog_normalizeEpochs(epochs, t, opts.baseline_time, 'subtractwithintrial');

    %-- Plot runs and electrodes
    if opts.doplots
        for kk = unique(runList(:)')
            runsidx   = runList==kk;
            %%-- Make trials & Set trial name 'PRF' & set specs
            trials          = [];
            trials.evoked   = permute(epochs(:,runsidx,:),[3,1,2]); % channels x t x events
            trials.time     = t;
            trials.channels = channels;
            trials.events   = events_prf(runsidx,:);
            
            %%-- Set specs
            specs = [];
            specs.dataTypes = {'evoked'};
            specs.plot      = opts.plot;

            if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
            else,         	whichElectrodes = trials.channels.name;
            end

            %%-- Plot averaged evoked potential for every runs
            figureName = sprintf('evoked-prf_%s_all', subject);

            ecog_plotGridTimecourses(trials, whichElectrodes, eventList, specs);
            hgexport(gcf, fullfile(opts.plotsavedir,'preselection', sprintf('%s-run%02d',figureName,kk)), hgexport('factorystyle'), 'Format', 'png'); close;
              if hasHDgrid
              [~,p]=ecog_plotGridTimecourses(trials, 'GB', eventList, specs);
              for pp=1:length(p)
                if length(p)==1,  hgexport(p(pp), fullfile(opts.plotsavedir,'preselection', sprintf('%s-GB-run%02d',figureName,kk)), hgexport('factorystyle'), 'Format', 'png');
                else,             hgexport(p(pp), fullfile(opts.plotsavedir,'preselection', sprintf('%s-GB-run%02d-%02d',figureName,kk,pp)), hgexport('factorystyle'), 'Format', 'png');
                end
              end
              close(p);
              end
        end
    end
    
%% STEP 1 Select epochs
    fprintf('[%s] Removing bad epochs...\n',mfilename);
    
    %-- Find bad epochs
    sessionList = grp2idx(events.session_name);
    nsession    = max(sessionList);
    outliers        = cell(nsession,1);
    max_powers      = cell(nsession,1);
    outlier_thresh  = cell(nsession,1);
    for mm = 1:nsession
        sesstrl = sessionList == mm;
        opts2 = copyfields(opts, [], {'baseline_time','stim_on_time','epoch_time','epoch_jump_thresh','epoch_dist_thresh','epoch_outlier_thresh','epoch_regress'});
        opts2.stims = stmlist(sesstrl);
        [~, outliers{mm}, max_powers{mm}, outlier_thresh{mm}] = ecog_selectEpochsStat(epochs(:,sesstrl,:), t, fsample, opts2);
    end
    
    %-- Plot bad epochs
    if opts.doplots
        for jj = 1:height(channels)
            for mm  = 1:nsession
              if any(outliers{mm}(:,jj))
                maxpnlspfig = 56;       % maximum subplots per figure
                outliers_list = find(outliers{mm}(:,jj)) + size(vertcat(outliers{1:(mm-1)}),1);
                nOutliersAll  = length(outliers_list);
                nfig = ceil(nOutliersAll./(maxpnlspfig-1));
                pnlspfig = ceil(nOutliersAll./nfig);
                for ll = 1:nfig
                    if nfig>1,  figpostfix = sprintf('_%d',ll);  pltpnls = [1:pnlspfig]+pnlspfig*(ll-1);
                    else        figpostfix = '';                 pltpnls = 1:(maxpnlspfig-1);
                    end
                    pltpnls(pltpnls > nOutliersAll) = [];
                    figureName = sprintf('outlierepochs_%s_%s_sess%02d%s', subject, channels.name{jj},mm,figpostfix);
                    figure('Name', figureName); hold on;
                    outliers_found = outliers_list(pltpnls);
                    nOutliers = length(outliers_found);
                    dim1 = round(sqrt(nOutliers+1)); dim2 = ceil((nOutliers+1)/dim1);
                    subplot(dim2,dim1,1); hold on; title(channels.name{jj});
                    histogram(max_powers{mm}(:,jj),100); line([outlier_thresh{mm}(jj) outlier_thresh{mm}(jj)], get(gca, 'YLim'), 'Color', 'r','LineStyle', ':', 'LineWidth', 2);
                    set(gca, 'fontsize', 14); xlabel('max powers'); ylabel('number of epochs');
                    for kk = 1:nOutliers
                        subplot(dim2,dim1,kk+1);
                        ecog_plotSingleTimeCourse(t, epochs(:,outliers_found(kk),jj), [], [], sprintf('epoch %d %s', outliers_found(kk), events.trial_name{outliers_found(kk)}));
                    end
                    set(gcf, 'Position', [150 100 300*dim1 300*dim2]);
                    hgexport(gcf, fullfile(opts.plotsavedir,'epochselection', figureName), hgexport('factorystyle'), 'Format', 'png'); close;
                end
              end
            end
        end
    end
    
    %-- concatenate outliers
    outliers = vertcat(outliers{:});
    
    %%% Exclude runs where 80% epochs were marked as bad
    for kk=1:length(runSet)
        outliers_runs = sum(outliers(runList==runSet(kk),:),1);
        outliers(runList==kk,outliers_runs>(trlprun(kk)*0.8)) = true;
    end
    
    %-- fill NaN for bad epochs
    epochs(:,outliers)    = nan;
    
%% STEP 2 Select electrodes   
    fprintf('[%s] Selecting electrodes...\n',mfilename);
    
    %-- Initialize selection to include all channels
    select_chs = ones(height(channels),3);
    
    %%% Exclude channels with many bad runs: (where at least 80% runs were rejected)
    thresh = 0.8;
    outliers_runs = zeros(nruns,height(channels));
    for kk=1:length(runSet)
        outliers_runs(kk,:) = sum(outliers(runList==runSet(kk),:),1);
    end
    select_chs(:,1) = sum(outliers_runs<reshape(trlprun,[],1),1)>(nruns.*(1-thresh));
    
    %%% Exclude channels with many bad epochs: (where at least 30% epochs were rejected in every runs: to reconstruct pRF time course)
    thresh = 0.3;
    select_chs(:,2) = any(outliers_runs<=reshape(trlprun,[],1)*thresh,1);
    
    %%% Exclude channels based on variance:
    switch opts.elec_selection_method
        case {'variance'}
        %%-- Compute variance across trials at each time point, and take average of the variance across times
            mean_resp_var   = mean(var(epochs,0,2,'omitnan'),1,'omitnan');
            select_chs(:,3) = mean_resp_var < median(mean_resp_var,'omitnan')+std(mean_resp_var,'omitnan')*opts.elec_thresh_var;
        otherwise
        %%-- ecog_selectElectrodes
            opts2 = copyfields(opts, [], {'elec_selection_method','elec_splitmethod','elec_max_thresh','elec_mean_thresh','elec_splithalf_thresh','elec_meanpredict_thresh'});
            opts2.stim_on   = opts.stim_on_time;
            switch opts.elec_selection_method
                case {'splithalf'}
            opts2.stimnames = setdiff(events_prf.trial_name,'BLANK','stable');
            [select_chs(:,3), R2, epochs_split] = ecog_selectElectrodes(epochs, channels, events_prf, t, opts2);
                otherwise      % case {'meanpredict'}
            opts2.stimnames = unique(events.trial_name,'stable');
            [select_chs(:,3), R2, epochs_split] = ecog_selectElectrodes(epochs, channels, events, t, opts2);
            end
            channels.noiseceilingR2 = round(R2,2);
    end

    %-- Combine criteria
    select_idx = all(select_chs,2);

    if opts.doplots
        nEl = height(channels); 
        
        %-- plot xval selection 
        if exist('epochs_split','var')&&~isempty(epochs_split)
            figureName = sprintf('halfepochs_%s_all', subject);
            figure('Name', figureName); plotDim1 = round(sqrt(nEl)); plotDim2 = ceil((nEl)/plotDim1);
            for el = 1:nEl
                subplot(plotDim1,plotDim2,el); hold on
                plot(squeeze(epochs_split(1,el,:)), 'r','LineWidth', 2);
                plot(squeeze(epochs_split(2,el,:)), 'b','LineWidth', 2);
                axis tight
                nSamp = size(epochs,1); nSampTot = nSamp * length(opts2.stimnames);
                set(gca, 'XTick', 1:nSamp:nSampTot, 'XTickLabel', opts2.stimnames);
                xtickangle(45)
                plotTitle = sprintf('%s %s %s R2 = %0.2f', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el}, R2(el));
                title(plotTitle);
            end
            set(gcf, 'Position', [150 100 1500 1250]);
            set(findall(gcf,'-property','FontSize'),'FontSize',14);
            hgexport(gcf, fullfile(opts.plotsavedir,'electrodeselection', figureName), hgexport('factorystyle'), 'Format', 'png'); close;
        end
            
            
        %-- Compute mean across all trials
        mean_resp = mean(epochs,2, 'omitnan');
        %-- Compute standard errors for plotting
        llim = (mean_resp - (std(epochs,0,2,'omitnan')./sqrt(sum(~isnan(epochs),2))));
        ulim = (mean_resp + (std(epochs,0,2,'omitnan')./sqrt(sum(~isnan(epochs),2))));
        mean_resp_sd = cat(2, llim, ulim);
        %-- Plot all channels
        figureName = sprintf('viselec_%s_all', subject);
        figure('Name', figureName); plotDim1 = round(sqrt(nEl)); plotDim2 = ceil((nEl)/plotDim1);
        for el = 1:nEl
            subplot(plotDim1,plotDim2,el); hold on
            plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        
            ecog_plotSingleTimeCourse(t, mean_resp(:,:,el), squeeze(mean_resp_sd(:,:,el)), [], plotTitle);
            %if el == 1; xlabel('Time (s)'); ylabel('Broadband signal change');end
        end
        set(gcf, 'Position', [150 100 1500 1250]);
        hgexport(gcf, fullfile(opts.plotsavedir,'electrodeselection', figureName), hgexport('factorystyle'), 'Format', 'png'); close;
        %-- Plot selected channels
        figureName = sprintf('viselec_%s_selected', subject);
        figure('Name', figureName); 
        for el = 1:nEl
            if select_idx(el)
                subplot(plotDim1,plotDim2,el); hold on
                plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        
                ecog_plotSingleTimeCourse(t, mean_resp(:,el), squeeze(mean_resp_sd(:,:,el)), [], plotTitle)    
                set(gcf, 'Position', [150 100 1500 1250]);
            end
        end
        hgexport(gcf, fullfile(opts.plotsavedir,'electrodeselection', figureName), hgexport('factorystyle'), 'Format', 'png'); close;
    end
    
    %-- Update datas
    epochs    = epochs(:,:,select_idx);
    channels  = channels(select_idx,:);
    outliers  = outliers(:,select_idx);
    
%% STEP 3 Select runs   
    fprintf('[%s] Selecting runs...\n',mfilename);
    
    %-- Initialize selection to include all epochs
    select_runs = ones(nruns,1);
    
    %%% Exclude runs with too much bad channels: (where at least 30% channels were rejected)
    thresh = 0.3;
    outliers_runs = zeros(nruns,height(channels));
    for kk=1:length(runSet)
        outliers_runs(kk,:) = sum(outliers(runList==runSet(kk),:),1);
    end
    select_runs(:,1) = sum(outliers_runs<reshape(trlprun,[],1),2) > (height(channels)*(1-thresh));
    
    %%% Combine all selection
    select_runs = all(select_runs,2);
    
    %-- Plot runs and electrodes
    if opts.doplots
        for kk = reshape(find(~select_runs),1,[])
            runsidx   = runList==runSet(kk);
            %%-- Make trials & Set trial name 'PRF' & set specs
            trials          = [];
            trials.evoked   = permute(epochs(:,runsidx,:),[3,1,2]); % channels x t x events
            trials.time     = t;
            trials.channels = channels;
            trials.events   = events_prf(runsidx,:);
            
            %%-- Set specs
            specs = [];
            specs.dataTypes = {'evoked'};
            specs.plot      = opts.plot;

            if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
            else,         	whichElectrodes = trials.channels.name;
            end

            %%-- Plot averaged evoked potential for every runs
            figureName = sprintf('evoked-prf_%s_outlierruns', subject);

            ecog_plotGridTimecourses(trials, whichElectrodes, eventList, specs);
            hgexport(gcf, fullfile(opts.plotsavedir,'runselection',sprintf('%s-run%02d',figureName,runSet(kk))), hgexport('factorystyle'), 'Format', 'png'); close;
              if hasHDgrid
              [~,p]=ecog_plotGridTimecourses(trials, 'GB', eventList, specs);
              for pp=1:length(p)
                if length(p)==1,  hgexport(p(pp), fullfile(opts.plotsavedir,'runselection', sprintf('%s-GB-run%02d',figureName,kk)), hgexport('factorystyle'), 'Format', 'png');
                else,             hgexport(p(pp), fullfile(opts.plotsavedir,'runselection', sprintf('%s-GB-run%02d-%02d',figureName,kk,pp)), hgexport('factorystyle'), 'Format', 'png');
                end
              end
              close(p);
              end
        end
    end
    
    %-- Update datas
    select_idx  = ismember(runList,runSet(select_runs));
    epochs      = epochs(:,select_idx,:);
    events      = events(select_idx,:);
    events_prf  = events_prf(select_idx,:);
    outliers    = outliers(select_idx,:);
    runList     = runList(select_idx);
    
%% STEP 4 Finalize processing
    %-- Collect selected information
    %-- Collect into an output struct
   idat.epochs   = epochs;
   idat.events   = events;
   idat.channels = channels;
      %-- Collect selected information
     idat.outliers        = outliers;
     idat.select_chs      = false(length(select_chs_pre),size(select_chs,2)+1);
     idat.select_chs(:,1) = select_chs_pre;
     idat.select_chs(select_chs_pre,2:end) = select_chs;
     idat.select_run      = false(length(select_runs_pre),size(select_runs,2)+1);
     idat.select_run(:,1) = select_runs_pre;
     idat.select_run(select_runs_pre,2:end) = select_runs;
    %-- Save out the data
    if isSave && ~isempty(epochs)
        fprintf('[%s] Saving selected data for subject %s to %s \n',mfilename, subject, opts.outputDir);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s.mat', subject, opts.fileid));
        saveauto(filename,'-struct','idat')
    end
    
    %-- Plot runs and electrodes
    if opts.doplots
        for kk = unique(runList(:)')
            runsidx   = runList==kk;
            %%-- Make trials & Set trial name 'PRF' & set specs
            trials          = [];
            trials.evoked   = permute(epochs(:,runsidx,:),[3,1,2]); % channels x t x events
            trials.time     = t;
            trials.channels = channels;
            trials.events   = events_prf(runsidx,:);
            
            %%-- Set specs
            specs = [];
            specs.dataTypes = {'evoked'};
            specs.plot      = opts.plot;

            if hasHDgrid,  	whichElectrodes = trials.channels.name(~contains(trials.channels.name,'GB'));
            else,         	whichElectrodes = trials.channels.name;
            end

            %%-- Plot averaged evoked potential for every runs
            figureName = sprintf('evoked-prf_%s_selected', subject);

            ecog_plotGridTimecourses(trials, whichElectrodes, eventList, specs);
            hgexport(gcf, fullfile(opts.plotsavedir, sprintf('%s-run%02d',figureName,kk)), hgexport('factorystyle'), 'Format', 'png'); close;
              if hasHDgrid
              [~,p]=ecog_plotGridTimecourses(trials, 'GB', eventList, specs);
              for pp=1:length(p)
                if length(p)==1,  hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB-run%02d',figureName,kk)), hgexport('factorystyle'), 'Format', 'png');
                else,             hgexport(p(pp), fullfile(opts.plotsavedir, sprintf('%s-GB-run%02d-%02d',figureName,kk,pp)), hgexport('factorystyle'), 'Format', 'png');
                end
              end
              close(p);
              end
        end
    end
    
    end  % skip everything when compute=false
    
    %-- Return data
    data{ii} = idat;
    
end
fprintf('[%s] Done! \n',mfilename);

end
