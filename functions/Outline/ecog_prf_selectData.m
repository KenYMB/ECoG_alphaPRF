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
%   - baseline_time     % time period across which to compute normalization baseline
%   - stim_on           % time period across which stimulus is presented
%   - 
%   - 
%   - 
%   - epoch_thresh_pow  % x-fold max magnitude above which epoch is ignored to estimate a distribution
%   - epoch_thresh_dist % percent tail of a distribution in which epoch will be labeled as outlier
%   - run_thresh_sn     % S/N below which run will be labeled as outlier
%   - elec_thresh_var   % x-fold variance above which electrode will be labeled as outlier
%   - elec_exclude_depth
%   - issave
%   - outputDir
%   - doplots
%   - plotsavedir
%   - plot
%   ---------------------
%   - fileid
% 
% Output
% - data

% To load selected data
%   use ecog_prf_getData as following:
% 
%   compute = false;
%   outputDir = fullfile(analysisRootPath, 'Data', 'Preprocessed');
%   [data] = ecog_prf_getData(compute, [], outputDir, subjectList);

% Dependency: ecog_regressout, SetDefault, cellstrfind, nanrms

% 20200127 - Yuasa: modified from tde_selectData in temporalECoG
% 20220215 - Yuasa: enable to set time range to check epochs
% 20220222 - Yuasa: load HDgrid information from subjectlist.tsv

%% Set options
%--Define inputs 
% <opts>
SetDefault('opts.baseline_time',[-0.2 0]);
SetDefault('opts.stim_on',[0 0.5]);
SetDefault('opts.issave',false);
SetDefault('opts.outputDir',fullfile(analysisRootPath, 'Data', 'Preprocessed'));
SetDefault('opts.epoch_thresh_pow',20);
SetDefault('opts.epoch_thresh_dist',1e-10);
SetDefault('opts.epoch_thresh_time',[]);
SetDefault('opts.run_thresh_sn',0.1);
SetDefault('opts.elec_thresh_var',10);
SetDefault('opts.elec_exclude_depth',false);
SetDefault('opts.doplots',false);
SetDefault('opts.plotsavedir',fullfile(analysisRootPath, 'Figures', 'dataselection'));
SetDefault('opts.plot.XLim',[-0.2 0.8]);
SetDefault('opts.plot.fontSize',16);
% <hidden opts>
SetDefault('opts.plot.RotGrid',false);
SetDefault('opts.fileid','data_visualelecs');

% <stimNames>
SetDefault('stimNames',{},'cell');
if isempty(stimNames) || isempty(stimNames{1}) || strcmpi(stimNames{1},'all')
    stimNames = {'*'};
end

%% Initialize
if ~exist(fullfile(opts.plotsavedir, 'preselection'), 'dir');       mkdir(fullfile(opts.plotsavedir, 'preselection'));       end
if ~exist(fullfile(opts.plotsavedir, 'epochselection'), 'dir');     mkdir(fullfile(opts.plotsavedir, 'epochselection'));     end
if ~exist(fullfile(opts.plotsavedir, 'electrodeselection'), 'dir'); mkdir(fullfile(opts.plotsavedir, 'electrodeselection')); end
if ~exist(fullfile(opts.plotsavedir, 'runselection'), 'dir');       mkdir(fullfile(opts.plotsavedir, 'runselection'));       end

%-- check inputs and outputs
assert(~isempty(data), 'Please provide the data struct');
if opts.issave && ~exist(opts.outputDir, 'dir'),  	mkdir(opts.outputDir); end
if opts.doplots && ~exist(opts.plotsavedir, 'dir'), mkdir(opts.plotsavedir); end

%-- check to save
isSave      = opts.issave;

%-- Correct Subject Information
subjectList_fname = 'subjectlist.tsv';
SbjInfo    = loadSbjInfo(subjectList_fname,'all');
hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');

%% Loop over subjects
for ii = 1:length(data) % Loop over subjects
    
    subject     = data{ii}.subject;
    epochs      = data{ii}.epochs;
	t           = data{ii}.t;
    channels    = data{ii}.channels;
    events      = data{ii}.events;
    fsample     = data{ii}.fsample;

    fprintf('[%s] Selecting data for subject %s \n',mfilename, subject);
    
    %-- HD grid flag
    if hasSbjInfo && istablefield(SbjInfo,'hasHDgrid')
        hasHDgrid = any(strcmpi(SbjInfo.hasHDgrid(ismember(SbjInfo.participant_id,subject)),'yes'));
    else
        hasHDgrid = false;
    end
          
    %-- Restrict selection to relevant stimuli only
    stimsForSelection = cellstrfind(events.trial_name, stimNames,1);
    epochs = epochs(:, stimsForSelection, :);
    events = events(stimsForSelection, :);
    
    %-- Get sessions x runs list
    runsList = grp2idx(cellfun(@(x,y) cat(2,x,y),events.session_name,events.run_name,'UniformOutput',false));
    nruns    = numel(unique(runsList));
    trlprun  = [];
    for kk = unique(runsList(:)')
        trlprun(kk)  = numel(find(runsList==kk));
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
          
%% STEP 0 Baseline correction
    fprintf('[%s] Apply baseline correction...\n',mfilename);
    stimon_idx = t > opts.stim_on(1) & t <= opts.stim_on(2);  
    bsl_idx  = t >= opts.baseline_time(1) & t < opts.baseline_time(2);
    epochs = bsxfun(@minus, epochs, mean(epochs(bsl_idx,:,:),1));

    %-- Plot runs and electrodes
    if opts.doplots
        for kk = unique(runsList(:)')
            runsidx   = runsList==kk;
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

            ecog_plotTimecourses(trials, whichElectrodes, eventList, specs);
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
    
    %-- Regress evoked responses
    regressed = ecog_regressout(epochs,bsl_idx,stmlist);
    
    %-- Find bad epochs
    [~, outliers, max_powers, outlier_thresh] = ecog_selectEpochsStat(epochs, regressed, t, opts.stim_on, stmlist~=0, opts.epoch_thresh_pow, opts.epoch_thresh_dist, opts.epoch_thresh_time);
    
    %-- Plot bad epochs
    if opts.doplots
        for jj = 1:height(channels)
            if any(outliers(:,jj))
                maxpnlspfig = 56;
                outliers_list = find(outliers(:,jj));
                nOutliersAll  = length(outliers_list);
                nfig = ceil(nOutliersAll./(maxpnlspfig-1));
                pnlspfig = ceil(nOutliersAll./nfig);
                for ll = 1:nfig
                    if nfig>1,  figpostfix = sprintf('_%d',ll);  pltpnls = [1:pnlspfig]+pnlspfig*(ll-1);
                    else        figpostfix = '';                 pltpnls = 1:(maxpnlspfig-1);
                    end
                    pltpnls(pltpnls > nOutliersAll) = [];
                    figureName = sprintf('outlierepochs_%s_%s%s', subject, channels.name{jj},figpostfix);
                    figure('Name', figureName); hold on;
                    outliers_found = outliers_list(pltpnls);
                    nOutliers = length(outliers_found);
                    dim1 = round(sqrt(nOutliers+1)); dim2 = ceil((nOutliers+1)/dim1);
                    subplot(dim2,dim1,1); hold on; title(channels.name{jj});
                    histogram(max_powers(:,jj),100); line([outlier_thresh(jj) outlier_thresh(jj)], get(gca, 'YLim'), 'Color', 'r','LineStyle', ':', 'LineWidth', 2);
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
    
    %%% Exclude other epochs in runs where 80% epochs were marked as bad
    for kk=unique(runsList(:)')
        outliers_runs = sum(outliers(runsList==kk,:),1);
        outliers(runsList==kk,outliers_runs>(trlprun(kk)*0.8)) = true;
    end
    
    %-- fill NaN for bad epochs
    epochs(:,outliers)    = nan;
    regressed(:,outliers) = nan;
    
%% STEP 2 Select electrodes   
    fprintf('[%s] Selecting electrodes...\n',mfilename);
    
    %-- Initialize selection to include all channels
    select_idx = ones(height(channels),4);
    
    %%% Exclude channels with bad runs: (where at least 2 runs were rejected)
    outliers_runs = zeros(nruns,height(channels));
    for kk=unique(runsList(:)')
        outliers_runs(kk,:) = sum(outliers(runsList==kk,:),1);
    end
    select_idx(:,1) = sum(outliers_runs<reshape(trlprun,[],1),1)>(nruns-2);
    
    %%% Exclude channels with many bad epochs: (where at least 20% epochs were rejected in all runs)
    select_idx(:,2) = any(outliers_runs<=reshape(trlprun,[],1)*0.2,1);
    
    %%% Exclude channels based on variance:
    %%-- Compute mean of variance across trials across times
    mean_resp_var = mean(var(epochs,0,2,'omitnan'),1,'omitnan');
    select_idx(:,3) = mean_resp_var < median(mean_resp_var,'omitnan')*opts.elec_thresh_var;

    %%% Exclude depth electrodes
    if opts.elec_exclude_depth, select_idx(:,4) = contains(lower(channels.type), 'ecog'); end  

    %-- Combine criteria
    select_idx = all(select_idx,2);

    if opts.doplots
        %-- Compute mean across all trials
        mean_resp = mean(epochs,2, 'omitnan');
        %-- Compute standard errors for plotting
        llim = (mean_resp - (std(epochs,0,2,'omitnan')./sqrt(sum(~isnan(epochs),2))));
        ulim = (mean_resp + (std(epochs,0,2,'omitnan')./sqrt(sum(~isnan(epochs),2))));
        mean_resp_sd = cat(2, llim, ulim);
        %-- Plot all channels
        nEl = size(mean_resp,3); 
        figureName = sprintf('viselec_%s_all', subject);
        figure('Name', figureName); plotDim1 = round(sqrt(nEl)); plotDim2 = ceil((nEl)/plotDim1);
        for el = 1:nEl
            subplot(plotDim1,plotDim2,el); hold on
            plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        
            ecog_plotSingleTimeCourse(t, mean_resp(:,:,el), squeeze(mean_resp_sd(:,:,el)), [], plotTitle);
            %if el == 1; xlabel('Time (s)'); ylabel('Broadband signal change');end
            set(gcf, 'Position', [150 100 1500 1250]);
        end
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
    regressed = regressed(:,:,select_idx);
    channels  = channels(select_idx,:);
    outliers  = outliers(:,select_idx);
    
%% STEP 3 Select runs   
    fprintf('[%s] Selecting runs...\n',mfilename);
    
    %-- Initialize selection to include all epochs
    select_runs = ones(nruns,3);
    
    %%% Exclude runs with too much bad channels: (where at least 30% channels were rejected)
    outliers_runs = zeros(nruns,height(channels));
    for kk=unique(runsList(:)')
        outliers_runs(kk,:) = sum(outliers(runsList==kk,:),1);
    end
    select_runs(:,1) = sum(outliers_runs<reshape(trlprun,[],1),2) > (height(channels)*0.7);
    
    %%% Exclude runs with too much noise: (only for pRF sessions)
    if all(events_prf.isprf)
    %%--Compute S/N (t-runs-stims-chs)
    residuals  = [];
    evoked = [];
    for kk=unique(runsList(:)')
        [residuals(:,kk,:,:)]  = sqrt(average_trials(regressed(:,runsList==kk,:).^2, events_prf(runsList==kk,:), stimNames));
        [evoked(:,kk,:,:)] = (average_trials(epochs(stimon_idx,stmlist==1&runsList==kk,:), events_prf(stmlist==1&runsList==kk,:), stimNames));
    end
    residuals = shiftdim(nanrms(nanrms(residuals,1),4),1);
    evoked    = shiftdim(nanrms(nanrms(evoked,1),4),1);
    
    select_runs(:,2) = max(evoked,[],2)./max(residuals,[],2) > opts.run_thresh_sn;
    end
    
    %%% Exclude runs based on experiment note
    switch subject
        case {'som674'}
            select_runs([1,2],3) = false;
    end
    
    %%% Combine all selection
    select_runs = all(select_runs,2);
    
    %-- Plot runs and electrodes
    if opts.doplots
        for kk = reshape(find(~select_runs),1,[])
            runsidx   = runsList==kk;
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

            ecog_plotTimecourses(trials, whichElectrodes, eventList, specs);
            hgexport(gcf, fullfile(opts.plotsavedir,'runselection',sprintf('%s-run%02d',figureName,kk)), hgexport('factorystyle'), 'Format', 'png'); close;
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
    select_idx  = ismember(runsList,find(select_runs));
    epochs      = epochs(:,select_idx,:);
    regressed   = regressed(:,select_idx,:);
    events      = events(select_idx,:);
    events_prf  = events_prf(select_idx,:);
    outliers    = outliers(select_idx,:);
    runsList    = runsList(select_idx);
    
%% STEP 4 Finalize processing
    %-- Collect into an output struct
    data{ii}.epochs   = epochs;
    data{ii}.events   = events;
    data{ii}.channels = channels;
    idat = data{ii};
    %-- Save out the data
    if isSave && ~isempty(epochs)
        fprintf('[%s] Saving data for subject %s to %s \n',mfilename, subject, opts.outputDir);
        filename    = fullfile(opts.outputDir, sprintf('%s_%s.mat', subject, opts.fileid));
        saveauto(filename,'-struct','idat')
    end
    
    %-- Plot runs and electrodes
    if opts.doplots
        for kk = unique(runsList(:)')
            runsidx   = runsList==kk;
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

            ecog_plotTimecourses(trials, whichElectrodes, eventList, specs);
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
    
end
fprintf('[%s] Done! \n',mfilename);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Trial averaging
function [epochs_averaged] = average_trials(epochs, events, stimNames)   
    
    %-- Average across trials within stimulus condition
    stimList = unique(events.trial_name(cellstrfind(events.trial_name, stimNames,1)),'stable');
    epochs_averaged = nan(size(epochs,1), length(stimList), size(epochs,3));
    for ii = 1:length(stimList)
        trial_idx = cellstrfind(events.trial_name, stimList{ii},1);
        epochs_averaged(:,ii,:) = mean(epochs(:,trial_idx,:),2, 'omitnan');
    end
end

