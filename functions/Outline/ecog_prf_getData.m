function [data] = ecog_prf_getData(compute, inputDir, outputDir, subjectList, sessionList, epochTime, fsample, tasks, stims, description, fileid)

% Description: 
%
% [data] = ecog_prf_getData(compute, [inputDir], [outputDir], [subjectList], [sessionList], [epochTime], [sampleRate], [taskList], [stimNames])
%
% Input (use varargin?)
% - compute         = [],0,1 % compute the data (1), load from disk (0), or
%                              select automatically for each subject
% - inputDir
% - outputDir
% - subjectslist    % subject name list (sub-OO)
% - sessionList     % session name list (ses-OO)
% - epochTime       = 1x2 vector
% - sampleRate      % target sampling frequency to resample [Hz]
% - taskList        % task name list (events.task_name)
% - stimNames       % stimuluss names (events.trial_name)
% 
% - inputFileID
% - outputFileID
%
% Output
% - data            % a cell array with for each cell a struct with the following fields:
%   - subject
%   - epochs        = t x events x channels
%   - t
%   - events
%   - channels
%   - fsample

% Dependency: ECoG_utils, cellstrfind, saveauto, istablefield

% 20200127 - Yuasa: modified from tde_getData in temporalECoG
% 20200317 - Yuasa: event correction for som674
% 20200322 - Yuasa: enable to load no-prf events
% 20201027 - Yuasa: update channel selection
%                   - include whole HD grids parts of which covered visual areas
%                   - include channels which have some FPM 
% 20201209 - Yuasa: update to load hemisphere of each electrode
% 20220222 - Yuasa: load HDgrid information from subjectlist.tsv

%% Define inputs 

% <compute>
if ~exist('compute', 'var') || isempty(compute)
    checkfile = true;
else
    checkfile = false;
end  

% <inputDir>
if ~exist('inputDir', 'var') || isempty(inputDir)
    inputDir = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR');
end  

% <outputDir>
if ~exist('outputDir', 'var') || isempty(outputDir)
	outputDir = fullfile(analysisRootPath, 'Data','Raw');
end 
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% <subjectList>
subjectList_fname = 'subjectlist.tsv';
if ~exist('subjectList', 'var') || isempty(subjectList)
    SbjInfo = loadSbjInfo(subjectList_fname);
    subjectList = SbjInfo.participant_id;
end

% <sessionList>
if ~exist('sessionList', 'var') || isempty(sessionList)
    sessionList = []; % default: all sessions
end

% <epochTime>
if ~exist('epochTime', 'var') || isempty(epochTime)
    epochTime = [-0.2 1];
end

% <fsample>
if ~exist('fsample', 'var')
    fsample = [];
elseif ~isempty(fsample)
    fsample = round(fsample);
end

% <tasks> (Use default for pRF project)
if ~exist('tasks', 'var') || isempty(tasks)
    tasks = {'prf'};
elseif ~iscell(tasks)
    tasks = {tasks};
end

% <stims> (Use default for pRF project)
if ~exist('stims', 'var') || isempty(stims)
    stims = {'*'};
elseif ~iscell(stims)
    stims = {stims};
end

% <preprocessing data type> (Use default for pRF project)
if ~exist('description', 'var') || isempty(description)
    description = 'reref';
end

% <save fille name>  (Use default for pRF project)
if ~exist('fileid', 'var') || isempty(fileid)
    fileid = 'data_visualelecs';
end

%% Loop across subjects
%-- Correct Subject Information
SbjInfo    = loadSbjInfo(subjectList_fname,'all');
hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');

%-- Main
data = cell(length(subjectList),1);
for ii = 1 : length(subjectList)
    
    subject = subjectList{ii};
    
    %-- Determine if we're loading or computing the data
    filename = fullfile(outputDir, sprintf('%s_%s.mat', subject, fileid));
    if checkfile
        compute = ~exist(filename,'file');
    end
    
    %-- HD grid flag
    if hasSbjInfo && istablefield(SbjInfo,'hasHDgrid')
        hasHDgrid = any(strcmpi(SbjInfo.hasHDgrid(ismember(SbjInfo.participant_id,subject)),'yes'));
    else
        hasHDgrid = false;
    end
    
    if ~compute
        
        %-- load from outputDir   
        fprintf('[%s] Loading data for subject %s from %s ',mfilename, subject, outputDir); 
        data{ii} = load(filename);
        fprintf('\n');
    
    else
        
        fprintf('[%s] Computing data for subject %s \n',mfilename, subject);

        %% STEP 0: GET visual area matches for this subject
        fprintf('[%s] Step 0: Computing matches with visual atlases...\n',mfilename);
        
        specs = [];
        specs.pID           = subject; 
        specs.plotmesh      = 'none';
        specs.atlasNames    =  {'wang2015_atlas','wang15_mplbl','wang15_fplbl', 'benson14_varea', 'benson14_eccen', 'benson14_angle', 'benson14_sigma','benson20_mplbl','benson20_fplbl'};
        BIDSformatted       = 1;
        visualelectrodes    = electrode_to_nearest_node(specs, BIDSformatted);
 
        %% STEP 1: GET THE DATA
        fprintf('[%s] Step 1: Loading data...\n',mfilename);

        %-- Different setting recordings and prepare variables (set for each subject)
        if ismember(subject,{'beilen'})
            ngroups     = 2;
            runnumsset  = {{{'1','2'}},{{'3','4'}}};     % @2048Hz, 512Hz
        else
            ngroups     = 1;
            runnumsset  = {[]};
        end
        subdata     = cell(1,ngroups);
        channels    = cell(1,ngroups);
        events      = cell(1,ngroups);
        datfsample  = cell(1,ngroups);
        
        if isempty(sessionList)
            [sessions] = bidsSpecifySessions(inputDir, subject);
        else
            sessions = sessionList;
        end
        
        for gg = 1:ngroups
            %-- Load data and exception handling (for bailen-prf)
            %%% In bidsEcogGetPreprocData.m,
            %%% if one or more electrodes are bad in one session but not another,
            %%% for now, label the electrode as good across all sessions.
            if ismember(subject,{'beilen'})
                subdata1 = [];  events1  = table([]);
                subdata2 = [];  events2  = table([]);
                if any(ismember(tasks,'prf'))
                    runnums = runnumsset{gg};
                    [subdata1, channels{gg}, events1, datfsample{gg}] = bidsEcogGetPreprocData(inputDir, subject, sessions, {'prf'}, runnums, description);
                    events1.trial_name   = num2cell(events1.trial_name);
                end
                if any(~ismember(tasks,'prf'))
                    runnums = repmat(runnumsset{gg},1,numel(find(~ismember(tasks,'prf'))));
                    [subdata2, channels{gg}, events2, datfsample{gg}] = bidsEcogGetPreprocData(inputDir, subject, sessions, tasks(~ismember(tasks,'prf')), runnums, description);
                end

                if length(tasks)>1 && any(ismember(tasks,'prf'))
                    events2.event_sample = events2.event_sample + size(subdata1,2);
                    events2.onset        = events2.onset + size(subdata1,2)./datfsample{gg};                
                    events{gg}  = [events1; events2];
                elseif any(ismember(tasks,'prf'))
                    events{gg}  = events1;
                else
                    events{gg}  = events2;
                end

                subdata{gg}     = cat(2,subdata1,subdata2);

            else
                runnums = repmat(runnumsset{gg},1,length(tasks));
                [subdata{gg}, channels{gg}, events{gg}, datfsample{gg}] = bidsEcogGetPreprocData(inputDir, subject, sessions, tasks, runnums, description);
            end
            
            %-- Correct channel names (for som800)
            if ismember(subject,{'som800'})
                [eleccat, elecnum] = strtok(channels{gg}.name,int2str(0:9));
                for el = 1:length(channels{gg}.name)
                    channels{gg}.name{el} = sprintf('%s%02d',eleccat{el},str2double(elecnum{el}));
                end
            end
            
            %-- Read in electrode cooordinates from session 1
            try  % skip if electrodes.tsv does not exist
                projectDir = dir(fullfile(inputDir,'..','..'));
                projectDir = projectDir(1).folder;
                [electrode_table] = bidsEcogReadElectrodeFile(projectDir, subject, sessions{1});
                if ismember('hemisphere',electrode_table.Properties.VariableNames)
                    [check_channels, check_electrodes] = ismember(channels{gg}.name,electrode_table.name);
                    channels{gg}.hemisphere(check_channels) = electrode_table.hemisphere(check_electrodes);
                end
            end

            %-- Correct wrong event onsets (for som674-prf)
            if ismember(subject,{'som674'}) && any(ismember(tasks,'prf'))
                warning('[%s] Event onset is wrong for %s. Trying to correct...',mfilename,subject)
                runlist = categorical(strcat(events{gg}.task_name,'-',events{gg}.session_name,'-',events{gg}.run_name));
                runcate = categories(runlist);
                runcate(~startsWith(runcate,'prf-')) = [];
                for jj=1:length(runcate)
                    currun      = find(runlist==runcate{jj});
                    ntrl        = numel(currun);
                    idealsoa    = 0.85;             % 1 trial = 0.50s(stimulus)+0.35s(ISI)
                    avgdly      = -4e-4;            % 849.6ms is the mode and the median of difference of onset trigger 
                    dlythresh   = 0.02;             % maxmum lag in the other subjects is 0.0175s
                    reconsets   = events{gg}.onset(currun);
                    idealonsets = reconsets(1)...
                                   + (0:idealsoa:((ntrl-1)*idealsoa))...
                                   + linspace(0,avgdly*(ntrl-1),ntrl);
                    badonsetidx = abs(reconsets(:)-idealonsets(:)) > dlythresh;
                    reconsets(badonsetidx)       = idealonsets(badonsetidx);
                    events{gg}.onset(currun)         = reconsets;
                    events{gg}.event_sample(currun)  = round(reconsets.*datfsample{gg});
                    if numel(find(badonsetidx))>0
                    fprintf('[%s] %d trials are corrected in %s\n',mfilename,numel(find(badonsetidx)),runcate{jj});
                    end
                end
            end

            %-- Correct events [specific to pRF analysis] (for bailen-prf,som674-prf)
            if any(ismember(tasks,'prf')) && (~iscellstr(events{gg}.trial_name) || ~any(ismember(events{gg}.trial_name,'BLANK')))
                warning('[%s] Event list is incomplete for %s. Trying to correct...',mfilename,subject)
                eventfile = fullfile(analysisRootPath, 'Data', 'eventlist.mat');
                if exist(eventfile,'file')
                    events0 = load(eventfile);   events0 = events0.events;
                    prftrial   = ismember(events{gg}.task_name,'prf');
                    events{gg}.trial_name(prftrial) = events0.trial_name(events{gg}.stim_file_index(prftrial));
                    events{gg}.trial_type(prftrial) = events0.trial_type(events{gg}.stim_file_index(prftrial));
                else
                    warning('[%s] templete event file is not exist',mfilename);
                    fprintf(2,'[%s] Failed to correct event list for %s\n',mfilename,subject);
                end
            end

            %-- select trials for analysis
            selevent    = cellstrfind(events{gg}.trial_name, stims,1);
            events{gg}      = events{gg}(selevent,:);
        end
        
        %% STEP 2: SELECT A SUBSET OF CHANNELS 

        for gg = 1:ngroups
            %-- Add visual area names (W and B) ecc, angle, sigma to channels table
            [channels{gg}] = bair_addVisualAtlasNamesToChannelTable(channels{gg},visualelectrodes);

            %-- Make selection on visual only, index into data + channels
            fprintf('[%s] Step 2: Selecting channels with visual matches \n',mfilename);
            if ismember('matchednode',channels{gg}.Properties.VariableNames)
                chan_idx = find(~isnan(channels{gg}.matchednode) & contains(channels{gg}.status, 'good'));
            else
                chan_idx1 = [];
                if ismember('wangarea',channels{gg}.Properties.VariableNames)
                  chan_idx1 = [chan_idx1;find(~contains(channels{gg}.wangarea, 'none') & contains(channels{gg}.status, 'good'))];
                end
                if ismember('bensonarea',channels{gg}.Properties.VariableNames)
                  chan_idx1 = [chan_idx1;find(~contains(channels{gg}.bensonarea, 'none') & contains(channels{gg}.status, 'good'))];
                end
                if ismember('hcparea',channels{gg}.Properties.VariableNames)
                  chan_idx1 = [chan_idx1;find(~contains(channels{gg}.hcparea, 'none') & contains(channels{gg}.status, 'good'))];
                end
                chan_idx2 = [];
                chan_idx2 = [chan_idx2; find(sum(channels{gg}(:,contains(channels{gg}.Properties.VariableNames,'wangprob_')).Variables,2)>0 ...
                                            & contains(channels{gg}.status, 'good'))];
                chan_idx2 = [chan_idx2; find(sum(channels{gg}(:,contains(channels{gg}.Properties.VariableNames,'hpcprob_')).Variables,2)>0 ...
                                            & contains(channels{gg}.status, 'good'))];
            
                chan_idx = unique([chan_idx1; chan_idx2]);  
            end
              %%-- add grid channels (GA*,GB*,G*) which covered visual area
              chan_idx3 = [];
              if hasHDgrid
                  chan_grid = find(contains(channels{gg}.name, 'GA') & contains(channels{gg}.status, 'good'));
                  if length(intersect(chan_idx,chan_grid))>(length(chan_grid)*0.35)
                    chan_idx3 = [chan_idx3; chan_grid];
                  end
                  chan_grid = find(contains(channels{gg}.name, 'GB') & contains(channels{gg}.status, 'good'));
                  if length(intersect(chan_idx,chan_grid))>(length(chan_grid)*0.35)
                    chan_idx3 = [chan_idx3; chan_grid];
                  end
                  chan_grid = find(startsWith(channels{gg}.name, 'G') & contains(channels{gg}.status, 'good')...
                                   & ~contains(channels{gg}.name, 'GA') & ~contains(channels{gg}.name, 'GB'));
                  if length(intersect(chan_idx,chan_grid))>(length(chan_grid)*0.35)
                    chan_idx3 = [chan_idx3; chan_grid];
                  end
              end
            chan_idx = unique([chan_idx; chan_idx3]);

            fprintf('[%s] Step 2: Found %d channels with visual matches out of %d ecog channels \n', ...
                mfilename, length(chan_idx), length(find(contains(lower(channels{gg}.type), {'ecog', 'seeg'}))));

            subdata{gg} = subdata{gg}(chan_idx,:);
            channels{gg} = channels{gg}(chan_idx,:);
        end
        
        %% STEP 3: CHECK SAMPLE RATES, CONCATENATE ANDS SHIFT DATA
        
        for gg = 1:ngroups
            %-- RESAMPLING DATA
            fprintf('[%s] Step 3: Checking sample rates...\n',mfilename);

            if ~isfield(events{gg},'event_sample')
                events{gg}.event_sample = round(events{gg}.onset * datfsample{gg});
            end
            datfsample{gg} = round(datfsample{gg});
            if ~isempty(fsample) && (datfsample{gg} ~= fsample)
                fprintf('[%s] Step 3: Sample rate does not match requested sample rate. Resampling \n',mfilename);
                subdata{gg} = resample(subdata{gg}', fsample, datfsample{gg})';
                events{gg}.event_sample = round(events{gg}.event_sample/(datfsample{gg}/fsample));
                channels{gg}.sampling_frequency(:) = fsample;
                datfsample{gg} = fsample;
            end
        end
        
        %-- CONCATENATING DATA
        for gg = 2:ngroups
            assert(datfsample{1}==datfsample{gg},'Detect various sampling frequencies. Please specify target frequency to resample data.');
            if isfield(summary(channels{1}), 'status') && isfield(summary(channels{gg}), 'status') && ~isequal(channels{1}.status, channels{gg}.status)
                %%% This means one or more electrodes are bad in one session
                %%% but not another. For now, label the electrode as good
                %%% across all sessions:
                channels{1}.status(~strcmp(channels{1}.status, channels{gg}.status)) = {'good'};
            end
            events{gg}.event_sample = events{gg}.event_sample + size(subdata{1},2);
            events{gg}.onset = events{gg}.onset + size(subdata{1},2)./datfsample{1};
            events{1}  = [events{1};events{gg}];
            subdata{1} = cat(2,subdata{1}, subdata{gg});
        end
        subdata     = subdata{1};
        channels    = channels{1};
        events      = events{1};
        datfsample  = datfsample{1};
        
        %-- SHIFT the UMCU data 
        if contains(subject, {'chaam', 'beilen'}) % SHOULD BE READ FROM participants.tsv, if site column = umcu
            fprintf('[%s] Step 3: This is a umcu patient. Shifting onsets \n',mfilename);
                       
            %-- Shift onsets
            shiftInSeconds = 0.072; % 72 ms; determined through cross correlation, see s_determineOnsetShiftUMCUvsNYU.m
            shiftInSamples = round(shiftInSeconds*datfsample); 
            events.onset = events.onset + shiftInSeconds;
            events.event_sample = events.event_sample + shiftInSamples; 
            
            %%% Remove electrodes identified as epileptic % TEMPORARY UNTIL GIO
            %%% FIXES BIDS FORMATTING FOR THIS PATIENT AFTER WHICH THESE
            %%% ELECTRODES SHOULD BE INDICATED AS "BAD" IN THE ELECTRODES FILES
            if contains(subject,'chaam')
                chan_idx = find(~contains(channels.name, {'Oc12', 'Oc13', 'Oc21', 'Oc22'}));
                subdata = subdata(chan_idx,:);
                channels = channels(chan_idx,:);
            end
        end
        
        %% STEP 4: EPOCH THE DATA 
        
        fprintf('[%s] Step 4: Epoching data \n',mfilename);
        
        [epochs, t] = ecog_makeEpochs(subdata, events.event_sample, epochTime, datfsample);  
        
        fprintf('[%s] Step 4: Found %d epochs across %d runs and %d sessions \n', ...
            mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));
        
        %% STEP 5: Save out a single preproc file for each subject 
        
        %-- Remove irrelevant/redundant columns from channels and events tables 
        if isfield(summary(events),'onset'), events = removevars(events,'onset');end
        if isfield(summary(events),'event_sample'), events = removevars(events,'event_sample');end
        if isfield(summary(events),'stim_file'), events = removevars(events,'stim_file');end
        if isfield(summary(channels),'notch'), channels = removevars(channels,'notch');end
        if isfield(summary(channels),'status'), channels = removevars(channels,'status');end
        if isfield(summary(channels),'description'), channels = removevars(channels,'description');end
        if isfield(summary(channels),'status_description'), channels = removevars(channels,'status_description');end
        
        %-- Add a subject index column to channels and events tables:
        events.subject_name = repmat({subject}, [height(events),1]);
        channels.subject_name = repmat({subject}, [height(channels),1]);
        
        %-- Collect into an output struct
        data{ii}.subject  = subject;
        data{ii}.epochs   = epochs;
        data{ii}.t        = t;
        data{ii}.events   = events;
        data{ii}.channels = channels;
        data{ii}.fsample  = datfsample;
        idat = data{ii};
        
        %-- Save out the data
        fprintf('[%s] Step 5: Saving data for subject %s to %s \n',mfilename, subject, outputDir);
        saveauto(filename,'-struct','idat')

    end   
end
fprintf('[%s] Done! \n',mfilename);
end
