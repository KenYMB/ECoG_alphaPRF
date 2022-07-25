function [data] = ecog_prf_getData(compute, subjects, sessions, tasks, stims, epochTime, fsample, inputDir, outputDir, fileid)

% Description: 
%
% [data] = ecog_prf_getData(compute, [subjectList, sessionList, taskList, stimNames, epochTime, sampleRate, inputDir, outputDir, fileid])
%
% Input (use varargin?)
% - compute         = [],0,1 % compute the data (1), load from disk (0), or
%                              select automatically for each subject
% - subjectslist    % subject name list (sub-OO)
% - sessionList     % session name list (ses-OO)
% - epochTime       = 1x2 vector
% - sampleRate      % target sampling frequency to resample [Hz]
% - taskList        % task name list (events.task_name)
% - stimNames       % stimuluss names (events.trial_name)
% - inputDir        % directory name to load ECoG data in BIDS derivatives
%                   % data will be loaded from 
%                   %   fullfile(bidsRootPath, 'derivatives', inputDir);
%                   % default = 'ECoGCAR'
% - outputDir       % directory name to save data
%                   % data will be saved in 
%                   %   fullfile(analysisRootPath, 'Data', outputDir);
%                   % default = 'Raw'
% 
% - fileid          % string to be added to filename for saved out data
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
% 20200322 - Yuasa: enable to load no-prf events
% 20201027 - Yuasa: update channel selection
%                   - include whole HD grids parts of which covered visual areas
%                   - include channels which have some FPM 
% 20201209 - Yuasa: update to load hemisphere of each electrode
% 20220222 - Yuasa: load HDgrid information from subjectlist.tsv
% 20220409 - Yuasa: major update not uot use legacy functions

%% Define inputs 

% <compute>
if ~exist('compute', 'var') || isempty(compute)
    checkfile = true;
else
    checkfile = false;
end  

% <inputDir> (passthrough if inputDir is path)
if ~exist('inputDir', 'var') || isempty(inputDir)
    inputDir = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR');
elseif numel(strsplit(inputDir,filesep)) == 1
    inputDir = fullfile(bidsRootPath, 'derivatives', inputDir);
end  

% <outputDir> (passthrough if outputDir is path)
if ~exist('outputDir', 'var') || isempty(outputDir)
	outputDir = fullfile(analysisRootPath, 'Data','Raw');
elseif numel(strsplit(outputDir,filesep)) == 1
	outputDir = fullfile(analysisRootPath, 'Data', outputDir);
end 
if ~exist(outputDir, 'dir'); mkdir(outputDir); end

% <subjectList>
subjectList_fname = 'subjectlist.tsv';
if ~exist('subjects', 'var') || isempty(subjects)
    SbjInfo = loadSbjInfo(subjectList_fname);
    subjects = SbjInfo.participant_id;
end

% <sessionList>
if ~exist('sessions', 'var') || isempty(sessions)
    sessions = []; % default: all sessions
elseif ~iscell(sessions)
    sessions = {sessions};
end

% <epochTime>
if ~exist('epochTime', 'var') || isempty(epochTime)
    epochTime = [-0.2 1];
end

% <fsample> (Use default for BAIR ECoG)
if ~exist('fsample', 'var')
    fsample = 512;
elseif ~isempty(fsample)
    fsample = round(fsample);
end

% <taskList> (Use default for pRF project)
if ~exist('tasks', 'var') || isempty(tasks)
    tasks = {'prf'};
elseif ~iscell(tasks)
    tasks = {tasks};
end

% <stimNames> (Use default for pRF project)
if ~exist('stims', 'var') || isempty(stims)
    stims = {'*'};
elseif ~iscell(stims)
    stims = {stims};
end

% <preprocessing data type> (Use default to load from ECoGCAR)
description = 'reref';

% <save fille name>  (Use default for pRF project)
if ~exist('fileid', 'var') || isempty(fileid)
    fileid = 'data_visualelecs';
end

%% Loop across subjects

%-- Main
data = cell(length(subjects),1);
for ii = 1 : length(subjects)
    
    subject = subjects{ii};
    
    %-- Determine if we're loading or computing the data
    filename = fullfile(outputDir, sprintf('%s_%s.mat', subject, fileid));
    if checkfile
        compute = ~exist(filename,'file');
    end
        
    if ~compute
        
        %-- load from outputDir   
        fprintf('[%s] Loading data for subject %s from %s ',mfilename, subject, outputDir); 
        data{ii} = load(filename);
        fprintf('\n');
    
    else
        
        fprintf('[%s] Computing data for subject %s \n',mfilename, subject);

        %% STEP 1: GET THE DATA
        fprintf('[%s] Step 1: Loading data...\n',mfilename);
        
        %-- Read in voltage data
        [subdata, channels, events, datfsample] = bidsEcogGetPreprocData(inputDir, subject, sessions, tasks, [], description, fsample);
        if isempty(subdata), warning('[%s] No voltage data found for subject %s!', mfilename, subject); continue; end
        
        %-- Read in electrode data and match to atlas
        atlasName = {'wang15_mplbl', 'wang15_fplbl','benson14_varea', 'benson14_eccen', 'benson14_angle', 'benson14_sigma'};
        [electrodes] = bidsEcogMatchElectrodesToAtlas(bidsRootPath, subject, sessions, atlasName, [], 0);
        [channels,chanidx] = bair_addVisualAtlasNamesToChannelTable(channels,electrodes);
        subdata = subdata(chanidx,:);
        
        %-- Correct data in personal
        %--- Correct wrong event onsets (for specific subject)
        cumulativecorrection = false;
        if ismember(subject,{'p05'}) && any(ismember(tasks,'prf'))
            warning('[%s] Event onset is wrong for %s. Trying to correct...',mfilename,subject)
            runlist = categorical(strcat(events.task_name,'-',events.session_name,'-',events.run_name));
            runcate = categories(runlist);
            runcate(~startsWith(runcate,'prf-')) = [];
            for jj=1:length(runcate)
                currun      = find(runlist==runcate{jj});
                ntrl        = numel(currun);
                idealsoa    = 0.85;             % 1 trial = 0.50s(stimulus)+0.35s(ISI)
                dlythresh   = 0.02;             % maxmum lag in the other subjects is 0.0175s
                reconsets   = events.onset(currun);
                avgdly      = median(diff(reconsets)-idealsoa);        % -4e-4 (849.6ms is the mode and the median of difference of onset trigger) 
                if cumulativecorrection
                idealonsets = zeros(size(reconsets));
                for kk=1:ntrl
                    idealonsets(kk,1) = reconsets(kk);
                    if kk > 1 && abs(diff(idealonsets([kk-1,kk]))-idealsoa) > dlythresh
                        idealonsets(kk,1) = idealonsets(kk-1)+idealsoa+avgdly;
                    end
                end
                else
                idealonsets = reconsets(1)...
                               + (0:idealsoa:((ntrl-1)*idealsoa))...
                               + linspace(0,avgdly*(ntrl-1),ntrl);
                end
                badonsetidx = abs(reconsets(:)-idealonsets(:)) > dlythresh;
                reconsets(badonsetidx)       = idealonsets(badonsetidx);
                events.onset(currun)         = reconsets;
                if isnumeric(events.event_sample)
                  events.event_sample(currun)  = round(reconsets.*datfsample);
                end
                if numel(find(badonsetidx))>0
                fprintf('[%s] %d trials are corrected in %s\n',mfilename,numel(find(badonsetidx)),runcate{jj});
                end
            end
        end

        %--- Correct events [correct trial names for pRF analysis]
        if any(ismember(tasks,'prf')) && (~iscellstr(events.trial_name) || ~any(ismember(events.trial_name,'BLANK')))
            warning('[%s] Event list is incomplete for %s. Trying to correct...',mfilename,subject)
            eventfile = 'eventlist.mat';
            if exist(eventfile,'file')
                fprintf(1,'[%s] Loading %s %s\n',mfilename,which(eventfile));
                events0 = load(eventfile);   events0 = events0.events;
                prftrial   = ismember(events.task_name,'prf');
                events.trial_name(prftrial) = events0.trial_name(events.stim_file_index(prftrial));
                events.trial_type(prftrial) = events0.trial_type(events.stim_file_index(prftrial));
            else
                warning('[%s] templete event file is not exist',mfilename);
                fprintf(2,'[%s] Failed to correct event list for %s\n',mfilename,subject);
            end
        end

        %-- select trials for analysis
        selevent    = cellstrfind(events.trial_name, stims,1);
        events      = events(selevent,:);
        
        
        %% STEP 2: SELECT A SUBSET OF CHANNELS 
        
        fprintf('[%s] Step 2: Selecting channels with visual matches \n',mfilename);

        %-- Make selection on visual only, index into data + channels
        [~, chan_idx1] = ecog_visualElectrodes(channels);
        
        %-- add HD grid channels (GA*,GB*,G*) which covered visual area (at least 10%)
        gridHDthresh = 64;
        visarearate = 0.1;
        [~, chan_idx2] = ecog_HDgridElectrodes(channels,gridHDthresh,visarearate,chan_idx1);        
        
        %-- Exclude bad channel
        goodchan = contains(channels.status, 'good');
        chan_idx = any([chan_idx1, chan_idx2],2);
        chan_idx = all([chan_idx,goodchan],2);

        fprintf('[%s] Step 2: Found %d channels with visual matches out of %d ecog channels \n', ...
            mfilename, sum(chan_idx), length(find(contains(lower(channels.type), {'ecog', 'seeg'}))));

        subdata = subdata(chan_idx,:);
        channels = channels(chan_idx,:);
        
        %% STEP 3: DEAL WITH UMCU DATA (shift onsets)
                
        %-- SHIFT the UMCU data 
        if contains(subject, {'p01','p02'}) % SHOULD BE READ FROM participants.tsv, if site column = umcu
            fprintf('[%s] Step 3: This is a umcu patient. Shifting onsets \n',mfilename);
                       
            %-- Shift onsets
            shiftInSeconds = 0.072; % 72 ms; determined through cross correlation, see s_determineOnsetShiftUMCUvsNYU.m
            events.onset = events.onset + shiftInSeconds;
            if isnumeric(events.event_sample)
              shiftInSamples = round(shiftInSeconds*datfsample); 
              events.event_sample = events.event_sample + shiftInSamples; 
            end
            
        end
        
        %% STEP 4: EPOCH THE DATA 
        
        fprintf('[%s] Step 4: Epoching data \n',mfilename);
        
        [epochs, t] = ecog_makeEpochs(subdata, events.onset, epochTime, datfsample);  
        
        fprintf('[%s] Step 4: Found %d epochs across %d runs and %d sessions \n', ...
            mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));
        
        %% STEP 5: Save out a single preproc file for each subject 
        
        %-- Remove irrelevant/redundant columns from events table
        if isfield(summary(events),'onset'), events = removevars(events,'onset');end
        if isfield(summary(events),'event_sample'), events = removevars(events,'event_sample');end
        if isfield(summary(events),'stim_file'), events = removevars(events,'stim_file');end

        %-- Remove irrelevant/redundant columns from channels table
        if isfield(summary(channels),'notch'), channels = removevars(channels,'notch');end
        if isfield(summary(channels),'status'), channels = removevars(channels,'status');end
        if isfield(summary(channels),'description'), channels = removevars(channels,'description');end
        if isfield(summary(channels),'status_description'), channels = removevars(channels,'status_description');end
        if isfield(summary(channels),'size'), channels = removevars(channels,'size');end
        if isfield(summary(channels),'material'), channels = removevars(channels,'material');end
        if isfield(summary(channels),'manufacturer'), channels = removevars(channels,'manufacturer');end
        
        %-- Add a subject index column to channels and events tables:
        if ~istablefield(events,'subject_name')
            events.subject_name(:) = {subject};
        end
        if ~istablefield(channels,'subject_name')
            channels.subject_name(:) = {subject};
        end
        
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
