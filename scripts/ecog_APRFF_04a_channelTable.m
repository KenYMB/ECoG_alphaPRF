% Save all channel information without bad channel exclusion

% 20211018 Yuasa

%% Define dataset
close all; clearvars;
%-- Set path
checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('rawPth',        'ECoGCAR');
SetDefault('chanPth',       'Channels');
SetDefault('description',   'reref');
else
rawPth         = 'ECoGCAR';
chanPth        = 'Channels';
description    = 'reref';
end
ouputPth = SetDefaultAnalysisPath('DAT',chanPth);
if ~exist(ouputPth,'dir'), mkdir(ouputPth);  end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);
 
%% Load channels
atlasName = {'wang15_mplbl', 'wang15_fplbl','benson14_varea', 'benson14_eccen', 'benson14_angle', 'benson14_sigma'};
fsample   = 512;
tasks     = {'prf'};
        
channelsC = cell(0);
fldnames = [];
for subject= reshape(string(subjectList),1,[])
    %-- Load channel information
    [~, channeltable] = bidsEcogGetPreprocData(SetDefaultAnalysisPath('BIDS',rawPth), subject, [], tasks, [], description, fsample);
    channeltable.subject_name(:) = cellstr(subject);
    %-- Load electrode information
    electrodes = bidsEcogMatchElectrodesToAtlas(bidsRootPath, subject, [], atlasName, [], 0);
    %-- Combine
    channels = bair_addVisualAtlasNamesToChannelTable(channeltable,electrodes);
    
    %-- Get common field names
    if isempty(fldnames),   fldnames = channels.Properties.VariableNames;
    else,                   fldnames = intersect(fldnames, channels.Properties.VariableNames,'stable');
    end
    
    %-- Exclude reference channels
    channels(~ismember(lower(channels.type),{'ecog','seeg'}),:) = [];
    
    %-- Rename group
    channels.group(~ismember(lower(channels.group),{'strip','grid','hdgrid','depth'}) & ismember(lower(channels.type),{'ecog'})) = {'strip'};
    channels.group(~ismember(lower(channels.group),{'strip','grid','hdgrid','depth'}) & ismember(lower(channels.type),{'seeg'})) = {'depth'};
    
    %-- Merge across subjects
    channelsC{end+1,1}  = channels;
end

for ii = 1:length(channelsC)
    %-- Exclude subject specific field
    channelsC{ii} = channelsC{ii}(:,fldnames);
    
    %-- Output
    channels  = channelsC{ii};
    subject   = channels.subject_name{1};
    save(fullfile(ouputPth,sprintf('%s-channels',subject)), 'channels');
end
