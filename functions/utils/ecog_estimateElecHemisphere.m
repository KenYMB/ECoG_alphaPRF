function    channels = ecog_estimateElecHemisphere(channels, BIDS_directory, subject_name)
% channels = ecog_estimateElecHemisphere(channels)
% channels = ecog_estimateElecHemisphere(channels, [], subject)
% channels = ecog_estimateElecHemisphere(channels, BIDS_directory)
% channels = ecog_estimateElecHemisphere(channels, BIDS_directory, subject)
%   add 'hemisphere' column to channel table.
%   The hemisphere is loaded from electrodes.tsv in the BIDS directory, 
%   or estimated based on subjectlist.tsv.
%   If channel table does not have subject_name field, you can specify the
%   subject name.

% 20220307 Yuasa

%-- Check available information
hashemisphere  = istablefield(channels,'hemisphere');
hasBIDSdir = exist('BIDS_directory','var') && ~isempty(BIDS_directory);
if istablefield(channels,'subject_name')
    hassubjectname = true;
    subjects = cellstr(channels.subject_name);
elseif exist('subject_name','var') && ~isempty(subject_name)
    hassubjectname = true;
    if size(subject_name,1) == 1
        subjects = repmat(cellstr(subject_name),height(channels),1);
    elseif size(subject_name,1) == height(channels)
        subjects = cellstr(subject_name);
    else
        error('''subject'' must be a charactor vector or has the same height to the channel table.');
    end
else    
    hassubjectname = false;
end

%-- Estimate electrodes located hemisphere
if hashemisphere
    %%% Pass through channels.hemisphere %%%
    
elseif ~hassubjectname
    %%% Error for unknown subject %%%
    error('Subject name is required to estimate electrode locations');
    
elseif hasBIDSdir
    %%% load from electrodes.tsv %%%
    assert(exist(BIDS_directory,'dir'),message('MATLAB:load:couldNotReadFile',BIDS_directory));
    subjectList = unique(subjects);
    %-- Outputs
    channels.hemisphere(:)         = {'n/a'};
    for isbj = 1:length(subjectList)
        subject = subjectList{isbj};
        %-- Estimate session directory
        [sessions] = bidsSpecifySessions(BIDS_directory, subject);
        idx = find(contains(lower(sessions), {'ecog', 'iemu'}));
        if isempty(idx), error('no ECOG sessions found for subject %s', subject); end
        session = sessions{idx(1)};
        %-- Read in electrode cooordinates
        [electrode_table] = bidsEcogReadElectrodeFile(BIDS_directory, subject, session);
        %-- Get channel index
        currchan = find(ismember(subjects,subject));
        [~,chanidx,tableidx]=intersect(channels.name(currchan),electrode_table.name,'stable');
        %-- Outputs
        channels.hemisphere(currchan(chanidx))         = electrode_table.hemisphere(tableidx);
    end
        
else
    %%% load from subjectlist.tsv %%%
    %-- Collect Subject Information
    subjectList_fname = 'subjectlist.tsv';
    SbjInfo    = loadSbjInfo(subjectList_fname,'all');
    hasSbjInfo = ~isempty(SbjInfo) && istablefield(SbjInfo,'participant_id');
    
    %-- Categorize subjects
    assert(hasSbjInfo,'subjectlist.tsv not found'); 
    subj_L = SbjInfo.participant_id(ismember(SbjInfo.hemi,'L'))';  % left hemisphere
    subj_R = SbjInfo.participant_id(ismember(SbjInfo.hemi,'R'))';  % right hemisphere
    subj_B = SbjInfo.participant_id(ismember(SbjInfo.hemi,'LR'))'; % bilateral
    assert(isempty(setdiff(subjects,[subj_L,subj_R,subj_B])),...
        'Data includes unknown subject. Please check ''subjectlist.tsv''');

    %-- Detect channels in left/right hemisphere
    islefths  = ismember(subjects,subj_L) | ...
                (ismember(subjects,subj_B) & startsWith(channels.name,'L'));
    isrighths = ismember(subjects,subj_R) | ...
                (ismember(subjects,subj_B) & startsWith(channels.name,'R'));

    %-- Outputs
    channels.hemisphere(:)         = {'n/a'};
    channels.hemisphere(islefths)  = {'L'};
    channels.hemisphere(isrighths) = {'R'};
    
end

end