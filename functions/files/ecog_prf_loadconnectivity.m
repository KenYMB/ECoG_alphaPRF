function coh_all = ecog_prf_loadconnectivity(subjectList,opts)

% Description: 
%   Load cross-spectral data during pRF experiment based on
%   ecog_prf_connectivity which are computed separately based on stimulus
%   categories, and concatenate them.
% 
% [coh_all] = ecog_prf_loadconnectivity(subjectList,opts)

% 20221109 Yuasa - segregate from ecog_prf_coherencedist

%--Define inputs 
SetDefault('subjectList',{},'cell');

%%% Subject Name
nsbj = length(subjectList);

%% Load Coherence
opts.compute    = false;
opts.issave     = false;
opts.stimNames      = 'HORIZONTAL*';
[coh_part1]     = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'VERTICAL*';
[coh_part2]     = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'DIAGONAL*';
[coh_part3]     = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'BLANK';
[coh_part0]     = ecog_prf_connectivity(subjectList, opts);

%% Concatenate and reorder events
coh_all = cell(size(coh_part0));
for isbj = 1:nsbj
%-- broadband
coh_all{isbj}  = coh_part0{isbj};
coh_all{isbj}.connectivity = cat(3,coh_all{isbj}.connectivity,coh_part1{isbj}.connectivity,coh_part2{isbj}.connectivity,coh_part3{isbj}.connectivity);
coh_all{isbj}.events       = cat(1,coh_all{isbj}.events,coh_part1{isbj}.events,coh_part2{isbj}.events,coh_part3{isbj}.events);

eventslist = sortrows([coh_all{isbj}.events.stim_file_index, [1:height(coh_all{isbj}.events)]']);

coh_all{isbj}.connectivity = coh_all{isbj}.connectivity(:,:,eventslist(:,2),:);
coh_all{isbj}.events       = coh_all{isbj}.events(eventslist(:,2),:);

end