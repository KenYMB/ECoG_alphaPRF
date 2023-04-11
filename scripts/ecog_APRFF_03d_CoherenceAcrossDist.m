% compute connectivities across electrodes for subjects with HDgrid
%   investigate how to analyze based on a patient across all electrodes

% 20220714 Yuasa - separate non-bootstrapping part
% 20221027 Yuasa - enable to change some variables from outside of the script

%%
close all; % clearvars;
% if isempty(gcp('nocreate')),  parpool([1 40]); end
% startupToolboxToolbox;

%% Define paths and dataset
checkPath;
outputDir      = 'xSpectrum';

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
HDsubjectList = SetSubjectsList(subjectList_fname, 'hasHDgrid','yes');

%%
SetDefault('cohmethod','mscoh');      % 'mscoh', 'imcoh'
SetDefault('disttype','norm');        % 'square','diamond','norm'
SetDefault('useChans','SELchs');      % 'pRFchs', 'SELchs', 'ALLchs'
    
for selsbj = 1:(length(HDsubjectList)+1)
%-- Dataset specs
if selsbj > length(HDsubjectList),  subjectList = HDsubjectList;
else,                               subjectList = HDsubjectList(selsbj);
end
%% load Coherence data
%%% Subject Name
nsbj = length(subjectList);
if nsbj==1
   subject = subjectList{1};
else
   subject = 'all';
end

%%% Coherence across distance file
filename = sprintf('%sdat_%s-%s-%s',cohmethod,subject,useChans,disttype);
filepath = fullfile(SetDefaultAnalysisPath('DAT',outputDir),filename);

if exist(filepath,'file')
    warning('Coherence over distance is already exist: %s',filepath);
    continue;
end

%%% Coherence across distance
opts = [];
opts.subjectname = subject;
opts.method      = cohmethod;
opts.arounddist  = 1:6;
[cohdatBBall,cohdatBBbsl,cohdatBBprf,cohdatBBout,...
             cohdatAall,cohdatAbsl,cohdatAprf,cohdatAout,...
             avgtsBBall,avgtsBBbsl,avgtsBBprf,avgtsBBout,...
             avgtsAall,avgtsAbsl,avgtsAprf,avgtsAout,...
             disttype,useChans,arounddist,channels,selch,...
             smoothingMode,smoothingN] ...
    = ecog_prf_coherencedist(subjectList,disttype,useChans,[],opts);    

%%% Save
saveauto(filepath,'cohdat*','avgts*','disttype','useChans','arounddist','channels','selch','smoothingMode','smoothingN');

end
