% compute connectivities across electrodes for subjects with HDgrid
%   investigate how to analyze based on a patient across all electrodes

% 20220714 Yuasa - separate non-bootstrapping part

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
disttype    = 'norm';            % 'square','diamond','norm'
% useChans = 'pRFchs';        % pRFchs, SELchs, ALLchs
useChans = 'SELchs';        % pRFchs, SELchs, ALLchs
% arounddist  = [1 2 3 6];
    arounddist  = 1:6;
    
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

%%% Coherence across distance
opts = [];
opts.subjectname = subject;
opts.outputDir   = outputDir;
opts.issave      = true;
[cohdatBBall,cohdatBBbsl,cohdatBBprf,cohdatBBout,...
             cohdatAall,cohdatAbsl,cohdatAprf,cohdatAout,...
             avgtsBBall,avgtsBBbsl,avgtsBBprf,avgtsBBout,...
             avgtsAall,avgtsAbsl,avgtsAprf,avgtsAout,...
             disttype,useChans,arounddist,channels,selch] ...
    = ecog_prf_coherencedist(subjectList,arounddist,disttype,useChans,[],opts);
    
end
