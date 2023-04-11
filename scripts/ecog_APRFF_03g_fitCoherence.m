% fit Coherence across distance

% 20220411 Yuasa
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

%% load coherence
%%% Subject Name
nsbj = length(subjectList);
if nsbj==1
   subject = subjectList{1};
else
   subject = 'all';
end

%%% Bootstrapping file
filename = sprintf('%sboot_%s-%s-%s.mat',cohmethod,subject,useChans,disttype);
filepath = fullfile(SetDefaultAnalysisPath('DAT',outputDir),filename);

load(filepath);

%%% Fitting file
filename = sprintf('%sbootfit_%s-%s-%s.mat',cohmethod,subject,useChans,disttype);
filepath = fullfile(SetDefaultAnalysisPath('DAT',outputDir),filename);

if exist(filepath,'file')
    warning('Fitted coherence is already exist: %s',filepath);
    continue;
end

%%
%% %%%%%%%%%%%%%
%% fit curve

cohfit = @(x,xdata)x(:,1).*exp(x(:,2).*xdata)+x(:,3);
nparams = 3;

%-- All
fitparamsBBall = zeros(nboot,nparams);
fitparamsAall = zeros(nboot,nparams);
for iboot=1:nboot
    y1 = cohbootBBall(iboot,:);
        %-- fitting curve
        params0 = [max(y1),-1,min(y1)];
        fitparamsBBall(iboot,:) = lsqcurvefit(cohfit,params0,arounddist,y1);
        
    y1 = cohbootAall(iboot,:);
        %-- fitting curve
        params0 = [max(y1),-1,min(y1)];
        fitparamsAall(iboot,:) = lsqcurvefit(cohfit,params0,arounddist,y1);
end

%-- in-pRF
fitparamsBBprf = zeros(nboot,nparams);
fitparamsAprf = zeros(nboot,nparams);
for iboot=1:nboot
    y1 = cohbootBBprf(iboot,:);
        %-- fitting curve
        params0 = [max(y1),-1,min(y1)];
        fitparamsBBprf(iboot,:) = lsqcurvefit(cohfit,params0,arounddist,y1);
        
    y1 = cohbootAprf(iboot,:);
        %-- fitting curve
        params0 = [max(y1),-1,min(y1)];
        fitparamsAprf(iboot,:) = lsqcurvefit(cohfit,params0,arounddist,y1);
end

%-- out-pRF
fitparamsBBout = zeros(nboot,nparams);
fitparamsAout = zeros(nboot,nparams);
for iboot=1:nboot
    y1 = cohbootBBout(iboot,:);
        %-- fitting curve
        params0 = [max(y1),-1,min(y1)];
        fitparamsBBout(iboot,:) = lsqcurvefit(cohfit,params0,arounddist,y1);
        
    y1 = cohbootAout(iboot,:);
        %-- fitting curve
        params0 = [max(y1),-1,min(y1)];
        fitparamsAout(iboot,:) = lsqcurvefit(cohfit,params0,arounddist,y1);
end

%-- BLANK
fitparamsBBbsl = zeros(nboot,nparams);
fitparamsAbsl = zeros(nboot,nparams);
for iboot=1:nboot
    y1 = cohbootBBbsl(iboot,:);
        %-- fitting curve
        params0 = [max(y1),-1,min(y1)];
        fitparamsBBbsl(iboot,:) = lsqcurvefit(cohfit,params0,arounddist,y1);
        
    y1 = cohbootAbsl(iboot,:);
        %-- fitting curve
        params0 = [max(y1),-1,min(y1)];
        fitparamsAbsl(iboot,:) = lsqcurvefit(cohfit,params0,arounddist,y1);
end

%-- save
saveauto(filepath,'fitparams*','cohfit','disttype','useChans');

end
