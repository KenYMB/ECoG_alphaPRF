% Save all channel information without bad channel exclusion

% 20211018 Yuasa

%% Define dataset
close all; clearvars;
%-- Set path
checkPath;

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

nsbj = length(subjectList);
 
%% Load data
%%% selected data
opts = [];
opts.doplots        = false;
opts.outputDir      = 'Preprocessed';
opts.issave         = false;
opts.compute        = false;
[data_post] = ecog_prf_selectData(subjectList,[],opts);

%% Show selected properties
% Select channels
% 1: many bad runs
% 2: many bad epochs
% 3: high variance

for isbj=1:nsbj
fprintf('Subject: p%02d\n',isbj+2*(isbj>2));
fprintf('Bad channels: %d %d %d %d\n',sum(~data_post{isbj}.select_chs(:,1),1),...
                sum(~data_post{isbj}.select_chs(:,2:end),1)-sum(~data_post{isbj}.select_chs(:,1),1));
fprintf('Bad epochs: %d / %d\n',[sum(data_post{isbj}.outliers,'all'),numel(data_post{isbj}.outliers)]);
fprintf('Bad runs:%s\n',sprintf(' %d',find(~all(data_post{isbj}.select_run,2))));
end

%% Show bad epochs in selected channels
bad_epochs = zeros(nsbj,1);
num_epochs = zeros(nsbj,1);
for isbj=1:nsbj
wangflds  = startsWith(data_post{isbj}.channels.Properties.VariableNames,'wangprob_')&~endsWith(data_post{isbj}.channels.Properties.VariableNames,'_FEF');
selectchs = (sum(data_post{isbj}.channels{:,wangflds},2)>0 | ismember(data_post{isbj}.channels.group,'HDgrid'));
    
fprintf('Subject: p%02d\n',isbj+2*(isbj>2));
bad_epochs(isbj) = sum(data_post{isbj}.outliers(:,selectchs),'all');
num_epochs(isbj) = numel(data_post{isbj}.outliers(:,selectchs));
fprintf('Bad epochs: %d / %d\n',bad_epochs(isbj),num_epochs(isbj));
end

fprintf('Total:\n');
fprintf('Bad epochs: %d / %d (%.2f%%)\n',sum(bad_epochs),sum(num_epochs),sum(bad_epochs)./sum(num_epochs).*100);
