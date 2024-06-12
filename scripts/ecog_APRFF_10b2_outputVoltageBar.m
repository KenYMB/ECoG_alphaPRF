% Show spectrum shapes for large Visul Filed stimuli
% for figure 2

% 20210119 Yuasa - update from ecog_APRF_10a_outputSpectrum
% 20210430 Yuasa - same figure from bar stimuli

%%
% close all; clear all;
% startupToolboxToolbox;
run_checkPath;
if contains(which('resample'),'Fieldtrip')
    rmpath(fileparts(which('resample')));
end

%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'spectra-representative');
SetDefault('preprocPth',     'Spectrum');
else
plotsavePth    = 'spectra-representative';
preprocPth     = 'Preprocessed';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'Bar');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end
%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%-- Plotting Setting
FontSize = 18;

%% Load data

  subjectList   = intersect(subjectList,{'chaam','p02'});  subject = subjectList{1};
  outputDir = SetDefaultAnalysisPath('DAT',preprocPth);
  
  %-- raw  
  fileid    = 'data_visualelecs';
  dataname  = fullfile(outputDir, sprintf('%s_%s.mat', subject,fileid));
  idat      = load(dataname);

  %-- set parameters
    t      	  = idat.t;
    channels  = idat.channels;
    events    = idat.events;
    epochs    = idat.epochs;
    fsample   = idat.fsample;

  %-- regressed  
  fileid    = 'data_regress';
  dataname  = fullfile(outputDir, sprintf('%s_%s.mat', subject,fileid));
  idat      = load(dataname);

  %-- set parameters
    regressed = idat.epochs;
  
  %% get selected channels & average across trial
  el = find(ismember(channels.name,{'Oc17','Oc18','Oc25'}))';
  % averagemethod = '2steps';  % '2steps', 'random'
  averagemethod = 'random';
    
  switch averagemethod
  %%% Two step average
      case '2steps'
        %-- take average across repeats (t x averaged events x channels)
        average        = 'stimuli';
        [data_epochs, events_avg] = ecog_averageEvents(epochs,events,average,@mean);
        [data_regress]            = ecog_averageEvents(regressed,events,average,@mean);

        data_epochs  = cat(2, mean(data_epochs(:,1:end-1,:),2,'omitnan'),  data_epochs(:,end,:));
        data_regress = cat(2, mean(data_regress(:,1:end-1,:),2,'omitnan'), data_epochs(:,end,:));

  %%% Random sampling
      case 'random'
        %-- take average across repeats (t x averaged events x channels)
        average        = 'trials';

        seltrl = randsample(find(~ismember(events.trial_name,'BLANK')),sum(~ismember(events.trial_name,'BLANK')),true);
        [data_epochs]             = mean(epochs(:,seltrl,:),2,'omitnan');
        [data_regress]            = mean(regressed(:,seltrl,:),2,'omitnan');
        seltrl = randsample(find(ismember(events.trial_name,'BLANK')),sum(ismember(events.trial_name,'BLANK')),true);
        data_epochs  = cat(2,data_epochs,mean(epochs(:,seltrl,:),2,'omitnan'));
        data_regress = cat(2,data_regress,mean(regressed(:,seltrl,:),2,'omitnan'));
  end
  
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot RAW spectrum in broadband and alpha
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if issaveplot
    pltel = 1:length(el);
else
    pltel = find(ismember(channels.name(el),{'Oc17'}));
end
hF = gobjects(0);
FntSiz = 21;
xshow = [-0.1 0.8];

for ee = pltel

hF(end+1)=figure('Position',[500 500 930 320],'MenuBar','none');
colororder(gcf,[hex2dec(["A0","64","40"])./255;0 0 0]);
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
%-- Plot evoked
nexttile;
h=plot(t,data_epochs(:,:,el(ee)),'LineWidth',2.5);
ha = gca;
ha.FontSize = FntSiz;
xlim(xshow);
if max(data_epochs(t>xshow(1)&t<xshow(2),1,el(ee)))>max(ylim)*.9, ylim(ylim.*[1 1.5]); end
if min(data_epochs(t>xshow(1)&t<xshow(2),1,el(ee)))<min(ylim)*.9, ylim(ylim.*[1.5 1]); end
yshow = ylim;
ha.XAxisLocation = 'origin';
ha.XTick = [];
box off;
title('Evoked Potential');

%-- Plot regressed
nexttile;
h=plot(t,data_regress(:,:,el(ee)),'LineWidth',2.5);
xlim(xshow);
% ylim(yshow.*ceil(range(ylim)*1.6/5)*5./range(yshow));
ylim(yshow);
ha = gca;
ha.FontSize = FntSiz;
ha.XAxisLocation = 'origin';
ha.XTick = [];
box off;
title('Regressed Signal');

figureName = sprintf('samplesignal_Bar-%s_%s_%s',averagemethod,subject,channels.name{el(ee)});
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName),'-vector');  end

end
%%
% close all;
