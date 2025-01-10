% Show spectrum shapes for large Visul Filed stimuli

% 20210119 Yuasa - update from ecog_APRF_10a_outputSpectrum
% 20210430 Yuasa - same figure from bar stimuli
% 20241118 Yuasa - for single trial

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
  el = find(ismember(channels.name,{'Oc17'}))';
  averagemethod = 'single';
    
% el = 12;  % Oc17
% reftrl = 8, 13, 14, 64, 144, 192, 207, 369
% bsltrl = 41, 44, 267, 272, 274, 328


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot RAW spectrum in broadband and alpha
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if issaveplot
    pltel = 1:length(el);
else
    pltel = find(ismember(channels.name(el),{'Oc17'}));
end
hF = gobjects(0);
FntSiz = 24;
xshow = [-0.1 0.6];

%-- Oc17
ee = pltel(1);
if issaveplot
    seltrl = [13,144,369; 44,267,328];
else
    seltrl = [13; 44];
end


%-- Grouping in BLANK/STIMULUS
if true
hF(end+1)=figure; tiledlayout('flow','Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[1,1,2.3,size(seltrl,2)]);
for ii=1:size(seltrl,2)
nexttile(2+(ii-1)*2); hold on; xlim(xshow);
plcol = 'A06440'; plcol = sscanf(plcol, '%2x')'./255;
    itrl = seltrl(1,ii);
    plot(t,epochs(:,itrl,el(ee)),'-.','Color',plcol,'LineWidth',2);
    plot(t,regressed(:,itrl,el(ee)),'-','Color',plcol,'LineWidth',2);
title('STIMULUS','FontSize',FntSiz);
ylin = ylim;
nexttile(1+(ii-1)*2); hold on; xlim(xshow); ylim(ylin);
plcol = [0,0,0];
    itrl = seltrl(2,ii);
    plot(t,epochs(:,itrl,el(ee)),'-.','Color',plcol,'LineWidth',2);
    plot(t,regressed(:,itrl,el(ee)),'-','Color',plcol,'LineWidth',2);
title('BLANK','FontSize',FntSiz);
if ii==1
legend({'Raw signal','Regressed out ERP'},'Location','southeast');
end
set(findall(gcf,'-property','FontSize'),'FontSize',FntSiz);
end
for ii=1:2:size(seltrl,2)
    nexttile(ii+1); ylin = ylim;
    nexttile(ii);   ylim(ylin);
end

figureName = sprintf('samplesignal_Bar-%s-blankVSstimulus_%s_%s',averagemethod,subject,channels.name{el(ee)});
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName),'-vector');  end

end


%-- Grouping in with/without ERP regression
if issaveplot
hF(end+1)=figure; tiledlayout('flow','Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[1,1,2.3,1]);
if size(seltrl,2)>1, adjrat = linspace(0,0.5,size(seltrl,2));
else,                adjrat = 0;
end
adjcol = [.8,.8,.8];
nexttile; hold on; xlim(xshow); 
plcol = [0,0,0];
for ii=1:size(seltrl,2)
    itrl = seltrl(2,ii);
    plot(t,epochs(:,itrl,el(ee)),'Color',plcol.*(1-adjrat(ii))+adjcol.*adjrat(ii),'LineWidth',1.5);
end
plcol = 'A06440'; plcol = sscanf(plcol, '%2x')'./255;
for ii=1:size(seltrl,2)
    itrl = seltrl(1,ii);
    plot(t,epochs(:,itrl,el(ee)),'Color',plcol.*(1-adjrat(ii))+adjcol.*adjrat(ii),'LineWidth',1.5);
end
title('Raw signal','FontSize',FntSiz);
ylin = ylim;
nexttile; hold on; xlim(xshow); ylim(ylin);
plcol = [0,0,0];
for ii=1:size(seltrl,2)
    itrl = seltrl(2,ii);
    plot(t,regressed(:,itrl,el(ee)),'Color',plcol.*(1-adjrat(ii))+adjcol.*adjrat(ii),'LineWidth',1.5);
end
plcol = 'A06440'; plcol = sscanf(plcol, '%2x')'./255;
for ii=1:size(seltrl,2)
    itrl = seltrl(1,ii);
    plot(t,regressed(:,itrl,el(ee)),'Color',plcol.*(1-adjrat(ii))+adjcol.*adjrat(ii),'LineWidth',1.5);
end
title('Regressed out ERP','FontSize',FntSiz);

figureName = sprintf('samplesignal_Bar-%s-ERPregression_%s_%s',averagemethod,subject,channels.name{el(ee)});
set(gcf,'Name',figureName);
if issaveplot,  savefigauto(gcf, fullfile(plotsavedir, figureName),'-vector');  end

end

%%
% close all;
