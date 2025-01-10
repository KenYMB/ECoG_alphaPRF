% Compare ERP regression
%   Folk from ecog_APRF_07b_checkmodel_fullts
 
% 20240903 Yuasa
 
%% prefix
% close all; clearvars;
% startupToolboxToolbox;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'pRFselections-representative');
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
else
plotsavePth    = 'pRFselections-representative';
prfPth         = 'pRFmodel';
prfstatPth     = 'pRFanalysis';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth), 'modelselection', 'regression', 'refch');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');
 
%% general parameter
clear alphaType broadbandType

average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';
selectchs      = 'wangprobchs';
    allowlag       = false;
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;
va_area = 'wangarea';

usefulltsR2   = false;
usefulltsxR2  = false;
usexvalparams = false;

%% load analyzePRF @DOG,OG,LFS
modelNames  = ["With ERP-Regression","Without ERP-Regression"];
modelFigNames = ["main","cntl"];

modeldata  = struct();
prf_params = struct();
model_all  = struct();
prf_all    = struct();
threshold  = struct();
for ii=1:length(modelNames)
if contains(modelNames(ii),'Without')
    modeldataID = 'freq_spectra-timeseries_noregress';
    prfID       = 'prf_noregress';
else
    modeldataID = 'freq_spectra-timeseries';
    prfID       = 'prf';
end
 
clear alphaType broadbandType
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITc_postfix;
if contains(modelNames(ii),'Without')
postfix = [postfix '_noregress'];
end
ecog_APRFF_INITd_threshold
 
%-- merge loaded data
modeldata(ii).a   = modeldata_a;
modeldata(ii).bb  = modeldata_bb;
prf_params(ii).a  = prf_params_a;
prf_params(ii).bb = prf_params_bb;
 
model_all(ii).a   = model_all_a;
model_all(ii).bb  = model_all_bb;
prf_all(ii).a     = prf_all_a;
prf_all(ii).bb    = prf_all_bb;

threshold(ii).a   = threshold_a;
threshold(ii).bb  = threshold_bb;
end

%% %%%%%%%%%%%%%%
%% Visualize
%% %%%%%%%%%%%%%%
%%% prepare visualize
res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;
 
alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');
hF = gobjects(0);

FntSiz = 24;

boundary = [0 40 56 96 112 152 168 208 224]+.5;         % for BLANK
bndfactor = 1/224*75;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FaCLb vs FaCLb-withERP
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelsets = [2 1];      % w/o vs w/ Regress
modelcolrat = [0.4 1.0];
targetBAND = ('noRegress');

% eclimit = 50;
idat = 1;
okch1bb = ~(prf_all(modelsets(idat)).bb.xval<=threshold(modelsets(idat)).bb | prf_all(modelsets(idat)).bb.ecc >= eclimit);  % use tilde for nan
okch1a  = ~(prf_all(modelsets(idat)).a.xval<=threshold(modelsets(idat)).a | prf_all(modelsets(idat)).a.ecc >= eclimit);  % use tilde for nan
idat = 2;
okch2bb = ~(prf_all(modelsets(idat)).bb.xval<=threshold(modelsets(idat)).bb | prf_all(modelsets(idat)).bb.ecc >= eclimit);  % use tilde for nan
okch2a  = ~(prf_all(modelsets(idat)).a.xval<=threshold(modelsets(idat)).a | prf_all(modelsets(idat)).a.ecc >= eclimit);  % use tilde for nan

okchORbb = okch1bb | okch2bb;
okchORa  = okch1a  | okch2a;

onlych1bb = okch1bb&~okch2bb;
onlych2bb = ~okch1bb&okch2bb;

% close all

%% Select electrodes
if issaveplot
    repelec = {'name',{'Oc17','Oc18','Oc25','GB103'},'subject_name',{'p02','p02','p02','p10'}};
else
    repelec = {'name',{'Oc17','Oc18','Oc25'},'subject_name',{'p02','p02','p02'}};
end

%% Time series
selchidx = findtable(channels,repelec{:});

opt = [];
opt.plot.fontSize   = FntSiz;
opt.plot.addWangToTitle   = 'no';
opt.plot.addBensonToTitle = 'no';
opt.plot.addSbjToTitle    = 'yes';
opt.plot.addR2ToTitle     = 'yes';
opt.skipprojection = 'yes';
opt.plot.reverseYaxis = 'no';
opt.plot.boundaries = boundary.*bndfactor;
opt.plot.options = {'XTick', boundary.*bndfactor,'XTickLabel',{}};
opt.plot.XLim    = [floor(boundary(1)) ceil(boundary(end))].*bndfactor;
opt.plot.nSubPlots = [1, length(selchidx)];

idat=1;
opt.plot.model_options    = {'Color',plcol(1,:).*modelcolrat(idat)};
p = ecog_plotGridPRFts(model_all(modelsets(idat)).bb.datats,model_all(modelsets(idat)).bb.stimulus,prf_all(modelsets(idat)).bb,selchidx, opt);
title(p{1}.Children,modelNames(modelsets(idat)),'FontSize',FntSiz);
p{1}.Position = p{1}.Position .* [1,1,2,2];
hF = [hF, p{1}];
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-ts_%s_%s','broadband',modelFigNames{idat})));  end
idat=2;
opt.plot.model_options    = {'Color',plcol(1,:).*modelcolrat(idat)};
p = ecog_plotGridPRFts(model_all(modelsets(idat)).bb.datats,model_all(modelsets(idat)).bb.stimulus,prf_all(modelsets(idat)).bb,selchidx, opt);
title(p{1}.Children,modelNames(modelsets(idat)),'FontSize',FntSiz);
p{1}.Position = p{1}.Position .* [1,1,2,2];
hF = [hF, p{1}];
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-ts_%s_%s','broadband',modelFigNames{idat})));  end

idat=1;
opt.plot.model_options    = {'Color',plcol(2,:).*modelcolrat(idat)};
p = ecog_plotGridPRFts(model_all(modelsets(idat)).a.datats,model_all(modelsets(idat)).a.stimulus,prf_all(modelsets(idat)).a,selchidx, opt);
title(p{1}.Children,modelNames(modelsets(idat)),'FontSize',FntSiz);
p{1}.Position = p{1}.Position .* [1,1,2,2];
hF = [hF, p{1}];
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-ts_%s_%s','alpha',modelFigNames{idat})));  end
idat=2;
opt.plot.model_options    = {'Color',plcol(2,:).*modelcolrat(idat)};
p = ecog_plotGridPRFts(model_all(modelsets(idat)).a.datats,model_all(modelsets(idat)).a.stimulus,prf_all(modelsets(idat)).a,selchidx, opt);
title(p{1}.Children,modelNames(modelsets(idat)),'FontSize',FntSiz);
p{1}.Position = p{1}.Position .* [1,1,2,2];
hF = [hF, p{1}];
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-ts_%s_%s','alpha',modelFigNames{idat})));  end


%% pRF locations
selchidx = findtable(channels,repelec{:});

opt = [];
opt.plot.pix2deg = cfactor;
opt.plot.XLim    = [-1 1].*12;
opt.plot.YLim    = [-1 1].*12;
opt.plot.addChsToTitle     = 'yes';
opt.plot.addSbjToTitle     = 'yes';
opt.plot.addBensonToTitle  = 'no';
opt.plot.addWangToTitle    = 'no';
opt.plot.fontSize          = FntSiz;
opt.plot.showaxis = 0;         % if show X & Y axis
opt.plot.nSubPlots = [1, length(selchidx)];
opt.plot.legend    = modelNames(modelsets);

opt.plot.colors    = plcol(1,:).*modelcolrat';
p = ecog_plotGridPRF(selchidx, opt, prf_all(modelsets(1)).bb,prf_all(modelsets(2)).bb);
p{1}.Position = p{1}.Position .* [1,1,2,2];
hF = [hF p{:}];
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-loc_%s','broadband')));  end

opt.plot.colors    = plcol(2,:).*modelcolrat';
p = ecog_plotGridPRF(selchidx, opt, prf_all(modelsets(1)).a,prf_all(modelsets(2)).a);
p{1}.Position = p{1}.Position .* [1,1,2,2];
hF = [hF p{:}];
if issaveplot,  savefigauto(gcf,fullfile(plotsavedir,sprintf('pRF-loc_%s','alpha')));  end

%%
% close all;