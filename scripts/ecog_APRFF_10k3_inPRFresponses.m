% compute ECoG power in/out pRFs

% 20211217 Yuasa

%% prefix
% close all; clearvars;
% startupToolboxToolbox;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'spectra-representative');
SetDefault('prfPth',        'pRFmodel');
else
plotsavePth    = 'spectra-representative';
prfPth         = 'pRFmodel';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'lowbb');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

%-- Plotting Setting
FntSiz = 22;
   
%% %%%%%%%%%%%%%%%%%%%%
%% test
%% %%%%%%%%%%%%%%%%%%%%

% tarBAND     = 'alpha';
%%%
try
%% load time series data
clear alphaType broadbandType
decN = 3;
% decN = 1;

average        ='runs';
prfmodel       = 'linear';
gaussianmode   = 'gs';
smoothingMode  ='decimate';
smoothingN     = decN;
selectchs      = 'wangprobchs';     % only use wangprobchs
    allowlag       = false;
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;
    noACorrect     = false;
    
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITd_threshold;

%% parameter for plot
   %%%
   if length(subjectList)==1
       subject = subjectList{1};
   else
       subject = 'all';
   end
   
   evlng = height(modeldata_bb{end}.events);
   chlng = height(model_all_bb.channels);
   
    decRate = round(evlng/decN)./evlng;
    boundaries = ([0 28 40 56 84 96 112 140 152 168 196 208 224]+0.5).*decRate;
    blankbnd   = boundaries([3,4;6,7;9,10;12,13]);
    medblank   = mean(blankbnd,2);
   evlng2 = round(evlng*decRate);

%% set pRF parameters
datbb = model_all_bb.datats;
data  = model_all_a.datats;

%-- BLANK
blankidx = false(evlng2,1);
for iset = 1:size(blankbnd,1)
    blankidx(ceil(blankbnd(iset,1)):fix(blankbnd(iset,2))) = true;
end
    %%% color for blank
    mkblcl = @(c) mean([c;ones(3,3)*0.6],1);
%     mkblcl = @(c) c;

%-- Hit pRF
    %%% Model computation
    hrf = 1;
    stimulus = modeldata_bb{1}.stimulus;
    res = size(stimulus{1},1,2);
    resmx = max(res);
    numruns = size(stimulus,2);
    %%-- Pre-compute cache for faster execution
    [~,xx,yy] = makegaussian2d(resmx,2,2,2,2);

    %%-- Prepare the stimuli for use in the model
    stimulusPP = repmat({},numruns,1);
    for pp=1:numruns
      stimulusPP{pp} = squish(stimulus{pp},2)';  % this flattens the image so that the dimensionality is now frames x pixels
      stimulusPP{pp} = [stimulusPP{pp} pp*ones(size(stimulusPP{pp},1),1)];  % this adds a dummy column to indicate run breaks
    end

    %%-- Set model
    switch gaussianmode
        case {'dog','gs','lfs'}
            modelfun = @(pp,dd) conv2run(modeldogcss(pp(1:5),pp(6:end),dd,res,xx,yy,0,0),hrf,dd(:,prod(res)+1));
        case {'og'}
            modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));
    end
    
    %%-- apply model
    modelbb = repmat({{}},chlng,1);  modela = repmat({{}},chlng,1);
    dattsbb = repmat({{}},chlng,1);  dattsa = repmat({{}},chlng,1);
    for ich=1:chlng
        for pp=1:numruns
            nanplace    = isnan(datbb{pp}(ich,:));
            dattsbb{ich}{pp}   = datbb{pp}(ich,~nanplace)';
            dattsa{ich}{pp}    = -data{pp}(ich,~nanplace)';
            modelbb{ich}{pp} = modelfun(prf_all_bb.params(1,:,ich),stimulusPP{pp}(~nanplace,:));
            modela{ich}{pp}  = -modelfun(prf_all_a.params(1,:,ich),stimulusPP{pp}(~nanplace,:));
        end
        dattsbb{ich} = cat(1,dattsbb{ich}{:});
        dattsa{ich}  = cat(1,dattsa{ich}{:});
        modelbb{ich} = cat(1,modelbb{ich}{:});
        modela{ich}  = cat(1,modela{ich}{:});
    end
    dattsbb = cat(2,dattsbb{:});
    dattsa  = cat(2,dattsa{:});
    modelbb = cat(2,modelbb{:});
    modela  = cat(2,modela{:});
        
%     gainbb  = prf_all_bb.gain';
%     gaina   = prf_all_a.gain';
    gainbb  = squeeze(prf_all_bb.params(:,4,:) .* (1-prf_all_bb.params(:,7,:)))';  % DoG correction
    gaina   = squeeze(prf_all_a.params(:,4,:) .* (1-prf_all_a.params(:,7,:)))';    % DoG correction
    
    %-- set threshold for pRF inside
    prfthreshbb = gainbb*0.05;
    prfthresha  = -gaina*0.05;
    prfidxbb = false(size(modelbb));  prfidxa = false(size(modela));
    for ich=1:chlng
        prfidxbb(:,ich) = (modelbb(:,ich) > prfthreshbb(ich)) & ~blankidx;
        prfidxa(:,ich)  = (modela(:,ich) < prfthresha(ich)) & ~blankidx;
    end
    
%% load model data (low broadband)
[modeldata_bbl, prf_params_bbl] = ecog_prf_loadprfs(subjectList,'bbL',prfPth,modeldataID,prfID,average,smoothingMode,smoothingN,prfmodel,gaussianmode,selectchs,selectch_exFEF,selectch_thresh,usefulltsR2,usefulltsxR2,skipsummarizeROIs);
[model_all_bbl, prf_all_bbl]    = ecog_prf_mergeprfs(modeldata_bbl,prf_params_bbl,va_area,usexvalparams,rearropt);

%% load model data (no broadband correction)
% [modeldata_a0, prf_params_a0] = ecog_prf_loadprfs(subjectList,'FaLb',prfPth,modeldataID,prfID,average,smoothingMode,smoothingN,prfmodel,gaussianmode,selectchs,selectch_exFEF,selectch_thresh,usefulltsR2,usefulltsxR2,skipsummarizeROIs);
% [model_all_a0, prf_all_a0]    = ecog_prf_mergeprfs(modeldata_bbl,prf_params_bbl,va_area,usexvalparams,rearropt);

catch ex
    rethrow(ex);
end

 %%
 nincl = 6;
 prfch = sum(prfidxbb,1) > nincl & sum(prfidxa,1) > nincl;
 acuch =  ~(prf_all_bb.xval<=threshold_bb | prf_all_a.xval<=threshold_a ...
            | prf_all_bb.ecc >= eclimit | prf_all_a.ecc >= eclimit)';  % use tilde for nan
%  prfch = sum(prfidxbb,1) > nincl;
%  acuch =  ~(prf_all_bb.xval<=threshold_bb | prf_all_bb.ecc >= eclimit)';  % use tilde for nan
%   prfch = sum(prfidxa,1) > nincl;
%   acuch =  ~(prf_all_a.xval<=threshold_a | prf_all_a.ecc >= eclimit)';  % use tilde for nan
 
 selch = prfch & acuch;
 
%%
% close all;
%%
%% %%%%%%%%%%%%%
%% Visualization in ROIs
useChans = 'pRFchs';        % pRFchs, SELchs, ALLchs

switch useChans
    case 'pRFchs',      selch = prfch & acuch;
    case 'SELchs',      selch = acuch;
    case 'ALLchs',      selch = true(size(acuch));
end

alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');
hF = gobjects(0);

%-- Set ROIs
% rois = {{'V1'},{'V2'},{'V3'},{'V3a','V3b','LO1','LO2','TO','IPS'}};
% roilabels = {'V1','V2','V3','Dorsolateral'};
if issaveplot
rois = {{'V1'},{'V2'},{'V3'},{'V1','V2','V3'},{'V3a','V3b','LO1','LO2','TO','IPS'}};
roilabels = {'V1','V2','V3','V1-V3','Dorsolateral'};
maxpow = 450;   maxlog = 5.5;   minlog = 0.3;
else
rois = {{'V1','V2','V3'},{'V3a','V3b','LO1','LO2','TO','IPS'}};
roilabels = {'V1-V3','Dorsolateral'};
maxpow = 300;   maxlog = 4.2;   minlog = 0.4;
end

nroi = length(rois);
ix = ceil(sqrt(nroi));    iy = floor(sqrt(nroi));

%% %%%%%%%%%%%%%%%%%%%%%%
%% Show broadband
%% %%%%%%%%%%%%%%%%%%%%%%
ylin = [-50 maxpow];
yticklog = [-100:100:maxpow];
    
%% Bar plot (low broadband VS broadband)
if issaveplot
datats = cat(2,model_all_bb.datats{:})'./100+1;
datats2 = cat(2,model_all_bbl.datats{:})'./100+1;
datats(datats<=0) = nan; datats2(datats2<=0) = nan;

meanfun = @(d) (geomean(d,'all','omitnan')-1)*100;
errfun  = @(d) (geosem(d,0,'all','omitnan'))*100;
errneg  = @(m,s) -(exp(log(m./100+1)-log(s./100))-1)*100 + m;
errpos  = @(m,s) (exp(log(m./100+1)+log(s./100))-1)*100 - m;

hF(end+1) = figure; tiledlayout(iy,ix,'Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[1./ix 1./iy 0.9.*ix 0.9.*iy]);
for iroi = 1:length(rois)
    roich = ismember(channels.wangarea,rois{iroi})';

inPRF_bb  = meanfun(datats(roich&selch&prfidxbb&~blankidx));
outPRF_bb = meanfun(datats(roich&selch&~prfidxbb&~blankidx));
E_inPRF_bb   = errfun(datats(roich&selch&prfidxbb&~blankidx));
E_outPRF_bb  = errfun(datats(roich&selch&~prfidxbb&~blankidx));

inPRF_bb2  = meanfun(datats2(roich&selch&prfidxbb&~blankidx));
outPRF_bb2 = meanfun(datats2(roich&selch&~prfidxbb&~blankidx));
E_inPRF_bb2   = errfun(datats2(roich&selch&prfidxbb&~blankidx));
E_outPRF_bb2 = errfun(datats2(roich&selch&~prfidxbb&~blankidx));

nexttile;
B=bar([inPRF_bb,outPRF_bb;inPRF_bb2,outPRF_bb2]');

hold on;
errorbar(B(1).XEndPoints,[inPRF_bb,outPRF_bb],...
            errneg([inPRF_bb,outPRF_bb],[E_inPRF_bb,E_outPRF_bb]),...
            errpos([inPRF_bb,outPRF_bb],[E_inPRF_bb,E_outPRF_bb]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');
errorbar(B(2).XEndPoints,[inPRF_bb2,outPRF_bb2],...
            errneg([inPRF_bb2,outPRF_bb2],[E_inPRF_bb2,E_outPRF_bb2]),...
            errpos([inPRF_bb2,outPRF_bb2],[E_inPRF_bb2,E_outPRF_bb2]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');

B(1).FaceColor = plcol(1,:);
B(2).FaceColor = plcol(5,:);
hold off;
ylim(ylin);
yticks(yticklog);
set(gca,'FontSize',FntSiz);
xticklabels({'in pRF','out pRF'});
% xticklabels({'In_{bb}','Out_{bb}','In_a','Out_a','Border'});
% xticklabels({'In_{bb}','Border','Out_a'});
% xticklabels({'Inside Broadband pRF','Outside Broadband pRF','Inside Alpja pRF','Outside Alpja pRF','Outside Broadband pRF & Inside Alpja pRF'});
title(roilabels{iroi});

end
legend({'Broadband (70–180 Hz)','Broadband (3–26 Hz)'});
set(gcf,'Name','ECoG Power');

figname = sprintf('Power-lowbb_%02d%%-%02d%%-ecc%02d_%s',threshold_bb,threshold_a,eclimit,useChans);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end

%% %%%%%%%%%%%%%%%%%%%%%%
%% Show alpha
%% %%%%%%%%%%%%%%%%%%%%%%
ylin = [0.5 maxlog];
yticklog = [-1:maxlog];

%% Bar plot (low broadband VS broadband VS alpha)
if issaveplot
datats = cat(2,model_all_bb.datats{:})'./100+1;
datats2 = cat(2,model_all_bbl.datats{:})'./100+1;
datats3 = exp(cat(2,model_all_a.datats{:})');
datats(datats<=0) = nan; datats2(datats2<=0) = nan;

meanfun = @(d) (geomean(d,'all','omitnan'));
errfun  = @(d) (geosem(d,0,'all','omitnan'));
errneg  = @(m,s) -m./s + m;   % -exp(log(m)-log(s)) + m
errpos  = @(m,s)  m.*s - m;   %  exp(log(m)+log(s)) - m

nchn      = size(datats,2);
prfinmat  = prfidxbb&~blankidx;
prfoutmat = ~prfidxbb&~blankidx;
datatsIn  = nan(1,nchn);   datatsOut  = nan(1,nchn);
datatsIn2  = nan(1,nchn);  datatsOut2  = nan(1,nchn);
datatsIn3  = nan(1,nchn);  datatsOut3  = nan(1,nchn);
for ich=1:size(datats,2)
  if sum(prfinmat(:,ich))>0
    datatsIn(ich) = meanfun(datats(prfinmat(:,ich),ich));
    datatsIn2(ich) = meanfun(datats2(prfinmat(:,ich),ich));
    datatsIn3(ich) = meanfun(datats3(prfinmat(:,ich),ich));
  end
  if sum(prfoutmat(:,ich))>0
    datatsOut(ich) = meanfun(datats(prfoutmat(:,ich),ich));
    datatsOut2(ich) = meanfun(datats2(prfoutmat(:,ich),ich));
    datatsOut3(ich) = meanfun(datats3(prfoutmat(:,ich),ich));
  end
end

hF(end+1) = figure; tiledlayout(iy,ix,'Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[1./ix 1./iy 0.9.*ix 0.9.*iy]);
for iroi = 1:length(rois)
    roich = ismember(channels.wangarea,rois{iroi})';

inPRF_bb  = meanfun(datatsIn(roich&selch));
outPRF_bb = meanfun(datatsOut(roich&selch));
E_inPRF_bb   = errfun(datatsIn(roich&selch));
E_outPRF_bb  = errfun(datatsOut(roich&selch));

inPRF_bb2  = meanfun(datatsIn2(roich&selch));
outPRF_bb2 = meanfun(datatsOut2(roich&selch));
E_inPRF_bb2   = errfun(datatsIn2(roich&selch));
E_outPRF_bb2 = errfun(datatsOut2(roich&selch));

inPRF_bb3  = meanfun(datatsIn3(roich&selch));
outPRF_bb3 = meanfun(datatsOut3(roich&selch));
E_inPRF_bb3   = errfun(datatsIn3(roich&selch));
E_outPRF_bb3 = errfun(datatsOut3(roich&selch));

nexttile;
B=bar([inPRF_bb,outPRF_bb;inPRF_bb2,outPRF_bb2;inPRF_bb3,outPRF_bb3]','BaseValue',1);

hold on;
errorbar(B(1).XEndPoints,[inPRF_bb,outPRF_bb],...
            errneg([inPRF_bb,outPRF_bb],[E_inPRF_bb,E_outPRF_bb]),...
            errpos([inPRF_bb,outPRF_bb],[E_inPRF_bb,E_outPRF_bb]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');
errorbar(B(2).XEndPoints,[inPRF_bb2,outPRF_bb2],...
            errneg([inPRF_bb2,outPRF_bb2],[E_inPRF_bb2,E_outPRF_bb2]),...
            errpos([inPRF_bb2,outPRF_bb2],[E_inPRF_bb2,E_outPRF_bb2]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');
errorbar(B(3).XEndPoints,[inPRF_bb3,outPRF_bb3],...
            errneg([inPRF_bb3,outPRF_bb3],[E_inPRF_bb3,E_outPRF_bb3]),...
            errpos([inPRF_bb3,outPRF_bb3],[E_inPRF_bb3,E_outPRF_bb3]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');

B(1).FaceColor = plcol(1,:);
B(2).FaceColor = plcol(5,:);
B(3).FaceColor = plcol(2,:);
hold off;
ylim(ylin);
yticks(yticklog);
set(gca,'FontSize',FntSiz);
xticklabels({'in pRF','out pRF'});
% xticklabels({'In_{bb}','Out_{bb}','In_a','Out_a','Border'});
% xticklabels({'In_{bb}','Border','Out_a'});
% xticklabels({'Inside Broadband pRF','Outside Broadband pRF','Inside Alpja pRF','Outside Alpja pRF','Outside Broadband pRF & Inside Alpja pRF'});
title(roilabels{iroi});

end
legend({'Broadband (70–180 Hz)','Broadband (3–26 Hz)','Alpha'});
set(gcf,'Name','ECoG Log Power');

figname = sprintf('Power-lowbb+alpha_%02d%%-%02d%%-ecc%02d_%s',threshold_bb,threshold_a,eclimit,useChans);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end

%% Bar plot (low broadband VS broadband VS alpha) in pRF only
if issaveplot
datats = cat(2,model_all_bb.datats{:})'./100+1;
datats2 = cat(2,model_all_bbl.datats{:})'./100+1;
datats3 = exp(cat(2,model_all_a.datats{:})');
datats(datats<=0) = nan; datats2(datats2<=0) = nan;

meanfun = @(d) (geomean(d,'all','omitnan'));
errfun  = @(d) (geosem(d,0,'all','omitnan'));
errneg  = @(m,s) -m./s + m;   % -exp(log(m)-log(s)) + m
errpos  = @(m,s)  m.*s - m;   %  exp(log(m)+log(s)) - m

nchn      = size(datats,2);
prfinmat  = prfidxbb&~blankidx;
prfoutmat = ~prfidxbb&~blankidx;
datatsIn  = nan(1,nchn);   datatsOut  = nan(1,nchn);
datatsIn2  = nan(1,nchn);  datatsOut2  = nan(1,nchn);
datatsIn3  = nan(1,nchn);  datatsOut3  = nan(1,nchn);
for ich=1:size(datats,2)
  if sum(prfinmat(:,ich))>0
    datatsIn(ich) = meanfun(datats(prfinmat(:,ich),ich));
    datatsIn2(ich) = meanfun(datats2(prfinmat(:,ich),ich));
    datatsIn3(ich) = meanfun(datats3(prfinmat(:,ich),ich));
  end
  if sum(prfoutmat(:,ich))>0
    datatsOut(ich) = meanfun(datats(prfoutmat(:,ich),ich));
    datatsOut2(ich) = meanfun(datats2(prfoutmat(:,ich),ich));
    datatsOut3(ich) = meanfun(datats3(prfoutmat(:,ich),ich));
  end
end

hF(end+1) = figure; tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[2./(nroi+1) 1 0.3.*nroi 1]);
inPRF_bb = []; E_inPRF_bb = []; inPRF_bb2 = []; E_inPRF_bb2 = []; inPRF_bb3 = []; E_inPRF_bb3 = [];
for iroi = 1:length(rois)
    roich = ismember(channels.wangarea,rois{iroi})';

inPRF_bb  = [inPRF_bb, meanfun(datatsIn(roich&selch))];
E_inPRF_bb   = [E_inPRF_bb, errfun(datatsIn(roich&selch))];

inPRF_bb2  = [inPRF_bb2, meanfun(datatsIn2(roich&selch))];
E_inPRF_bb2   = [E_inPRF_bb2, errfun(datatsIn2(roich&selch))];

inPRF_bb3  = [inPRF_bb3, meanfun(datatsIn3(roich&selch))];
E_inPRF_bb3   = [E_inPRF_bb3, errfun(datatsIn3(roich&selch))];
end

nexttile;
B=bar([inPRF_bb;inPRF_bb2;inPRF_bb3]','BaseValue',1);

hold on;
errorbar(B(1).XEndPoints,[inPRF_bb],...
            errneg([inPRF_bb],[E_inPRF_bb]),...
            errpos([inPRF_bb],[E_inPRF_bb]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');
errorbar(B(2).XEndPoints,[inPRF_bb2],...
            errneg([inPRF_bb2],[E_inPRF_bb2]),...
            errpos([inPRF_bb2],[E_inPRF_bb2]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');
errorbar(B(3).XEndPoints,[inPRF_bb3],...
            errneg([inPRF_bb3],[E_inPRF_bb3]),...
            errpos([inPRF_bb3],[E_inPRF_bb3]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');

B(1).FaceColor = plcol(1,:);
B(2).FaceColor = plcol(5,:);
B(3).FaceColor = plcol(2,:);
hold off;
ylim(ylin);
yticks(yticklog);
set(gca,'FontSize',FntSiz);
xticklabels(roilabels);
% xticklabels({'In_{bb}','Out_{bb}','In_a','Out_a','Border'});
% xticklabels({'In_{bb}','Border','Out_a'});
% xticklabels({'Inside Broadband pRF','Outside Broadband pRF','Inside Alpja pRF','Outside Alpja pRF','Outside Broadband pRF & Inside Alpja pRF'});
title('in pRF');

legend({['Broadband' newline '(70–180 Hz)'],['Broadband' newline '(3–26 Hz)'],'Alpha'});
set(gcf,'Name','ECoG Log Power');

figname = sprintf('Power-lowbb+alpha-inPRF_%02d%%-%02d%%-ecc%02d_%s',threshold_bb,threshold_a,eclimit,useChans);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end

%% %%%%%%%%%%%%%%%%%%%%%%
%% Show alpha in log axis
%% %%%%%%%%%%%%%%%%%%%%%%
scalefun = @(d) log2(d)+1;      % log axis
% scalefun = @(d) -(1./d)+2;      % reciprocal axis

ylin = [scalefun(minlog) maxlog];
yticklog = [scalefun([0.25 0.5]) 1:maxlog];
ytickloglabel = [[0.25 0.5] 1:maxlog];

%% Bar plot (low broadband VS broadband VS alpha) in pRF only <dual axis>
if issaveplot
datats = cat(2,model_all_bb.datats{:})'./100+1;
datats2 = cat(2,model_all_bbl.datats{:})'./100+1;
datats3 = 10.^(cat(2,model_all_a.datats{:})');
datats(datats<=0) = nan; datats2(datats2<=0) = nan;

meanfun = @(d) (geomean(d,'all','omitnan'));
errfun  = @(d) (geosem(d,0,'all','omitnan'));
errneg  = @(m,s) -m./s + m;   % -exp(log(m)-log(s)) + m
errpos  = @(m,s)  m.*s - m;   %  exp(log(m)+log(s)) - m
revfun  = @(d) scalefun(1./d);
revneg  = @(m,s) -revfun(m.*s) + revfun(m);  % -revfun(exp(log(m)+log(s))) + revfun(m)
revpos  = @(m,s)  revfun(m./s) - revfun(m);  %  revfun(exp(log(m)-log(s))) - revfun(m)

nchn      = size(datats,2);
prfinmat  = prfidxbb&~blankidx;
prfoutmat = ~prfidxbb&~blankidx;
datatsIn  = nan(1,nchn);   datatsOut  = nan(1,nchn);
datatsIn2  = nan(1,nchn);  datatsOut2  = nan(1,nchn);
datatsIn3  = nan(1,nchn);  datatsOut3  = nan(1,nchn);
for ich=1:size(datats,2)
  if sum(prfinmat(:,ich))>0
    datatsIn(ich) = meanfun(datats(prfinmat(:,ich),ich));
    datatsIn2(ich) = meanfun(datats2(prfinmat(:,ich),ich));
    datatsIn3(ich) = meanfun(datats3(prfinmat(:,ich),ich));
  end
  if sum(prfoutmat(:,ich))>0
    datatsOut(ich) = meanfun(datats(prfoutmat(:,ich),ich));
    datatsOut2(ich) = meanfun(datats2(prfoutmat(:,ich),ich));
    datatsOut3(ich) = meanfun(datats3(prfoutmat(:,ich),ich));
  end
end

hF(end+1) = figure; tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[2./(nroi+1) 1 0.3.*nroi 1]);
inPRF_bb = []; E_inPRF_bb = []; inPRF_bb2 = []; E_inPRF_bb2 = []; inPRF_bb3 = []; E_inPRF_bb3 = [];
for iroi = 1:length(rois)
    roich = ismember(channels.wangarea,rois{iroi})';

inPRF_bb  = [inPRF_bb, meanfun(datatsIn(roich&selch))];
E_inPRF_bb   = [E_inPRF_bb, errfun(datatsIn(roich&selch))];

inPRF_bb2  = [inPRF_bb2, meanfun(datatsIn2(roich&selch))];
E_inPRF_bb2   = [E_inPRF_bb2, errfun(datatsIn2(roich&selch))];

inPRF_bb3  = [inPRF_bb3, meanfun(datatsIn3(roich&selch))];
E_inPRF_bb3   = [E_inPRF_bb3, errfun(datatsIn3(roich&selch))];
end

nexttile;
B=bar([inPRF_bb;inPRF_bb2;revfun(inPRF_bb3)]','BaseValue',1);

hold on;
errorbar(B(1).XEndPoints,[inPRF_bb],...
            errneg([inPRF_bb],[E_inPRF_bb]),...
            errpos([inPRF_bb],[E_inPRF_bb]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');
errorbar(B(2).XEndPoints,[inPRF_bb2],...
            errneg([inPRF_bb2],[E_inPRF_bb2]),...
            errpos([inPRF_bb2],[E_inPRF_bb2]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');
errorbar(B(3).XEndPoints,revfun([inPRF_bb3]),...
            revneg([inPRF_bb3],[E_inPRF_bb3]),...
            revpos([inPRF_bb3],[E_inPRF_bb3]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');

B(1).FaceColor = plcol(1,:);
B(2).FaceColor = plcol(5,:);
B(3).FaceColor = plcol(2,:);
hold off;
ylim(ylin);
yticks(yticklog);
yticklabels(arrayfun(@num2str, ytickloglabel, 'UniformOutput', false));
    yticklabels(strrep(yticklabels,num2str(1/8),'1/8'));
    yticklabels(strrep(yticklabels,num2str(1/4),'1/4'));
    yticklabels(strrep(yticklabels,num2str(1/2),'1/2'));
set(gca,'FontSize',FntSiz);
xticklabels(roilabels);
% xticklabels({'In_{bb}','Out_{bb}','In_a','Out_a','Border'});
% xticklabels({'In_{bb}','Border','Out_a'});
% xticklabels({'Inside Broadband pRF','Outside Broadband pRF','Inside Alpja pRF','Outside Alpja pRF','Outside Broadband pRF & Inside Alpja pRF'});
title('in pRF');

legend({['Broadband' newline '(70–180 Hz)'],['Broadband' newline '(3–26 Hz)'],'Alpha'});
set(gcf,'Name','ECoG Log Power');

figname = sprintf('Power-lowbb+logalpha-inPRF_%02d%%-%02d%%-ecc%02d_%s',threshold_bb,threshold_a,eclimit,useChans);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
end

%% %%%%%%%%%%%%%%%%%%%%%%
%% Plot in log axis
%% %%%%%%%%%%%%%%%%%%%%%%

ylin = [minlog maxlog+0.6];
yticklog = [[0.25 0.5] 2.^(0:log2(maxlog))];
ytickloglabel = [[0.25 0.5] 2.^(0:log2(maxlog))];

%% Bar plot (low broadband VS broadband VS alpha) in pRF only <dual axis>
datats = cat(2,model_all_bb.datats{:})'./100+1;
datats2 = cat(2,model_all_bbl.datats{:})'./100+1;
datats3 = 10.^-(cat(2,model_all_a.datats{:})');
datats(datats<=0) = nan; datats2(datats2<=0) = nan;

meanfun = @(d) (geomean(d,'all','omitnan'));
errfun  = @(d) (geosem(d,0,'all','omitnan'));
errneg  = @(m,s) -m./s + m;   % -exp(log(m)-log(s)) + m
errpos  = @(m,s)  m.*s - m;   %  exp(log(m)+log(s)) - m

nchn      = size(datats,2);
prfinmat  = prfidxbb&~blankidx;
prfoutmat = ~prfidxbb&~blankidx;
datatsIn  = nan(1,nchn);   datatsOut  = nan(1,nchn);
datatsIn2  = nan(1,nchn);  datatsOut2  = nan(1,nchn);
datatsIn3  = nan(1,nchn);  datatsOut3  = nan(1,nchn);
for ich=1:size(datats,2)
  if sum(prfinmat(:,ich))>0
    datatsIn(ich) = meanfun(datats(prfinmat(:,ich),ich));
    datatsIn2(ich) = meanfun(datats2(prfinmat(:,ich),ich));
    datatsIn3(ich) = meanfun(datats3(prfinmat(:,ich),ich));
  end
  if sum(prfoutmat(:,ich))>0
    datatsOut(ich) = meanfun(datats(prfoutmat(:,ich),ich));
    datatsOut2(ich) = meanfun(datats2(prfoutmat(:,ich),ich));
    datatsOut3(ich) = meanfun(datats3(prfoutmat(:,ich),ich));
  end
end

hF(end+1) = figure; tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[2./(nroi+1) 1 0.3.*nroi 1]);
inPRF_bb = []; E_inPRF_bb = []; inPRF_bb2 = []; E_inPRF_bb2 = []; inPRF_bb3 = []; E_inPRF_bb3 = [];
for iroi = 1:length(rois)
    roich = ismember(channels.wangarea,rois{iroi})';

inPRF_bb  = [inPRF_bb, meanfun(datatsIn(roich&selch))];
E_inPRF_bb   = [E_inPRF_bb, errfun(datatsIn(roich&selch))];

inPRF_bb2  = [inPRF_bb2, meanfun(datatsIn2(roich&selch))];
E_inPRF_bb2   = [E_inPRF_bb2, errfun(datatsIn2(roich&selch))];

inPRF_bb3  = [inPRF_bb3, meanfun(datatsIn3(roich&selch))];
E_inPRF_bb3   = [E_inPRF_bb3, errfun(datatsIn3(roich&selch))];
end

nexttile;
B=bar([inPRF_bb;inPRF_bb2;inPRF_bb3]','BaseValue',1);

hold on;
errorbar(B(1).XEndPoints,[inPRF_bb],...
            errneg([inPRF_bb],[E_inPRF_bb]),...
            errpos([inPRF_bb],[E_inPRF_bb]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');
errorbar(B(2).XEndPoints,[inPRF_bb2],...
            errneg([inPRF_bb2],[E_inPRF_bb2]),...
            errpos([inPRF_bb2],[E_inPRF_bb2]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');
errorbar(B(3).XEndPoints,[inPRF_bb3],...
            errneg([inPRF_bb3],[E_inPRF_bb3]),...
            errpos([inPRF_bb3],[E_inPRF_bb3]),...
    '.','LineStyle','none',...
    'LineWidth',2.5,'CapSize',8,'Color','k');

B(1).FaceColor = plcol(1,:);
B(2).FaceColor = plcol(5,:);
B(3).FaceColor = plcol(2,:);
hold off;
ylim(ylin);
yticks(yticklog);
yticklabels(arrayfun(@num2str, ytickloglabel, 'UniformOutput', false));
    yticklabels(strrep(yticklabels,num2str(1/8),'1/8'));
    yticklabels(strrep(yticklabels,num2str(1/4),'1/4'));
    yticklabels(strrep(yticklabels,num2str(1/2),'1/2'));
set(gca,'FontSize',FntSiz,'YScale','log');
% set(gca,'FontSize',FntSiz);
xticklabels(roilabels);
% xticklabels({'In_{bb}','Out_{bb}','In_a','Out_a','Border'});
% xticklabels({'In_{bb}','Border','Out_a'});
% xticklabels({'Inside Broadband pRF','Outside Broadband pRF','Inside Alpja pRF','Outside Alpja pRF','Outside Broadband pRF & Inside Alpja pRF'});
title('in pRF');

legend({['Broadband' newline '(70–180 Hz)'],['Broadband' newline '(3–26 Hz)'],'Alpha'});
set(gcf,'Name','ECoG Log Power');

figname = sprintf('Power-loglowbb+logalpha-inPRF_%02d%%-%02d%%-ecc%02d_%s',threshold_bb,threshold_a,eclimit,useChans);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end
