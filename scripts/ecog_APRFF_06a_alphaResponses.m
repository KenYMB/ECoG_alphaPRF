% compute ECoG power in/out pRFs & compare alpha model

% 20211210 Yuasa

%% %%%%%%%%%%%%%%%%%%%%
%% test
%% %%%%%%%%%%%%%%%%%%%%
%%
 close all; clear all;
% if isempty(gcp('nocreate')),  parpool([1 40]); end
% startupToolboxToolbox;
%%
checkPath;
%-- Input & Output path
plotsavePth    = 'ConnectivityTS';
%-- Set save figure dirctory
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'Distance');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'all');

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
    
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITd_threshold;

%% parameter for plot
   FntSiz = 20;
   alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd
   plcol = get(groot,'defaultAxesColorOrder');
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
    
catch ex
    rethrow(ex);
end

 %%
 nincl = 6;
%  prfch = sum(prfidxbb,1) > nincl & sum(prfidxa,1) > nincl;
%  acuch =  ~(prf_all_bb.xval<=threshold_bb | prf_all_a.xval<=threshold_a ...
%             | prf_all_bb.ecc >= eclimit | prf_all_a.ecc >= eclimit)';  % use tilde for nan
 prfch = sum(prfidxbb,1) > nincl;
 acuch =  ~(prf_all_bb.xval<=threshold_bb | prf_all_bb.ecc >= eclimit)';  % use tilde for nan
 
 selch = prfch & acuch;
 
%% load model data (no broadband correction)
opts = [];
opts.fileid         = modeldataID;
opts.average        = average;
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;
opts.issave         = false;
opts.compute        = false;
opts.outputDir      = 'pRFmodel';

opts.targetBAND     = 'FaLb';
modeldata_a0  = ecog_prf_constructTimeSeries(subjectList, opts);
modeldata_a0  = ecog_summarizeROIs(modeldata_a0);
modeldata_a0  = ecog_channelSelection(modeldata_a0,selectchs,selectch_exFEF,selectch_thresh);
model_all_a0  = ecog_rearrangePRF(ecog_realignTimeSeries(modeldata_a0,1));
 

%%
%%
close all
%%
%% %%%%%%%%%%%%%
%% Visualization
useChans = 'pRFchs';        % pRFchs, SELchs, ALLchs

switch useChans
    case 'pRFchs',      selch = prfch & acuch;
    case 'SELchs',      selch = acuch;
    case 'ALLchs',      selch = true(size(acuch));
end

ylin = [log10(0.6) 0];
yticklog = log10([0.6:0.1:1]);

ylinD = [log10(0.2) log10(3)];
yticklogD = log10([0.25 0.5 1 2.5 5]);

ylinW = [log10(0.08) log10(10)];
yticklogW = log10([0.1 0.3 1 3 10]);

   FntSiz = 20;
   plcol = get(groot,'defaultAxesColorOrder');

%% Bar plot (Alpha Suppression)
datats_a = -cat(2,model_all_a.datats{:})';

inPRF_bb  = mean(datats_a(selch&prfidxbb&~blankidx),'all','omitnan');
outPRF_bb = mean(datats_a(selch&~prfidxbb&~blankidx),'all','omitnan');
inPRF_a   = mean(datats_a(selch&prfidxa&~blankidx),'all','omitnan');
outPRF_a  = mean(datats_a(selch&~prfidxa&~blankidx),'all','omitnan');
inoutPRF  = mean(datats_a(selch&prfidxa&~prfidxbb&~blankidx),'all','omitnan');
BLANK     = mean(datats_a(blankidx,selch),'all','omitnan');

figure; bar([inPRF_bb,outPRF_bb,inPRF_a,outPRF_a,inoutPRF],'FaceColor',plcol(2,:));
hold on;
plot(xlim,[1 1].*BLANK,'k--');       % BLANK
hold off;
ylim(ylin);
yticks(yticklog);
yticklabels(arrayfun(@(x) sprintf('%g',10.^x),yticks,'UniformOutput',false));
set(gca,'FontSize',FntSiz);
xticklabels({'In_{bb}','Out_{bb}','In_a','Out_a','Border'});
% xticklabels({'Inside Broadband pRF','Outside Broadband pRF','Inside Alpja pRF','Outside Alpja pRF','Outside Broadband pRF & Inside Alpja pRF'});

%% Bar plot (Alpha Change)
datats_a = -cat(2,model_all_a0.datats{:})';

inPRF_bb  = mean(datats_a(selch&prfidxbb&~blankidx),'all','omitnan');
outPRF_bb = mean(datats_a(selch&~prfidxbb&~blankidx),'all','omitnan');
inPRF_a   = mean(datats_a(selch&prfidxa&~blankidx),'all','omitnan');
outPRF_a  = mean(datats_a(selch&~prfidxa&~blankidx),'all','omitnan');
inoutPRF  = mean(datats_a(selch&prfidxa&~prfidxbb&~blankidx),'all','omitnan');
BLANK     = mean(datats_a(blankidx,selch),'all','omitnan');

figure; bar([inPRF_bb,outPRF_bb,inPRF_a,outPRF_a,inoutPRF],'FaceColor',plcol(2,:));
hold on;
plot(xlim,[1 1].*BLANK,'k--');       % BLANK
hold off;
ylim(ylin);
yticks(yticklog);
yticklabels(arrayfun(@(x) sprintf('%g',10.^x),yticks,'UniformOutput',false));
set(gca,'FontSize',FntSiz);
xticklabels({'In_{bb}','Out_{bb}','In_a','Out_a','Border'});

%%
%%

%% Box plot (Alpha Suppression)
datats_a = -cat(2,model_all_a.datats{:})';
inPRF_bb = datats_a; outPRF_bb = datats_a; inPRF_a = datats_a; outPRF_a = datats_a; inoutPRF = datats_a; BLANK = datats_a;
inPRF_bb(~(selch&prfidxbb&~blankidx))           = nan;   inPRF_bb(:,~selch)  = [];
outPRF_bb(~(selch&~prfidxbb&~blankidx))         = nan;   outPRF_bb(:,~selch) = [];
inPRF_a(~(selch&prfidxa&~blankidx))             = nan;   inPRF_a(:,~selch)   = [];
outPRF_a(~(selch&~prfidxa&~blankidx))           = nan;   outPRF_a(:,~selch)  = [];
inoutPRF(~(selch&prfidxa&~prfidxbb&~blankidx))  = nan;   inoutPRF(:,~selch)  = [];
BLANK(~blankidx,:) = nan;                                BLANK(:,~selch)  = [];

figure; boxplot([inPRF_bb(:),outPRF_bb(:),inPRF_a(:),outPRF_a(:),inoutPRF(:),BLANK(:)],...
                'Labels',{'In_{bb}','Out_{bb}','In_a','Out_a','Border','Blank'});
            
hold on;
hl = plot(xlim,[1 1].*mean(BLANK,'all','omitnan'),'k--');       % BLANK
uistack(hl,'bottom');
hold off;
ylim(ylinW);
yticks(yticklogW);
yticklabels(arrayfun(@(x) sprintf('%g',10.^x),yticks,'UniformOutput',false));
set(gca,'FontSize',FntSiz);
set(get(gca,'XAxis'),'TickLabelInterpreter','tex');
set(gcf,'Position',get(gcf,'Position').*[1 1 1.05 1]);

%% Box plot (Alpha Change)
datats_a = -cat(2,model_all_a0.datats{:})';
inPRF_bb = datats_a; outPRF_bb = datats_a; inPRF_a = datats_a; outPRF_a = datats_a; inoutPRF = datats_a; BLANK = datats_a;
inPRF_bb(~(selch&prfidxbb&~blankidx))           = nan;   inPRF_bb(:,~selch)  = [];
outPRF_bb(~(selch&~prfidxbb&~blankidx))         = nan;   outPRF_bb(:,~selch) = [];
inPRF_a(~(selch&prfidxa&~blankidx))             = nan;   inPRF_a(:,~selch)   = [];
outPRF_a(~(selch&~prfidxa&~blankidx))           = nan;   outPRF_a(:,~selch)  = [];
inoutPRF(~(selch&prfidxa&~prfidxbb&~blankidx))  = nan;   inoutPRF(:,~selch)  = [];
BLANK(~blankidx,:) = nan;                                BLANK(:,~selch)  = [];

figure; boxplot([inPRF_bb(:),outPRF_bb(:),inPRF_a(:),outPRF_a(:),inoutPRF(:),BLANK(:)],...
                'Labels',{'In_{bb}','Out_{bb}','In_a','Out_a','Border','Blank'});
            
hold on;
hl = plot(xlim,[1 1].*mean(BLANK,'all','omitnan'),'k--');       % BLANK
uistack(hl,'bottom');
hold off;
ylim(ylinW);
yticks(yticklogW);
yticklabels(arrayfun(@(x) sprintf('%g',10.^x),yticks,'UniformOutput',false));
set(gca,'FontSize',FntSiz);
set(get(gca,'XAxis'),'TickLabelInterpreter','tex');
set(gcf,'Position',get(gcf,'Position').*[1 1 1.05 1]);

%%
%%

%% Histogram (Alpha Suppression)
datats_a = -cat(2,model_all_a.datats{:})';
inPRF_bb = datats_a; outPRF_bb = datats_a; inPRF_a = datats_a; outPRF_a = datats_a; inoutPRF = datats_a; BLANK = datats_a;
inPRF_bb(~(selch&prfidxbb&~blankidx))           = nan;   inPRF_bb(:,~selch)  = [];
outPRF_bb(~(selch&~prfidxbb&~blankidx))         = nan;   outPRF_bb(:,~selch) = [];
inPRF_a(~(selch&prfidxa&~blankidx))             = nan;   inPRF_a(:,~selch)   = [];
outPRF_a(~(selch&~prfidxa&~blankidx))           = nan;   outPRF_a(:,~selch)  = [];
inoutPRF(~(selch&prfidxa&~prfidxbb&~blankidx))  = nan;   inoutPRF(:,~selch)  = [];
BLANK(~blankidx,:) = nan;                                BLANK(:,~selch)  = [];

Labels = {'In_{bb}','Out_{bb}','In_a','Out_a','Border','Blank'};

figure; ht=tiledlayout(1,6,'Padding','compact','TileSpacing','none');
ii=1;
for idat = {inPRF_bb,outPRF_bb,inPRF_a,outPRF_a,inoutPRF,BLANK}
nexttile;
    H=histogram(idat{:},'BinWidth',0.02,'FaceColor',plcol(2,:),'EdgeColor','none');
    set(gca,'CameraUpVector',[1 0 0]);
    xlim(ylinD);
    xticks(yticklogD);
    xticklabels({});
    yticks([]);
    set(gca,'FontSize',FntSiz);
    
    %-- Median
    hold on;
    plot([1 1].*median(idat{:},'all','omitnan'),ylim,'k-.');
    plot([1 1].*median(BLANK,'all','omitnan'),ylim,'k--');
    hold off;
    
    %-- Title
    tl=title(Labels{ii},'VerticalAlignment','top');
    tl.Position(1) = max(xlim) - diff(xlim)*.01;
    ii=ii+1;
end
xticklabels(arrayfun(@(x) sprintf('%g',10.^x),xticks,'UniformOutput',false));
            
%% Histogram (Alpha Change)
datats_a = -cat(2,model_all_a0.datats{:})';
inPRF_bb = datats_a; outPRF_bb = datats_a; inPRF_a = datats_a; outPRF_a = datats_a; inoutPRF = datats_a; BLANK = datats_a;
inPRF_bb(~(selch&prfidxbb&~blankidx))           = nan;   inPRF_bb(:,~selch)  = [];
outPRF_bb(~(selch&~prfidxbb&~blankidx))         = nan;   outPRF_bb(:,~selch) = [];
inPRF_a(~(selch&prfidxa&~blankidx))             = nan;   inPRF_a(:,~selch)   = [];
outPRF_a(~(selch&~prfidxa&~blankidx))           = nan;   outPRF_a(:,~selch)  = [];
inoutPRF(~(selch&prfidxa&~prfidxbb&~blankidx))  = nan;   inoutPRF(:,~selch)  = [];
BLANK(~blankidx,:) = nan;                                BLANK(:,~selch)  = [];

Labels = {'In_{bb}','Out_{bb}','In_a','Out_a','Border','Blank'};

figure; ht=tiledlayout(1,6,'Padding','compact','TileSpacing','none');
ii=1;
for idat = {inPRF_bb,outPRF_bb,inPRF_a,outPRF_a,inoutPRF,BLANK}
nexttile;
    H=histogram(idat{:},'BinWidth',0.02,'FaceColor',plcol(2,:),'EdgeColor','none');
    set(gca,'CameraUpVector',[1 0 0]);
    xlim(ylinD);
    xticks(yticklogD);
    xticklabels({});
    yticks([]);
    set(gca,'FontSize',FntSiz);
    
    %-- Median
    hold on;
    plot([1 1].*median(idat{:},'all','omitnan'),ylim,'k-.');
    plot([1 1].*median(BLANK,'all','omitnan'),ylim,'k--');
    hold off;
    
    %-- Title
    tl=title(Labels{ii},'VerticalAlignment','top');
    tl.Position(1) = max(xlim) - diff(xlim)*.01;
    ii=ii+1;
end
xticklabels(arrayfun(@(x) sprintf('%g',10.^x),xticks,'UniformOutput',false));
            

%%
%%

%% Check data distribution (Alpha Suppression)
datats_a = -cat(2,model_all_a.datats{:})';
%-- all time points
figure;
boxplot(datats_a(:,selch));
ylim([-.45 1.05]);
%-- Blanks
figure;
boxplot(datats_a(blankidx,selch));
ylim([-.45 1.05]);

%-- Out pRF
datats_aS = datats_a;
datats_aS(~(selch&~prfidxa&~blankidx)) = nan;
figure; boxplot(datats_aS(:,selch));
ylim([-.45 1.05]);
%-- In pRF
datats_aS = datats_a;
datats_aS(~(selch&prfidxa&~blankidx)) = nan;
figure; boxplot(datats_aS(:,selch));
ylim([-.45 1.05]);

%% Check data distribution (Alpha Change)
datats_a = -cat(2,model_all_a0.datats{:})';
%-- all time points
figure;
boxplot(datats_a(:,selch));
ylim([-.45 1.05]);
%-- Blanks
figure;
boxplot(datats_a(blankidx,selch));
ylim([-.45 1.05]);

%-- Out pRF
datats_aS = datats_a;
datats_aS(~(selch&~prfidxa&~blankidx)) = nan;
figure; boxplot(datats_aS(:,selch));
ylim([-.45 1.05]);
%-- In pRF
datats_aS = datats_a;
datats_aS(~(selch&prfidxa&~blankidx)) = nan;
figure; boxplot(datats_aS(:,selch));
ylim([-.45 1.05]);


%%
%%
close all
%%
%% %%%%%%%%%%%%%
%% Visualization in ROIs
useChans = 'pRFchs';        % pRFchs, SELchs, ALLchs

switch useChans
    case 'pRFchs',      selch = prfch & acuch;
    case 'SELchs',      selch = acuch;
    case 'ALLchs',      selch = true(size(acuch));
end

ylin = [log10(0.45) log10(1.8)];
yticklog = log10([0.5 0.7 1 1.4 2]);

   FntSiz = 20;
   plcol = get(groot,'defaultAxesColorOrder');

rois = {{'V1'},{'V2'},{'V3'},{'V3a','V3b','LO1','LO2','TO','IPS'}};
roilabels = {'V1','V2','V3','Dorsolateral'};

%% Bar plot (Alpha Suppression)
datats_a = -cat(2,model_all_a.datats{:})';
meanfun = @(d) mean(d,'all','omitnan');

figure; tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[.5 .5 1.8 1.8]);
for iroi = 1:length(rois)
    roich = ismember(channels.wangarea,rois{iroi})';

inPRF_bb  = meanfun(datats_a(roich&selch&prfidxbb&~blankidx));
outPRF_bb = meanfun(datats_a(roich&selch&~prfidxbb&~blankidx));
inPRF_a   = meanfun(datats_a(roich&selch&prfidxa&~blankidx));
outPRF_a  = meanfun(datats_a(roich&selch&~prfidxa&~blankidx));
inoutPRF  = meanfun(datats_a(roich&selch&prfidxa&~prfidxbb&~blankidx));
BLANK     = meanfun(datats_a(blankidx,roich&selch));

nexttile;
bar([inPRF_bb,outPRF_bb,inPRF_a,outPRF_a,inoutPRF],'FaceColor',plcol(2,:));
hold on;
plot(xlim,[1 1].*BLANK,'k--');       % BLANK
hold off;
ylim(ylin);
yticks(yticklog);
yticklabels(arrayfun(@(x) sprintf('%g',10.^x),yticks,'UniformOutput',false));
set(gca,'FontSize',FntSiz);
xticklabels({'In_{bb}','Out_{bb}','In_a','Out_a','Border'});
% xticklabels({'Inside Broadband pRF','Outside Broadband pRF','Inside Alpja pRF','Outside Alpja pRF','Outside Broadband pRF & Inside Alpja pRF'});
title(roilabels{iroi});

end

%% Bar plot (Alpha change)
datats_a = -cat(2,model_all_a0.datats{:})';
meanfun = @(d) mean(d,'all','omitnan');

figure; tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[.5 .5 1.8 1.8]);
for iroi = 1:length(rois)
    roich = ismember(channels.wangarea,rois{iroi})';

inPRF_bb  = meanfun(datats_a(roich&selch&prfidxbb&~blankidx));
outPRF_bb = meanfun(datats_a(roich&selch&~prfidxbb&~blankidx));
inPRF_a   = meanfun(datats_a(roich&selch&prfidxa&~blankidx));
outPRF_a  = meanfun(datats_a(roich&selch&~prfidxa&~blankidx));
inoutPRF  = meanfun(datats_a(roich&selch&prfidxa&~prfidxbb&~blankidx));
BLANK     = meanfun(datats_a(blankidx,roich&selch));

nexttile;
bar([inPRF_bb,outPRF_bb,inPRF_a,outPRF_a,inoutPRF],'FaceColor',plcol(2,:));
hold on;
plot(xlim,[1 1].*BLANK,'k--');       % BLANK
hold off;
ylim(ylin);
yticks(yticklog);
yticklabels(arrayfun(@(x) sprintf('%g',10.^x),yticks,'UniformOutput',false));
set(gca,'FontSize',FntSiz);
xticklabels({'In_{bb}','Out_{bb}','In_a','Out_a','Border'});
% xticklabels({'Inside Broadband pRF','Outside Broadband pRF','Inside Alpja pRF','Outside Alpja pRF','Outside Broadband pRF & Inside Alpja pRF'});
title(roilabels{iroi});

end


%% Bar plot (alpha suppression VS alpha power change)
datats = -cat(2,model_all_a.datats{:})';
datats2 = -cat(2,model_all_a0.datats{:})';
meanfun = @(d) mean(d,'all','omitnan');

figure; tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[.5 .5 1.8 1.8]);
for iroi = 1:length(rois)
    roich = ismember(channels.wangarea,rois{iroi})';

inPRF_bb  = meanfun(datats(roich&selch&prfidxbb&~blankidx));
outPRF_bb = meanfun(datats(roich&selch&~prfidxbb&~blankidx));
inPRF_a   = meanfun(datats(roich&selch&prfidxa&~blankidx));
outPRF_a  = meanfun(datats(roich&selch&~prfidxa&~blankidx));
inoutPRF  = meanfun(datats(roich&selch&prfidxa&~prfidxbb&~blankidx));
BLANK     = meanfun(datats(blankidx,roich&selch));

inPRF_bb2  = meanfun(datats2(roich&selch&prfidxbb&~blankidx));
outPRF_bb2 = meanfun(datats2(roich&selch&~prfidxbb&~blankidx));
inPRF_a2   = meanfun(datats2(roich&selch&prfidxa&~blankidx));
outPRF_a2  = meanfun(datats2(roich&selch&~prfidxa&~blankidx));
inoutPRF2  = meanfun(datats2(roich&selch&prfidxa&~prfidxbb&~blankidx));
BLANK2     = meanfun(datats2(blankidx,roich&selch));

nexttile;
% bar([inPRF_bb,outPRF_bb,inPRF_a,outPRF_a,inoutPRF],'FaceColor',plcol(2,:));
B=bar([inPRF_bb,inoutPRF,outPRF_a;inPRF_bb2,inoutPRF2,outPRF_a2]');
B(1).FaceColor = plcol(2,:);
B(2).FaceColor = plcol(2,:).*0.6;
hold on;
% plot(xlim,[1 1].*BLANK,'k--');        % BLANK
hold off;
ylim(ylin);
yticks(yticklog);
yticklabels(arrayfun(@(x) sprintf('%g',10.^x),yticks,'UniformOutput',false));
set(gca,'FontSize',FntSiz);
% xticklabels({'In_{bb}','Out_{bb}','In_a','Out_a','Border'});
xticklabels({'In_{bb}','Border','Out_a'});
% xticklabels({'Inside Broadband pRF','Outside Broadband pRF','Inside Alpja pRF','Outside Alpja pRF','Outside Broadband pRF & Inside Alpja pRF'});
title(roilabels{iroi});

end
legend({'Model-based Alpha Suppression','Power Change in Alpha'});
