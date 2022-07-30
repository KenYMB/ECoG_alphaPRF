% compute connectivities across electrodes for subjects with HDgrid
%   investigate how to analyze across all electrodes

% 20210601 Yuasa - compute time series (test)
%                  for alpha

%%
close all; % clearvars;
% if isempty(gcp('nocreate')),  parpool([1 40]); end
% startupToolboxToolbox;

%% Define paths and dataset
checkPath;
%-- Input & Output path
plotsavePth    = 'ConnectivityTS';
xspctrmPth     = 'xSpectrum';
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'Distance');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'hasHDgrid','yes');
% subjectList = subjectList(1);
% subjectList = subjectList(2);

tarBAND     = 'alpha';
disttype    = 'norm';            % 'square','diamond','norm'
% arounddist  = [1 2 3 6];
    arounddist  = 1:8;
%%
try
%-- Dataset specs
if ~iscell(subjectList),       subjectList = {subjectList};
elseif iscell(subjectList{1}), subjectList = subjectList{:};
end

%% Cross-spectra
%% Load Coherence
opts = [];
opts.compute    = false;
opts.issave     = false;
opts.average    = 'runs';
opts.isavgfreq  = false;
opts.method     = 'mscoh';
opts.targetBAND     = tarBAND;
opts.stimNames      = 'HORIZONTAL*';
[coh_a1]    = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'VERTICAL*';
[coh_a2]    = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'DIAGONAL*';
[coh_a3]    = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'BLANK';
[coh_a0]    = ecog_prf_connectivity(subjectList, opts);

%% load time series data
clear alphaType broadbandType
decN = 3;
% decN = 1;

average        ='runs';
prfmodel       = 'linear';
gaussianmode   = 'gs';
smoothingMode  ='decimate';
smoothingN     = decN;
selectchs      = 'wangprobchs';     % only use wangprobchs as seed but use all grid channels for averaged coherence
    allowlag       = false;
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;
    
ecog_APRFF_INITa_loaddata;

%% load iaf
opts = [];
opts.average        = 'runs';
opts.compute         = false;
spcrm_params = ecog_prf_fitalpha(subjectList, opts);

%% rearrange events 
coh_aAll = cell(size(coh_a0));
avgdat   = cell(size(coh_a0));
avgbsl   = cell(size(coh_a0));
iafidx   = cell(size(coh_a0));
for isbj = 1:length(subjectList)
coh_aAll{isbj}  = coh_a0{isbj};
coh_aAll{isbj}.connectivity = cat(3,coh_aAll{isbj}.connectivity,coh_a1{isbj}.connectivity,coh_a2{isbj}.connectivity,coh_a3{isbj}.connectivity);
coh_aAll{isbj}.events       = cat(1,coh_aAll{isbj}.events,coh_a1{isbj}.events,coh_a2{isbj}.events,coh_a3{isbj}.events);

eventslist = sortrows([coh_aAll{isbj}.events.stim_file_index, [1:height(coh_aAll{isbj}.events)]']);

coh_aAll{isbj}.connectivity = coh_aAll{isbj}.connectivity(:,:,eventslist(:,2),:);
coh_aAll{isbj}.events       = coh_aAll{isbj}.events(eventslist(:,2),:);

%% compute coherence for around channels
whichHDgrid = 'GB';

iconn     = coh_aAll{isbj};
connbase  = coh_a0{isbj};
pltdat  = iconn.connectivity;
pltbase = mean(connbase.connectivity,3,'omitnan');

datdim = size(pltdat);
datdim(2) = length(arounddist);

   avgdat{isbj} = nan(datdim);
   datdim(3) = 1;
   avgbsl{isbj} = nan(datdim);
   dd = 0;
   for arddist = arounddist
       dd = dd + 1;
       for ich=1:height(iconn.channels)
           nRows = 16; nCols = 8;
           tarchnum  = str2double(strtok(iconn.channels.name{ich},whichHDgrid));
           tarchx = mod(tarchnum-1,nRows)+1;
           tarchy = ceil(tarchnum./nRows);
           [dGridx,dGridy]= meshgrid((1:nRows)-tarchx,(nCols:-1:1)-tarchy);
           switch disttype
               case 'square'
           dGrid = (abs(dGridx)==arddist & abs(dGridy)<=arddist) |...
                    (abs(dGridy)==arddist & abs(dGridx)<=arddist);
               case 'diamond'
           dGrid = abs(dGridx)+abs(dGridy);
           dGrid = dGrid == arddist;
               case 'norm'
           dGrid = sqrt(dGridx.^2+dGridy.^2);
           dGrid = (dGrid >= arddist) & (dGrid < (arddist+1));
           end
           avgchnum = find(flipud(dGrid)');
           avgchname = unique(...
                         [arrayfun(@(E) sprintf('%s%d',whichHDgrid,E),avgchnum,'UniformOutput',false),...
                          arrayfun(@(E) sprintf('%s%02d',whichHDgrid,E),avgchnum,'UniformOutput',false),...
                          arrayfun(@(E) sprintf('%s%03d',whichHDgrid,E),avgchnum,'UniformOutput',false)],...
                       'stable');
           avgchidx  = ismember(iconn.channels.name,avgchname);
           avgdat{isbj}(ich,dd,:,:) = mean(pltdat(avgchidx,ich,:,:),1,'omitnan');
           avgbsl{isbj}(ich,dd,:,:) = mean(pltbase(avgchidx,ich,:,:),1,'omitnan');
       end
   end

%% align channels
nchan = height(modeldata_a{isbj}.channels.name);
[~,ichconn,ichdat]=intersect(iconn.channels.name,modeldata_a{isbj}.channels.name,'stable');

iconn.connectivity = iconn.connectivity(ichconn,ichconn,:,:);
modeldata_a{isbj}.datats  = cellfun(@(x) x(ichdat,:),modeldata_a{isbj}.datats,'UniformOutput',false);
modeldata_bb{isbj}.datats = cellfun(@(x) x(ichdat,:),modeldata_bb{isbj}.datats,'UniformOutput',false);
subfld = fieldnames(prf_params_a{isbj});
for ifld = reshape(subfld,1,[])
  if ~ismember(ifld,'channels')
    if size(prf_params_a{isbj}.(ifld{:}),1) == nchan
        prf_params_a{isbj}.(ifld{:})  = prf_params_a{isbj}.(ifld{:})(ichdat,:,:);
        prf_params_bb{isbj}.(ifld{:}) = prf_params_bb{isbj}.(ifld{:})(ichdat,:,:);
    elseif size(prf_params_a{isbj}.(ifld{:}),3) == nchan
        prf_params_a{isbj}.(ifld{:})  = prf_params_a{isbj}.(ifld{:})(:,:,ichdat,:);
        prf_params_bb{isbj}.(ifld{:}) = prf_params_bb{isbj}.(ifld{:})(:,:,ichdat,:);
    end
  end
end

iconn.channels  = iconn.channels(ichconn,:);
avgdat{isbj}          = avgdat{isbj}(ichconn,:,:,:);
avgbsl{isbj}          = avgbsl{isbj}(ichconn,:,:,:);
modeldata_a{isbj}.channels   = modeldata_a{isbj}.channels(ichdat,:);
modeldata_bb{isbj}.channels  = modeldata_bb{isbj}.channels(ichdat,:);
prf_params_a{isbj}.channels  = prf_params_a{isbj}.channels(ichdat,:);
prf_params_bb{isbj}.channels = prf_params_bb{isbj}.channels(ichdat,:);


[~,ichconn,ichparams]=intersect(iconn.channels.name,spcrm_params{isbj}.channels.name,'stable');
assert(isequal(ichconn',1:numel(ichconn)),'Something wrong');
spcrm_params{isbj}.channels     = spcrm_params{isbj}.channels(ichparams,:);
spcrm_params{isbj}.resamp_parms = spcrm_params{isbj}.resamp_parms(ichparams,:);

%% iaf
iaf = cellfun(@(x) mean(10.^x(:,10)),spcrm_params{isbj}.resamp_parms);
iafidx{isbj} = arrayfun(@(x) find(abs(iconn.f - x) == min(abs(iconn.f - x))), iaf, 'UniformOutput', true);

end

%% update channel information
% % iconn = ecog_updatePRFdata(iconn);
% iconn   = ecog_summarizeROIs(iconn);
% [~,iconn]   = ecog_rearrangePRF(iconn,'wang','norm',.05);

ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITd_threshold;

avgdat = cat(1,avgdat{:});
avgbsl = cat(1,avgbsl{:});
iafidx = cat(1,iafidx{:});

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
   
   evlng = height(iconn.events);
   flng  = length(iconn.f);
   chlng = height(prf_all_a.channels);
   ndist = length(arounddist);
   
    decRate = round(evlng/decN)./evlng;
    boundaries = ([0 28 40 56 84 96 112 140 152 168 196 208 224]+0.5).*decRate;
    blankbnd   = boundaries([3,4;6,7;9,10;12,13]);
    medblank   = mean(blankbnd,2);

catch ex
    rethrow(ex);
end

%-- pick iaf & apply decimate
coha = zeros(size(avgdat,[3 1 2]));
bsla = zeros(size(avgbsl,[3 1 2]));
for ich = 1:chlng
    coha(:,ich,:)  = permute(avgdat(ich,:,:,iafidx(ich)),[3 1 2 4]);
    bsla(:,ich,:)  = permute(avgbsl(ich,:,:,iafidx(ich)),[3 1 2 4]);
end
coha  = decimate_col(coha,decN,[],1);

evlng2 = size(coha,1);

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
    
    %-- set threshold for pRF inside  (1sd~60.7%, 2sd~13.5%, 2.5sd~4.4%, 3sd~1.1%)
    prfthreshbb = gainbb*0.05;
    prfthresha  = -gaina*0.05;
    prfidxbb = false(size(modelbb));  prfidxa = false(size(modela));
    for ich=1:chlng
        prfidxbb(:,ich) = (modelbb(:,ich) > prfthreshbb(ich)) & ~blankidx;
        prfidxa(:,ich)  = (modela(:,ich) < prfthresha(ich)) & ~blankidx;
    end
    

 %%
 nincl = 6;
 prfch = sum(prfidxbb,1) > nincl & sum(prfidxa,1) > nincl;
 acuch =  ~(prf_all_bb.xval<=threshold_bb | prf_all_a.xval<=threshold_a ...
            | prf_all_bb.ecc >= eclimit | prf_all_a.ecc >= eclimit)';  % use tilde for nan
 
 selch = prfch & acuch;
 
%%
%%
%%
%% Chance level
% filepath = fullfile(datPth, 'xSpectrum',sprintf('coh-shuffle-%s-%s-%s',subject,tarBAND,disttype));
% 
% if exist([filepath '.mat'],'file')
%     load(filepath);
% else
%     error('File does not exist');
% end
% 
% chcavgtsAdat = mean(chcavgtsAdat,'all','omitnan');
% chcavgtsAbsl = mean(chcavgtsAbsl,'all','omitnan');

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
switch tarBAND
    case 'alpha'
        prfidx = prfidxa;  inPRFmode = 'A';
    case 'broadband'
        prfidx = prfidxbb; inPRFmode = 'B';
end

ylindiff  = [-1 1].*0.081;
ylinraw   = [0.25 0.7];

%% Plot coherence for distance

%%% Segregate in inPRF & outPRF
cohdatAbsl = mat2cell(coha,size(coha,1),ones(size(coha,2),1),size(coha,3));
cohdatAprf = cohdatAbsl;
cohdatAout = cohdatAbsl;
for ich=1:chlng
    cohdatAbsl{ich} = cohdatAbsl{ich}(blankidx,:);
    cohdatAprf{ich} = cohdatAprf{ich}(prfidx(:,ich),:);
    cohdatAout{ich} = cohdatAout{ich}(~blankidx&~prfidx(:,ich),:);
end

%%% Take average in each elec
avgtsAbsl = cellfun(@mean,cohdatAbsl(selch),repmat({1},1,sum(selch)),'UniformOutput',false)';
avgtsAprf = cellfun(@mean,cohdatAprf(selch),repmat({1},1,sum(selch)),'UniformOutput',false)';
avgtsAout = cellfun(@mean,cohdatAout(selch),repmat({1},1,sum(selch)),'UniformOutput',false)';
    avgtsAbsl = cat(1,avgtsAbsl{:});
    avgtsAprf = cat(1,avgtsAprf{:});
    avgtsAout = cat(1,avgtsAout{:});
    
    
%%% Plot Difference
figure('MenuBar','none');
ys = avgtsAprf-avgtsAbsl;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl = plot(arounddist,y1,'Color',plcol(1,:),'LineWidth',1.5);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(1,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl,'top');
ys = avgtsAout-avgtsAbsl;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(2) = plot(arounddist,y1,'Color',plcol(2,:),'LineWidth',1.5);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(2,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(2),'top');
%     hb = plot(xlim,[1,1].*(chcavgtsAdat-chcavgtsAbsl),'k--','LineWidth',1.2);
%     uistack(hb,'bottom');
set(gca,'FontSize',FntSiz);
legend(hl,{'in-pRF - BLANK','out-pRF - BLANK'},'Location','southeast');
xlabel('Distance'); ylabel('Coherence');
xlim(minmax(arounddist));
ylim(ylindiff);

   figname = sprintf('CoherenceTS_%s-%s-%s-%s_Distance%s-diff',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   savefigauto(gcf,fullfile(plotsavedir,figname),{'png','eps'})

figure('MenuBar','none');
ys = avgtsAprf-avgtsAout;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl = plot(arounddist,y1,'Color',plcol(3,:),'LineWidth',1.5);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(3,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl,'top');
    hb = plot(xlim,[0,0],'k--','LineWidth',1.2);
    uistack(hb,'bottom');
set(gca,'FontSize',FntSiz);
legend(hl,{'in-pRF - out-pRF'},'Location','southeast');
xlabel('Distance'); ylabel('Coherence');
xlim(minmax(arounddist));
ylim(ylindiff);

   figname = sprintf('CoherenceTS_%s-%s-%s-%s_Distance%s-In-Out',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   savefigauto(gcf,fullfile(plotsavedir,figname),{'png','eps'})
   
   
%%% Plot Raw
figure('MenuBar','none');
ys = avgtsAbsl;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl = plot(arounddist,y1,'Color',plcol(4,:),'LineWidth',1.5);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(4,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl,'top');
ys = avgtsAprf;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(2) = plot(arounddist,y1,'Color',plcol(5,:),'LineWidth',1.5);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(5,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(2),'top');
ys = avgtsAout;
%     y1 = median(ys,1,'omitnan');
    y1 = mean(ys,1,'omitnan');
    y2 = prctile(ys,(alpha/2)*100,1);
    y3 = prctile(ys,(1-alpha/2)*100,1);
  hl(3) = plot(arounddist,y1,'Color',plcol(6,:),'LineWidth',1.5);
  hold on;
  fill([arounddist, fliplr(arounddist)], [y2, fliplr(y3)],plcol(6,:),'EdgeColor','none','FaceAlpha',0.3);
  uistack(hl(3),'top');
%     hb = plot(xlim,[1,1].*chcavgtsAdat,'k--','LineWidth',1.2);
%     uistack(hb,'bottom');
set(gca,'FontSize',FntSiz);
legend(hl,{'BLANK','in-pRF','out-pRF'},'Location','northeast');
xlabel('Distance'); ylabel('Coherence');
xlim(minmax(arounddist));
ylim(ylinraw);

   figname = sprintf('CoherenceTS_%s-%s-%s-%s_Distance%s',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   savefigauto(gcf,fullfile(plotsavedir,figname),{'png','eps'})

%% anova
valch  = ~any(isnan([avgtsAbsl avgtsAprf avgtsAout]),2);
avgtsA = [avgtsAbsl(valch,:);avgtsAprf(valch,:);avgtsAout(valch,:)];
ps = anova2(avgtsA,sum(valch));


%%
%%
%%
%%
%% %%%%%%%%%%%%%
%% Visualization for distance #1
irnd = 1;       % index of arounddist

ylindiff  = [-1 1].*0.12;
ylinraw   = [0.2 0.8];
%% plot correlation

hglblank = true;        % use gray color for BLANK
hglinprf = false;       % use bold line for in-pRF
showelps = true;

dcoha = coha(:,:,irnd) - bsla(:,:,irnd);

figure('Menubar','none');
yyaxis left
scatter(dcoha(:,selch),dattsbb(:,selch),'LineWidth',1);
hold on
yy1 = ylim;
 if hglblank,  scatter(dcoha(blankidx,selch),dattsbb(blankidx,selch),[],mkblcl(plcol(1,:)),'LineWidth',1);  end
 if hglinprf,  scatter(dcoha(prfidxbb&selch),dattsbb(prfidxbb&selch),[],plcol(1,:),'LineWidth',2.0);  end
yyaxis right
scatter(dcoha(:,selch),dattsa(:,selch),'LineWidth',1);
hold on
yy2 = ylim;
 if hglblank,  scatter(dcoha(blankidx,selch),dattsa(blankidx,selch),[],mkblcl(plcol(2,:)),'LineWidth',1);  end
 if hglinprf,  scatter(dcoha(prfidxa&selch),dattsa(prfidxa&selch),[],plcol(2,:),'LineWidth',2.0);  end
hA = get(gca,'YAxis');
    hA(1).Limits = [-1 1] .* max(abs(yy1(2)));
    hA(2).Limits = [-1 1] .* max(abs(yy2(1)));
    
%-- zero line
yyaxis left
    xwidth = xlim;
    h0 = plot(xwidth,[0 0],'k:','LineWidth',1.8);
    xlim(xwidth);
    hA(1).Parent.Children = cat(1,hA(1).Parent.Children(2:end),h0);
    
%-- blank line
    h0 = plot([0 0],ylim,'k:','LineWidth',1.8);
    hA(1).Parent.Children = cat(1,hA(1).Parent.Children(2:end),h0);
    
%-- show ellipse
if showelps
  yyaxis left
  h=error_ellipse(cov(reshape(dcoha(:,selch),[],1),reshape(dattsbb(:,selch),[],1)),...
                [median(reshape(dcoha(:,selch),[],1)),median(reshape(dattsbb(:,selch),[],1))],'conf',1-alpha);
  h.LineWidth = 2.6;  h.LineStyle = '-.'; h.Marker='none';
  h.Color = [0 0 1];
  yyaxis right
  h=error_ellipse(cov(reshape(dcoha(:,selch),[],1),reshape(dattsa(:,selch),[],1)),...
                [median(reshape(dcoha(:,selch),[],1)),median(reshape(dattsa(:,selch),[],1))],'conf',1-alpha);
  h.LineWidth = 2.6;  h.LineStyle = '-.'; h.Marker='none';
  h.Color = [1 0 0];
end
    
%-- parameters
set(gca,'FontSize',FntSiz);
hA(1).TickLabel = arrayfun(@(x) sprintf('%g',x./100+1),hA(1).TickValues,'UniformOutput',false);
if max(abs(hA(2).TickValues)) > 1
    hA(2).TickValues = log10([0.03 0.1 0.3 1 3 10 30]);
else
    hA(2).TickValues = log10([1/8 1/4 1/2 1 2 4 8]);
end
hA(2).TickLabel = arrayfun(@(x) sprintf('%g',10.^x),hA(2).TickValues,'UniformOutput',false);

yyaxis left
ylabel('Broadband Power Ratio');
yyaxis right
A=ylabel('Alpha Power Ratio','Rotation',270,'VerticalAlignment','bottom');
xlabel(sprintf('Difference of Coherence from BLANK')); %  - Peak Alpha Frequency
% title(sprintf('%s',useChans));


   figname = sprintf('CoherenceTS_%s-%s-%s-%s_2Dhist',subject,useChans,tarBAND,disttype);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   savefigauto(gcf,fullfile(plotsavedir,figname),{'png','eps'})
   
  
  %% separate time points into in- & out-pRFs 
switch tarBAND
    case 'alpha'
        datts = dattsa;
    case 'broadband'
        datts = dattsbb;
end
  
figure('Menubar','none');
h1=scatter(dcoha(~prfidx&~blankidx&selch),datts(~prfidx&~blankidx&selch),[],plcol(5,:),'LineWidth',1.0);
hold on
h2=scatter(dcoha(prfidx&selch),datts(prfidx&selch),[],plcol(7,:),'LineWidth',1.0);
yy1=ylim;

hA = get(gca,'YAxis');
%-- zero line
    xwidth = xlim;
    h0 = plot(xwidth,[0 0],'k:','LineWidth',1.8);
    xlim(xwidth);
    hA(1).Parent.Children = cat(1,hA(1).Parent.Children(2:end),h0);
    
%-- blank line
    h0 = plot([0 0],yy1,'k:','LineWidth',1.8);
    hA(1).Parent.Children = cat(1,hA(1).Parent.Children(2:end),h0);
    ylim(yy1);
  
%-- show ellipse
if showelps
  h=error_ellipse(cov(reshape(dcoha(~prfidx&~blankidx&selch),[],1),reshape(datts(~prfidx&~blankidx&selch),[],1)),...
                [median(reshape(dcoha(~prfidx&~blankidx&selch),[],1)),median(reshape(datts(~prfidx&~blankidx&selch),[],1))],'conf',1-alpha);
  h.LineWidth = 2.6;  h.LineStyle = '--'; h.Marker='none';
  h.Color = [0 1 0];
  
  h=error_ellipse(cov(reshape(dcoha(prfidx&selch),[],1),reshape(datts(prfidx&selch),[],1)),...
                [median(reshape(dcoha(prfidx&selch),[],1)),median(reshape(datts(prfidx&selch),[],1))],'conf',1-alpha);
  h.LineWidth = 2.6;  h.LineStyle = '-.'; h.Marker='none';
  h.Color = [1 0 0];
  
end
    
%-- parameters
set(gca,'FontSize',FntSiz);
switch tarBAND
    case 'alpha'
        if max(abs(hA(1).TickValues)) > 1
            hA(1).TickValues = log10([0.03 0.1 0.3 1 3 10 30]);
        else
            hA(1).TickValues = log10([1/8 1/4 1/2 1 2 4 8]);
        end
        hA(1).TickLabel = arrayfun(@(x) sprintf('%g',10.^x),hA(1).TickValues,'UniformOutput',false);
    case 'broadband'
        hA(1).TickLabel = arrayfun(@(x) sprintf('%g',x./100+1),hA(1).TickValues,'UniformOutput',false);
end

legend([h2,h1],{'in pRF','out pRF'},'Location','best');

A=ylabel('Alpha Power Ratio');
xlabel(sprintf('Difference of Coherence from BLANK')); %  - Peak Alpha Frequency
% title(sprintf('%s',useChans));


   figname = sprintf('CoherenceTS_%s-%s-%s-%s_2Dhist%s_In-Out',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   savefigauto(gcf,fullfile(plotsavedir,figname),{'png','eps'})
  
   
%%
%% based on pRF model

%%% box plot of single chs
cohdatbsl = mat2cell(coha(:,:,irnd),size(coha,1),ones(size(coha,2),1));
cohdatprf = cohdatbsl;
cohdatout = cohdatbsl;
for ich=1:chlng
    cohdatbsl{ich} = cohdatbsl{ich}(blankidx);
    cohdatprf{ich} = cohdatprf{ich}(prfidx(:,ich));
    cohdatout{ich} = cohdatout{ich}(~blankidx&~prfidx(:,ich));
end

close all
%% %%%%%%%%%%%%%%%%%%%
%% Distribution of coherence

%%% box plot of all elec
avgtsbsl = cellfun(@mean,cohdatbsl(selch))';
avgtsprf = cellfun(@mean,cohdatprf(selch))';
avgtsout = cellfun(@mean,cohdatout(selch))';

figure('Menubar','none');
hb=boxplot([avgtsbsl;avgtsprf;avgtsout],...
        [repmat({'BLANK'},numel(avgtsbsl),1);...
         repmat({'in pRF'},numel(avgtsprf),1);...
         repmat({'out pRF'},numel(avgtsout),1)]);
ylim(ylinraw);
set(hb,{'linew'},{1.6});
set(gca,'FontSize',FntSiz);
title(sprintf('%s','pRF channels'));
ylabel(sprintf('Coherence - Peak Alpha Frequency'));
   
   figname = sprintf('CoherenceTS_%s-%s-%s-%s_BoxDist%s',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   savefigauto(gcf,fullfile(plotsavedir,figname),{'png','eps'})
   
%%

%%% difference box plot of all elec

figure('Menubar','none');
set(gcf,'Position',get(gcf,'Position').*[1 1 1.2 1]);
hb=boxplot([avgtsprf-avgtsbsl;avgtsout-avgtsbsl;avgtsprf-avgtsout],...
        [repmat({'in pRF - BLANK'},numel(avgtsprf),1);...
         repmat({'out pRF - BLANK'},numel(avgtsout),1);...
         repmat({'in pRF - out pRF'},numel(avgtsout),1)]);
ylim(ylindiff);
hold on; hl = plot(xlim,[0 0],'k:','LineWidth',1.2);  uistack(hl,'bottom');
set(hb,{'linew'},{1.6});
set(gca,'FontSize',FntSiz);
title(sprintf('%s','pRF channels'));
ylabel(sprintf('Coherence - Peak Alpha Frequency'));
   
   figname = sprintf('CoherenceTS_%s-%s-%s-%s_BoxDist%s-diff',subject,useChans,tarBAND,disttype,inPRFmode);
   if decN>1, figname = sprintf('%s-decimate%d',figname,decN); end
   set(gcf,'Name',figname);
   savefigauto(gcf,fullfile(plotsavedir,figname),{'png','eps'})

%%
%%
%%
%%
% end
