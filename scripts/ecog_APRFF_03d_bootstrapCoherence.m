% compute connectivities across electrodes for subjects with HDgrid
%   investigate how to analyze based on a patient across all electrodes

% 20220207 Yuasa - bootstrapping for coherence
% 20220317 Yuasa - save all trials condition

%% %%%%%%%%%%%%%%%%%%%%
%% test
%% %%%%%%%%%%%%%%%%%%%%
%%
 close all; clear all;
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

nboot = 5000;
    
for selsbj = 1:(length(HDsubjectList)+1)
%-- Dataset specs
if selsbj > length(HDsubjectList),  subjectList = HDsubjectList;
else,                               subjectList = HDsubjectList(selsbj);
end
%% load time series data
clear alphaType broadbandType
decN = 3;
% decN = 1;

average        ='runs';
prfmodel       = 'linear';
gaussianmode   = 'gs';
smoothingMode  ='decimate';
smoothingN     = decN;
% selectchs      = 'wangprobchs';     % only use wangprobchs as seed but use all grid channels for averaged coherence
selectchs      = 'GB*';             % use all grid channels
    allowlag       = false;
    allowbeta      = true;
    allowwide      = true;
    allowmixbeta   = true;
    
ecog_APRFF_INITa_loaddata;

%%% Subject Name
nsbj = length(subjectList);
if nsbj==1
   subject = subjectList{1};
else
   subject = 'all';
end

%%% Bootstrapping file
filename = sprintf('cohboot-%s-%s-%s.mat',subject,useChans,disttype);
filepath = fullfile(SetDefaultAnalysisPath('DAT',outputDir),filename);

%%% Coherence
try
%% Cross-spectra
%% Load Coherence for alpha
tarBAND     = 'alpha';
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

%% Load Coherence for broadband
tarBAND     = 'broadband';
opts = [];
opts.compute    = false;
opts.issave     = false;
opts.average    = 'runs';
opts.isavgfreq  = true;
opts.method     = 'mscoh';
opts.targetBAND     = tarBAND;
opts.stimNames      = 'HORIZONTAL*';
[coh_bb1]    = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'VERTICAL*';
[coh_bb2]    = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'DIAGONAL*';
[coh_bb3]    = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'BLANK';
[coh_bb0]    = ecog_prf_connectivity(subjectList, opts);

%% load iaf
opts = [];
opts.average        = 'runs';
opts.compute         = false;
spcrm_params = ecog_prf_fitalpha(subjectList, opts);

%% rearrange events 
coh_aAll  = cell(size(coh_a0));
coh_bbAll = cell(size(coh_bb0));
avgdat_a  = coh_aAll;
avgdat_bb = coh_bbAll;
iafidx   = {size(coh_a0)};
for isbj = 1:length(subjectList)
%-- alpha
coh_aAll{isbj}  = coh_a0{isbj};
coh_aAll{isbj}.connectivity = cat(3,coh_aAll{isbj}.connectivity,coh_a1{isbj}.connectivity,coh_a2{isbj}.connectivity,coh_a3{isbj}.connectivity);
coh_aAll{isbj}.events       = cat(1,coh_aAll{isbj}.events,coh_a1{isbj}.events,coh_a2{isbj}.events,coh_a3{isbj}.events);

eventslist = sortrows([coh_aAll{isbj}.events.stim_file_index, [1:height(coh_aAll{isbj}.events)]']);

coh_aAll{isbj}.connectivity = coh_aAll{isbj}.connectivity(:,:,eventslist(:,2),:);
coh_aAll{isbj}.events       = coh_aAll{isbj}.events(eventslist(:,2),:);

%-- broadband
coh_bbAll{isbj}  = coh_bb0{isbj};
coh_bbAll{isbj}.connectivity = cat(3,coh_bbAll{isbj}.connectivity,coh_bb1{isbj}.connectivity,coh_bb2{isbj}.connectivity,coh_bb3{isbj}.connectivity);
coh_bbAll{isbj}.events       = cat(1,coh_bbAll{isbj}.events,coh_bb1{isbj}.events,coh_bb2{isbj}.events,coh_bb3{isbj}.events);

eventslist = sortrows([coh_bbAll{isbj}.events.stim_file_index, [1:height(coh_bbAll{isbj}.events)]']);

coh_bbAll{isbj}.connectivity = coh_bbAll{isbj}.connectivity(:,:,eventslist(:,2),:);
coh_bbAll{isbj}.events       = coh_bbAll{isbj}.events(eventslist(:,2),:);

%% align channels
whichHDgrid = 'GB';

iconn_a     = coh_aAll{isbj};
iconn_bb    = coh_bbAll{isbj};
nchan = height(modeldata_a{isbj}.channels.name);
[~,ichconn,ichdat]=intersect(iconn_a.channels.name,modeldata_a{isbj}.channels.name,'stable');

iconn_a.channels  = iconn_a.channels(ichconn,:);
iconn_a.connectivity = iconn_a.connectivity(ichconn,ichconn,:,:);
iconn_bb.channels  = iconn_bb.channels(ichconn,:);
iconn_bb.connectivity = iconn_bb.connectivity(ichconn,ichconn,:,:);
modeldata_a{isbj}.channels   = modeldata_a{isbj}.channels(ichdat,:);
modeldata_bb{isbj}.channels  = modeldata_bb{isbj}.channels(ichdat,:);
prf_params_a{isbj}.channels  = prf_params_a{isbj}.channels(ichdat,:);
prf_params_bb{isbj}.channels = prf_params_bb{isbj}.channels(ichdat,:);
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

[~,ichconn,ichparams]=intersect(iconn_a.channels.name,spcrm_params{isbj}.channels.name,'stable');
assert(isequal(ichconn',1:numel(ichconn)),'Something wrong');
spcrm_params{isbj}.channels     = spcrm_params{isbj}.channels(ichparams,:);
spcrm_params{isbj}.resamp_parms = spcrm_params{isbj}.resamp_parms(ichparams,:);

%% iaf
iaf = cellfun(@(x) mean(10.^x(:,10)),spcrm_params{isbj}.resamp_parms);
iafidx{isbj} = arrayfun(@(x) find(abs(iconn_a.f - x) == min(abs(iconn_a.f - x))), iaf, 'UniformOutput', true);

end

%% update channel information
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITd_threshold;

%% set pRF parameters
datbb = model_all_bb.datats;
data  = model_all_a.datats;

%-- BLANK
   evlng = height(iconn_a.events);
   flng  = length(iconn_a.f);
   chlng = height(prf_all_a.channels);
   ndist = length(arounddist);
    evlng_dec = round(evlng/decN);
    decRate   = evlng_dec./evlng;
    boundaries = ([0 28 40 56 84 96 112 140 152 168 196 208 224]+0.5).*decRate;
    blankbnd   = boundaries([3,4;6,7;9,10;12,13]);
    medblank   = mean(blankbnd,2);
    
blankidx = false(evlng_dec,1);
for iset = 1:size(blankbnd,1)
    blankidx(ceil(blankbnd(iset,1)):fix(blankbnd(iset,2))) = true;
end

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
    
 %%
 nincl = 6;
 prfch = sum(prfidxbb,1) > nincl & sum(prfidxa,1) > nincl;
 acuch =  ~(prf_all_bb.xval<=threshold_bb | prf_all_a.xval<=threshold_a ...
            | prf_all_bb.ecc >= eclimit | prf_all_a.ecc >= eclimit)';  % use tilde for nan
 
 switch useChans
     case 'pRFchs',      selch = prfch & acuch;
     case 'SELchs',      selch = acuch;
     case 'ALLchs',      selch = true(size(acuch));
 end
 
 selchlng = sum(selch);
 prfidxa  = prfidxa(:,selch);
 prfidxbb = prfidxbb(:,selch);

%% compute coherence for around channels
for isbj = 1:length(subjectList)

% pltdat_a  = iconn_a.connectivity;
% pltdat_bb = iconn_bb.connectivity;

sbjselch     = selch(ismember(channels.subject_name,subjectList{isbj}));
sbjselchlng  = sum(sbjselch);
seedchannels = iconn_a.channels(sbjselch,:);
pairchannels = iconn_a.channels;

%%% apply decimation

pltdat_a   = decimate_col(iconn_a.connectivity(sbjselch,:,:,:),decN,[],3,'omitnan');
pltdat_bb  = decimate_col(iconn_bb.connectivity(sbjselch,:,:,:),decN,[],3,'omitnan');

datdim = size(pltdat_a);
datdim(2) = length(arounddist);
datdim(4) = datdim(3);
datdim(3) = nboot;

   dd = 0;
   avgdat_a{isbj}  = nan(datdim);
   avgdat_bb{isbj} = nan(datdim);
   for arddist = arounddist
       dd = dd + 1;
       for ich=1:sbjselchlng
           nRows = 16; nCols = 8;
           tarchnum  = str2double(strtok(seedchannels.name{ich},whichHDgrid));
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
           avgchidx  = find(ismember(pairchannels.name,avgchname));
           bootidx = randi(length(avgchidx),length(avgchidx),nboot);
           avgdat_a{isbj}(ich,dd,:,:) = mean( reshape(...
                               pltdat_a(ich,avgchidx(bootidx),:,iafidx{isbj}(ich)),...
                               [],nboot,datdim(4)) ,1,'omitnan');
           avgdat_bb{isbj}(ich,dd,:,:) = mean( reshape(...
                               pltdat_bb(ich,avgchidx(bootidx),:,1),...
                               [],nboot,datdim(4)) ,1,'omitnan');
       end
   end

end

avgdat_a  = cat(1,avgdat_a{:});
avgdat_bb = cat(1,avgdat_bb{:});
   
catch ex
    rethrow(ex);
end

%-- permute
coha   = permute(avgdat_a,[4 2 3 1]);
cohbb  = permute(avgdat_bb,[4 2 3 1]);

%% Grouping coherence into in-pRF, out-pRF, BLANK

%%% Segregate in inPRF & outPRF (alpha)
cohdatAall = squeeze(num2cell(coha,1:3));
cohdatAbsl = cohdatAall;
cohdatAprf = cohdatAall;
cohdatAout = cohdatAall;
for ich=1:selchlng
    cohdatAbsl{ich} = cohdatAbsl{ich}(blankidx,:,:,:);
    cohdatAprf{ich} = cohdatAprf{ich}(prfidxa(:,ich),:,:,:);
    cohdatAout{ich} = cohdatAout{ich}(~blankidx&~prfidxa(:,ich),:,:,:);
end

%%% Take average in each elec (alpha)
avgtsAall = cellfun(@mean,cohdatAall,repmat({1},selchlng,1),'UniformOutput',false);
avgtsAbsl = cellfun(@mean,cohdatAbsl,repmat({1},selchlng,1),'UniformOutput',false);
avgtsAprf = cellfun(@mean,cohdatAprf,repmat({1},selchlng,1),'UniformOutput',false);
avgtsAout = cellfun(@mean,cohdatAout,repmat({1},selchlng,1),'UniformOutput',false);
    avgtsAall = cat(1,avgtsAall{:});
    avgtsAbsl = cat(1,avgtsAbsl{:});
    avgtsAprf = cat(1,avgtsAprf{:});
    avgtsAout = cat(1,avgtsAout{:});
    
%%% Segregate in inPRF & outPRF (broadband)
cohdatBBall = squeeze(num2cell(cohbb,1:3));
cohdatBBbsl = cohdatBBall;
cohdatBBprf = cohdatBBall;
cohdatBBout = cohdatBBall;
for ich=1:selchlng
    cohdatBBbsl{ich} = cohdatBBbsl{ich}(blankidx,:,:,:);
    cohdatBBprf{ich} = cohdatBBprf{ich}(prfidxbb(:,ich),:,:,:);
    cohdatBBout{ich} = cohdatBBout{ich}(~blankidx&~prfidxbb(:,ich),:,:,:);
end

%%% Take average in each elec (broadband)
avgtsBBall = cellfun(@mean,cohdatBBall,repmat({1},selchlng,1),'UniformOutput',false);
avgtsBBbsl = cellfun(@mean,cohdatBBbsl,repmat({1},selchlng,1),'UniformOutput',false);
avgtsBBprf = cellfun(@mean,cohdatBBprf,repmat({1},selchlng,1),'UniformOutput',false);
avgtsBBout = cellfun(@mean,cohdatBBout,repmat({1},selchlng,1),'UniformOutput',false);
    avgtsBBall = cat(1,avgtsBBall{:});
    avgtsBBbsl = cat(1,avgtsBBbsl{:});
    avgtsBBprf = cat(1,avgtsBBprf{:});
    avgtsBBout = cat(1,avgtsBBout{:});
    
%%% bootstrapping
bootch = randi(selchlng,selchlng,nboot);
cohbootAall  = nan(nboot,ndist);
cohbootAbsl  = cohbootAall;
cohbootAprf  = cohbootAall;
cohbootAout  = cohbootAall;
cohbootBBall = nan(nboot,ndist);
cohbootBBbsl = cohbootBBall;
cohbootBBprf = cohbootBBall;
cohbootBBout = cohbootBBall;
for iboot = 1:nboot
    cohbootAall(iboot,:)  = squeeze(mean(avgtsAall(bootch(:,iboot),:,iboot),1,'omitnan'));
    cohbootAbsl(iboot,:)  = squeeze(mean(avgtsAbsl(bootch(:,iboot),:,iboot),1,'omitnan'));
    cohbootAprf(iboot,:)  = squeeze(mean(avgtsAprf(bootch(:,iboot),:,iboot),1,'omitnan'));
    cohbootAout(iboot,:)  = squeeze(mean(avgtsAout(bootch(:,iboot),:,iboot),1,'omitnan'));
    cohbootBBall(iboot,:) = squeeze(mean(avgtsBBall(bootch(:,iboot),:,iboot),1,'omitnan'));
    cohbootBBbsl(iboot,:) = squeeze(mean(avgtsBBbsl(bootch(:,iboot),:,iboot),1,'omitnan'));
    cohbootBBprf(iboot,:) = squeeze(mean(avgtsBBprf(bootch(:,iboot),:,iboot),1,'omitnan'));
    cohbootBBout(iboot,:) = squeeze(mean(avgtsBBout(bootch(:,iboot),:,iboot),1,'omitnan'));
end

saveauto(filepath,'cohboot*','nboot','disttype','useChans','arounddist','channels','selch');

end
