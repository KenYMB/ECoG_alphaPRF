function    [coh_bbAll,coh_aAll,prfidxbb,prfidxa,blankidx,...
             useChans,seedchannels,pairchannels,selch,...
             iaf,smoothingMode,smoothingN] ...
    = ecog_prf_getprfcoherence(subjectList,useChans,opts)

% Description: 
%
% [coh_bbAll,coh_aAll,prfidxbb,prfidxa,blankidx,...
%  useChans,seedchannels,pairchannels,selch,...
%  iaf,smoothingMode,smoothingN] ...
%    = ecog_prf_getprfcoherence(subjectList,useChans,opts)

%
% Input
%   subjectList     cell-array of subject names
%   useChans:       thresholding for seed channels: 'pRFchs','SELchs'(default),'ALLchs'
%   opts:    
% 
% Output
%   coh_All:        {sbj} (seedchs x pairchs x stims x 1): cell-array of coherence
%   prfidx:         index indicating if stimulus is inside the pRF or not
%   blankidx:       index indicating blank trials
%   iaf:            peak alpha frequnecy for each channel

% 20221110 Yuasa - segregate from ecog_prf_coherencedist
% 20230308 Yuasa - update for low-broadband

%%
narginchk(2,inf);

%--Define inputs 
SetDefault('subjectList',{},'cell');
SetDefault('useChans','SELchs');            % pRFchs, SELchs, ALLchs
% <opts>
SetDefault('opts.bandname_a','alpha');
SetDefault('opts.bandname_bb','broadband');
% <method>
SetDefault('opts.average','runs');
SetDefault('opts.prfmodel','linear');
SetDefault('opts.gaussianmode','gs');
SetDefault('opts.smoothingMode','decimate');
SetDefault('opts.smoothingN',3);
SetDefault('opts.method','mscoh');
% <hidden opts>
SetDefault('opts.allowlag',false);
SetDefault('opts.allowbeta',true);
SetDefault('opts.allowwide',true);
SetDefault('opts.allowmixbeta',true);
SetDefault('opts.cohfileid',[]);
SetDefault('opts.doloadiaf',false);

%% load time series data
average        = opts.average;
prfmodel       = opts.prfmodel;
gaussianmode   = opts.gaussianmode;
smoothingMode  = opts.smoothingMode;
smoothingN     = opts.smoothingN;
selectchs      = 'GB*';             % use all grid channels
% selectchs      = 'wangprobchs';     % only use wangprobchs as seed but use all grid channels for averaged coherence
    allowlag       = opts.allowlag;
    allowbeta      = opts.allowbeta;
    allowwide      = opts.allowwide;
    allowmixbeta   = opts.allowmixbeta;
    
ecog_APRFF_INITa_loaddata;

%% set up options
%%% Coherence parameters
cohmethod = opts.method;

%%% Subjects
nsbj = length(subjectList);

%% Load Coherence for alpha
opt = [];
opt.compute    = false;
opt.issave     = false;
opt.average    = average;
opt.isavgfreq  = false;
opt.method     = cohmethod;
opt.targetBAND = opts.bandname_a;
opt.fileid     = opts.cohfileid;
coh_aAll    = ecog_prf_loadconnectivity(subjectList, opt);

%% Load Coherence for broadband
opt = [];
opt.compute    = false;
opt.issave     = false;
opt.average    = average;
opt.isavgfreq  = true;
opt.method     = cohmethod;
opt.targetBAND = opts.bandname_bb;
opt.fileid     = opts.cohfileid;
coh_bbAll   = ecog_prf_loadconnectivity(subjectList, opt);

%% load iaf
opt = [];
opt.average    = average;
opt.compute         = false;
spcrm_params = ecog_prf_fitalpha(subjectList, opt);

%% rearrange events 
iaf      = cell(size(coh_aAll));
iafidx   = cell(size(coh_aAll));
for isbj = 1:nsbj

%% align channels

nchan = height(modeldata_a{isbj}.channels.name);
[~,ichconn,ichdat]=intersect(coh_aAll{isbj}.channels.name,modeldata_a{isbj}.channels.name,'stable');

coh_aAll{isbj}.channels  = coh_aAll{isbj}.channels(ichconn,:);
coh_aAll{isbj}.connectivity = coh_aAll{isbj}.connectivity(ichconn,ichconn,:,:);
coh_bbAll{isbj}.channels  = coh_bbAll{isbj}.channels(ichconn,:);
coh_bbAll{isbj}.connectivity = coh_bbAll{isbj}.connectivity(ichconn,ichconn,:,:);
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

[~,ichconn,ichparams]=intersect(coh_aAll{isbj}.channels.name,spcrm_params{isbj}.channels.name,'stable');
assert(isequal(ichconn',1:numel(ichconn)),'Something wrong');
spcrm_params{isbj}.channels     = spcrm_params{isbj}.channels(ichparams,:);
spcrm_params{isbj}.resamp_parms = spcrm_params{isbj}.resamp_parms(ichparams,:);

%% iaf
if opts.doloadiaf
iaf{isbj} = cellfun(@(x) mean(10.^x(:,10)),spcrm_params{isbj}.resamp_parms);
iafidx{isbj} = arrayfun(@(x) find(abs(coh_aAll{isbj}.f - x) == min(abs(coh_aAll{isbj}.f - x))), iaf{isbj}, 'UniformOutput', true);
end

end

%% update channel information
ecog_APRFF_INITb_mergedata;
ecog_APRFF_INITd_threshold;

%% set pRF parameters
datbb = model_all_bb.datats;
% data  = model_all_a.datats;

%-- BLANK
   evlng = height(coh_aAll{isbj}.events);
%    flng  = length(coh_aAll{isbj}.f);
   chlng = height(prf_all_a.channels);
%    ndist = length(arounddist);
    evlng_dec = round(evlng/smoothingN);
    decRate   = evlng_dec./evlng;
    boundaries = ([0 28 40 56 84 96 112 140 152 168 196 208 224]+0.5).*decRate;
    blankbnd   = boundaries([3,4;6,7;9,10;12,13]);
%     medblank   = mean(blankbnd,2);
    
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
%     dattsbb = repmat({{}},chlng,1);  dattsa = repmat({{}},chlng,1);
    for ich=1:chlng
        for pp=1:numruns
            nanplace    = isnan(datbb{pp}(ich,:));
%             dattsbb{ich}{pp}   = datbb{pp}(ich,~nanplace)';
%             dattsa{ich}{pp}    = -data{pp}(ich,~nanplace)';
            modelbb{ich}{pp} = modelfun(prf_all_bb.params(1,:,ich),stimulusPP{pp}(~nanplace,:));
            modela{ich}{pp}  = -modelfun(prf_all_a.params(1,:,ich),stimulusPP{pp}(~nanplace,:));
        end
%         dattsbb{ich} = cat(1,dattsbb{ich}{:});
%         dattsa{ich}  = cat(1,dattsa{ich}{:});
        modelbb{ich} = cat(1,modelbb{ich}{:});
        modela{ich}  = cat(1,modela{ich}{:});
    end
%     dattsbb = cat(2,dattsbb{:});
%     dattsa  = cat(2,dattsa{:});
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
    
 %% Select channels
 nincl = 6;
 prfch = sum(prfidxbb,1) > nincl & sum(prfidxa,1) > nincl;
 acuch =  ~(prf_all_bb.xval<=threshold_bb | prf_all_a.xval<=threshold_a ...
            | prf_all_bb.ecc >= eclimit | prf_all_a.ecc >= eclimit)';  % use tilde for nan
 
 switch useChans
     case 'pRFchs',      selch = prfch & acuch;
     case 'SELchs',      selch = acuch;
     case 'ALLchs',      selch = true(size(acuch));
 end
 
 %% Apply selection
 %  selchlng = sum(selch);
 prfidxa  = prfidxa(:,selch);
 prfidxbb = prfidxbb(:,selch);

 seedchannels = channels(selch,:);
 pairchannels = channels;

 for isbj = 1:nsbj
     %% seed channel selection 
     sbjselch   = selch(ismember(channels.subject_name,subjectList{isbj}));
     cohdat_a   = coh_aAll{isbj}.connectivity(sbjselch,:,:,:);
     cohdat_bb  = coh_bbAll{isbj}.connectivity(sbjselch,:,:,:);

     %% apply decimation
     switch smoothingMode
         case {'smooth'}
             stimtypelist = grp2idx(cellfun(@(x) strtok(x,'-'),prf_all_bb.events.trial_name,'UniformOutput',false));
             boundaries = [0 reshape(find(diff(stimtypelist)~=0),1,[]) height(prf_all_bb.events)];
             for jj=1:length(boundaries)-1
                 cohdat_a(:,:,(boundaries(jj)+1):boundaries(jj+1),:) ...
                     = smoothdata(cohdat_a(:,:,(boundaries(jj)+1):boundaries(jj+1),:),3,'movmean',smoothingN);
                 cohdat_bb(:,:,(boundaries(jj)+1):boundaries(jj+1),:) ...
                     = smoothdata(cohdat_bb(:,:,(boundaries(jj)+1):boundaries(jj+1),:),3,'movmean',smoothingN);
             end
         case {'decimate'}
             cohdat_a   = decimate_col(cohdat_a,smoothingN,[],3,'omitnan');
             cohdat_bb  = decimate_col(cohdat_bb,smoothingN,[],3,'omitnan');
     end

     %% return
     coh_aAll{isbj}.subchannels     = coh_aAll{isbj}.channels;
     coh_bbAll{isbj}.subchannels    = coh_bbAll{isbj}.channels;
     coh_aAll{isbj}.channels        = coh_aAll{isbj}.channels(sbjselch,:);
     coh_bbAll{isbj}.channels       = coh_bbAll{isbj}.channels(sbjselch,:);
     coh_aAll{isbj}.connectivity    = cohdat_a;
     coh_bbAll{isbj}.connectivity   = cohdat_bb;
     
     %% iaf
     if opts.doloadiaf
     sbjiafidx                   = iafidx{isbj}(sbjselch);
     else
     cohAvg = squeeze(mean(coh_aAll{isbj}.connectivity,[2,3],'omitnan'))';
     regF   = reshape(coh_aAll{isbj}.f,[],1);  regF   = [regF,regF];  regF(:,2) = 1;
     cohSlp = regF\cohAvg;
     cohAvg = cohAvg - cohSlp(1,:).*regF(:,1);          % regress out linear trend
     cohDif = diff(cohAvg,1,1);   cohDif = cat(1,ones(1,size(cohDif,2)),cohDif(1:(end-1),:).*cohDif(2:end,:),ones(1,size(cohDif,2)));
     cohDif = cohDif<=0;  cohDif(:,sum(cohDif,1)==0)=1; % detect peaks OR keep all values if no-peak
     cohAvg(~cohDif) = nan;
     [~,sbjiafidx]=max(cohAvg,[],1,'omitnan');          % find max at peaks
     sbjiafidx = reshape(sbjiafidx,[],1);
     end
     coh_aAll{isbj}.f            = coh_aAll{isbj}.f(sbjiafidx);
     outmatdim  = size(coh_aAll{isbj}.connectivity);  outmatdim(end) = [];
     iafidxmat  = reshape(1:prod(outmatdim),outmatdim) + (sbjiafidx-1)*prod(outmatdim);
     coh_aAll{isbj}.connectivity = coh_aAll{isbj}.connectivity(iafidxmat);

 end

end
