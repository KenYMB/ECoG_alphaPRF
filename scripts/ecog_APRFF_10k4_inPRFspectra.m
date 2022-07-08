% compute ECoG spectra in/out pRFs

% 20211230 Yuasa

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
SbjInfo    = loadSbjInfo(subjectList_fname,'all');

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
    
%% load spectra
opts = [];
opts.compute        = false;
opts.doplots        = false;
opts.allowlag       = allowlag;

[freq] = ecog_prf_spectra(subjectList, opts);
[freq]   = ecog_summarizeROIs(freq);
[freq]   = ecog_channelSelection(freq,selectchs,selectch_exFEF,selectch_thresh);
%--- exclude power artifact & concatenate across subjects
for isbj=1:length(freq)
    f = freq{isbj}.f;
    RECsite = SbjInfo.site(ismember(SbjInfo.participant_id,freq{isbj}.subject));
    RECsite = RECsite{1};
    switch upper(RECsite)
        case {'NYU'}    % NYU subjects: 60Hz line noise
%         f_noise  = f((f>52 & f<70) | (f>115 & f<126) | (f>175 & f<186));
        f_noise  = f((f>52 & f<70) | (f>111 & f<130) | (f>171 & f<190));
        case {'UMCU'}   % Utrecht subjects: 50Hz line noise
        f_noise  = f((f>42 & f<60) | (f>95 & f<106) | (f>145 & f<156));
        otherwise
        f_noise = [];
    end
    freq{isbj}.spectra(:,:,f_noise)     = nan;
    freq{isbj}.spectra_off(:,:,f_noise) = nan;
end
[freq_all]    = ecog_rearrangePRF(ecog_averageTimeSeries(freq,'runs','realign'));
%-- apply decimation
freq_all.spectra = decimate_col(freq_all.spectra,3,[],2,'omitnan');
freq_all.spectra_rate = freq_all.spectra./freq_all.spectra_off;

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
else
rois = {{'V1','V2','V3'},{'V3a','V3b','LO1','LO2','TO','IPS'}};
roilabels = {'V1-V3','Dorsolateral'};
end

nroi = length(rois);

%% %%%%%%%%%%%%%%%%%%%%%%
%% Show figure
%% %%%%%%%%%%%%%%%%%%%%%%
lowF   = [3  26];
highF  = [70 180];
alphaF = [6 17];
%% Power Spectra
datats  = permute(freq_all.spectra,[2,1,3]);
datatsR = 1./permute(freq_all.spectra_rate,[2,1,3]);
datats(datats<=0) = nan;  datatsR(datatsR<=0) = nan;
f      = freq_all.f;

meanfun = @(d) (geomean(d,1,'omitnan'));
% errfun  = @(d) (geosem(d,0,1,'omitnan'));
% errneg  = @(m,s) -exp(log(m)-log(s)) + m;
% errpos  = @(m,s) exp(log(m)+log(s)) - m;

nchn      = size(datats,2);
nfreq     = size(datats,3);
prfinmat  = prfidxbb&~blankidx;
prfoutmat = ~prfidxbb&~blankidx;
blankmat  = blankidx;
datatsIn  = nan(nchn,nfreq);   datatsOut  = nan(nchn,nfreq);
datatsBlank = nan(nchn,nfreq); datatsRate = nan(nchn,nfreq);
for ich=1:size(datats,2)
  if sum(prfinmat(:,ich))>0
    datatsIn(ich,:)   = squeeze(meanfun(datats(prfinmat(:,ich),ich,:)));
    datatsRate(ich,:) = squeeze(meanfun(datatsR(prfinmat(:,ich),ich,:)));
  end
  if sum(prfoutmat(:,ich))>0
    datatsOut(ich,:)  = squeeze(meanfun(datats(prfoutmat(:,ich),ich,:)));
  end
  if sum(blankmat)>0
    datatsBlank(ich,:) = squeeze(meanfun(datats(blankmat,ich,:)));
  end
end

hF(end+1) = figure; hT=tiledlayout(1,nroi,'Padding','compact','TileSpacing','compact');
set(gcf,'Position',get(gcf,'Position').*[3./(nroi+2) .5 .2+.4.*nroi 1.1]);
for iroi = 1:length(rois)
    roich = ismember(channels.wangarea,rois{iroi})';

inPRF_bb  = meanfun(datatsIn(roich&selch,:));
outPRF_bb = meanfun(datatsOut(roich&selch,:));
blnak_bb = meanfun(datatsBlank(roich&selch,:));
% E_inPRF_bb   = errfun(datatsIn(roich&selch,:));
% E_outPRF_bb  = errfun(datatsOut(roich&selch,:));

nexttile;
B=loglog(f,inPRF_bb,f,blnak_bb,'LineWidth',2.0);
B(1).Color = '#A06440';
B(2).Color = 'k';
set(gca,'FontSize',FntSiz);

%-- visual
if iroi==1
    yticklabels(cellfun(@(x) sprintf('%g',x),num2cell(yticks),'UniformOutput',false));
else
    yticklabels({});
end
xlim(minmax([lowF,highF]));
xticks(unique([lowF,highF]))
xticklabels(cellfun(@(x) sprintf('%g',x),num2cell(xticks),'UniformOutput',false));
title(roilabels{iroi});

end
legend(B,{'In pRF','BLANK'},'AutoUpdate','off');

ylin  = ylim;
hax   = get(hT,'Children');
faidx = f(f>=alphaF(1) & f<=alphaF(2));
for iroi = 1:length(rois)
    roich = ismember(channels.wangarea,rois{iroi})';
    
hA = hax(end-iroi+1);
hold(hA,'on');
%-- show broadband Box
pb1=plotbox(hA,highF,ylim,plcol(1,:),'FaceAlpha',0.3);
pb2=plotbox(hA,lowF,ylim,plcol(5,:),'FaceAlpha',0.3);
uistack(pb1,'bottom');
uistack(pb2,'bottom');

%-- show alpha
ratePRF_bb  = meanfun(datatsRate(roich&selch,:));
[~,f_alpha] = findpeaks(ratePRF_bb(faidx),f(faidx),'NPeaks',1);
plot(hA,[f_alpha f_alpha],ylim,'--','Color',plcol(2,:),'LineWidth',2.0);

ylim(hA,ylin);
hold(hA,'off');
end

xlabel(hT,'Frequency (Hz)','FontSize',FntSiz,'FontWeight','normal');
ylabel(hT,'Power (\muV^2/Hz)','FontSize',FntSiz,'FontWeight','normal');
set(gcf,'Name','ECoG Spectra');

figname = sprintf('Power-spectral_%02d%%-%02d%%-ecc%02d_%s',threshold_bb,threshold_a,eclimit,useChans);
if issaveplot,   savefigauto(gcf,fullfile(plotsavedir, figname));   end

