function    [cohdatBBall,cohdatBBbsl,cohdatBBprf,cohdatBBout,...
             cohdatAall,cohdatAbsl,cohdatAprf,cohdatAout,...
             avgtsBBall,avgtsBBbsl,avgtsBBprf,avgtsBBout,...
             avgtsAall,avgtsAbsl,avgtsAprf,avgtsAout,...
             disttype,useChans,arounddist,channels,selch,...
             smoothingMode,smoothingN] ...
    = ecog_prf_coherencedist(subjectList,disttype,useChans,nboot,opts)

% Description: 
%
% [cohdatBBall,cohdatBBbsl,cohdatBBprf,cohdatBBout,
%  cohdatAall,cohdatAbsl,cohdatAprf,cohdatAout,
%  avgtsBBall,avgtsBBbsl,avgtsBBprf,avgtsBBout,
%  avgtsAall,avgtsAbsl,avgtsAprf,avgtsAout,
%  disttype,useChans,arounddist,channels,selch,
%  smoothingMode,smoothingN]
%    = ecog_prf_coherencedist(subjectList,disttype,useChans,nboot,opts)

%
% Input
%   subjectList     cell-array of subject names
%   disttype:       how to compute the distance: 'square','diamond','norm'(default)
%   useChans:       thresholding for seed channels: 'pRFchs','SELchs'(default),'ALLchs'
%   nboot:          number of bootstrapping, no bootstrapping if 1(default)
%                   shuffling paired electrodes for each distance for each seed channel
%   opts:    
%       distwin:    distance window to average coherence in this range
% 
%       diststep:   distance step to compute coherence
%       distmax:    maximum distance to output
%           or
%       arounddist: distance between electrodes to compute coherence (ex. [1:6])
% 
% Output
%   cohdat:         {sbj_chs} (stims x dist x boot): cell-array of coherence for each distance
%   avgts:          sbj_chs x dist x boot: coherence-matrix averaging across stimulus locations
%                   If opts.isavgchs is true, averaged coherence across electrodes and stimulus locations
%                   as 1 x dist, or boot-chs x dist
%   arounddist:     distance between electrodes to compute coherence

% 20220715 Yuasa
% 20221102 Yuasa - enable to compute across precomputed electrode distances
% 20230308 Yuasa - update for low-broadband

%%
narginchk(2,inf);

%--Define inputs 
SetDefault('subjectList',{},'cell');
SetDefault('arounddist',1:10); 
SetDefault('disttype','norm');              % square, diamond, norm, <filetag>
SetDefault('useChans','SELchs');            % pRFchs, SELchs, ALLchs
SetDefault('nboot',[]);                     % 1: no bootstrapping
if isempty(nboot) || nboot < 1,    nboot = 1;
else,                              nboot = round(nboot(1));
end    
% <opts>
SetDefaultAnalysisPath('DATA','Channels','opts.eldistDir');
SetDefault('opts.isavgchs',false);       % if true, avgts is averaged across channels
SetDefault('opts.bandname_a','alpha');
SetDefault('opts.bandname_bb','broadband');
% <method>
SetDefault('opts.average','runs');
SetDefault('opts.prfmodel','linear');
SetDefault('opts.gaussianmode','gs');
SetDefault('opts.smoothingMode','decimate');
SetDefault('opts.smoothingN',3);
SetDefault('opts.method','mscoh');
SetDefault('opts.distwin',[]);
SetDefault('opts.diststep',[]);
SetDefault('opts.distmax',[]);
SetDefault('opts.arounddist','all');
% <hidden opts>
SetDefault('opts.allowlag',false);
SetDefault('opts.allowbeta',true);
SetDefault('opts.allowwide',true);
SetDefault('opts.allowmixbeta',true);
SetDefault('opts.subjectname',[]);
SetDefault('opts.cohfileid',[]);
SetDefault('opts.doloadiaf',false);

%% set up options

%%% Coherence parameters
cohmethod = opts.method;
isavgchs  = opts.isavgchs;

%%% Subject Name
nsbj = length(subjectList);
if isempty(opts.subjectname)
    if nsbj==1
       subject = subjectList{1};
    else
       subject = 'all';
    end
else
    subject = opts.subjectname;
end

%%% Distance options
needeldist = ~ismember(disttype, {'square', 'diamond', 'norm'});
arounddist = opts.arounddist;
diststep = opts.diststep;
distmax  = opts.distmax;
distwin  = opts.distwin;

%% %%%%%%  Cross-spectra  %%%%%%%%
%% get selected coherence & prf information
[coh_bbAll,coh_aAll,prfidxbb,prfidxa,blankidx,...
             useChans,~,channels,selch,...
             ~,smoothingMode,smoothingN] ...
    = ecog_prf_getprfcoherence(subjectList,useChans,opts);

%% load electrode distance
nRows = 16; nCols = 8;
if needeldist
    elecdist = cell(1,nsbj);
    for isbj = 1:nsbj
        elecdist{isbj} = load(fullfile(opts.eldistDir,sprintf('%s-elecdistance',subjectList{isbj})),'channels',disttype);
    end
end
if ~isnumeric(arounddist)
    if ismember(arounddist,{'all'})
        if needeldist
            arounddist = cellfun(@(C) C.(disttype)(:)',elecdist,'UniformOutput',false);
            arounddist = unique([arounddist{:}]);
        else
            switch disttype
                case 'square',      arounddist = 1:max(nRows-1,nCols-1);
                case 'diamond',     arounddist = 1:(nRows-1 + nCols-1);
                case 'norm',        arounddist = unique(sqrt((0:(nRows-1)).^2+(0:(nCols-1))'.^2))';
            end
        end
        arounddist(arounddist==0) = [];
        if ~isempty(distmax), arounddist(arounddist>distmax) = [];  end
        if ~isempty(distmax), arounddist = diststep:diststep:max(arounddist);   end
        if isempty(distwin), distwin = 0;  end
    else
        error('''arounddist'' is invalid');
    end
elseif isempty(distwin)     % isnumeric(arounddist)
    distwin = median(diff(arounddist),'omitnan');
    if isempty(distwin), distwin = 1;  end
end

%% compute coherence for around channels
avgdat_a  = cell(size(coh_aAll));
avgdat_bb = cell(size(coh_bbAll));
whichHDgrid = 'GB';

for isbj = 1:nsbj

sbjseedchan  = coh_aAll{isbj}.channels;
sbjpairchan  = coh_aAll{isbj}.subchannels;
sbjselchlng  = height(sbjseedchan);

%%% get distances
pltdat_a   = coh_aAll{isbj}.connectivity;
pltdat_bb  = coh_bbAll{isbj}.connectivity;

datdim = size(pltdat_a);
datdim(2) = length(arounddist);
datdim(4) = datdim(3);
datdim(3) = nboot;

   dd = 0;
   avgdat_a{isbj}  = nan(datdim);
   avgdat_bb{isbj} = nan(datdim);
   for arddist = arounddist
       dd = dd + 1;
       wineps =  eps.*arddist;
       for ich=1:sbjselchlng
           if needeldist
               seedchidx = ismember(elecdist{isbj}.channels.name,sbjseedchan.name{ich});
               dGrid     = elecdist{isbj}.(disttype)(seedchidx,:);
               dGrid     = (dGrid >= (arddist-distwin/2)) & (dGrid < (arddist+distwin/2+wineps)) & dGrid~=0;
               avgchname = elecdist{isbj}.channels.name(dGrid);
           else
               tarchnum  = str2double(strtok(sbjseedchan.name{ich},whichHDgrid));
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
               dGrid = (dGrid >= arddist) & (dGrid < (arddist+distwin+wineps));
               end
               avgchnum = find(flipud(dGrid)');
               avgchname = unique(...
                             [arrayfun(@(E) sprintf('%s%d',whichHDgrid,E),avgchnum,'UniformOutput',false),...
                              arrayfun(@(E) sprintf('%s%02d',whichHDgrid,E),avgchnum,'UniformOutput',false),...
                              arrayfun(@(E) sprintf('%s%03d',whichHDgrid,E),avgchnum,'UniformOutput',false)],...
                           'stable');
           end
           avgchidx  = find(ismember(sbjpairchan.name,avgchname));
           if nboot > 1
               if isempty(avgchidx)
               bootidx = zeros(0,nboot);
               else
               bootidx = randi(length(avgchidx),length(avgchidx),nboot);
               end
           else
               bootidx = 1:length(avgchidx);
           end
           avgdat_a{isbj}(ich,dd,:,:) = mean( reshape(...
                               pltdat_a(ich,avgchidx(bootidx),:,1),...
                               [],nboot,datdim(4)) ,1,'omitnan');
           avgdat_bb{isbj}(ich,dd,:,:) = mean( reshape(...
                               pltdat_bb(ich,avgchidx(bootidx),:,1),...
                               [],nboot,datdim(4)) ,1,'omitnan');
       end
   end

end

avgdat_a  = cat(1,avgdat_a{:});
avgdat_bb = cat(1,avgdat_bb{:});
   
%-- permute
coha   = permute(avgdat_a,[4 2 3 1]);
cohbb  = permute(avgdat_bb,[4 2 3 1]);

%% Grouping coherence into in-pRF, out-pRF, BLANK
selchlng = sum(selch);

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
avgtsAall = cell2mat(cellfun(@(COH) mean(COH,1,'omitnan'),cohdatAall,'UniformOutput',false));
avgtsAbsl = cell2mat(cellfun(@(COH) mean(COH,1,'omitnan'),cohdatAbsl,'UniformOutput',false));
avgtsAprf = cell2mat(cellfun(@(COH) mean(COH,1,'omitnan'),cohdatAprf,'UniformOutput',false));
avgtsAout = cell2mat(cellfun(@(COH) mean(COH,1,'omitnan'),cohdatAout,'UniformOutput',false));
    
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
avgtsBBall = cell2mat(cellfun(@(COH) mean(COH,1,'omitnan'),cohdatBBall,'UniformOutput',false));
avgtsBBbsl = cell2mat(cellfun(@(COH) mean(COH,1,'omitnan'),cohdatBBbsl,'UniformOutput',false));
avgtsBBprf = cell2mat(cellfun(@(COH) mean(COH,1,'omitnan'),cohdatBBprf,'UniformOutput',false));
avgtsBBout = cell2mat(cellfun(@(COH) mean(COH,1,'omitnan'),cohdatBBout,'UniformOutput',false));

%%% Take average across electrode
if isavgchs
    if nboot > 1
        %-- bootstrapping seed electrode
        bootch = randi(selchlng,selchlng,nboot);
        bootch = bootch + (0:selchlng:selchlng*(nboot-1));
        %-- alpha
        avgtsAall = permute(avgtsAall,[2,1,3]);
        avgtsAall = shiftdim(mean(reshape(avgtsAall(:,bootch),length(arounddist),selchlng,nboot),2,'omitnan'),2);
        avgtsAbsl = permute(avgtsAbsl,[2,1,3]);
        avgtsAbsl = shiftdim(mean(reshape(avgtsAbsl(:,bootch),length(arounddist),selchlng,nboot),2,'omitnan'),2);
        avgtsAprf = permute(avgtsAprf,[2,1,3]);
        avgtsAprf = shiftdim(mean(reshape(avgtsAprf(:,bootch),length(arounddist),selchlng,nboot),2,'omitnan'),2);
        avgtsAout = permute(avgtsAout,[2,1,3]);
        avgtsAout = shiftdim(mean(reshape(avgtsAout(:,bootch),length(arounddist),selchlng,nboot),2,'omitnan'),2);
        %-- broadband
        avgtsBBall = permute(avgtsBBall,[2,1,3]);
        avgtsBBall = shiftdim(mean(reshape(avgtsBBall(:,bootch),length(arounddist),selchlng,nboot),2,'omitnan'),2);
        avgtsBBbsl = permute(avgtsBBbsl,[2,1,3]);
        avgtsBBbsl = shiftdim(mean(reshape(avgtsBBbsl(:,bootch),length(arounddist),selchlng,nboot),2,'omitnan'),2);
        avgtsBBprf = permute(avgtsBBprf,[2,1,3]);
        avgtsBBprf = shiftdim(mean(reshape(avgtsBBprf(:,bootch),length(arounddist),selchlng,nboot),2,'omitnan'),2);
        avgtsBBout = permute(avgtsBBout,[2,1,3]);
        avgtsBBout = shiftdim(mean(reshape(avgtsBBout(:,bootch),length(arounddist),selchlng,nboot),2,'omitnan'),2);

    else
        %-- alpha
        avgtsAall = mean(avgtsAall,1,'omitnan');
        avgtsAbsl = mean(avgtsAbsl,1,'omitnan');
        avgtsAprf = mean(avgtsAprf,1,'omitnan');
        avgtsAout = mean(avgtsAout,1,'omitnan');
        %-- broadband
        avgtsBBall = mean(avgtsBBall,1,'omitnan');
        avgtsBBbsl = mean(avgtsBBbsl,1,'omitnan');
        avgtsBBprf = mean(avgtsBBprf,1,'omitnan');
        avgtsBBout = mean(avgtsBBout,1,'omitnan');
    end
end

end

