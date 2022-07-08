% visualize connectivities across electrodes for subjects with HDgrid

% 20201115 Yuasa
% 20210516 Yuasa - modify 

%%
% close all; clear all;
% if isempty(gcp('nocreate')),  parpool([1 40]); end
% startupToolboxToolbox;

%% Define paths and dataset
checkPath;
%-- Input & Output path
plotsavePth    = 'Connectivity';
%-- Set save figure dirctory
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth));
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'hasHDgrid','yes');

%% %%%%%%%%%%%%%%%%%%%%%%
average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';
va_area        = 'wangarea';
connmethod     ='mscoh';

%% Load Data
%%% Load Reference
ecog_APRFF_INITa_loaddata;
ecog_APRFF_INITb_mergedata;

%%% Set threshold
ecog_APRFF_INITd_threshold;

%%% Load Coherence

opts = [];
opts.compute    = false;
opts.issave     = false;
opts.average    = 'stimuli';
opts.method     = connmethod;
opts.targetBAND     = 'broadband';
opts.stimNames      = 'HORIZONTAL*';
[coh_bb1]   = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'VERTICAL*';
[coh_bb2]   = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'DIAGONAL*';
[coh_bb3]   = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'BLANK';
[coh_bb0]   = ecog_prf_connectivity(subjectList, opts);
opts.targetBAND     = 'alpha';
opts.stimNames      = 'HORIZONTAL*';
[coh_a1]    = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'VERTICAL*';
[coh_a2]    = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'DIAGONAL*';
[coh_a3]    = ecog_prf_connectivity(subjectList, opts);
opts.stimNames      = 'BLANK';
[coh_a0]    = ecog_prf_connectivity(subjectList, opts);


%% %%%%%%%%%%%%%%%%%%%%%
%% Plot parameters
res       = [100 100];
resmx     = max(res);
cfactor   = 16.6./resmx;
close all;

alpha = 0.32;   % 0.05 for 2sd, 0.32 for 1sd

plcol = get(groot,'defaultAxesColorOrder');

plotsavedir    = fullfile(plotsavedir, connmethod);
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

%%% Grid plot
whichHDgrid = 'GB';

%% plot Coherence (& average across around electrodes)
close all;

isbj = 1;

band = 'Alpha';
stim = 'HORIZONTAL';

setmask   = false;
basecorr  = true;
avgaround = false;

for isbj = 1:2

for band = ["Alpha","Broadband"]
for stim = ["HORIZONTAL","VERTICAL","DIAGONAL"]

for basecorr  = [true, false]
for avgaround = [false, true]
if avgaround
    arounddist = [1 2 3];
else
    arounddist = 1;
end
for arddist = arounddist

if isbj==1
seedch = 'GB102';
else
seedch = 'GB068';
end

%-- main
tarBAND = band;
switch lower(tarBAND)
    case {'bb','broadband'}
        thresh  = threshold_bb;
        ifreq  = 90;
        connbase  = coh_bb0{isbj};
        switch stim
            case 'HORIZONTAL',  iconn  = coh_bb1{isbj};
            case 'VERTICAL',    iconn  = coh_bb2{isbj};
            case 'DIAGONAL',    iconn  = coh_bb3{isbj};
        end
    case {'a','alpha'}
        thresh  = threshold_a;
        ifreq  = 13;
        connbase  = coh_a0{isbj};
        switch stim
            case 'HORIZONTAL',  iconn  = coh_a1{isbj};
            case 'VERTICAL',    iconn  = coh_a2{isbj};
            case 'DIAGONAL',    iconn  = coh_a3{isbj};
        end
end
allconn  = cat(4,cat(3,coh_bb1{isbj}.connectivity,coh_bb2{isbj}.connectivity,coh_bb3{isbj}.connectivity),...
                 cat(3,coh_a1{isbj}.connectivity,coh_a2{isbj}.connectivity,coh_a3{isbj}.connectivity));
allbase  = cat(4,coh_bb0{isbj}.connectivity,coh_a0{isbj}.connectivity);
subject  = iconn.subject;

if setmask
    [~,selchan] = intersect(prf_params_bb{isbj}.channels.name,iconn.channels.name,'stable');
    elec_ok = ~(prf_params_bb{isbj}.xval<=threshold_bb | prf_params_a{isbj}.xval<=threshold_a ...
         | prf_params_bb{isbj}.ecc >= eclimit | prf_params_a{isbj}.ecc >= eclimit); 
    elec_ok = elec_ok(selchan);
else
    elec_ok = true(height(iconn.channels),1);
end

%% plot Coherence with interactive UI
specs = [];
specs.channels          = iconn.channels;
specs.plot.nSubPlots    = [];
specs.plot.RotGrid      = true;
specs.plot.fontSize     = 16;

specs.seed = seedch;
specs.plot.ShowSeed = true;
specs.plot.Time = round(height(iconn.events)/2)+1;

switch connmethod
    case {'coh','mscoh'}
        %-- Varianve Explained
        specs.plot.FigName      = sprintf(sprintf('Coherence (%s)',tarBAND));
        specs.plot.colorMap     = colormapskew('hot',0.8,[0 0.95]);
        specs.plot.AlphaData = elec_ok;
        if basecorr
%             maxconn = max(abs(iconn.connectivit) -  abs(connbase.connectivity),[],'all');
            maxconn = max(abs(abs(allconn)-abs(allbase)),[],'all');
            specs.plot.YLim      = [-1 1]* round(maxconn,2);
            pltdat      = abs(iconn.connectivity) - ...
                abs(connbase.connectivity);
            poscol = colormapskew('hot',0.8,[0 0.95]);
            negcol = flipud(poscol(:,[3,2,1]));
            specs.plot.colorMap     = cat(1,negcol,poscol);
        else
            specs.plot.YLim      = [0 1.0];
            pltdat      = abs(iconn.connectivity);
            specs.plot.colorMap     = colormapskew('hot',0.8,[0 0.95]);
        end
end
if ~avgaround 
   specs.timelabel      = iconn.events.trial_name;
else
   specs.timelabel      = cellfun(@(C) sprintf('Around %d - %s',arddist,C),iconn.events.trial_name,'UniformOutput',false);
   for ich=1:height(iconn.channels)
       nRows = 16;
       tarchnum  = str2double(strtok(iconn.channels.name{ich},whichHDgrid));
       aroundch  = reshape(unique([[nRows -nRows]'*arddist+[-arddist:arddist]; [(-nRows*arddist):nRows:(nRows*arddist)]+[-arddist;arddist]]),1,[]);
       avgchnum = tarchnum+aroundch;
       if mod(tarchnum-1-arddist,nRows)+1 > mod(tarchnum-1,nRows)+1   % when leftside
           avgchnum(mod(avgchnum-1,nRows)+1 >= mod(tarchnum-1-arddist,nRows)+1) = [];
       end
       if mod(tarchnum-1+arddist,nRows)+1 < mod(tarchnum-1,nRows)+1   % when rightside
           avgchnum(mod(avgchnum-1,nRows)+1 <= mod(tarchnum-1+arddist,nRows)+1) = [];
       end
       avgchnum(avgchnum<1 | avgchnum>height(iconn.channels)) = [];   % when top or bottom
       avgchname = unique(...
                     [arrayfun(@(E) sprintf('%s%d',whichHDgrid,E),avgchnum,'UniformOutput',false),...
                      arrayfun(@(E) sprintf('%s%02d',whichHDgrid,E),avgchnum,'UniformOutput',false),...
                      arrayfun(@(E) sprintf('%s%03d',whichHDgrid,E),avgchnum,'UniformOutput',false)],...
                   'stable');
       avgchidx  = ismember(iconn.channels.name,avgchname);
       pltdat(ich,1,:) = mean(pltdat(avgchidx,ich,:),1,'omitnan');
   end
   pltdat(:,2:end,:) = [];
   seedch = 'AvgRound';
end

hF = ecog_plotGridXSCts(pltdat, whichHDgrid, specs);

%% save
if setmask
    fignamebase = fullfile(plotsavedir, sprintf('%sGrid_%s-%02d%%-%02d%%-ecc%02d_%s-%s',connmethod,subject,threshold_bb,threshold_a,eclimit,tarBAND,stim));
else
    fignamebase = fullfile(plotsavedir, sprintf('%sGrid_%s_%s-%s',connmethod,subject,tarBAND,stim));
end
if basecorr
    fignamebase = sprintf('%s-basecorr',fignamebase);
end
if avgaround
    fignamebase = sprintf('%s_AvgRound%d',fignamebase,arddist);
    figname     = fignamebase;
else
    figname     = sprintf('%s_ref-%s',fignamebase,seedch);
end
figname = sprintf('%s_%s',figname,iconn.events.trial_name{specs.plot.Time});
if numel(hF)==1
    saveas(hF{1},fignamebase);
    hgexport(hF{1},figname,hgexport('factorystyle'),'Format','png');
else
    for ii=1:numel(hF)
    saveas(hF{1},sprintf('%s-%d',fignamebase,ii));
    hgexport(hF{ii},sprintf('%s-%d',figname,ii),hgexport('factorystyle'),'Format','png');
    end
end
close(hF{:});
end
end
end
end
end
end
%%

%% %%%%%%%%%%%%%%%%%%%%%%%%
%% plot Coherence for BLANK
close all;

isbj = 1;

band = 'Alpha';
stim = 'BLANK';

setmask   = false;
basecorr  = false;
avgaround = false;

for isbj = 1:2

for band = ["Alpha","Broadband"]
for avgaround = [false, true]
if avgaround
    arounddist = [1 2 3];
else
    arounddist = 1;
end
for arddist = arounddist

if isbj==1
seedch = 'GB103';
% seedch = 'GB102';
else
seedch = 'GB082';
% seedch = 'GB068';
end

%-- main
tarBAND = band;
switch lower(tarBAND)
    case {'bb','broadband'}
        thresh  = threshold_bb;
        ifreq  = 90;
        connbase  = coh_bb0{isbj};
        switch stim
            case 'BLANK',  iconn  = coh_bb0{isbj};
        end
    case {'a','alpha'}
        thresh  = threshold_a;
        ifreq  = 13;
        connbase  = coh_a0{isbj};
        switch stim
            case 'BLANK',  iconn  = coh_a0{isbj};
        end
end
allconn  = cat(4,cat(3,coh_bb1{isbj}.connectivity,coh_bb2{isbj}.connectivity,coh_bb3{isbj}.connectivity),...
                 cat(3,coh_a1{isbj}.connectivity,coh_a2{isbj}.connectivity,coh_a3{isbj}.connectivity));
allbase  = cat(4,coh_bb0{isbj}.connectivity,coh_a0{isbj}.connectivity);
subject  = iconn.subject;

if setmask
    [~,selchan] = intersect(prf_params_bb{isbj}.channels.name,iconn.channels.name,'stable');
    elec_ok = ~(prf_params_bb{isbj}.xval<=threshold_bb | prf_params_a{isbj}.xval<=threshold_a ...
         | prf_params_bb{isbj}.ecc >= eclimit | prf_params_a{isbj}.ecc >= eclimit); 
    elec_ok = elec_ok(selchan);
else
    elec_ok = true(height(iconn.channels),1);
end

%% plot Coherence with interactive UI
specs = [];
specs.channels          = iconn.channels;
specs.plot.nSubPlots    = [];
specs.plot.RotGrid      = true;
specs.plot.fontSize     = 16;

specs.seed = seedch;
specs.plot.ShowSeed = true;
specs.plot.Time = round(height(iconn.events)/2);

switch connmethod
    case {'coh','mscoh'}
        %-- Varianve Explained
        specs.plot.FigName      = sprintf(sprintf('Coherence (%s)',tarBAND));
        specs.plot.colorMap     = colormapskew('hot',0.8,[0 0.95]);
        specs.plot.AlphaData = elec_ok;
        if basecorr
%             maxconn = max(abs(iconn.connectivit) -  abs(connbase.connectivity),[],'all');
            maxconn = max(abs(abs(allconn)-abs(allbase)),[],'all');
            specs.plot.YLim      = [-1 1]* round(maxconn,2);
            pltdat      = abs(iconn.connectivity) - ...
                abs(connbase.connectivity);
            poscol = colormapskew('hot',0.8,[0 0.95]);
            negcol = flipud(poscol(:,[3,2,1]));
            specs.plot.colorMap     = cat(1,negcol,poscol);
        else
            specs.plot.YLim      = [0 1.0];
            pltdat      = abs(iconn.connectivity);
            specs.plot.colorMap     = colormapskew('hot',0.8,[0 0.95]);
        end
end
if ~avgaround 
   specs.timelabel      = iconn.events.trial_name;
else
   specs.timelabel      = cellfun(@(C) sprintf('Around %d - %s',arddist,C),iconn.events.trial_name,'UniformOutput',false);
   for ich=1:height(iconn.channels)
       nRows = 16;
       tarchnum  = str2double(strtok(iconn.channels.name{ich},whichHDgrid));
       aroundch  = reshape(unique([[nRows -nRows]'*arddist+[-arddist:arddist]; [(-nRows*arddist):nRows:(nRows*arddist)]+[-arddist;arddist]]),1,[]);
       avgchnum = tarchnum+aroundch;
       if mod(tarchnum-1-arddist,nRows)+1 > mod(tarchnum-1,nRows)+1   % when leftside
           avgchnum(mod(avgchnum-1,nRows)+1 >= mod(tarchnum-1-arddist,nRows)+1) = [];
       end
       if mod(tarchnum-1+arddist,nRows)+1 < mod(tarchnum-1,nRows)+1   % when rightside
           avgchnum(mod(avgchnum-1,nRows)+1 <= mod(tarchnum-1+arddist,nRows)+1) = [];
       end
       avgchnum(avgchnum<1 | avgchnum>height(iconn.channels)) = [];   % when top or bottom
       avgchname = unique(...
                     [arrayfun(@(E) sprintf('%s%d',whichHDgrid,E),avgchnum,'UniformOutput',false),...
                      arrayfun(@(E) sprintf('%s%02d',whichHDgrid,E),avgchnum,'UniformOutput',false),...
                      arrayfun(@(E) sprintf('%s%03d',whichHDgrid,E),avgchnum,'UniformOutput',false)],...
                   'stable');
       avgchidx  = ismember(iconn.channels.name,avgchname);
       pltdat(ich,1,:) = mean(pltdat(avgchidx,ich,:),1,'omitnan');
   end
   pltdat(:,2:end,:) = [];
   seedch = 'AvgRound';
end

hF = ecog_plotGridXSCts(pltdat, whichHDgrid, specs);

%% save
if setmask
    fignamebase = sprintf('%sGrid_%s-%02d%%-%02d%%-ecc%02d_%s-%s',connmethod,subject,threshold_bb,threshold_a,eclimit,tarBAND,stim);
else
    fignamebase = sprintf('%sGrid_%s_%s-%s',connmethod,subject,tarBAND,stim);
end
if basecorr
    fignamebase = sprintf('%s-basecorr',fignamebase);
end
if avgaround
    fignamebase = sprintf('%s_AvgRound%d',fignamebase,arddist);
    figname     = fignamebase;
else
    figname     = sprintf('%s_ref-%s',fignamebase,seedch);
end
% figname = sprintf('%s_%s',figname,iconn.events.trial_name{specs.plot.Time});
if numel(hF)==1
    saveas(hF{1},fullfile(plotsavedir,fignamebase));
    hgexport(hF{1},fullfile(plotsavedir,figname),hgexport('factorystyle'),'Format','png');
else
    for ii=1:numel(hF)
    saveas(hF{ii},fullfile(plotsavedir,sprintf('%s-%d',fignamebase,ii)));
    hgexport(hF{ii},fullfile(plotsavedir,sprintf('%s-%d',figname,ii)),hgexport('factorystyle'),'Format','png');
    end
end
close(hF{:});
end
end
end
end
%%


%% %%%%%%%%%%%%%%%%%%%%%%%%
%% Timecourse


isbj = 1;

band = 'Alpha';
stim = 'HORIZONTAL';

setmask   = false;
basecorr  = false;
avgaround = true;

    arounddist = [1 2 3 6];

for isbj = 1:2

for band = ["Alpha","Broadband"]
for stim = ["HORIZONTAL","VERTICAL","DIAGONAL"]

if isbj==1
seedch = 'GB102';
else
seedch = 'GB068';
end

%-- main
tarBAND = band;
switch lower(tarBAND)
    case {'bb','broadband'}
        thresh  = threshold_bb;
        ifreq  = 90;
        connbase  = coh_bb0{isbj};
        switch stim
            case 'HORIZONTAL',  iconn  = coh_bb1{isbj};
            case 'VERTICAL',    iconn  = coh_bb2{isbj};
            case 'DIAGONAL',    iconn  = coh_bb3{isbj};
            case 'BLANK',       iconn  = coh_bb0{isbj};
        end
    case {'a','alpha'}
        thresh  = threshold_a;
        ifreq  = 13;
        connbase  = coh_a0{isbj};
        switch stim
            case 'HORIZONTAL',  iconn  = coh_a1{isbj};
            case 'VERTICAL',    iconn  = coh_a2{isbj};
            case 'DIAGONAL',    iconn  = coh_a3{isbj};
            case 'BLANK',       iconn  = coh_a0{isbj};
        end
end
allconn  = cat(4,cat(3,coh_bb1{isbj}.connectivity,coh_bb2{isbj}.connectivity,coh_bb3{isbj}.connectivity),...
                 cat(3,coh_a1{isbj}.connectivity,coh_a2{isbj}.connectivity,coh_a3{isbj}.connectivity));
allbase  = cat(4,coh_bb0{isbj}.connectivity,coh_a0{isbj}.connectivity);
subject  = iconn.subject;

if setmask
    [~,selchan] = intersect(prf_params_bb{isbj}.channels.name,iconn.channels.name,'stable');
    elec_ok = ~(prf_params_bb{isbj}.xval<=threshold_bb | prf_params_a{isbj}.xval<=threshold_a ...
         | prf_params_bb{isbj}.ecc >= eclimit | prf_params_a{isbj}.ecc >= eclimit); 
    elec_ok = elec_ok(selchan);
else
    elec_ok = true(height(iconn.channels),1);
end

%%% plot Coherence with interactive UI
specs = [];
specs.channels          = iconn.channels;
specs.plot.nSubPlots    = [];
specs.plot.RotGrid      = true;
specs.plot.fontSize     = 16;

specs.seed = seedch;
specs.plot.ShowSeed = true;
specs.plot.Time = round(height(iconn.events)/2);

switch connmethod
    case {'coh','mscoh'}
        %-- Varianve Explained
        specs.plot.FigName      = sprintf(sprintf('Coherence (%s)',tarBAND));
        specs.plot.colorMap     = colormapskew('hot',0.8,[0 0.95]);
        specs.plot.AlphaData = elec_ok;
        if basecorr
%             maxconn = max(abs(iconn.connectivit) -  abs(connbase.connectivity),[],'all');
            maxconn = max(abs(abs(allconn)-abs(allbase)),[],'all');
            specs.plot.YLim      = [-1 1]* round(maxconn,2);
            pltdat      = abs(iconn.connectivity) - ...
                abs(connbase.connectivity);
            poscol = colormapskew('hot',0.8,[0 0.95]);
            negcol = flipud(poscol(:,[3,2,1]));
            specs.plot.colorMap     = cat(1,negcol,poscol);
        else
            specs.plot.YLim      = [0 1.0];
            pltdat      = abs(iconn.connectivity);
            specs.plot.colorMap     = colormapskew('hot',0.8,[0 0.95]);
        end
        pltbsl          = abs(connbase.connectivity);
end
if ~avgaround 
   avgdat               = pltdat;
   avgbsl               = pltbsl;
else
   specs.timelabel      = iconn.events.trial_name;
   avgdat = nan(height(iconn.channels),length(arounddist),height(iconn.events));
   avgbsl = nan(height(iconn.channels),length(arounddist),1);
   dd = 0;
   for arddist = arounddist
       dd = dd + 1;
       for ich=1:height(iconn.channels)
           nRows = 16;
           tarchnum  = str2double(strtok(iconn.channels.name{ich},whichHDgrid));
           aroundch  = reshape(unique([[nRows -nRows]'*arddist+[-arddist:arddist]; [(-nRows*arddist):nRows:(nRows*arddist)]+[-arddist;arddist]]),1,[]);
           avgchnum = tarchnum+aroundch;
           if mod(tarchnum-1-arddist,nRows)+1 > mod(tarchnum-1,nRows)+1   % when leftside
               avgchnum(mod(avgchnum-1,nRows)+1 >= mod(tarchnum-1-arddist,nRows)+1) = [];
           end
           if mod(tarchnum-1+arddist,nRows)+1 < mod(tarchnum-1,nRows)+1   % when rightside
               avgchnum(mod(avgchnum-1,nRows)+1 <= mod(tarchnum-1+arddist,nRows)+1) = [];
           end
           avgchnum(avgchnum<1 | avgchnum>height(iconn.channels)) = [];   % when top or bottom
           avgchname = unique(...
                         [arrayfun(@(E) sprintf('%s%d',whichHDgrid,E),avgchnum,'UniformOutput',false),...
                          arrayfun(@(E) sprintf('%s%02d',whichHDgrid,E),avgchnum,'UniformOutput',false),...
                          arrayfun(@(E) sprintf('%s%03d',whichHDgrid,E),avgchnum,'UniformOutput',false)],...
                       'stable');
           avgchidx  = ismember(iconn.channels.name,avgchname);
           avgdat(ich,dd,:) = mean(pltdat(avgchidx,ich,:),1,'omitnan');
           avgbsl(ich,dd,:) = mean(pltbsl(avgchidx,ich,:),1,'omitnan');
       end
   end
end

%% plot
% close all;

seedidx = find(ismember(iconn.channels.name,seedch),1);
figure;
h  = plot(permute(avgdat(seedidx,:,:),[3,2,1]),'LineWidth',2);
set(gca,'FontSize',20);
hold on;
hb = plot(repmat(permute(avgbsl(seedidx,:,:),[3,2,1]),height(iconn.events),1),':','LineWidth',2);
for ih=1:length(hb)
    hb(ih).Color = plcol(ih,:);
end

xlim([1 height(iconn.events)]);
xticks([]);
xlabel(stim);
ylim(ylim+[0 diff(ylim)*0.15]);
legend(h,arrayfun(@(x) sprintf('distance %d',x),arounddist,'UniformOutput',false));
title(specs.plot.FigName)


%% save
if setmask
    fignamebase = sprintf('%sGrid-ts_%s-%02d%%-%02d%%-ecc%02d_%s-%s',connmethod,subject,threshold_bb,threshold_a,eclimit,tarBAND,stim);
else
    fignamebase = sprintf('%sGrid-ts_%s_%s-%s',connmethod,subject,tarBAND,stim);
end
if basecorr
    fignamebase = sprintf('%s-basecorr',fignamebase);
end
if avgaround
    fignamebase = sprintf('%s_AvgRound',fignamebase);
end
figname     = sprintf('%s_ref-%s',fignamebase,seedch);
    
% figname = sprintf('%s_%s',figname,iconn.events.trial_name{specs.plot.Time});
saveas(gcf,fullfile(plotsavedir, fignamebase));
hgexport(gcf,fullfile(plotsavedir, figname),hgexport('factorystyle'),'Format','png');

end
end
end



%% plot Coherence average across around electrodes
close all;

flds    = {'mscoh'};
bands   = {'broadband','alpha'};
setmask = false;

isbj = 1;

if setmask
    elec_ok = ~(prf_params_bb{isbj}.xval<=threshold_bb | prf_params_a{isbj}.xval<=threshold_a ...
         | prf_params_bb{isbj}.ecc >= eclimit | prf_params_a{isbj}.ecc >= eclimit); 
end

if isbj==1
seedch = 'GB102';
else
seedch = 'GB066';
end
istim  = 1;

%-- main
iprf     = conn{isbj};
subject  = iprf.subject;
seedidx  = find(strcmp(iprf.channels.name,seedch),1);
stimname = iprf.events.trial_name(istim);
if ~exist(fullfile(plotsavedir,stimname),'dir'), mkdir(fullfile(plotsavedir,stimname)); end

for tarBAND = bands

switch tarBAND{:}
    case {'bb','broadband'}
        thresh  = threshold_bb;
        ifreq  = 90;
    case {'a','alpha'}
        thresh  = threshold_a;
        ifreq  = 13;
end

for sbfld = flds

specs = [];
specs.channels          = iprf.channels;
specs.plot.nSubPlots    = [];
specs.plot.RotGrid      = true;
specs.plot.fontSize     = 16;
specs.plot.FigName      = sprintf('pRF property');

switch sbfld{:}
    case {'coh'}
        %-- Varianve Explained
        specs.plot.YLim         = [0 1.0];
        specs.plot.colorMap     = colormapskew('hot',0.8,[0 0.95]);
        specs.timelabel         = {sprintf('Coherence (%s)',tarBAND{:})};
        pltdat      = abs(iprf.connectivity(:,seedidx,istim,ifreq));
end
if setmask
    specs.plot.AlphaData = elec_ok;
end

hF = ecog_plotGridSC(pltdat, whichHDgrid, specs);

if setmask
    figname = sprintf('prfGrid_%s-%02d%%-%02d%%-ecc%02d_%s-%s_ref-%s',subject,threshold_bb,threshold_a,eclimit,tarBAND{:},sbfld{:},seedch);
else
    figname = sprintf('prfGrid_%s_%s-%s_ref-%s',subject,tarBAND{:},sbfld{:},seedch);
end
if numel(hF)==1
    hgexport(hF{1},fullfile(plotsavedir,stimname, figname),hgexport('factorystyle'),'Format','png');
else
    for ii=1:numel(hF)
    hgexport(hF{ii},fullfile(plotsavedir,stimname, sprintf('%s-%d',figname,ii)),hgexport('factorystyle'),'Format','png');
    end
end

end
end






%% %%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%

%% Load pRF parameters
average        ='runs';
smoothingMode  ='decimate';
smoothingN     = 3;
prfmodel       ='linear';
gaussianmode   ='gs';

opts = [];
opts.average        = average;
opts.smoothingMode  = smoothingMode;
opts.smoothingN     = smoothingN;
opts.prfmodel       = prfmodel;
opts.gaussianmode   = gaussianmode;
opts.issave         = false;
opts.compute        = false;

opts.targetBAND     ='FaCLb';
prf_params_a = ecog_prf_analyzePRF(subjectList, opts);
opts.targetBAND     ='bbS';
prf_params_bb = ecog_prf_analyzePRF(subjectList, opts);

%% plot xR2 distribution for specific stimulus

isbj = 1;

thresh_bb = 31;
thresh_a  = 15;
eclimit   = 50;
angrange  = 10;
angtarg   = 45;     % [0, 45, 90, 135]

refparam = prf_params_bb{isbj};
pickchs_bb = (mod(refparam.ang,360) < (angrange/2+angtarg) & mod(refparam.ang,360) > (-angrange/2+angtarg)) | ...
             (mod(refparam.ang,360) < (360+angrange/2+angtarg) & mod(refparam.ang,360) > (360-angrange/2+angtarg)) | ...
             (mod(refparam.ang,360) < (180+angrange/2+angtarg) & mod(refparam.ang,360) > (180-angrange/2+angtarg));
okchs_bb = refparam.xval > thresh_bb & refparam.ecc < eclimit;

refparam = prf_params_a{isbj};
pickchs_a = (mod(refparam.ang,360) < (angrange/2+angtarg) & mod(refparam.ang,360) > (-angrange/2+angtarg)) | ...
             (mod(refparam.ang,360) < (360+angrange/2+angtarg) & mod(refparam.ang,360) > (360-angrange/2+angtarg)) | ...
             (mod(refparam.ang,360) < (180+angrange/2+angtarg) & mod(refparam.ang,360) > (180-angrange/2+angtarg));
okchs_a = refparam.xval > thresh_a & refparam.ecc < eclimit;

find(pickchs_bb & okchs_bb)
find(pickchs_a & okchs_a)
find(pickchs_bb & okchs_bb & pickchs_a & okchs_a)

        
        
