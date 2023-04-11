% Plot coherence between electrodes in HD grid across distance
%   for all frequency range

% 20230329 Yuasa

%% %%%%%%%%%%%%%%%%%%%%
%% test
%% %%%%%%%%%%%%%%%%%%%%
%% prefix
% close all; clear all;
% startupToolboxToolbox;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true); 
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'Connectivity-representative');
SetDefault('prfPth',        'pRFmodel');
SetDefault('prfstatPth',    'pRFanalysis');
SetDefault('xspctrmPth',    'xSpectrum');
else
plotsavePth    = 'Connectivity-representative';
prfPth         = 'pRFmodel';
prfstatPth     = 'pRFanalysis';
xspctrmPth     = 'xSpectrum';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth));
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
HDsubjectList = SetSubjectsList(subjectList_fname, 'hasHDgrid','yes');

%-- Plotting Setting
FntSiz = 22;

%%
cohmethod = 'mscoh';                 % 'mscoh', 'imcoh'             % needed for loading channels selection
disttype  = 'norm';                  % square','diamond','norm'     % needed for loading channels selection
SetDefault('useChans','SELchs');     % 'pRFchs', 'SELchs', 'ALLchs'
SetDefault('iswideWin',false);       % if use wide window
arounddist  = 1:6;

if issaveplot,  plsbjs = 1:(length(HDsubjectList)+1);
else,           plsbjs = 1:length(HDsubjectList);
end
    
hF = gobjects(0);
for selsbj = plsbjs
%-- Dataset specs
if selsbj > length(HDsubjectList),  subjectList = HDsubjectList;
else,                               subjectList = HDsubjectList(selsbj);
end

%% load coherence
%%% Subject Name
nsbj = length(subjectList);
if nsbj==1
   subject = subjectList{1};
else
   subject = 'all';
end

if iswideWin,   dattypeid = '-Wide';
else,           dattypeid = '';
end

%%% Coherence file
opts = [];
opts.outputDir      = 'xSpectrum';
opts.compute        = false;
opts.targetBAND     = 'all';
opts.calcmode       = 'coh';
if iswideWin
opts.fileid         = 'freq_coherenceWide';
end
opts.stimNames      = 'HORIZONTAL*';
[coh_all1] = ecog_prf_crossspectra(subjectList, opts);
opts.stimNames      = 'VERTICAL*';
[coh_all2] = ecog_prf_crossspectra(subjectList, opts);
opts.stimNames      = 'DIAGONAL*';
[coh_all3] = ecog_prf_crossspectra(subjectList, opts);
opts.stimNames      = 'BLANK';
[coh_all0] = ecog_prf_crossspectra(subjectList, opts);
opts.stimNames      = 'SHUFFLE';
[coh_sfl] = ecog_prf_crossspectra(subjectList, opts);

%%% channel selection
filename = sprintf('%sdat%s_%s-%s-%s.mat',cohmethod,dattypeid,subject,useChans,disttype);
filepath = fullfile(SetDefaultAnalysisPath('DAT',xspctrmPth),filename);
assert(exist(filepath,'file'),'Need to run %s in advance.','ecog_APRFF_03d_CoherenceAcrossDist');
load(filepath,'channels','selch');

%%% concatenate coherence
chancmb = cell(1,length(selsbj));
chandist = cell(1,length(selsbj));
coh_avg_blnk = cell(1,length(selsbj));
coh_avg_stim = cell(1,length(selsbj));
coh_avg_all  = cell(1,length(selsbj));
coh_avg_shfl = cell(1,length(selsbj));
for isbj = 1:length(coh_all0)

    %-- select frequency
    f_coh    = coh_all0{isbj}.f; % f_show = true(size(f_coh));
    f_show   = f_coh>=1 & f_coh<=200;
    f_coh    = f_coh(f_show);
    
    %-- select channels (select seed-channel satisfying threshold & exclude the same electrode-pair)
    sbjch      = ismember(channels.subject_name,coh_all0{isbj}.subject);
    selchidx   = find(selch(sbjch));
    selchancmb = cat(1,find(ismember(coh_all0{isbj}.chancmb(:,1),selchidx)),...
                       find(ismember(coh_all0{isbj}.chancmb(:,2),selchidx)));
    selchancmb(diff(coh_all0{isbj}.chancmb(selchancmb,:),1,2)==0) = []; % exclude auto-coherence
    chancmb{isbj}  = coh_all0{isbj}.chancmb(selchancmb,:);

    %-- get channel pair and distance
    chanintvl  = 3;     % 3mm
    chanpos = coh_all0{isbj}.channels(:,{'name','x','y','z'});
    chanpos{:,["grid_x","grid_y"]} = ...
            [mod(str2double(strrep(chanpos.name,'GB',''))-1,16)+1,...
             ceil(str2double(strrep(chanpos.name,'GB',''))/16)].*chanintvl;
%     chandist{isbj} = vecnorm(squeeze(diff(reshape(...
%             chanpos{permute(chancmb{isbj},[2,3,1]),["x","y","z"]}...
%             ,2,[],3),1,1)),2,2);
    chandist{isbj} = vecnorm(squeeze(diff(reshape(...
            chanpos{permute(chancmb{isbj},[2,3,1]),["grid_x","grid_y"]}...
            ,2,[],2),1,1)),2,2);
    
    %-- get averaged coherence
    coh_avg0 = mean(coh_all0{isbj}.crossspectra(selchancmb,:,f_show),2,'omitnan'); 
        n_avg0 = sum(~isnan(coh_all0{isbj}.crossspectra(selchancmb,:,f_show)),2);
    coh_avg1 = mean(coh_all1{isbj}.crossspectra(selchancmb,:,f_show),2,'omitnan'); 
        n_avg1 = sum(~isnan(coh_all1{isbj}.crossspectra(selchancmb,:,f_show)),2);
    coh_avg2 = mean(coh_all2{isbj}.crossspectra(selchancmb,:,f_show),2,'omitnan'); 
        n_avg2 = sum(~isnan(coh_all2{isbj}.crossspectra(selchancmb,:,f_show)),2);
    coh_avg3 = mean(coh_all3{isbj}.crossspectra(selchancmb,:,f_show),2,'omitnan'); 
        n_avg3 = sum(~isnan(coh_all3{isbj}.crossspectra(selchancmb,:,f_show)),2);
        n_avgS = n_avg1 + n_avg2 + n_avg3;
    coh_avgS = sum(cat(2,coh_avg1.*n_avg1,coh_avg2.*n_avg2,coh_avg3.*n_avg3),2,'omitnan') ./ n_avgS;
        n_avgA = n_avg0 + n_avgS;
    coh_avgA = sum(cat(2,coh_avg0.*n_avg0,coh_avgS.*n_avgS),2,'omitnan') ./ n_avgA;
    coh_avgR = mean(coh_sfl{isbj}.crossspectra(selchancmb,:,f_show),2,'omitnan'); 
        
    coh_avg_blnk{isbj} = squeeze(coh_avg0);
    coh_avg_stim{isbj} = squeeze(coh_avgS);
    coh_avg_all{isbj}  = squeeze(coh_avgA);
    coh_avg_shfl{isbj} = squeeze(coh_avgR);
    
end
chancmb = cat(1,chancmb{:});
chandist = cat(1,chandist{:});
coh_avg_blnk = cat(1,coh_avg_blnk{:});
coh_avg_stim = cat(1,coh_avg_stim{:});
coh_avg_all  = cat(1,coh_avg_all{:});
coh_avg_shfl = cat(1,coh_avg_shfl{:});

%% %%%%%%%%%%%%%
%% Visualize
arounddist = 3:3:12;
distrange  = 3;
fardist    = 20;

if iswideWin,   figprefix = 'CoherenseW-allfreq';
else,           figprefix = 'Coherense-allfreq';
end

%% Plot coherence <Based on long-range>
n_dist = length(arounddist)+1;
hF(end+1) = figure;  hold on; 
colororder(flipud(copper(n_dist)));  xlim(minmax(f_coh'));

%-- set long-range coherence as reference
ch_far   = chandist>=fardist;
coh_ref  = mean(coh_avg_all(ch_far,:),1,'omitnan');

%-- plot coherence for each distance
for idist = arounddist
    ch_dist = chandist>=(idist-distrange/2) & chandist<(idist+distrange/2);
    plot(f_coh, (mean(coh_avg_all(ch_dist,:),1,'omitnan')-coh_ref));
end
hold off;
lgndnames = strsplit(sprintf('%dmm ',arounddist));
lgndnames(end) = [];
legend(lgndnames,'AutoUpdate','off');
title("difference from long-range (>2cm)")

%-- set figure properties
set(findobj(gca,'Type','Line'),'LineWidth',4);
set(findobj(gca,'-property','FontSize'),'FontSize',FntSiz);
set(gca,'XScale','log'); xticks(unique([xticks xlim])); xticklabels(xticks);
ylim([-0.02 max(ylim)]);

%-- add baseline (0)
hold on; hb = plot(xlim,[0 0],'k:','LineWidth',2); uistack(hb,'bottom'); hold off;

if issaveplot
savefigauto(gcf,fullfile(plotsavedir,sprintf('%s-DiffDistNoRatio-avg_%s',figprefix,subject)));
end


%% Plot coherence <Based on Shuffle>
if issaveplot
n_dist = length(arounddist)+2;
hF(end+1) = figure;  hold on; 
colororder(flipud(copper(n_dist)));  xlim(minmax(f_coh'));

%-- set shuffled coherence as reference
ch_far   = chandist>=fardist;
coh_ref  = mean(coh_avg_shfl,1,'omitnan');

%-- plot coherence for each distance
for idist = arounddist
    ch_dist = chandist>=(idist-distrange/2) & chandist<(idist+distrange/2);
    plot(f_coh, (mean(coh_avg_all(ch_dist,:),1,'omitnan')-coh_ref));
end
%-- plot coherence for long-range
plot(f_coh, (mean(coh_avg_all(ch_far,:),1,'omitnan')-coh_ref),'k');

hold off;
lgndnames = strsplit(sprintf('%dmm ',arounddist));
lgndnames{end} = sprintf('>%dmm',fardist);
legend(lgndnames,'AutoUpdate','off');
title("difference from Shuffled")

%-- set figure properties
set(findobj(gca,'Type','Line'),'LineWidth',4);
set(findobj(gca,'-property','FontSize'),'FontSize',FntSiz);
set(gca,'XScale','log'); xticks(unique([xticks xlim])); xticklabels(xticks);
ylim([-0.02 max(ylim)]);

%-- add baseline (0)
hold on; hb = plot(xlim,[0 0],'k:','LineWidth',2); uistack(hb,'bottom'); hold off;

if issaveplot
savefigauto(gcf,fullfile(plotsavedir,sprintf('%s-DiffShflNoRatio-avg_%s',figprefix,subject)));
end
end

end
