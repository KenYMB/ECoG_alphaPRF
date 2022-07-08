% visualize pRF parameters on HD-Grid computed in 06c (equivalant 03c2 + 03d2)
%   Applicable xR2 with full time-series (no decimation)

% 20201111 Yuasa - update from test_Output_VisualAreas
% 20220621 Yuasa - Update for new environment

%% prefix
% close all; clearvars;
% startupToolboxToolbox;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'pRF_Grid-representative');
SetDefault('prfPth',        'pRFmodel');
else
plotsavePth    = 'pRF_Grid-representative';
prfPth         = 'pRFmodel';
end
%-- Set save figure dirctory
if issaveplot
plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),'Atlas');
if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
subjectList = SetSubjectsList(subjectList_fname, 'hasHDgrid', 'yes');

%-- Plotting Setting
FntSiz  = 18;

%%
%-- Channels and Grid Setting
projectDir      = bidsRootPath;
whichHDgrid     = 'GB';
gridSiz         = [16 8];
if issaveplot
atlasNames      = {'smry_wang15_fplbl_norm_5','wang15_fplbl_norm_5','smry_wang15_mplbl','smry_benson14_varea'};
else
atlasNames      = {'smry_wang15_fplbl_norm_5'};
end

%-- Plot
hF= gobjects(0);
for isbj = 1:length(subjectList)
    subject     = subjectList{isbj};
    
    for ii = 1:length(atlasNames)
        atlasName   = atlasNames{ii};
        
        [channels] = ...
            bidsEcogMatchElectrodesToAtlas(projectDir, subject, [], atlasName);
        [area_name, area_cmap, units] = getAtlasLabels(atlasName);
        
        chname = categorical(channels.(atlasName));
        chcate = categories(chname);
        
        %-- add none & get idx
        area_name = ['none',area_name];
        area_cmap = [0.8 0.8 0.8; area_cmap];
        [~,chcolr] = match_str(chname,area_name);
        
        %%
        if strcmpi(units,'area')
            
            %%% Plot Area Color
            figname = sprintf('VisualArea_%s_HDgrid-%s-%s',subject,whichHDgrid,atlasName);
            
            % close all;
            chidx = startsWith(channels.name,whichHDgrid);
            fillIndex = intersect(1:prod(gridSiz),  str2double(strtok(channels.name(chidx),whichHDgrid)));
            %-- Area Cluster
            AC = nan(gridSiz);
            AC(fillIndex) = chcolr(chidx);
            iflefthemi = mean(ismember(channels.hemisphere(chidx),'L')) > 0.7;
            if iflefthemi                   % left hemisphere
                rotAC = rot90(AC,-1);
            else                            % right hemisphere
                rotAC = rot90(AC);
            end
            hF(end+1) = figure('Name',figname,'Position',[200 200 900 420],'MenuBar','none');
            hA = imagesc(rotAC,'AlphaData',~isnan(rotAC));
            daspect([1 1 1]);
            set(gca, 'FontSize', FntSiz, 'XTick', [], 'YTick',[]);
            colormap(area_cmap); clim([1 length(area_name)+1]);
            hcb = colorbar;
            hcb.Ticks = (1:length(area_name))+0.5;
            hcb.TickLabels = area_name;
            
            if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');
            else,           hcb.Visible = 'off';
            end
            
            %%% Plot Border
            figname = sprintf('VisualAreaBorder_%s_HDgrid-%s-%s',subject,whichHDgrid,atlasName);
            
            %-- set boundary image
            imAC = imresize(rotAC,20,'nearest');
            imAC(isnan(imAC))=0;
            sizAC = size(imAC);
            bwAC = false(sizAC);
            for ii=1:numel(imAC)
                %-- get surround pixel index
                rndpix = ii + [[-1;0;1]+[-sizAC(1),0,sizAC(1)]];
                rndpix(5) = [];
                rndpix(rndpix<1 | rndpix>numel(imAC)) = [];
                if mod(ii,sizAC(1))==1
                    rndpix(mod(rndpix,sizAC(1))==0) = [];
                elseif mod(ii,sizAC(1))==0
                    rndpix(mod(rndpix,sizAC(1))==1) = [];
                end
                
                %-- detect edge & return
                bwAC(ii) = ~all(imAC(rndpix) == imAC(ii));
            end
            %-- resize boundary image to make the line bold
            bwAC = imresize(bwAC,5,'nearest');
            
            hF(end+1) = figure('Name',figname,'Position',[200 200 900 420],'MenuBar','none');
            hI = imshow(~bwAC);
            daspect([1 1 1]);
            set(gca, 'FontSize', FntSiz, 'XTick', [], 'YTick',[]);
            hF(end).Position = hF(end-1).Position;
            hI.Parent.Position = hA.Parent.Position;
            
            if issaveplot,  savefigauto(gcf,fullfile(plotsavedir, figname),'-vector');  end
            
            
        end
    end
end
