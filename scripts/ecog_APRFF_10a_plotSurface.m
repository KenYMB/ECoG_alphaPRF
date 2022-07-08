% plot mesh surface

% 20210203 Yuasa 
% 20220601 Yuasa - major update with bidsEcogPlotElectrodesOnMesh
%%
% close all; clearvars;
run_checkPath;
%-- Input & Output path
SetDefault('issaveplot',true);
SetDefault('closefig',true);
if exist('KEEPCURRENTPATH','var')&&KEEPCURRENTPATH
SetDefault('plotsavePth',   'Surface');
else
plotsavePth    = 'Surface';
end

%-- Subjects
SetDefault('subjectList_fname','subjectlist.tsv');
SetDefault('selsbj','all');
subjectList = SetSubjectsList(subjectList_fname, selsbj);

%% Plot surface with ECoG electrodes

projectDir      = bidsRootPath;
if issaveplot
atlasNames      =  {'smry_wang15_fplbl_norm_5','wang15_fplbl_norm_5','smry_benson14_varea','smry_wang15_mplbl'};
SetDefault('plotelecrad',[1 2]);
iscolorbar      = 'yes';
else
atlasNames      =  {'smry_wang15_fplbl_norm_5'};
SetDefault('plotelecrad',1);
closefig        = false;
iscolorbar      = 'no';
end

hF= gobjects(0);
for elecrad = reshape(plotelecrad,1,[])
for subject = string(reshape(subjectList,1,[]))

%% plot surface

specs = [];
specs.plotelecs         = 'yes';
specs.plotlabel         = 'no';
specs.plotmatchednodes  = 'no';
specs.plotelecrad       = elecrad;
specs.plotcbar          = iscolorbar;
specs.adjustLRdistance  = false;
hFc = bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasNames, [], specs);
hF = [hF hFc];

    %% save Figures
    
if issaveplot

    %-- Set save figure dirctory
    plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),sprintf('elecrad%d',elecrad));
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

    for ii = reshape(find(ishandle(hFc)),1,[])
        figname = strrep(get(hFc(ii),'Name'),' ','_');
        saveas(hFc(ii),fullfile(plotsavedir,figname),'fig');

        %-- front
        view(hFc(ii).Children(end),[0 0]);
        savefigauto(hFc(ii),fullfile(plotsavedir,[figname '-front']),{'png','eps'},'Renderer','auto');

        %-- right
        view(hFc(ii).Children(end),[30 0]);
        savefigauto(hFc(ii),fullfile(plotsavedir,[figname '-right']),{'png','eps'},'Renderer','auto');

        %-- left
        view(hFc(ii).Children(end),[-30 0]);
        savefigauto(hFc(ii),fullfile(plotsavedir,[figname '-left']),{'png','eps'},'Renderer','auto');

    end

    if closefig,  close(hFc);   end

end


%% plot w/ label

if issaveplot

specs = [];
specs.plotelecs         = 'yes';
specs.plotlabel         = 'yes';
specs.plotmatchednodes  = 'no';
specs.plotelecrad       = elecrad;
specs.plotcbar          = iscolorbar;
specs.adjustLRdistance  = false;
specs.labelsize         = 14;
hFc = bidsEcogPlotElectrodesOnMesh(projectDir, subject, [], atlasNames, [], specs);
hF = [hF hFc];

    %% save Figures

    %-- Set save figure dirctory
    plotsavedir    = fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),sprintf('elecrad%d-label',elecrad));
    if ~exist(plotsavedir,'dir'), mkdir(plotsavedir); end

    for ii = 1:reshape(find(ishandle(hFc)),1,[])
        %-- label
        figname = strrep(get(hFc(ii),'Name'),' ','_');
        savefigauto(hFc(ii),fullfile(plotsavedir,[figname '-label']),'Renderer','auto');

    end    

    if closefig,  close(hFc);   end

end

end
end


%% %%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%
%% Color Labels
%% %%%%%%%%%%%%%%%%%%%%%%%

if issaveplot
    
    
for ii = 1:length(atlasNames)
    [atlasLabel,~,~,issmry] = interpretAtlasNames(atlasNames{ii});
    [area_name, ColorOrder, units] = getAtlasLabels(atlasNames{ii});
    
    if strcmpi(units,'area')
        %-- Set Name
        atlasType = strtok(atlasLabel,'_');
        if issmry,      postfix = '-summary';
        else,           postfix = '';
        end
        figname = sprintf('ColorLabel_%s%s',atlasType,postfix);
        
        %-- Plot colorlabel
        figure('Name',figname,'MenuBar','none');
        meshres = 300;
        nroi    = length(area_name);
        [X,Y] = meshgrid(0:meshres,0:meshres/nroi);
        Z = X;
        surf(X,Y,Z,'EdgeColor','none');
        axis equal;
        colormap(ColorOrder);
        set(gcf, 'InvertHardcopy', 'off', 'Color', [1 1 1], 'Position', [200 100 100 (nroi+5)*17]);
        view([90 0]); box on;
        set(gca,'XTick',[],'YTick',[]);
        set(gca,'ZTick',linspace(0,meshres,size(ColorOrder,1)+1)+meshres/size(ColorOrder,1)/2,...
            'ZTickLabel',area_name,...
            'FontSize',18);
        Zax = get(gca,'ZRuler');  Zax.FirstCrossoverValue = Inf;
        
        set(gca,'Position',get(gca,'Position').*[-0.5 .3 1 1.15]);
        Zax.TickDirection = 'in';  Zax.TickLength = [0 0];

        savefigauto(gcf,fullfile(SetDefaultAnalysisPath('FIG',plotsavePth),[figname]));
    end
end


end