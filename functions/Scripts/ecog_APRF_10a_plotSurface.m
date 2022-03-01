% plot mesh surface

% 20210203 Yuasa 
%%
% close all;
checkPath;

fsDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/freesurfer/';

%% run script for a specific patient (with plotmesh on) to look at electrodes on brain
close all;

% atlasNames      =  {'wang15_fplbl_norm_5','wang15_fplbl_norm','benson14_varea','wang15_mplbl'};
atlasNames      =  {'wang15_fplbl_norm_5'};
summarizedAtlas = 'yes';

plotsavePATH = fullfile(analysisRootPath,'Figures','Surface_representative');

% subject = 'chaam';
% subject = 'som708';
% subject = 'som726';
% subject = 'som748';
for subject= ["beilen" "chaam" "som674" "som692" "som708" "som718" "som723" "som726" "som748"]
% for subject= ["chaam" "som708" "som726" "som748"]

switch subject
    case {'chaam','som661','som674','som718','som726'} % right hemisphere
        elecmesh = 'right';
    case {'beilen','som648','som692','som708','som748','som799','som800'} % left hemisphere
        elecmesh = 'left';
    case{'som723'} % bilateral
        elecmesh = 'both';
    otherwise
        error('''%s'' is unknown subject', subject);
end

if strcmpi(summarizedAtlas,'yes')
    plotsaveDir = fullfile(plotsavePATH,'summarize');
else
    plotsaveDir = plotsavePATH;
end
if ~exist(plotsaveDir,'dir'), mkdir(plotsaveDir); end

%% plot surface

specs = [];

specs.pID           = subject; % patient ID 
specs.atlasNames    = atlasNames;
specs.addRECatlas   = 'no';
specs.plotmesh      = elecmesh;
specs.plotelecs     = 'yes';
specs.plotlabel     = 'no';
specs.plotmatchednodes = 'no';
specs.plotelecrad    = 1;
specs.fsDir         = fsDir;
specs.smrylabel     = summarizedAtlas;
visualelectrodes    = electrode_to_nearest_node(specs,1);

%% save Figures

hF = get(groot,'Children');

for ii = 1:(length(hF)-1)
    figname = strrep(get(hF(ii),'Name'),' ','_');
    saveas(hF(ii),fullfile(plotsaveDir,figname),'fig');
    
    %-- front
    view(hF(ii).Children(end),[0 0])
    savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-front']),{'png','eps'},'Renderer','auto');
    
    %-- right
    view(hF(ii).Children(end),[30 0])
    savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-right']),{'png','eps'},'Renderer','auto');
    
    %-- left
    view(hF(ii).Children(end),[-30 0])
    savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-left']),{'png','eps'},'Renderer','auto');
    
    switch subject
        case 'chaam'
        %-- left2
        view(hF(ii).Children(end),[-45 0])
        savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-left2']),{'png','eps'},'Renderer','auto');
        
        %-- left3
        view(hF(ii).Children(end),[-90 0])
        savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-left3']),{'png','eps'},'Renderer','auto');
    end
            
end    

close(hF)


%% plot label

specs = [];

specs.pID           = subject; % patient ID 
specs.atlasNames    = atlasNames;
specs.addRECatlas   = 'no';
specs.plotmesh      = elecmesh;
specs.plotelecs     = 'yes';
specs.plotlabel     = 'yes';
specs.plotmatchednodes = 'no';
specs.plotelecrad    = 1;
specs.labelsize     = 14;
specs.fsDir         = fsDir;
specs.smrylabel     = summarizedAtlas;
visualelectrodes    = electrode_to_nearest_node(specs,1);

%% save Figures

hF = get(groot,'Children');

for ii = 1:(length(hF)-1)
    %-- label
    figname = strrep(get(hF(ii),'Name'),' ','_');
    savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-label']));
    
end    

close(hF)

end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% large dot
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run script for a specific patient (with plotmesh on) to look at electrodes on brain
close all;

% subject = 'chaam';
% subject = 'som708';
% subject = 'som726';
% subject = 'som748';
for subject= ["beilen" "chaam" "som674" "som692" "som708" "som718" "som723" "som726" "som748"]
% for subject= ["chaam" "som708" "som726" "som748"]

switch subject
    case {'chaam','som661','som674','som718','som726'} % right hemisphere
        elecmesh = 'right';
    case {'beilen','som648','som692','som708','som748','som799','som800'} % left hemisphere
        elecmesh = 'left';
    case{'som723'} % bilateral
        elecmesh = 'both';
    otherwise
        error('''%s'' is unknown subject', subject);
end

plotsaveDir = fullfile(analysisRootPath,'Figures','Surface_representative','summarize-largedot');
if ~exist(plotsaveDir,'dir'), mkdir(plotsaveDir); end

%% plot surface

specs = [];

specs.pID           = subject; % patient ID 
specs.atlasNames    = atlasNames;
specs.addRECatlas   = 'no';
specs.plotmesh      = elecmesh;
specs.plotelecs     = 'yes';
specs.plotlabel     = 'no';
specs.plotmatchednodes = 'no';
specs.plotelecrad    = 2;
specs.fsDir         = fsDir;
specs.smrylabel     = summarizedAtlas;
visualelectrodes    = electrode_to_nearest_node(specs,1);

%% save Figures

hF = get(groot,'Children');

for ii = 1:(length(hF)-1)
    figname = strrep(get(hF(ii),'Name'),' ','_');
    saveas(hF(ii),fullfile(plotsaveDir,figname),'fig');
    
    %-- front
    view(hF(ii).Children(end),[0 0])
    savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-front']),{'png','eps'},'Renderer','auto');
    
    %-- right
    view(hF(ii).Children(end),[30 0])
    savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-right']),{'png','eps'},'Renderer','auto');
    
    %-- left
    view(hF(ii).Children(end),[-30 0])
    savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-left']),{'png','eps'},'Renderer','auto');
    
    switch subject
        case 'chaam'
        %-- left2
        view(hF(ii).Children(end),[-45 0])
        savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-left2']),{'png','eps'},'Renderer','auto');
        
        %-- left3
        view(hF(ii).Children(end),[-90 0])
        savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-left3']),{'png','eps'},'Renderer','auto');
    end
            
end    

close(hF)


%% plot label

specs = [];

specs.pID           = subject; % patient ID 
specs.atlasNames    = atlasNames;
specs.addRECatlas   = 'no';
specs.plotmesh      = elecmesh;
specs.plotelecs     = 'yes';
specs.plotlabel     = 'yes';
specs.plotmatchednodes = 'no';
specs.plotelecrad    = 2;
specs.labelsize     = 14;
specs.fsDir         = fsDir;
specs.smrylabel     = summarizedAtlas;
visualelectrodes    = electrode_to_nearest_node(specs,1);

%% save Figures

hF = get(groot,'Children');

for ii = 1:(length(hF)-1)
    %-- label
    figname = strrep(get(hF(ii),'Name'),' ','_');
    savefigauto(hF(ii),fullfile(plotsaveDir,[figname '-label']));
    
end    

close(hF)

end

%% %%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%
%% Color Label for Wang
%% %%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(summarizedAtlas,'yes')
area_name =    {'V1','V2','V3', ...
                'hV4','VO','PHC', ...
                'TO','LO','V3AB', ...
                'IPS','SPL','FEF'}; 
ColorOrder  = [255   0   0; 255 128   0; 255 255   0;
                 0 128 255;   0 102  51; 153 153   0;
               153  76   0; 102   0 204;   0 255   0;
               255 153 153; 255 153 255; 255 178 102]./255;
else
area_labels =  {'V1v','V1d', ...
                'V2v','V2d', ...
                'V3v','V3d', ...
                'hV4', ...
                'VO1','VO2', ...
                'PHC1','PHC2', ...
                'TO2','TO1', ...
                'LO2','LO1', ...
                'V3b','V3a', ...
                'IPS0','IPS1','IPS2','IPS3','IPS4','IPS5', ...
                'SPL1','FEF'}; 
ColorOrder  = [102  51   0; 255   0   0; 
               102   0  51; 255 128   0;
                51   0 102; 255 255   0;
                 0 128 255;
                 0 102  51; 153 153   0;
               204 102   0; 204   0   0;
                 0  76 153; 153  76   0;
               255  51 153; 102   0 204;
                 0 255 255;   0 255   0;
               255 153 153; 255 204 153; 255 255 153; 153 255 153; 153 255 255; 153 153 255; 
               255 153 255; 255 178 102]./255;
end

%%
if contains(atlasNames,'wang15')
figure;
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

    figname = 'ColorLabel_wang15';
    if strcmpi(summarizedAtlas,'yes')
        figname = [figname '-summary'];
    end
    savefigauto(gcf,fullfile(plotsavePATH,[figname]));
end
    
%% %%%%%%%%%%%%%%%%%%%%%%%
%% Color Label for Benson
%% %%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(summarizedAtlas,'yes')
area_name =    {'V1','V2','V3', ...
                'hV4','VO', ...
                'LO','TO','V3AB'}; 
ColorOrder  = [255   0   0; 255 128   0; 255 255   0;
                 0 128 255;   0 102  51;
               102   0 204; 153  76   0;   0 255   0]./255;
else
area_labels = {'V1', ...
               'V2', ...
               'V3', ...
               'hV4', ...
               'VO1', 'VO2', ...
               'LO1', 'LO2', ...
               'TO1', 'TO2', ...
               'V3b', 'V3a'}; 
ColorOrder  = [255   0   0; 
               255 128   0; 
               255 255   0; 
                 0 128 255;
                 0  76 153; 153  76   0;
                 0 102  51; 153 153   0;  
               255  51 153; 102   0 204;
                 0 255 255;   0 255   0]./255;
end
%%
if contains(atlasNames,'benson14')
figure;
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

    figname = 'ColorLabel_benson14';
    if strcmpi(summarizedAtlas,'yes')
        figname = [figname '-summary'];
    end
    savefigauto(gcf,fullfile(plotsavePATH,[figname]));
end
    