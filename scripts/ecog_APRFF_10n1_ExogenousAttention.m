%% Spatial profile of Exogenous attention

% 20230109
% 20230215 - update for makeFigure

% close all; clear all;
SetDefault('plotavg',true);
SetDefault('plotovl',true);
if ~plotavg,  plotovl = false;  end

hF = gobjects(0);

%% Downing 1985
% Data from
% Downing, C. and S. Pinker (1985). "The spatial structure of visual
% attention." Mechanisms of attention: Attention and performance XI:
% 171-187.

% Data from figure 8.2
% Rows are x-values from (1,10) to (10,1)
% Columns are data series from cue at (5 or 6) to cue at (1 or 10)
data = [...
 -32.5806 -23.5484 -15.8065 7.4194 48.7097
 -20.6452 -15.8065 -21.2903 -10.3226 40.0000
 5.4839 -20.0000 -17.4194 -17.4194 19.0323
 20.6452 4.5161 -10.3226 -21.2903 1.9355
 25.4839 25.4839 15.1613 6.7742 -30.3226
 41.9355 34.1935 41.2903 50.9677 29.3548
 42.9032 45.8065 42.5806 45.4839 45.4839
 49.6774 42.2581 42.5806 47.4194 44.1935
 49.3548 43.8710 50.6452 52.2581 60.6452
 46.1290 40.3226 30.3226 48.0645 60.6452];

% cue locations, in degrees
cue = 11.25:-2.5:1.25;  % str = cellfun(@num2str, num2cell(cue'), 'UniformOutput', false);

% target locations, in degrees
target = 11.25:-2.5:-11.25;

[cues, targets] = meshgrid(cue, target);

% flip cue for ascending
data = fliplr(data);
cues = fliplr(cues);
targets = fliplr(targets);

decval    = 5;
upscale   = 140;
targetsNormed = targets./cues;
if plotavg
xlin = linspace(min(targetsNormed(:)), max(targetsNormed(:)), decval*upscale);
xlog = logspace(log10(min(abs(targetsNormed(:)))), log10(max(targetsNormed(:))), round(decval*upscale./4));
xswc = xlog(find(diff(xlog)>=mean(diff(xlin)),1));
x = cat(2,xlin(xlin<xswc),xlog(xlog>=xswc));
targetsInterp = NaN(length(x), size(targetsNormed,2));
for ii = 1:size(targetsNormed,2)

targetsInterp(:,ii) = interp1(targetsNormed(:,ii), data(:,ii), x, 'linear', NaN);

end
mn = mean(targetsInterp,2, 'omitnan');
n = sum(isfinite(targetsInterp), 2);
se = std(targetsInterp,[], 2, 'omitnan')./sqrt(n);
end

%%%% Figure %%%
hF(end+1) = figure('Color', 'w');

if plotavg
%%% plot average
    %-- Set axis
    if ~plotovl
    set(gcf, 'Position', round(get(gcf,'Position').*[1/3,1,1,1]));
    end
    colororder([0,0,0]);
    hl = gobjects(0);
    
    %-- Plot
    xdec  = decimate(x,decval);
    mndec = decimate(mn,decval);
    
    hold on; box off;
    hl(end+1) = plot(xdec, smooth(mndec,9,'sgolay'), '-', 'LineWidth', 5);
    set(gca, 'FontSize', 20); %, 'XTick', 0:9)
    ylabel('Response time cost (ms)', 'FontSize',24)
    xlabel('Target eccentricity (normalized)', 'FontSize',24)
    
    xlim([0 5.5])
    % ylim([-25 30])
    plot(xlim, [0 0], '--', 'LineWidth',2);
    plot([1 1], ylim, ':', 'LineWidth',2);
    text(1.1, prctile(ylim,65), 'Cue location', 'FontSize', 20)
    
    ylinl = ylim;
    yminl = min(get(hl(end),'YData'));
    yrngl = ylinl-yminl;
else
%%% plot each
    ipsitarget = target>0;
    %-- Set axis
    set(gcf, 'Position', round(get(gcf,'Position').*[1/3,1/2,1,7/4]));
    hT = tiledlayout('flow','TileSpacing','compact','Padding','compact');
    colororder(fliplr(copper(size(cues,2)+0)));
    
    %-- Plot
    nexttile; hold on; box off;
        plot(targets(ipsitarget,:), data(ipsitarget,:), '-o', 'LineWidth', 3, 'MarkerSize', 6)
        xlim([0 max(targets,[],'all')+1])
        xticks(unique(targets));
        
        set(gca, 'FontSize', 20); %, 'XTick', 0:9)
        ylabel('Response time cost (ms)', 'FontSize',24)
        xlabel('Target eccentricity (deg)', 'FontSize',24)
        
        plot(xlim, [0 0], 'k--', 'LineWidth',2);
    
    nexttile; hold on; box off;
        plot(targetsNormed(ipsitarget,:), data(ipsitarget,:), '-o', 'LineWidth', 3, 'MarkerSize', 6)
        xlim([0 max(targetsNormed,[],'all')+0.5])
        
        set(gca, 'FontSize', 20); %, 'XTick', 0:9)
        ylabel('Response time cost (ms)', 'FontSize',24)
        xlabel('Target eccentricity (normalized)', 'FontSize',24)
        hl = legend(string(cues(1,:)),'AutoUpdate','off','Location','southeast');
        hl.Title.String = 'Cue (deg)';
        
        plot(xlim, [0 0], 'k--', 'LineWidth',2);
        plot([1 1], ylim, 'k:', 'LineWidth',2);
        text(1.1, prctile(ylim,65), 'Cue location', 'FontSize', 20)
    
    hT.Title.String = "Downing and Pinker (1985)";
    hT.Title.FontSize = 24;
end

%% Shulman 1986
% Data from
% Shulman, G. L., et al. (1986). "Gradients of spatial attention." Acta
% Psychol (Amst) 61(2): 167-181.

% Data from Table 2
% Rows are x-values from 0.5 to 24.5
% Columns are data series from cue at 0.5 to 24.5
data = [...
 0	19	23	29	35
 15	0	5	13	20
 29	12	0	2	5
 34	24	4	0	-3
 40	30	16	4	0]; % Mean cost (ms)
% % Another data from Table 1
% data0 = [...
%  299 318 322 328 334
%  319 304 309 317 324
%  343 325 313 315 318
%  369 359 340 335 332
%  363 353 338 326 323]; % Mean reaction time (ms)
% data = data0 - diag(data0); % Mean cost (ms)
 
% cue locations, in degrees
cue = (0:6:24)+0.5;  % str = cellfun(@num2str, num2cell(cue'), 'UniformOutput', false);

% target locations, in degrees
target = (0:6:24)+0.5;

[cues, targets] = meshgrid(cue, target);

% exclude cue or target location at 0.5 (too foveal)
data(:,1) = []; %data(1,:) = [];
cues(:,1) = []; %cues(1,:) = [];
targets(:,1) = []; %targets(1,:) = [];

decval    = 5;
upscale   = 150;
targetsNormed = targets./cues;
% x = linspace(min(targetsNormed(:)), max(targetsNormed(:)), decval*upscale);
xlin = linspace(min(targetsNormed(:)), max(targetsNormed(:)), decval*upscale);
xlog = logspace(log10(min(abs(targetsNormed(:)))), log10(max(targetsNormed(:))), round(decval*upscale./4));
xswc = xlog(find(diff(xlog)>=mean(diff(xlin)),1));
x = cat(2,xlin(xlin<xswc),xlog(xlog>=xswc));
targetsInterp = NaN(length(x), size(targetsNormed,2));
for ii = 1:size(targetsNormed,2)

targetsInterp(:,ii) = interp1(targetsNormed(:,ii), data(:,ii), x, 'linear', NaN);

end
mn = mean(targetsInterp,2, 'omitnan');
n = sum(isfinite(targetsInterp), 2);
se = std(targetsInterp,[], 2, 'omitnan')./sqrt(n);

%%%% Figure %%%
if ~(plotavg&&plotovl)
hF(end+1) = figure('Color', 'w');
end

if plotavg
%%% plot average
    %-- Set axis
    if plotovl
    colororder([0,0,0;0.65,0.50,0.05]);
    yyaxis right
    else
    set(gcf, 'Position', round(get(gcf,'Position').*[0,1,1,1]) + (get(hF(end-1),'Position')*[1,0,1,0]').*[1,0,0,0]);
    colororder([0,0,0]);
    end
    
    %-- Plot
    xdec  = decimate(x,decval);
    mndec = decimate(mn,decval);
    
    hold on;  box off;
    hl(end+1) = plot(xdec, smooth(mndec,9,'sgolay'), '.-.', 'LineWidth', 5);
    set(gca, 'FontSize', 20); %, 'XTick', 0:9)
    ylabel('Response time cost (ms)', 'FontSize',24)
    
    if plotovl
    ylinr = ylim;
    yminr = min(get(hl(end),'YData'));
    yrngr = ylinr-yminr;
    ylim(yrngl+yminr);
    else
    xlabel('Target eccentricity (normalized)', 'FontSize',24)
    xlim([0 4.0])
    plot([1 1], ylim, ':', 'LineWidth',2);
    text(1.1, prctile(ylim,65), 'Cue location', 'FontSize', 20);
    end


else
%%% plot each
    ipsitarget = target>0;
    %-- Set axis
    set(gcf, 'Position', round(get(gcf,'Position').*[0,1/2,1,7/4]) + (get(hF(end-1),'Position')*[1,0,1,0]').*[1,0,0,0]);
    hT = tiledlayout('flow','TileSpacing','compact','Padding','compact');
    if size(data,1)==size(data,2)
        colororder(fliplr(copper(size(cues,2)+0)));
    else
        colororder(circshift(fliplr(copper(size(cues,2)+1)),-1,1));
    end
    
    %-- Plot
    nexttile; hold on;
        plot(targets(ipsitarget,:), data(ipsitarget,:), '-o', 'LineWidth', 3, 'MarkerSize', 6)
        xlim([0 max(targets,[],'all')+1]);
        xticks(unique(targets));
        
        set(gca, 'FontSize', 20); %, 'XTick', 0:9)
        ylabel('Response time cost (ms)', 'FontSize',24)
        xlabel('Target eccentricity (deg)', 'FontSize',24)
        
    nexttile; hold on;
        plot(targetsNormed(ipsitarget,:), data(ipsitarget,:), '-o', 'LineWidth', 3, 'MarkerSize', 6)
        if min(cues,[],'all')<1
            xlim([0 4.5]);
        else
            xlim([0 max(targetsNormed,[],'all')+0.5])
        end
        
        set(gca, 'FontSize', 20); %, 'XTick', 0:9)
        ylabel('Response time cost (ms)', 'FontSize',24)
        xlabel('Target eccentricity (normalized)', 'FontSize',24)
        hl = legend(string(cues(1,:)),'AutoUpdate','off','Location','southeast');
        hl.Title.String = 'Cue (deg)';
        
        plot([1 1], ylim, 'k:', 'LineWidth',2);
        text(1.1, prctile(ylim,65), 'Cue location', 'FontSize', 20)
        
    hT.Title.String = "Shulman et al. (1986)";
    hT.Title.FontSize = 24;
end

%%% add legend
if plotavg&&plotovl
legend(hl,["Downing and Pinker (1985)","Shulman et al. (1986)"],'Location','best','AutoUpdate','off');
end
